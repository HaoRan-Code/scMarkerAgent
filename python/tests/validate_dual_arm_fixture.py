#!/usr/bin/env python3
"""Compare a Python-arm and an R-arm run of the matched synthetic fixture.

SCOPE LIMIT, read before running: this script asserts DIGIT-FOR-DIGIT equality of the
two arms. That was achievable only while the R arm could reuse the Python arm's frozen
normalization and PCA matrices. Those were removed deliberately -- the project requires
the two arms to run the SAME ALGORITHM, not to reproduce each other bit for bit -- so
every assertion downstream of clustering (cluster ids, per-cell candidates, cluster
summary) can now legitimately differ.

What still holds, and what this script is useful for, is the format and gene-universe
contract: identical cells, identical measured genes, identical candidate menu and panel
keys, and identical cluster_summary COLUMNS. Treat a failure of a clustering-dependent
key as information, not as a regression, and compare partitions with ARI instead.

Each arm is read in its own native serialization: the Python caches as pickles, the R
caches as .rds through a one-shot Rscript export. There is no interchange bridge --
that is the point of the two-arm design, not an inconvenience to it.
"""

from __future__ import annotations

import argparse
import json
import pickle
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from scmarkeragent.config import OUTPUT_SCHEMA

NUMERIC_TOLERANCE = {
    "auc": 1e-10,
    "avgExpr": 1e-8,
    "logFC": 1e-8,
    "pct_in": 1e-10,
    "pct_out": 1e-10,
    "pval": 1e-8,
    "padj": 1e-8,
}

ANNOTATION_FIELDS = [
    "annotation",
    "subtype",
    "state",
    "resolution_status",
    "annotation_source",
    "llm_status",
    # Two arms that disagree on how a label left quality control are two different
    # algorithms, whatever the label itself says.
    "annotation_qc",
]

R_EXPORT = r"""
arguments <- commandArgs(TRUE)
cache <- arguments[1]
tag <- arguments[2]
out <- arguments[3]
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
artifact <- function(kind) readRDS(file.path(cache, sprintf("%s_%s.rds", tag, kind)))
dm <- artifact("de_meta")
scoring <- artifact("candidate_scoring")
annotations <- artifact("annotations")
fields <- c(
  "annotation", "subtype", "state", "resolution_status",
  "annotation_source", "llm_status", "annotation_qc"
)
results <- lapply(annotations$results, function(record) {
  values <- lapply(fields, function(field) as.character(record[[field]] %||% ""))
  names(values) <- fields
  values
})
de_columns <- intersect(
  c("feature", "group", "auc", "avgExpr", "logFC", "pct_in", "pct_out", "pval", "padj"),
  names(dm$de)
)
panels <- scoring$measured_panels
panel_keys <- unlist(lapply(names(panels), function(candidate) {
  vapply(
    as.character(unlist(panels[[candidate]])),
    function(gene) paste(candidate, gene, sep = "\x1f"),
    ""
  )
}), use.names = FALSE)
payload <- list(
  cells = as.list(as.character(dm$meta$cell)),
  clusters = as.list(as.character(dm$meta$cluster)),
  menu_genes = as.list(as.character(dm$menu_genes)),
  de = as.data.frame(dm$de)[, de_columns, drop = FALSE],
  candidates = as.list(sort(unique(as.character(scoring$scored$candidate)))),
  panel_keys = as.list(sort(panel_keys)),
  results = results
)
jsonlite::write_json(payload, out, digits = NA, auto_unbox = TRUE)
"""


def load_pickle(path: Path):
    with path.open("rb") as handle:
        return pickle.load(handle)


def load_r_arm(cache: Path, tag: str) -> dict:
    """The R arm's artifacts, exported by R itself rather than re-parsed by guesswork."""
    with tempfile.TemporaryDirectory() as scratch:
        script = Path(scratch) / "export.R"
        script.write_text(R_EXPORT, encoding="utf-8")
        payload = Path(scratch) / "r_arm.json"
        subprocess.run(
            ["Rscript", str(script), str(cache), tag, str(payload)],
            check=True,
        )
        return json.loads(payload.read_text(encoding="utf-8"))


def dataframe_difference(
    left: pd.DataFrame,
    right: pd.DataFrame,
    keys: list[str],
    values: list[str],
) -> dict:
    merged = left.merge(
        right, on=keys, suffixes=("_python", "_r"), how="outer", indicator=True
    )
    if not merged["_merge"].eq("both").all():
        return {"row_membership_equal": False}
    differences = {}
    for value in values:
        a = pd.to_numeric(merged[f"{value}_python"], errors="coerce")
        b = pd.to_numeric(merged[f"{value}_r"], errors="coerce")
        valid = a.notna() | b.notna()
        delta = (a - b).abs()
        differences[value] = float(delta[valid].max()) if valid.any() else 0.0
    return {"row_membership_equal": True, "max_abs_difference": differences}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--python-cache", required=True)
    parser.add_argument("--r-cache", required=True)
    parser.add_argument("--python-tag", required=True)
    parser.add_argument("--r-tag", required=True)
    parser.add_argument("--python-results", required=True)
    parser.add_argument("--r-results", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    pc = Path(args.python_cache)

    py_meta = load_pickle(pc / f"{args.python_tag}_de_meta.pkl")
    py_scoring = load_pickle(pc / f"{args.python_tag}_candidate_scoring.pkl")
    py_annotations = load_pickle(pc / f"{args.python_tag}_annotations.pkl")
    r_arm = load_r_arm(Path(args.r_cache), args.r_tag)

    result = {
        "schema_version": "dual-arm-fixture-parity-v2",
        "tolerances": NUMERIC_TOLERANCE,
        "cells_equal": (
            py_meta["meta"]["cell"].astype(str).tolist() == list(r_arm["cells"])
        ),
        "clusters_equal": (
            py_meta["meta"]["cluster"].astype(str).tolist() == list(r_arm["clusters"])
        ),
        "genes_equal": (
            [str(gene) for gene in py_meta["menu_genes"]] == list(r_arm["menu_genes"])
        ),
    }
    r_de = pd.DataFrame(r_arm["de"])
    de_values = [
        column
        for column in NUMERIC_TOLERANCE
        if column in py_meta["de"].columns and column in r_de.columns
    ]
    result["differential_expression"] = dataframe_difference(
        py_meta["de"], r_de, ["feature", "group"], de_values
    )
    result["candidate_menu_equal"] = sorted(
        set(py_scoring["scored"]["candidate"].astype(str))
    ) == list(r_arm["candidates"])
    py_panel_keys = sorted(
        f"{candidate}\x1f{gene}"
        for candidate, genes in py_scoring["measured_panels"].items()
        for gene in genes
    )
    result["candidate_panel_keys_equal"] = py_panel_keys == list(r_arm["panel_keys"])

    annotation_diffs = []
    for cluster in sorted(py_annotations["results"]):
        left = py_annotations["results"][cluster]
        right = r_arm["results"].get(str(cluster), {})
        for field in ANNOTATION_FIELDS:
            if str(left.get(field) or "") != str(right.get(field) or ""):
                annotation_diffs.append(
                    {
                        "cluster": cluster,
                        "field": field,
                        "python": left.get(field),
                        "r": right.get(field),
                    }
                )
    result["annotation_differences"] = annotation_diffs

    # Read from the schema rather than restating it: a second copy of the column list
    # goes stale silently.
    table_fields = list(OUTPUT_SCHEMA["cluster_summary"]["columns"])
    py_table = pd.read_csv(
        Path(args.python_results) / "cluster_summary.csv",
        dtype=str,
        keep_default_na=False,
    )
    r_table = pd.read_csv(
        Path(args.r_results) / "cluster_summary.csv", dtype=str, keep_default_na=False
    )
    result["cluster_summary_columns_equal"] = (
        py_table.columns.tolist() == table_fields
        and r_table.columns.tolist() == table_fields
    )
    result["cluster_summary_equal"] = result[
        "cluster_summary_columns_equal"
    ] and py_table.to_dict("records") == r_table.to_dict("records")

    failures = []
    for key in (
        "cells_equal",
        "clusters_equal",
        "genes_equal",
        "candidate_menu_equal",
        "candidate_panel_keys_equal",
        "cluster_summary_columns_equal",
        "cluster_summary_equal",
    ):
        if not result[key]:
            failures.append(key)
    if not result["differential_expression"]["row_membership_equal"]:
        failures.append("differential_expression_rows")
    for field, delta in (
        result["differential_expression"].get("max_abs_difference", {}).items()
    ):
        if delta > NUMERIC_TOLERANCE[field]:
            failures.append(f"differential_expression_{field}")
    if result["annotation_differences"]:
        failures.append("cluster_annotation")
    result["failures"] = failures
    result["status"] = "PASS" if not failures else "FAIL"
    Path(args.output).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
