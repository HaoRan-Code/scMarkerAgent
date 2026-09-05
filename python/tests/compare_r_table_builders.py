#!/usr/bin/env python3
"""Feed both arms the same summary-table inputs and require the same table out.

The frozen golden runs the report layer end to end, which is the right acceptance test and
a weak branch test: `validation/synthetic_golden` has four clusters, one candidate each,
three markers at most and no co-occurring identity, so `build_table_a` agreeing on it says
nothing about the twelve-slot cap, the fair share between identities, the de-duplication by
gene, or the co-occurring columns. Those are most of the function.

This harness skips the pipeline entirely and hands both implementations the same
`marker_evidence` and `clustermap` tables, built to enter each of those branches. It is
cheap enough to run on every change and it fails on the thing an end-to-end fixture cannot
see: one arm taking a branch the other does not.

  Usage: PYTHONPATH=src python tests/compare_r_table_builders.py
"""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from scmarkeragent.reporting import TABLE_A_COLS, build_table_a  # noqa: E402

MARKER_COLS = [
    "dataset", "cluster_id", "candidate_annotation", "is_selected_annotation", "gene",
    "marker_polarity", "marker_role", "marker_provenance", "cluster_detection_fraction",
    "out_of_cluster_detection_fraction", "average_log2_fold_change", "adjusted_p_value",
    "cross_cluster_percentile", "marker_specificity_weight", "publication_support_count",
    "evidence_tier", "pmid", "pmcid", "source_sentence",
]


def marker(cluster, identity, selected, gene, *, pub, pct_in, pct_out=0.05,
           polarity="positive", pmcid="PMC1"):
    return {
        "dataset": "fixture",
        "cluster_id": cluster,
        "candidate_annotation": identity,
        "is_selected_annotation": "True" if selected else "False",
        "gene": gene,
        "marker_polarity": polarity,
        "marker_role": "key",
        "marker_provenance": "curated",
        "cluster_detection_fraction": f"{pct_in:.4f}",
        "out_of_cluster_detection_fraction": f"{pct_out:.4f}",
        "average_log2_fold_change": "1.0000",
        "adjusted_p_value": "0.001",
        "cross_cluster_percentile": "0.9",
        "marker_specificity_weight": "0.5",
        "publication_support_count": str(pub),
        "evidence_tier": "A",
        "pmid": "1",
        "pmcid": pmcid,
        "source_sentence": "a sentence with, a comma and \"quotes\"",
    }


def cluster_row(cluster, order, *, alternatives="N/A", primary="alpha cell"):
    return {
        "cluster": cluster,
        "cluster_order": order,
        "n_cells": 100 + order,
        "cell_type_annotation": primary,
        "cell_subtype_annotation": "N/A",
        "primary_annotation": primary,
        "cell_state": "N/A",
        "cell_ontology": "CL:0000171",
        "annotation_confidence": "high",
        "annotation_rationale": "because the panel says so",
        "resolution_status": "resolved",
        "resolution_detail": "agent_selected",
        "annotation_source": "cluster_annotation",
        "llm_status": "annotated",
        "annotation_qc": "passed",
        "cluster_celltype": f"{cluster}: {primary}",
        "alternative_candidates": alternatives,
        "claim_warnings": "N/A",
        "unlisted_identity": "N/A",
        "identity_groups_json": "[]",
    }


def build_fixture() -> tuple[pd.DataFrame, pd.DataFrame]:
    markers: list[dict] = []
    clusters: list[dict] = []

    # 1. More than twelve candidates for one identity: exercises the cap, and the
    #    publication-support ordering that decides which twelve survive it.
    for index in range(18):
        markers.append(
            marker("1", "alpha cell", True, f"GENE{index:02d}", pub=100 - index,
                   pct_in=0.9 - index / 100)
        )
    clusters.append(cluster_row("1", 0))

    # 2. Two identities, one co-occurring: exercises the fair share, the cooc columns and
    #    the identity tagging. GENE_SHARED is curated for both, so it also exercises the
    #    one-slot-per-gene rule -- the identity holding the earlier slot keeps it.
    for index in range(9):
        markers.append(
            marker("2", "beta cell", True, f"BGENE{index:02d}", pub=90 - index,
                   pct_in=0.8, pmcid=f"PMC{index}")
        )
    markers.append(marker("2", "beta cell", True, "GENE_SHARED", pub=50, pct_in=0.7))
    for index in range(9):
        markers.append(
            marker("2", "delta cell", False, f"DGENE{index:02d}", pub=80 - index,
                   pct_in=0.6, pmcid=f"PMC{index}")
        )
    markers.append(marker("2", "delta cell", False, "GENE_SHARED", pub=95, pct_in=0.7))
    clusters.append(cluster_row("2", 1, alternatives="delta cell", primary="beta cell"))

    # 3. Negative polarity and a missing publication count: the first must never reach
    #    key_markers, the second must sort as zero rather than crash.
    markers.append(marker("3", "gamma cell", True, "POSGENE", pub=5, pct_in=0.5))
    markers.append(
        marker("3", "gamma cell", True, "NEGGENE", pub=999, pct_in=0.9,
               polarity="negative")
    )
    blank = marker("3", "gamma cell", True, "NOPUB", pub=0, pct_in=0.4)
    blank["publication_support_count"] = ""
    blank["pmcid"] = "N/A"
    markers.append(blank)
    clusters.append(cluster_row("3", 2, primary="gamma cell"))

    # 4. No marker rows at all, and a UTF-8 / punctuation name: the empty branch must give
    #    the not-applicable token, not a blank cell.
    clusters.append(
        cluster_row("4", 3, primary="Mac-2 (Atf3) \u03b1\u03b2 cell \u2014 \u4e2d\u6587")
    )

    # 5. Only non-selected rows: `build_table_a` falls back to the whole group.
    markers.append(marker("5", "epsilon cell", False, "EGENE", pub=7, pct_in=0.3))
    clusters.append(cluster_row("5", 4, primary="epsilon cell"))

    # 6. TWO selected identities, together holding more than the cap, sharing a gene.
    #    Nothing above reaches the fair share or the de-duplication: a group with one
    #    identity gets the whole cap, and a gene curated for two identities never lands in
    #    the same group. Both branches were green with the rule deleted until this cluster
    #    existed, which is the only reason it is here.
    for index in range(10):
        markers.append(
            marker("6", "parent cell", True, f"PGENE{index:02d}", pub=200 - index,
                   pct_in=0.95)
        )
    for index in range(10):
        markers.append(
            marker("6", "parent cell subtype", True, f"SGENE{index:02d}", pub=20 - index,
                   pct_in=0.55)
        )
    # Curated for both, best-published on the parent: one slot, and the parent's.
    markers.append(marker("6", "parent cell", True, "BOTH_GENE", pub=150, pct_in=0.9))
    markers.append(
        marker("6", "parent cell subtype", True, "BOTH_GENE", pub=15, pct_in=0.9)
    )
    clusters.append(cluster_row("6", 5, primary="parent cell subtype"))

    # 7. The same two pressures on the CO-OCCURRING column, which allocates separately.
    markers.append(marker("7", "host cell", True, "HGENE", pub=60, pct_in=0.7))
    for index in range(8):
        markers.append(
            marker("7", "delta cell", False, f"D2GENE{index:02d}", pub=70 - index,
                   pct_in=0.5)
        )
    for index in range(8):
        markers.append(
            marker("7", "PP cell", False, f"P2GENE{index:02d}", pub=65 - index,
                   pct_in=0.45)
        )
    markers.append(marker("7", "delta cell", False, "COOC_SHARED", pub=99, pct_in=0.6))
    markers.append(marker("7", "PP cell", False, "COOC_SHARED", pub=98, pct_in=0.6))
    clusters.append(
        cluster_row("7", 6, alternatives="delta cell; PP cell", primary="host cell")
    )

    return pd.DataFrame(markers, columns=MARKER_COLS), pd.DataFrame(clusters)


def main() -> int:
    marker_evidence, clustermap = build_fixture()
    with tempfile.TemporaryDirectory() as work_name:
        work = Path(work_name)
        marker_evidence.to_csv(work / "marker_evidence.csv", index=False)
        clustermap.to_csv(work / "clustermap.csv", index=False)

        expected = build_table_a("fixture", clustermap, marker_evidence)
        expected.to_csv(work / "expected.csv", index=False)

        script = f"""
suppressPackageStartupMessages({{library(jsonlite); library(data.table); library(Matrix)}})
Sys.setenv(SCMA_OFFLINE = "1")
pkg <- "{ROOT / 'src' / 'scmarkeragent'}"
source(file.path(pkg, "rflow", "config.R"))
source(file.path(pkg, "rflow", "reporting.R"))
me <- fread("{work}/marker_evidence.csv", colClasses = "character", na.strings = NULL)
cm <- fread("{work}/clustermap.csv", colClasses = "character", na.strings = NULL)
ta <- scma_build_table_a("fixture", cm, me)
fwrite(ta, "{work}/actual.csv")
"""
        run = subprocess.run(
            ["Rscript", "-e", script], capture_output=True, text=True, timeout=600
        )
        if run.returncode != 0:
            print("R side failed:\n" + (run.stderr or run.stdout)[-2000:])
            return 1

        def read(path: Path) -> list[dict]:
            with path.open(newline="", encoding="utf-8") as handle:
                return list(csv.DictReader(handle))

        want, got = read(work / "expected.csv"), read(work / "actual.csv")

    if len(want) != len(got):
        print(f"row count: {len(want)} expected, {len(got)} found")
        return 1
    failures = []
    for index, (left, right) in enumerate(zip(want, got)):
        if list(left) != list(right):
            failures.append(f"row {index}: column set differs")
            continue
        for column in TABLE_A_COLS:
            if str(left[column]) != str(right[column]):
                failures.append(
                    f"row {index} cluster {left['cluster_id']}, column '{column}'\n"
                    f"    python: {left[column]!r}\n    r     : {right[column]!r}"
                )
    print(f"compared {len(want)} rows x {len(TABLE_A_COLS)} columns")
    if failures:
        print(f"FAILED with {len(failures)} difference(s):")
        for failure in failures[:12]:
            print(f"  - {failure}")
        return 1
    print("PASSED: both arms build the same Table A on every branch")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
