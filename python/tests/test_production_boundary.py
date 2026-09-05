"""The installable package must not be able to see author reference labels.

Two products are built from one tree: the package under `software/`, which a user
runs on their own dataset and which therefore has no reference labels to compare
against, and the benchmark workspace under `benchmark/`, which joins author labels
by cell id and computes reference-derived metrics. The separation is only worth
claiming if it is enforced mechanically, so these tests fail on any reference
vocabulary reaching the package, on any benchmark-only input living inside it, and
on any reference-derived field appearing in the published output schema.
"""

from __future__ import annotations

import json
import re
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent
PACKAGE = ROOT / "src" / "scmarkeragent"
BENCHMARK = PROJECT / "benchmark"

# Reference-derived vocabulary. Assembled from fragments so this file can name what
# it forbids without matching itself.
FORBIDDEN = tuple(
    re.compile(pattern)
    for pattern in (
        "adjusted" + "_rand",
        "reference" + "_celltype",
        "reference" + "_manifest",
        "author" + "_label",
        "ground" + "_truth",
        r"\b" + "reference" + "_cl_id" + r"\b",
        r"\b" + "puri" + "ty" + r"\b",
    )
)

# The viewer states in user-facing prose that no reference metric is shown. Naming
# the absent quantity is the point of that line, so it is the one allowed mention.
ALLOWED_MENTIONS = {("viewer.py", "puri" + "ty")}

SCANNED_SUFFIXES = {".py", ".R", ".json", ".txt", ".toml", ".cfg", ".md"}


def _package_files() -> list[Path]:
    return [
        path
        for path in sorted(PACKAGE.rglob("*"))
        if path.is_file()
        and path.suffix in SCANNED_SUFFIXES
        and "__pycache__" not in path.parts
        # the vendored ontologies are third-party data, not our source
        and "resources" not in path.parts
    ]


def test_package_source_is_free_of_reference_vocabulary():
    failures = []
    for path in _package_files():
        text = path.read_text(encoding="utf-8", errors="ignore")
        for pattern in FORBIDDEN:
            for match in pattern.finditer(text):
                if (path.name, match.group(0)) in ALLOWED_MENTIONS:
                    continue
                line = text.count("\n", 0, match.start()) + 1
                failures.append(
                    f"{path.relative_to(ROOT)}:{line}: {match.group(0)}"
                )
    assert not failures, "reference vocabulary inside the installable package: " + "; ".join(
        failures
    )


def test_benchmark_inputs_live_outside_the_installable_package():
    # The author-label column mapping is the benchmark's only reference input; it
    # must not be reachable from the package or from its packaged data.
    #
    # Found by search rather than at a fixed path. The benchmark tree gets reorganised
    # (the file moved into scripts/ on 2026-08-07), and a stale path here fails as "the
    # reference input is missing" -- which reads like the opposite of what this guards,
    # and invites someone to satisfy it by putting a copy back in the wrong place.
    #
    # The published tree ships without the benchmark workspace at all. The boundary is
    # then held by construction, and the package-side assertion still binds.
    if BENCHMARK.is_dir():
        assert list(BENCHMARK.rglob("reference_manifest.csv"))
    assert not list(ROOT.rglob("reference_manifest.csv"))


def test_packaging_cannot_ship_the_benchmark_workspace():
    config = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    # Everything shipped is resolved from `src/`, and the benchmark workspace is a
    # sibling of `software/`, so no packaging glob can reach it.
    assert config["tool"]["setuptools"]["package-dir"] == {"": "src"}
    assert config["tool"]["setuptools"]["packages"]["find"]["where"] == ["src"]
    assert BENCHMARK.resolve() not in ROOT.resolve().parents
    assert not (ROOT / "src" / "benchmark").exists()
    for glob in config["tool"]["setuptools"]["package-data"]["scmarkeragent"]:
        assert not glob.startswith(".."), glob
        # every declared glob must actually match something, or the wheel silently
        # ships without a resource the package resolves at import time
        assert list(PACKAGE.glob(glob)), f"package-data glob matches nothing: {glob}"


def test_published_schema_declares_no_reference_derived_field():
    schema = json.loads(
        (PACKAGE / "schemas" / "output_schema.json").read_text(encoding="utf-8")
    )
    declared = []
    for section in ("cluster_summary", "marker_evidence"):
        declared.extend(schema[section]["columns"])
    declared.extend(schema["cluster_evidence"]["required"])
    for name in declared:
        for pattern in FORBIDDEN:
            assert not pattern.search(name), f"{name} is reference-derived"
    for name in ("reference_composition", "reference_purity", "secondary_annotation"):
        assert name not in declared
