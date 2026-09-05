"""Point the suite at the in-tree curated bundle before scmarkeragent is imported.

The package ships only the ontology files under `resources/static`; the marker
and source exports are far too large to vendor into the package and live in
`datasets/resources/current` instead. Production supplies that location through
`SCMA_RESOURCE_DIR`, so the tests declare the same precondition rather than
depending on the caller's environment.
"""

import os
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[2]

os.environ.setdefault(
    "SCMA_RESOURCE_DIR", str(PROJECT / "datasets" / "resources" / "current")
)

import pandas as pd  # noqa: E402
import pytest  # noqa: E402

# One synthetic dataset carrying the two readings an unraised marker can have, because
# every check in the suite turns on telling them apart.
#
#   enterocyte     -- a dominant population. Its best-corroborated markers are detected
#                     in most of the cluster AND in most cells outside it, so the
#                     one-versus-rest fold change is negative and they read as unraised.
#   enteroendocrine -- an identity the cells do not carry. Its best-corroborated markers
#                     are near zero inside the cluster, and also unraised.
#   tuft cell      -- a straightforward positive: its top marker is raised.
_DE_ROWS = [
    # gene,      cluster, pct_in, pct_out, avg_log2FC
    ("FABP2", "1", 92.4, 68.5, -0.10),
    ("VIL1", "1", 70.8, 72.1, -0.22),
    ("SI", "1", 86.5, 50.1, 0.01),
    ("RBP2", "1", 95.8, 45.3, 1.43),
    ("CHGA", "1", 2.3, 4.1, -2.41),
    ("NEUROD1", "1", 0.0, 2.2, -0.22),
    ("PROX1", "1", 77.7, 21.7, 1.19),
    ("POU2F3", "1", 1.0, 40.0, -3.00),
    ("PTPRC", "1", 0.0, 1.2, -0.04),
    ("FABP2", "2", 10.0, 80.0, -2.00),
    ("VIL1", "2", 12.0, 74.0, -1.80),
    ("SI", "2", 9.0, 60.0, -1.20),
    ("RBP2", "2", 8.0, 55.0, -1.40),
    ("CHGA", "2", 91.0, 1.8, 7.41),
    ("NEUROD1", "2", 80.5, 0.2, 2.87),
    ("PROX1", "2", 30.0, 25.0, 0.20),
    ("POU2F3", "2", 89.3, 4.6, 3.10),
    ("PTPRC", "2", 0.0, 1.1, -0.03),
]

# (candidate, gene, polarity, n_pub, tier). Publication support is what orders a panel,
# so the fixture gives the defining marker of each identity the largest count.
#
# `M cell` is here and NOT in `_CANDIDATES`, because that is the shape of the real payload:
# these rows are every cell type eligible in the context, while a cluster is handed the
# retrieved subset -- 15 of a median 67 in the shipped run. It claims PROX1, the gene this
# cluster raises most cleanly, with more publications than either candidate that claims it,
# so the fixture carries the one case retrieval can get wrong and nothing downstream used
# to be able to see.
_PANEL_ROWS = [
    ("enterocyte", "VIL1", "positive", 53, "high"),
    ("enterocyte", "SI", "positive", 34, "high"),
    ("enterocyte", "FABP2", "positive", 24, "high"),
    ("enterocyte", "RBP2", "positive", 3, "low"),
    ("enterocyte", "PTPRC", "negative", 1, "low"),
    ("enteroendocrine cell", "CHGA", "positive", 160, "high"),
    ("enteroendocrine cell", "NEUROD1", "positive", 15, "high"),
    ("enteroendocrine cell", "PROX1", "positive", 4, "high"),
    ("enteroendocrine cell", "MISSING1", "positive", 2, "medium"),
    ("tuft cell", "POU2F3", "positive", 40, "high"),
    ("tuft cell", "PROX1", "positive", 1, "low"),
    ("M cell", "SPIB", "positive", 12, "high"),
    ("M cell", "PROX1", "positive", 8, "high"),
]

_CANDIDATES = ["enterocyte", "enteroendocrine cell", "tuft cell"]


def make_scoring_and_prep():
    """A scoring payload and DE table shaped exactly like the cached production ones."""
    de = pd.DataFrame(
        [
            {
                "feature": gene,
                "group": cluster,
                "pct_in": pct_in,
                "pct_out": pct_out,
                "avg_log2FC": lfc,
                "padj": 1e-9,
                "auc": 0.9,
                "avgExpr": 1.0,
            }
            for gene, cluster, pct_in, pct_out, lfc in _DE_ROWS
        ]
    )
    panel_records = pd.DataFrame(
        [
            {
                "candidate": candidate,
                "gene": gene.capitalize(),
                "gene_key": gene,
                "marker_polarity": polarity,
                "n_pub": n_pub,
                "tier": tier,
            }
            for candidate, gene, polarity, n_pub, tier in _PANEL_ROWS
        ]
    )
    measured = {gene for gene, *_ in _DE_ROWS}
    positive_panels: dict[str, list[str]] = {}
    negative_panels: dict[str, list[str]] = {}
    for candidate, gene, polarity, *_ in _PANEL_ROWS:
        if gene not in measured:
            continue
        target = positive_panels if polarity == "positive" else negative_panels
        target.setdefault(candidate, []).append(gene)

    scored = pd.DataFrame(
        [
            {
                "cluster": cluster,
                "candidate": candidate,
                "retrieval_rank": rank,
                "hits": 3,
                "panel_size": 4,
                "program_in_median": 0.6,
                "program_out_median": 0.4,
            }
            for cluster in ("1", "2")
            for rank, candidate in enumerate(_CANDIDATES, start=1)
        ]
    )
    scoring = {
        "context": {
            "species": "Human",
            "tissue": "intestine",
            "disease": ["Normal"],
            "development_stage": "",
            "n_clusters": 2,
        },
        "clusters": {
            "1": {"n_cells": 288, "status": "pool", "candidates": list(_CANDIDATES)},
            "2": {"n_cells": 133, "status": "pool", "candidates": list(_CANDIDATES)},
        },
        "scored": scored,
        "measured_panels": positive_panels,
        "negative_panels": negative_panels,
        "panel_records": panel_records,
        "marker_specificity": {gene: 0.5 for gene in measured},
        "top_candidates": len(_CANDIDATES),
    }
    return scoring, {"de": de}


@pytest.fixture
def pool():
    from scmarkeragent.annotator_pool import build_pool

    scoring, prep = make_scoring_and_prep()
    return build_pool(scoring, prep)


class FakeSources:
    """Stands in for the curated source index, returning one sentence per marker."""

    def __init__(self):
        self.asked = []

    def records_for_marker(self, cell_type, gene, k=3):
        self.asked.append((cell_type, str(gene).upper()))
        return [{"pmid": "1", "pmcid": "PMC1", "sentence": f"{gene} is absent from {cell_type}."}]


@pytest.fixture
def pool_with_detected_exclusion():
    """The enterocyte panel's negative marker raised above the provenance threshold.

    The stock fixture holds PTPRC at 0.0, which is the case where nothing should ship;
    this one is the mirror image, so both sides of the threshold are covered.
    """
    from scmarkeragent.annotator_pool import build_pool

    scoring, prep = make_scoring_and_prep()
    de = prep["de"]
    detected = (de["feature"] == "PTPRC") & (de["group"] == "1")
    de.loc[detected, ["pct_in", "pct_out", "avg_log2FC"]] = [66.1, 72.7, -0.54]
    sources = FakeSources()
    return build_pool(scoring, {"de": de}, sources), sources
