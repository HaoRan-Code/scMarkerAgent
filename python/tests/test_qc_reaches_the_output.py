"""Quality control is only worth paying for if its verdict reaches the reader.

The failure this guards against has already happened once on this pipeline: the R arm
computed an exclusion audit, paid a model turn for it, and dropped it on the way to the
record because one interchange field was missing. Nothing in the logs said so. So the
chain is asserted hop by hop rather than at either end -- the record, the summary column
and the sidecar are three different places a verdict can go missing between.
"""

from __future__ import annotations

import pandas as pd

from scmarkeragent import cluster_annotation as agents
from scmarkeragent import reporting

TAG = "fixture"


def _review(**overrides):
    review = {
        "checked": True,
        "rounds": 2,
        "selected_round": 2,
        "judged_label": "enterocyte",
        "demoted_subtype": "",
        "passed": True,
        "rank": 3,
        "subtype_verdict": None,
        "annotation_verdict": {
            "verdict": "supported",
            "supporting_markers": ["SI"],
            "conflicting_markers": [],
            "reason": "SI is tied to enterocytes by the supplied sentence.",
        },
    }
    review.update(overrides)
    return review


def _frame(record):
    scored = pd.DataFrame(
        [
            {
                "cluster": "1",
                "candidate": "enterocyte",
                "retrieval_rank": 1,
                "retrieval_score": 1.0,
                "marker_level": 1.0,
                "cluster_level": 1.0,
                "single_cell_level": 1.0,
                "hits": 2,
                "panel_size": 4,
            }
        ]
    )
    return {
        "res": {"1": record},
        "scoring": {
            "clusters": {"1": {"n_cells": 288}},
            "candidate_cl_id": {"enterocyte": "CL:0000584"},
            "scored": scored,
            "top_candidates": 1,
        },
    }


def _marker_evidence():
    return pd.DataFrame(
        [
            {
                "dataset": TAG,
                "cluster_id": "1",
                "candidate_annotation": "enterocyte",
                "is_selected_annotation": "True",
                "gene": "SI",
                "cluster_detection_fraction": "0.865",
                "out_of_cluster_detection_fraction": "0.685",
                "publication_support_count": 4,
                "pmcid": "PMC1",
            }
        ]
    )


def _record(pool, review):
    outcome = {
        "final": {
            "action": "final",
            "schema_version": agents.ANNOTATOR_SCHEMA,
            "selected": "enterocyte",
            "subtype": "",
            "state": "",
            "co_occurring_identities": [],
            "claim_evidence": [
                {
                    "identity": "enterocyte",
                    "decisive_gene": "FABP2",
                    "pct_in": 92.4,
                    "pct_out": 68.5,
                }
            ],
            "identity_groups": [
                {"identity": "absorptive", "candidates": ["enterocyte"]},
                {
                    "identity": "secretory",
                    "candidates": ["enteroendocrine cell", "tuft cell"],
                },
            ],
            "support_markers": ["FABP2"],
            "confidence": "high",
            "reason": "FABP2 at 92.4 against 68.5.",
        },
        "turns": 2,
        "transcript": [],
        "quality_control": review,
    }
    return agents._result("1", pool, {}, outcome, pd.DataFrame())


def test_the_verdict_travels_from_the_record_to_both_delivered_files(pool):
    record = _record(pool, _review())
    # Hop 1: the annotation stage's own record.
    assert record["annotation_qc"] == agents.QC_REVISED
    assert record["quality_control"]["annotation_verdict"]["verdict"] == "supported"

    frame = _frame(record)
    marker_evidence = _marker_evidence()

    # Hop 2: the reader-facing summary carries the one word, in the declared column set.
    clustermap = reporting.build_clustermap(frame)
    table_a = reporting.build_table_a(TAG, clustermap, marker_evidence)
    assert "annotation_qc" in table_a.columns
    assert list(table_a.columns) == reporting.TABLE_A_COLS
    assert table_a.loc[0, "annotation_qc"] == agents.QC_REVISED

    # Hop 3: the sidecar carries the verdicts themselves, beside the other audit.
    evidence = reporting.build_cluster_evidence(TAG, frame, marker_evidence)
    audit = evidence[0]["audit"]
    assert audit["judge"]["annotation_verdict"]["verdict"] == "supported"
    assert audit["judge"]["selected_round"] == 2
    # ... and the reason, which is the part a reader actually needs when it did not pass.
    assert "SI is tied to enterocytes" in audit["judge"]["annotation_verdict"]["reason"]


def test_a_failed_check_is_delivered_as_failed_rather_than_dropped(pool):
    """A label that never passed still ships, and says so. Withholding it would turn a
    quality flag into an abstention, which is a different product decision."""
    review = _review(
        passed=False,
        rank=0,
        annotation_verdict={
            "verdict": "contradicted",
            "conflicting_markers": ["PTPRC"],
            "reason": "PTPRC excludes this identity.",
        },
    )
    record = _record(pool, review)
    assert record["annotation_qc"] == agents.QC_FAILED
    # The label, its rationale and its markers are all still there.
    assert record["annotation"] == "enterocyte"
    assert record["rationale"]
    assert record["support_markers"] == ["FABP2"]
    table_a = reporting.build_table_a(
        TAG, reporting.build_clustermap(_frame(record)), _marker_evidence()
    )
    assert table_a.loc[0, "annotation_qc"] == agents.QC_FAILED
    assert table_a.loc[0, "primary_annotation"] == "enterocyte"


def test_the_summary_shows_the_strongest_markers_rather_than_the_alphabet(pool):
    """`key_markers` is capped at twelve, so which twelve it keeps is the whole column.

    Table B is written gene-ascending and Table A used to take its first twelve rows. On
    `human_bladder` cluster 1 that dropped PECAM1 and VWF -- 2,640 and 1,169 publications,
    both detected in over 91% -- and kept a one-publication gene with a log2FC of 0.18,
    because P and V sort after K. The panel here is built the same way: the two best-
    corroborated genes sort last by name.
    """
    genes = [
        # gene, n_pub, pct_in, pct_out
        ("AAA1", 1, "0.100", "0.090"),
        ("AAA2", 2, "0.110", "0.100"),
        ("AAA3", 3, "0.120", "0.110"),
        ("AAA4", 4, "0.130", "0.120"),
        ("AAA5", 5, "0.140", "0.130"),
        ("AAA6", 6, "0.150", "0.140"),
        ("AAA7", 7, "0.160", "0.150"),
        ("AAA8", 8, "0.170", "0.160"),
        ("AAA9", 9, "0.180", "0.170"),
        ("BBB1", 10, "0.190", "0.180"),
        ("BBB2", 11, "0.200", "0.190"),
        ("BBB3", 12, "0.210", "0.200"),
        ("ZPECAM1", 2640, "0.919", "0.086"),
        ("ZVWF", 1169, "0.913", "0.076"),
    ]
    marker_evidence = pd.DataFrame(
        [
            {
                "dataset": TAG,
                "cluster_id": "1",
                "candidate_annotation": "enterocyte",
                "is_selected_annotation": "True",
                "gene": gene,
                "cluster_detection_fraction": pct_in,
                "out_of_cluster_detection_fraction": pct_out,
                "publication_support_count": n_pub,
                "pmcid": f"PMC{n_pub}",
            }
            for gene, n_pub, pct_in, pct_out in genes
        ]
    )
    record = _record(pool, _review())
    table_a = reporting.build_table_a(
        TAG, reporting.build_clustermap(_frame(record)), marker_evidence
    )
    shown = table_a.loc[0, "key_markers"]
    assert shown.startswith("ZPECAM1(0.919); ZVWF(0.913);")
    assert "AAA1(" not in shown          # one publication, last by evidence, first by name
    assert shown.count(";") == 11        # still capped at twelve
    assert table_a.loc[0, "cooccurring_markers"] == reporting.NOT_AVAILABLE
    assert table_a.loc[0, "cooccurring_pmcid"] == reporting.NOT_AVAILABLE


def test_summary_splits_cooccurring_markers_from_the_assigned_identity(pool):
    """Assigned and co-occurring evidence are separate columns, each capped at twelve.

    A mixed cluster used to show only is_selected rows in key_markers, so the reader had
    to open marker_evidence.csv to see what carried the co-occurring identity. The
    co-occurring shortlist is tagged with its identity because a cluster can carry more
    than one, and it is fair-shared so one well-published co-occurring panel cannot hide
    another.
    """
    rows = []
    for gene, n_pub, pct_in in [
        ("Cd3e", 100, "0.670"),
        ("Cxcr6", 80, "0.710"),
        ("Tbx21", 60, "0.030"),
    ]:
        rows.append(
            {
                "dataset": TAG,
                "cluster_id": "1",
                "candidate_annotation": "type I NK T cell",
                "is_selected_annotation": "True",
                "gene": gene,
                "marker_polarity": "positive",
                "cluster_detection_fraction": pct_in,
                "out_of_cluster_detection_fraction": "0.010",
                "publication_support_count": n_pub,
                "pmcid": f"PMC{n_pub}",
            }
        )
    for gene, n_pub, pct_in in [
        ("Gata3", 50, "0.480"),
        ("Il1rl1", 40, "0.370"),
        ("Il7r", 30, "0.830"),
    ]:
        rows.append(
            {
                "dataset": TAG,
                "cluster_id": "1",
                "candidate_annotation": "group 2 innate lymphoid cell",
                "is_selected_annotation": "False",
                "gene": gene,
                "marker_polarity": "positive",
                "cluster_detection_fraction": pct_in,
                "out_of_cluster_detection_fraction": "0.020",
                "publication_support_count": n_pub,
                "pmcid": f"PMC{n_pub + 1000}",
            }
        )
    marker_evidence = pd.DataFrame(rows)
    record = _record(pool, _review())
    record["annotation"] = "type I NK T cell"
    record["co_occurring_identities"] = ["group 2 innate lymphoid cell"]
    record["resolution_status"] = "mixed"
    frame = _frame(record)
    frame["scoring"]["candidate_cl_id"] = {"type I NK T cell": "CL:0000000"}
    table_a = reporting.build_table_a(
        TAG, reporting.build_clustermap(frame), marker_evidence
    )
    assert "Cd3e(0.670)" in table_a.loc[0, "key_markers"]
    assert "Gata3|" not in table_a.loc[0, "key_markers"]
    assert (
        "Gata3|group 2 innate lymphoid cell|0.480"
        in table_a.loc[0, "cooccurring_markers"]
    )
    assert (
        "Il7r|group 2 innate lymphoid cell|0.830"
        in table_a.loc[0, "cooccurring_markers"]
    )
    assert "PMC1050" in table_a.loc[0, "cooccurring_pmcid"]
    assert "PMC100" in table_a.loc[0, "pmcid"]

    # A rejected_subtype row must not appear as co-occurring when alternative_candidates
    # does not list it — marker_evidence keeps those rows for audit, not as a mix claim.
    rejected = {
        "dataset": TAG,
        "cluster_id": "1",
        "candidate_annotation": "proliferative β-cell",
        "is_selected_annotation": "False",
        "gene": "MKI67",
        "marker_polarity": "positive",
        "cluster_detection_fraction": "0.900",
        "out_of_cluster_detection_fraction": "0.010",
        "publication_support_count": 99,
        "pmcid": "PMC9999",
    }
    marker_evidence = pd.concat(
        [marker_evidence, pd.DataFrame([rejected])], ignore_index=True
    )
    table_a = reporting.build_table_a(
        TAG, reporting.build_clustermap(frame), marker_evidence
    )
    assert "MKI67|" not in table_a.loc[0, "cooccurring_markers"]
    assert "proliferative" not in table_a.loc[0, "cooccurring_markers"]


def test_a_gene_curated_for_two_identities_takes_one_slot(pool):
    """A shared gene is one measurement, so it gets one slot and is printed once.

    Both a parent and its subtype routinely curate the same gene. Printed twice it reads
    as two independent pieces of evidence, and it spends two of the twelve slots on one
    measurement -- on the freeze that hit 48 clusters, `INS` and `MAFA` among them.
    """
    rows = []
    for identity, genes in [
        ("type B pancreatic cell", ["INS", "PDX1", "MAFA"]),
        ("mature β-cell", ["INS", "UCN3"]),
    ]:
        for n_pub, gene in enumerate(genes, start=1):
            rows.append(
                {
                    "dataset": TAG,
                    "cluster_id": "1",
                    "candidate_annotation": identity,
                    "is_selected_annotation": "True",
                    "gene": gene,
                    "marker_polarity": "positive",
                    "cluster_detection_fraction": "0.900",
                    "out_of_cluster_detection_fraction": "0.100",
                    "publication_support_count": 10 - n_pub,
                    "pmcid": f"PMC{gene}",
                }
            )
    record = _record(pool, _review())
    table_a = reporting.build_table_a(
        TAG, reporting.build_clustermap(_frame(record)), pd.DataFrame(rows)
    )
    shown = table_a.loc[0, "key_markers"].split("; ")
    genes = [item.split("(")[0] for item in shown]
    assert genes.count("INS") == 1
    assert sorted(genes) == ["INS", "MAFA", "PDX1", "UCN3"]


def test_every_reported_status_is_one_the_schema_documents():
    schema_text = (
        reporting.OUTPUT_SCHEMA["cluster_summary"]["fields"]["annotation_qc"]
    )
    for status in agents.QC_VALUES:
        assert status in schema_text
