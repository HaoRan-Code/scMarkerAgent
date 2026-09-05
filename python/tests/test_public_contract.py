from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

from scmarkeragent import config

ROOT = Path(__file__).resolve().parents[1]


def test_public_output_schema_is_concise_and_auditable():
    schema = json.loads(
        (ROOT / "src/scmarkeragent/schemas/output_schema.json").read_text()
    )
    assert schema["schema_version"] == "scmarkeragent-output-v1"
    assert schema["cluster_summary"]["columns"] == [
        "dataset",
        "cluster_id",
        "n_cells",
        # The label to read, plot and score against, derived in code from the two identity
        # columns that follow it rather than typed by a model.
        "primary_annotation",
        "cell_type_annotation",
        "cell_subtype_annotation",
        "cell_state",
        "cell_ontology",
        "annotation_confidence",
        "annotation_rationale",
        "resolution_status",
        "annotation_source",
        "llm_status",
        # One word for how the label left quality control. The verdicts behind it stay in
        # the sidecar: a reader scanning the summary wants to know whether to look, not
        # to read two judgments per row.
        "annotation_qc",
        "key_markers",
        "pmcid",
        # Markers/PMCIDs for alternative_candidates, split out of key_markers so a mixed
        # cluster's co-occurring evidence is visible without joining marker_evidence.csv.
        "cooccurring_markers",
        "cooccurring_pmcid",
        "alternative_candidates",
        # An audit flag, not a correction: which claimed identities have curated markers
        # this cluster does not raise. The annotation stays exactly what the agent
        # decided, and a reader can see which of its claims to check.
        "claim_warnings",
        # What the agent says retrieval did not offer. A column and not a sentence in the
        # rationale, because "how often was the answer not on the menu" is a rate somebody
        # has to be able to count, and retrieval hands over 15 of a median 67 admitted
        # cell types.
        "unlisted_identity",
    ]
    assert not {
        "subtype_refinement",
        "review_evidence",
        "refinement_status",
        "secondary_annotation",
        "reference_purity",
        "reference_composition",
        "candidate_fractions",
    }.intersection(schema["cluster_summary"]["columns"])
    assert schema["cluster_evidence"]["file"] == "cluster_evidence.jsonl"
    assert schema["cluster_evidence"]["format"] == "jsonl"
    marker_columns = schema["marker_evidence"]["columns"]
    for column in (
        "out_of_cluster_detection_fraction",
        "average_log2_fold_change",
        "cross_cluster_percentile",
        "marker_polarity",
        "marker_role",
        "pmid",
        "pmcid",
    ):
        assert column in marker_columns
    assert "log_fold_change" not in marker_columns


def test_not_applicable_token_is_one_shared_value():
    """Every table has to write the same token where a field does not apply."""
    schema = json.loads(
        (ROOT / "src/scmarkeragent/schemas/output_schema.json").read_text()
    )
    assert config.NOT_AVAILABLE == "N/A"
    assert schema["not_applicable_token"] == config.NOT_AVAILABLE
    assert config.na_display(None) == config.NOT_AVAILABLE
    assert config.na_display("") == config.NOT_AVAILABLE
    assert config.na_display(float("nan")) == config.NOT_AVAILABLE
    assert config.na_display("CL:0000236") == "CL:0000236"


def test_ontology_is_display_only():
    """The CL id may be shown. It may not be resolved, expanded or compared."""
    assert "cell_ontology" not in config.DEFAULTS["resources"]
    package = ROOT / "src/scmarkeragent"
    assert not (package / "cl_ontology.py").exists()
    assert not (package / "rflow/cl_ontology.R").exists()
    forbidden = re.compile(r"\.(ancestors|lca|is_desc|same_axis|ndesc|depth_of)\(")
    offenders = []
    for path in sorted(package.rglob("*.py")) + sorted(package.rglob("*.R")):
        if "uberon" in path.name:
            continue  # UBERON drives tissue subsumption and stays
        if forbidden.search(path.read_text(encoding="utf-8")):
            offenders.append(str(path.relative_to(package)))
    assert not offenders


def test_negative_markers_never_enter_a_score():
    """Negative markers are evidence for the model, never a term in a number."""
    source = (ROOT / "src/scmarkeragent/candidate_scoring.py").read_text(
        encoding="utf-8"
    )
    assert "negative_panels" in source
    for line in source.splitlines():
        if "negative" not in line.lower():
            continue
        assert "marker_level" not in line
        assert "cluster_level" not in line
        assert "single_cell_level" not in line
        assert "retrieval_score" not in line
    assert "negative" not in config.DEFAULTS["retrieval"]


def test_hits_gate_never_empties_a_pool_and_relaxes_one_step_at_a_time():
    """A thin curated menu must not be turned into a confident Unknown.

    The gate exists to keep one- and two-hit retrieval false positives out of the
    shortlist, but three of the benchmark's contexts hold under 110 eligible cell types,
    and there the gate would leave a handful of candidates or none at all. So it is
    checked on both counts: it must drop the weak candidates when there are enough strong
    ones, and it must hand back a non-empty pool with the applied threshold recorded when
    there are not.
    """
    import pandas as pd

    from scmarkeragent.candidate_scoring import MIN_HITS, MIN_POOL_FLOOR, admit_by_hits

    rich = pd.DataFrame({"candidate": [f"c{i}" for i in range(40)], "hits": list(range(1, 41))})
    admitted, threshold = admit_by_hits(rich)
    assert threshold == MIN_HITS
    assert int(admitted["hits"].min()) >= MIN_HITS

    # Only one- and two-hit candidates exist: the gate has to come down to reach them.
    thin = pd.DataFrame({"candidate": ["a", "b", "c"], "hits": [1, 2, 2]})
    admitted, threshold = admit_by_hits(thin)
    assert threshold == 1
    assert len(admitted) == len(thin)

    # Just under the floor at the strict threshold, plenty one step below it.
    borderline = pd.DataFrame(
        {
            "candidate": [f"c{i}" for i in range(MIN_POOL_FLOOR + 4)],
            "hits": [MIN_HITS] * (MIN_POOL_FLOOR - 1) + [MIN_HITS - 1] * 5,
        }
    )
    admitted, threshold = admit_by_hits(borderline)
    assert threshold == MIN_HITS - 1
    assert len(admitted) >= MIN_POOL_FLOOR

    for frame in (rich, thin, borderline):
        admitted, _ = admit_by_hits(frame)
        assert len(admitted) > 0


def test_user_facing_files_exclude_internal_terms():
    targets = [
        ROOT / "README.md",
        ROOT / ".env.example",
        ROOT / "src/scmarkeragent/config/defaults.json",
        ROOT / "src/scmarkeragent/schemas/output_schema.json",
        ROOT / "src/scmarkeragent/cli.py",
        ROOT / "src/scmarkeragent/r_cli.py",
        *sorted((ROOT / "src/scmarkeragent/prompts").glob("*.txt")),
        *sorted((ROOT / "examples").glob("*")),
        *sorted((ROOT / "tests").glob("*")),
    ]
    patterns = [
        re.compile(r"\b" + "B" + "0" + r"\b"),
        re.compile("Route" + "-A"),
        re.compile(r"\b" + "roll" + "up" + r"\b", re.IGNORECASE),
    ]
    failures = []
    for path in targets:
        if not path.is_file() or path.suffix in {".h5ad", ".rds", ".mtx"}:
            continue
        text = path.read_text(encoding="utf-8")
        for pattern in patterns:
            if pattern.search(text):
                failures.append(f"{path.relative_to(ROOT)}:{pattern.pattern}")
    assert not failures


def test_cli_help_uses_public_stage_names():
    result = subprocess.run(
        [sys.executable, "-m", "scmarkeragent.cli", "annotate", "--help"],
        check=True,
        capture_output=True,
        text=True,
    )
    # argparse rewraps help text, so compare on collapsed whitespace.
    helptext = " ".join(result.stdout.split())
    assert "--no-cluster-annotation" in helptext
    assert "top of the retrieval order" in helptext
