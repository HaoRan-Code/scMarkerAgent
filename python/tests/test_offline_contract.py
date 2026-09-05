from __future__ import annotations

from pathlib import Path
from unittest import mock

from scmarkeragent import cluster_annotation, llm_client

ROOT = Path(__file__).resolve().parents[1]


def test_offline_mode_hides_credentials(monkeypatch):
    monkeypatch.setenv("SCMA_OFFLINE", "1")
    monkeypatch.setenv("OPENAI_API_KEY", "must-not-be-used")
    monkeypatch.setenv("OPENAI_BASE_URL", "https://example.invalid/v1")
    assert llm_client.resolve_api() == ("", "")


def test_missing_credentials_never_reach_the_network(tmp_path, monkeypatch):
    monkeypatch.setenv("SCMA_LLM_CACHE", str(tmp_path / "responses.json"))
    monkeypatch.setenv("SCMA_OFFLINE", "1")
    with mock.patch.object(
        llm_client,
        "call_llm",
        side_effect=AssertionError("network call attempted"),
    ):
        response, error, replayed = llm_client.cached_call_llm(
            "packet", "", "", reasoning_effort="high"
        )
    assert response is None
    assert error == "cache_miss_no_credentials"
    assert replayed is False


def test_cached_response_replays_without_credentials(tmp_path, monkeypatch):
    monkeypatch.setenv("SCMA_LLM_CACHE", str(tmp_path / "responses.json"))
    prompt = "packet"
    llm_client._RESPONSE_CACHE.path = None
    llm_client._RESPONSE_CACHE.data = None
    cache = llm_client._load_vcache()
    cache[llm_client._vkey(prompt, "high")] = '{"selected":"Unknown"}'
    with mock.patch.object(
        llm_client,
        "call_llm",
        side_effect=AssertionError("network call attempted"),
    ):
        response, error, replayed = llm_client.cached_call_llm(
            prompt, "", "", reasoning_effort="high"
        )
    assert response == '{"selected":"Unknown"}'
    assert error is None
    assert replayed is True


def test_fallback_is_never_presented_as_a_model_judgement():
    """A label produced without a model must say so in every field that records how."""
    import pandas as pd

    frame = pd.DataFrame(
        [{"candidate": "hepatocyte", "hits": 12, "panel_size": 20}]
    )
    record = cluster_annotation._fallback(
        "3", frame, "skipped_no_credentials", "relative_score_fallback_top1"
    )
    assert record["annotation"] == "hepatocyte"
    assert record["annotation_source"] == "relative_score_fallback"
    assert record["resolution_status"] == cluster_annotation.RESOLVED
    assert record["resolution_detail"] == "relative_score_fallback_top1"
    assert record["llm_status"] == "skipped_no_credentials"
    assert record["confidence"] == cluster_annotation.NOT_AVAILABLE
    assert "not a confidence" in record["rationale"]
    assert record["claim_warnings"] == []


def test_empty_candidate_pool_abstains_rather_than_guessing():
    import pandas as pd

    record = cluster_annotation._fallback(
        "7", pd.DataFrame(), "skipped_no_credentials", "empty_candidate_pool"
    )
    assert record["annotation"] == cluster_annotation.UNKNOWN
    assert record["resolution_status"] == cluster_annotation.UNRESOLVED
    assert record["resolution_detail"] == "empty_candidate_pool"
    assert record["annotation_source"] == "no_candidate"


def test_resolution_status_has_exactly_three_reader_facing_values():
    """A reader should not have to learn five words for three outcomes."""
    assert {
        cluster_annotation.RESOLVED,
        cluster_annotation.MIXED,
        cluster_annotation.UNRESOLVED,
    } == {"resolved", "mixed", "unresolved"}


def test_cell_type_column_carries_only_a_cell_type_or_unknown():
    """The reader-facing label must not require a vocabulary to interpret.

    A `Mixed: A + B` form and an `Uncertain` token were removed: neither told a reader
    which cell they were looking at. A cluster carrying several identities keeps the
    majority one as its label, and a cluster establishing none is Unknown.
    """
    from scmarkeragent import config

    agent_config = config.DEFAULTS["cluster_annotation"]
    assert "mixed_prefix" not in agent_config
    assert "uncertain_token" not in agent_config
    assert agent_config["unknown_token"] == "Unknown"
    source = (ROOT / "src/scmarkeragent/cluster_annotation.py").read_text(
        encoding="utf-8"
    )
    for removed in ("MIXED_PREFIX", "UNCERTAIN"):
        assert removed not in source
    prompt = (
        ROOT / "src/scmarkeragent/prompts/cluster_annotator.txt"
    ).read_text(encoding="utf-8")
    assert "Mixed:" not in prompt
    assert "Uncertain" not in prompt


def test_a_mixture_is_recorded_beside_the_label_not_inside_it(pool):
    def answer(**overrides):
        value = {
            "action": "final",
            "schema_version": cluster_annotation.ANNOTATOR_SCHEMA,
            "selected": "enterocyte",
            "subtype": "",
            "state": "",
            "co_occurring_identities": ["tuft cell"],
            "claim_evidence": [
                {
                    "identity": "enterocyte",
                    "decisive_gene": "FABP2",
                    "pct_in": 92.4,
                    "pct_out": 68.5,
                },
                {
                    "identity": "tuft cell",
                    "decisive_gene": "POU2F3",
                    "pct_in": 1.0,
                    "pct_out": 40.0,
                },
            ],
            "identity_groups": [
                {"identity": "absorptive", "candidates": ["enterocyte"]},
                {"identity": "endocrine", "candidates": ["enteroendocrine cell"]},
                {"identity": "tuft", "candidates": ["tuft cell"]},
            ],
            "support_markers": ["FABP2"],
            "confidence": "high",
            "reason": "FABP2 defines the majority; POU2F3 marks a separate population.",
        }
        value.update(overrides)
        return value

    assert cluster_annotation.valid_final(answer(), pool, "1")

    # A co-occurring identity has to be a different identity, not a synonym: a name from
    # the selected candidate's own group is the same cell read twice.
    assert not cluster_annotation.valid_final(
        answer(
            identity_groups=[
                {
                    "identity": "absorptive",
                    "candidates": ["enterocyte", "tuft cell"],
                },
                {"identity": "endocrine", "candidates": ["enteroendocrine cell"]},
            ]
        ),
        pool,
        "1",
    )
    # Unknown means nothing was established, so nothing can be reported beside it.
    assert not cluster_annotation.valid_final(
        answer(selected=cluster_annotation.UNKNOWN, support_markers=[]), pool, "1"
    )
    assert cluster_annotation.valid_final(
        answer(
            selected=cluster_annotation.UNKNOWN,
            co_occurring_identities=[],
            claim_evidence=[],
            support_markers=[],
        ),
        pool,
        "1",
    )
