"""The conversation: asking for evidence, being answered, and terminating."""

from __future__ import annotations

import json

from scmarkeragent import cluster_annotation as agents


def _audit_genes(outcome, identity="enterocyte"):
    """The genes returned for one identity's block of the audit."""
    for block in outcome["exclusion_audit"]:
        if block["identity"] == identity:
            return [row["gene"] for row in block["detected_exclusions"]]
    return None


def _final_text(**overrides):
    value = {
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
        "reason": "FABP2 at 92.4 against 68.5 with a high pct_out is a dominant population.",
    }
    value.update(overrides)
    return json.dumps(value)


class _Script:
    """A model that returns the scripted replies in order, recording each prompt."""

    def __init__(self, replies):
        self.replies = list(replies)
        self.prompts = []

    def __call__(self, prompt, api_key, api_url, trace_id, turn):
        self.prompts.append(prompt)
        reply = self.replies[min(turn - 1, len(self.replies) - 1)]
        return reply, None, False


class _Judges:
    """Quality control that answers from a script instead of a model.

    Both judges pass unless a test says otherwise, so a test about the conversation is
    not also a test about the verdicts. Every packet is kept: what a judge was shown is
    as much a part of the contract as what it replied, because the one thing it must
    never be shown is the candidate list.
    """

    def __init__(self, subtype="separable", annotation="supported"):
        self.replies = {"subtype_judge": subtype, "annotation_judge": annotation}
        self.packets = []

    def __call__(
        self, name, packet, pool, cluster, sources, api_key, api_url, trace_id, turn
    ):
        self.packets.append((name, packet))
        verdict = self.replies[name]
        if verdict is None:
            return None
        if isinstance(verdict, list):
            verdict = verdict[min(len(self.packets) - 1, len(verdict) - 1)]
        return {
            "verdict": verdict,
            "conflicting_markers": ["PTPRC"],
            "reason": f"scripted {name} verdict",
        }


def _run(monkeypatch, pool, replies, sources=None, judges=None, cluster="1"):
    script = _Script(replies)
    judges = judges if judges is not None else _Judges()
    monkeypatch.setattr(agents, "_turn", script)
    monkeypatch.setattr(agents, "_judge", judges)
    outcome = agents.annotate_cluster(
        pool, cluster, agents._prompt("cluster_annotator"), "k", "u", sources
    )
    return outcome, script


def test_an_exclusion_binding_the_pick_comes_back_before_the_answer_is_accepted(
    monkeypatch, pool_with_detected_exclusion
):
    """The audit runs while the agent can still answer it, not after the answer is fixed.

    PTPRC is curated negative for enterocyte and detected in 66.1% of these cells. The
    first answer selects enterocyte, so the row is returned and the agent answers again --
    here unchanged, which is a legitimate outcome: the check asks for the reading, not for
    a different call.
    """
    pool, _sources = pool_with_detected_exclusion
    silent = _final_text()
    addressed = _final_text(
        reason="FABP2 at 92.4 against 68.5 is a dominant population; PTPRC at 66.1 "
               "against 72.7 is detected across the tissue and does not refute it."
    )
    outcome, script = _run(monkeypatch, pool, [silent, addressed])

    assert outcome["final"]["selected"] == "enterocyte"
    assert outcome["turns"] == 2
    assert _audit_genes(outcome) == ["PTPRC"]
    assert outcome["exclusion_audit"][0]["role"] == "selected"
    assert outcome["exclusion_audit"][0]["detected_exclusions"][0]["pct_in"] == 66.1
    # Returned as an observation on the running prompt, and only the second turn sees it.
    # The key in its minified JSON form, since the prompt text names the field too.
    assert '"exclusion_audit":[' not in script.prompts[0]
    assert '"exclusion_audit":[' in script.prompts[1]
    # Not a tool call: the agent never asked for it, so it does not enter the transcript.
    assert outcome["transcript"] == []


def test_the_audit_is_raised_once_so_a_cluster_costs_at_most_one_extra_turn(
    monkeypatch, pool_with_detected_exclusion
):
    """A second silent answer is accepted rather than looping on the same observation."""
    pool, _sources = pool_with_detected_exclusion
    outcome, script = _run(monkeypatch, pool, [_final_text(), _final_text()])

    assert outcome["turns"] == 2
    assert outcome["final"]["selected"] == "enterocyte"
    assert _audit_genes(outcome) == ["PTPRC"]


def test_naming_the_gene_does_not_by_itself_skip_the_audit(
    monkeypatch, pool_with_detected_exclusion
):
    """An exclusion quoted against a RIVAL is the failure this check exists for.

    Whether the reason mentions the gene is deliberately not consulted, because a text
    match cannot separate weighing PTPRC against the chosen identity from citing it to
    reject a different one.
    """
    pool, _sources = pool_with_detected_exclusion
    against_rival = _final_text(
        reason="Tuft cell is rejected because PTPRC at 66.1 contradicts it."
    )
    outcome, _script = _run(monkeypatch, pool, [against_rival, _final_text()])

    assert outcome["turns"] == 2
    assert _audit_genes(outcome) == ["PTPRC"]


def test_a_demoted_identity_carrying_no_exclusion_says_so_explicitly(
    monkeypatch, pool_with_detected_exclusion
):
    """The asymmetry is the point, so the empty side has to be on the page.

    A pick with four detected exclusions reads as ordinary until the identity it demoted
    is shown carrying none. Reporting only the selected side leaves that second fact to
    be inferred from silence, and silence is also what a missing measurement looks like.
    """
    pool, _sources = pool_with_detected_exclusion
    with_rival = _final_text(
        co_occurring_identities=["tuft cell"],
        claim_evidence=[
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
    )
    outcome, script = _run(monkeypatch, pool, [with_rival, with_rival])

    assert outcome["turns"] == 2
    assert [block["role"] for block in outcome["exclusion_audit"]] == [
        "selected",
        "co_occurring",
    ]
    assert _audit_genes(outcome, "enterocyte") == ["PTPRC"]
    assert _audit_genes(outcome, "tuft cell") == []
    # An empty list, present in the observation rather than left out of it.
    assert '"identity":"tuft cell","role":"co_occurring"' in script.prompts[1]
    assert '"detected_exclusions":[]' in script.prompts[1]


def test_an_undetected_exclusion_is_never_raised(monkeypatch, pool):
    """The stock fixture holds PTPRC at 0.0: nothing to answer, so no extra turn."""
    outcome, _script = _run(monkeypatch, pool, [_final_text()])

    assert outcome["turns"] == 1
    assert outcome["exclusion_audit"] == []


def test_a_query_is_executed_and_only_its_result_is_appended(monkeypatch, pool):
    """The opening packet stays byte-identical from turn to turn.

    Everything before the appended observation is the cacheable prefix, and it is also
    the whole of the evidence the agent was given; rebuilding it each turn would both
    destroy the prefix and let the packet drift mid-conversation.
    """
    query = json.dumps(
        {"action": "query", "tool": "gene_across_clusters", "args": {"gene": "FABP2"}}
    )
    outcome, script = _run(monkeypatch, pool, [query, _final_text()])

    assert outcome["turns"] == 2
    assert outcome["transcript"] == [
        {
            "turn": 1,
            "call": {"tool": "gene_across_clusters", "args": {"gene": "FABP2"}},
            "duplicate": False,
        }
    ]
    first, second = script.prompts
    assert second.startswith(first)
    appended = second[len(first) :]
    assert "# Observation 1" in appended
    assert '"pct_in":92.4' in appended
    assert "cells_in_cluster" not in appended


def test_a_repeated_call_returns_the_same_observation_and_then_asks_for_an_answer(
    monkeypatch, pool
):
    """Spinning is the failure a turn budget was supposed to prevent, so it is the brake.

    A budget cuts off a cluster still making progress. Re-asking an identical question
    cannot make progress by construction, so that is what the loop counts: the repeat is
    answered from what was already fetched, and after a few in a row the agent is told
    to answer.
    """
    query = json.dumps(
        {"action": "query", "tool": "candidates_with_gene", "args": {"gene": "PROX1"}}
    )
    outcome, script = _run(
        monkeypatch,
        pool,
        [query, query, query, query, query, _final_text()],
    )

    assert "final" in outcome
    assert outcome["turn_budget_exhausted"] is True
    assert [item["duplicate"] for item in outcome["transcript"]] == [
        False,
        True,
        True,
        True,
    ]
    assert '"duplicate_query":true' in script.prompts[-1]
    assert "Return the final answer now" in script.prompts[-1]


def test_a_query_after_the_agent_was_told_to_answer_does_not_run_a_tool(
    monkeypatch, pool
):
    query = json.dumps(
        {"action": "query", "tool": "candidates_with_gene", "args": {"gene": "PROX1"}}
    )
    outcome, _script = _run(
        monkeypatch, pool, [query, query, query, query, query, query, _final_text()]
    )
    # Four calls were recorded before the directive; the fifth and sixth ran no tool.
    assert len(outcome["transcript"]) == 4


def test_malformed_output_is_retried_a_bounded_number_of_times(monkeypatch, pool):
    outcome, script = _run(monkeypatch, pool, ["not json at all"])
    assert outcome["error"] == "schema violation"
    assert len(script.prompts) == agents.SCHEMA_RETRIES + 1
    assert "# Retry 1" in script.prompts[1]


def test_an_answer_that_misquotes_a_measurement_is_rejected_like_malformed_output(
    monkeypatch, pool
):
    """A number the agent typed rather than read is not an answer this cluster can keep."""
    bad = _final_text(
        claim_evidence=[
            {
                "identity": "enterocyte",
                "decisive_gene": "FABP2",
                "pct_in": 12.0,
                "pct_out": 68.5,
            }
        ]
    )
    outcome, _script = _run(monkeypatch, pool, [bad])
    assert outcome["error"] == "schema violation"


def test_a_recovered_answer_after_one_bad_turn_still_succeeds(monkeypatch, pool):
    outcome, _script = _run(monkeypatch, pool, ["{}", _final_text()])
    assert "final" in outcome
    assert outcome["turns"] == 2


def test_missing_credentials_stop_the_conversation_immediately(monkeypatch, pool):
    def no_credentials(prompt, api_key, api_url, trace_id, turn):
        return None, "cache_miss_no_credentials", False

    monkeypatch.setattr(agents, "_turn", no_credentials)
    outcome = agents.annotate_cluster(
        pool, "1", agents._prompt("cluster_annotator"), "", ""
    )
    assert outcome["error"] == "cache_miss_no_credentials"
    assert outcome["turns"] == 1


def _record(pool, cluster, final):
    import pandas as pd

    outcome = {"final": final, "turns": 1, "transcript": []}
    return agents._result(cluster, pool, {}, outcome, pd.DataFrame())


def _answer(**overrides):
    value = json.loads(_final_text())
    value.update(overrides)
    return value


def test_a_subtype_without_its_own_markers_falls_back_to_the_parent(pool):
    """A finer name on the parent's own genes is the parent, so the parent is reported.

    tuft cell raises POU2F3 and PROX1 in this cluster, but PROX1 is curated for the
    parent too, leaving one exclusive marker where two are required. The reported
    identity is the level the measurements separate.
    """
    record = _record(
        pool,
        "2",
        _answer(
            selected="enteroendocrine cell",
            subtype="tuft cell",
            claim_evidence=[
                {
                    "identity": "enteroendocrine cell",
                    "decisive_gene": "CHGA",
                    "pct_in": 91.0,
                    "pct_out": 1.8,
                },
                {
                    "identity": "tuft cell",
                    "decisive_gene": "POU2F3",
                    "pct_in": 89.3,
                    "pct_out": 4.6,
                },
            ],
            support_markers=[],
        ),
    )
    assert record["annotation"] == "enteroendocrine cell"
    assert record["subtype"] == ""
    # The rejected name keeps its role and its evidence: why the cluster is NOT the finer
    # type has to be as checkable as why it is the parent.
    roles = {
        entry["cell_type"]: entry["claim_role"] for entry in record["candidate_entries"]
    }
    assert roles["tuft cell"] == "rejected_subtype"
    assert record["claim_warnings"] == [
        "rejected_subtype tuft cell: exclusive_defining_markers 1/2"
        " | raised panel is shared with enteroendocrine cell | exclusive: POU2F3"
    ]


def test_a_subtype_with_its_own_markers_is_kept(pool):
    """The gate is not a retreat from finer calls: enough exclusive markers and it stands."""
    record = _record(
        pool,
        "2",
        _answer(
            selected="enterocyte",
            subtype="enteroendocrine cell",
            claim_evidence=[
                {
                    "identity": "enterocyte",
                    "decisive_gene": "FABP2",
                    "pct_in": 10.0,
                    "pct_out": 80.0,
                },
                {
                    "identity": "enteroendocrine cell",
                    "decisive_gene": "CHGA",
                    "pct_in": 91.0,
                    "pct_out": 1.8,
                },
            ],
            support_markers=[],
        ),
    )
    assert record["subtype"] == "enteroendocrine cell"
    roles = {
        entry["cell_type"]: entry["claim_role"] for entry in record["candidate_entries"]
    }
    assert roles["enteroendocrine cell"] == "subtype"
    assert not any("rejected_subtype" in line for line in record["claim_warnings"])


def test_the_record_keeps_the_whole_pool_and_flags_the_claim(monkeypatch, pool):
    """Every candidate stays on the record, and the warning names what was not raised.

    Nothing filters the pool now, so the reason a cluster was not called something is
    the candidate's own panel sitting beside the one that was chosen.
    """
    import pandas as pd

    outcome, _script = _run(monkeypatch, pool, [_final_text()])
    record = agents._result("1", pool, {}, outcome, pd.DataFrame())

    assert record["annotation"] == "enterocyte"
    assert record["annotation_source"] == "cluster_annotation"
    assert record["resolution_status"] == agents.RESOLVED
    assert [entry["cell_type"] for entry in record["candidate_entries"]] == [
        "enterocyte",
        "enteroendocrine cell",
        "tuft cell",
    ]
    roles = {
        entry["cell_type"]: entry["claim_role"] for entry in record["candidate_entries"]
    }
    assert roles == {
        "enterocyte": "selected",
        "enteroendocrine cell": "",
        "tuft cell": "",
    }
    assert record["claim_warnings"] == [
        "selected enterocyte: not_raised 2/4 | top_n_pub: "
        "VIL1(70.8/72.1),FABP2(92.4/68.5) | +0 more"
    ]
    # The published tables read fractions; the agent read percent. The scale changes
    # once, here, and the marker rows are the raised ones.
    measurements = record["candidate_entries"][0]["decisive_marker_measurements"]
    assert [row["gene"] for row in measurements] == ["SI", "RBP2"]
    assert measurements[0]["detection_fraction_in"] == 0.865


def _final_with_subtype():
    """A valid answer naming a refinement, for cluster 2.

    A subtype has to sit in the SAME identity group as its parent, which is what makes
    it a refinement rather than a second population -- so the pair is the secretory one,
    and cluster 2 is where both of them are actually raised.
    """
    return _final_text(
        selected="enteroendocrine cell",
        subtype="tuft cell",
        support_markers=["CHGA"],
        claim_evidence=[
            {
                "identity": "enteroendocrine cell",
                "decisive_gene": "CHGA",
                "pct_in": 91.0,
                "pct_out": 1.8,
            },
            {
                "identity": "tuft cell",
                "decisive_gene": "POU2F3",
                "pct_in": 89.3,
                "pct_out": 4.6,
            },
        ],
        reason="CHGA at 91.0 against 1.8, with POU2F3 raised on top of it.",
    )


# ---------------------------------------------------------------- quality control ----
def _qc_of(outcome):
    return outcome.get("quality_control") or {}


def test_a_failing_verdict_is_carried_into_the_conversation_rather_than_only_recorded(
    monkeypatch, pool
):
    """An error found and not passed on is an error nobody acted on.

    This is the whole reason quality control is inside the loop instead of beside it, so
    it is asserted on the prompt the agent actually received: the verdict, the reason and
    the genes have to be in it, and the packet must NOT be in it a second time.
    """
    judges = _Judges(annotation=["too_specific", "supported"])
    outcome, script = _run(monkeypatch, pool, [_final_text()], judges=judges)
    assert len(script.prompts) == 2
    fed_back = script.prompts[1][len(script.prompts[0]):]
    assert "quality_control" in fed_back
    assert "too_specific" in fed_back
    assert "scripted annotation_judge verdict" in fed_back
    assert "PTPRC" in fed_back
    # Tightened, not replayed: the candidates were already sent once and are not resent.
    assert "unmeasured_curated_genes" not in fed_back
    assert _qc_of(outcome)["passed"] is True
    assert _qc_of(outcome)["selected_round"] == 2


def test_the_rounds_are_bounded_so_a_judge_checks_rather_than_negotiates(
    monkeypatch, pool
):
    judges = _Judges(annotation="too_specific")
    outcome, script = _run(monkeypatch, pool, [_final_text()], judges=judges)
    assert _qc_of(outcome)["rounds"] == agents.JUDGE_MAX_ROUNDS
    assert len(script.prompts) == agents.JUDGE_MAX_ROUNDS


def test_the_delivered_answer_is_the_best_rated_round_not_the_last(monkeypatch, pool):
    """Two rounds, the second worse than the first: the first is what ships.

    Answering again is allowed to fail. Delivering whichever answer happened to come
    last would make a second opinion a risk rather than a check.
    """
    judges = _Judges(annotation=["too_specific", "contradicted"])
    outcome, _script = _run(monkeypatch, pool, [_final_text()], judges=judges)
    assert _qc_of(outcome)["selected_round"] == 1
    assert _qc_of(outcome)["rank"] == 2


def test_no_reachable_judge_leaves_the_answer_alone_rather_than_coarsening_it(
    monkeypatch, pool
):
    """An offline replay must not demote every refinement it could not check."""
    judges = _Judges(subtype=None, annotation=None)
    outcome, script = _run(
        monkeypatch, pool, [_final_with_subtype()], judges=judges, cluster="2"
    )
    assert len(script.prompts) == 1
    assert _qc_of(outcome)["checked"] is False
    assert _qc_of(outcome)["demoted_subtype"] == ""


def test_the_judges_are_never_shown_the_candidates_together(monkeypatch, pool):
    """A judge able to see the pool could prefer another candidate, which is the one
    thing it is not for. It is handed one label, its own panel and its own sentences."""
    judges = _Judges()
    _run(monkeypatch, pool, [_final_with_subtype()], judges=judges, cluster="2")
    assert [name for name, _packet in judges.packets] == [
        "subtype_judge",
        "annotation_judge",
    ]
    for _name, packet in judges.packets:
        assert packet["label_under_test"]
        assert set(packet) <= {
            "query",
            "label_under_test",
            "panel",
            "unmeasured_curated_genes",
            "sources",
            "parent_label",
            "parent_panel",
            "raised_and_absent_from_parent_panel",
        }
    subtype_packet = judges.packets[0][1]
    assert subtype_packet["label_under_test"] == "tuft cell"
    assert subtype_packet["parent_label"] == "enteroendocrine cell"
    assert subtype_packet["raised_and_absent_from_parent_panel"] == ["POU2F3"]
    # The refinement passed, so the label that goes on to be checked is the finer one.
    assert judges.packets[1][1]["label_under_test"] == "tuft cell"


class _JudgeModel:
    """A judge model returning scripted replies, recording every prompt it was sent."""

    def __init__(self, replies):
        self.replies = list(replies)
        self.prompts = []

    def __call__(self, prompt, api_url, api_key, **_kwargs):
        self.prompts.append(prompt)
        reply = self.replies[min(len(self.prompts) - 1, len(self.replies) - 1)]
        return reply, None, False


def test_a_judge_can_ask_for_a_sentence_the_packet_did_not_carry(monkeypatch, pool):
    """Code decided what evidence a judge had, and it decided wrong.

    A panel cut at the fourteen best-published rows left 10% of the genes the answers
    rested on off the judge's page entirely, and it had no way to ask. It now reads the
    whole panel and asks for the sentences past the head -- one growing prompt, so the
    packet is never re-sent and the head a provider caches stays byte-identical.
    """
    from conftest import FakeSources

    model = _JudgeModel(
        [
            json.dumps({"action": "query", "tool": "sources", "args": {"genes": ["SI"]}}),
            json.dumps({"verdict": "supported", "reason": "SI is defining and raised."}),
        ]
    )
    monkeypatch.setattr(agents.llm, "cached_call_llm", model)
    packet = agents.pool_api.judge_packet(pool, "1", "enterocyte", FakeSources())

    verdict = agents._judge(
        "annotation_judge", packet, pool, "1", FakeSources(), "k", "u", "t", 1
    )

    assert verdict["verdict"] == "supported"
    assert len(model.prompts) == 2
    # The observation is appended and the head is untouched, so the second prompt extends
    # the first rather than replacing it.
    assert model.prompts[1].startswith(model.prompts[0])
    assert '"tool":"sources"' in model.prompts[1]
    # Answered for the label under test, so the sentence it gets back is that label's.
    assert '"sentence":"SI is absent from enterocyte."' in model.prompts[1]


def test_a_reply_carrying_a_question_and_its_own_answer_is_read_as_the_question(
    monkeypatch, pool
):
    """Two objects in one reply used to parse as neither, and cost the cluster its check.

    Observed on human_liver c24: the judge asked for KRT19's sentence and then, in the same
    reply, answered without it. Matching from the first brace to the last swallowed both,
    parsed as nothing, and the label was filed `not_checked` -- which reads as "no judge was
    reachable" when a judge had plainly replied. The first object wins, so the question is
    honoured and the answer that skipped the evidence it asked for is dropped.
    """
    from conftest import FakeSources

    both = (
        '{"action":"query","tool":"sources","args":{"genes":["SI"]}}\n'
        '{"verdict":"supported","reason":"answered without waiting"}'
    )
    model = _JudgeModel(
        [both, json.dumps({"verdict": "supported", "reason": "SI is defining and raised."})]
    )
    monkeypatch.setattr(agents.llm, "cached_call_llm", model)
    packet = agents.pool_api.judge_packet(pool, "1", "enterocyte", FakeSources())

    verdict = agents._judge(
        "annotation_judge", packet, pool, "1", FakeSources(), "k", "u", "t", 1
    )

    assert verdict["reason"] == "SI is defining and raised."
    assert len(model.prompts) == 2
    # The tool ran, so the second answer was made with the sentence rather than without it.
    assert '"sentence":"SI is absent from enterocyte."' in model.prompts[1]


def test_a_judge_reply_that_is_neither_is_asked_again_not_filed_as_unchecked(
    monkeypatch, pool
):
    """`not_checked` has to keep meaning "no judge was reachable".

    A label reported as never checked because one reply was malformed tells a reader the
    check did not run, which is a different fact about the run.
    """
    from conftest import FakeSources

    model = _JudgeModel(
        ["not json at all", json.dumps({"verdict": "supported", "reason": "second try."})]
    )
    monkeypatch.setattr(agents.llm, "cached_call_llm", model)
    packet = agents.pool_api.judge_packet(pool, "1", "enterocyte", FakeSources())

    verdict = agents._judge(
        "annotation_judge", packet, pool, "1", FakeSources(), "k", "u", "t", 1
    )
    assert verdict["reason"] == "second try."
    assert "# Retry 1" in model.prompts[1]

    # And it still gives up rather than retrying forever.
    model = _JudgeModel(["not json at all"])
    monkeypatch.setattr(agents.llm, "cached_call_llm", model)
    assert (
        agents._judge(
            "annotation_judge", packet, pool, "1", FakeSources(), "k", "u", "t", 1
        )
        is None
    )
    assert len(model.prompts) == agents.SCHEMA_RETRIES + 1


def test_a_judges_questions_are_bounded_so_it_checks_rather_than_stalls(monkeypatch, pool):
    """A check that could ask forever would be negotiating. It is told to answer, once."""
    from conftest import FakeSources

    query = json.dumps({"action": "query", "tool": "sources", "args": {"genes": ["SI"]}})
    model = _JudgeModel([query, query, query, query])
    monkeypatch.setattr(agents.llm, "cached_call_llm", model)
    packet = agents.pool_api.judge_packet(pool, "1", "enterocyte", FakeSources())

    verdict = agents._judge(
        "annotation_judge", packet, pool, "1", FakeSources(), "k", "u", "t", 1
    )

    assert verdict is None
    # Exactly the budgeted number of questions were answered, and no more.
    assert model.prompts[-1].count("# Observation") == agents.JUDGE_QUERY_TURNS
    assert "No further queries" in model.prompts[-1]
    # Asking `sources` again is a request for the sentences it has not seen, not a repeat
    # of the ones it has, so it is answered rather than marked duplicate. The query budget
    # is what bounds it.
    assert '"duplicate_query":true' not in model.prompts[-1]
    # Past the budget it is nudged for a verdict a bounded number of times, then dropped.
    assert len(model.prompts) == agents.JUDGE_QUERY_TURNS + 1 + agents.SCHEMA_RETRIES


def test_an_identity_retrieval_never_offered_is_reported_and_never_answered(
    monkeypatch, pool
):
    """The menu is 15 of a median 67 admitted cell types, so a miss has to be sayable.

    Sayable, and only that: a name with no supplied panel has no measurement to quote, so
    it cannot stand in an identity column beside labels that were verified. Naming a
    supplied candidate here is rejected -- that one could just be answered with.
    """
    outcome, _script = _run(
        monkeypatch, pool, [_final_text(unlisted_identity="M cell")]
    )
    assert outcome["final"]["selected"] == "enterocyte"
    assert outcome["final"]["unlisted_identity"] == "M cell"

    # A candidate belongs in `selected`, not here: the answer is rejected and the agent is
    # asked again rather than the contradiction being filed.
    assert not agents.valid_final(
        json.loads(_final_text(unlisted_identity="tuft cell")), pool, "1"
    )
    assert agents.valid_final(
        json.loads(_final_text(unlisted_identity="M cell")), pool, "1"
    )
    # Absent is the ordinary case, and an old answer that never carried the field is still
    # a valid answer.
    assert agents.valid_final(json.loads(_final_text()), pool, "1")

    # And it reaches the record as its own field, because "how often was the answer not on
    # the menu" is a rate somebody has to be able to count.
    import pandas as pd

    record = agents._result("1", pool, {}, outcome, pd.DataFrame())
    assert record["unlisted_identity"] == "M cell"
    assert record["annotation"] == "enterocyte"


def _demotion_outcome():
    return {
        "final": _answer(
            selected="enterocyte",
            subtype="tuft cell",
            claim_evidence=[
                {
                    "identity": "enterocyte",
                    "decisive_gene": "FABP2",
                    "pct_in": 92.4,
                    "pct_out": 68.5,
                },
                {
                    "identity": "tuft cell",
                    "decisive_gene": "PROX1",
                    "pct_in": 77.7,
                    "pct_out": 21.7,
                },
            ],
        ),
        "turns": 1,
        "transcript": [],
        "quality_control": {
            "checked": True,
            "selected_round": 1,
            "judged_label": "enterocyte",
            "demoted_subtype": "tuft cell",
            "subtype_verdict": {
                "verdict": "not_separable",
                "reason": "PROX1 is curated for the parent as well.",
            },
            "annotation_verdict": {"verdict": "supported"},
        },
    }


def test_a_refinement_the_judge_cannot_separate_is_reported_as_the_parent(pool):
    """The demotion has to reach the record, not just the verdict log."""
    import pandas as pd

    record = agents._result("1", pool, {}, _demotion_outcome(), pd.DataFrame())
    assert record["subtype"] == ""
    assert record["annotation"] == "enterocyte"
    assert record["annotation_qc"] == agents.QC_DEMOTED
    assert any(
        line.startswith("rejected_subtype tuft cell") for line in record["claim_warnings"]
    )
    roles = {
        entry["cell_type"]: entry["claim_role"] for entry in record["candidate_entries"]
    }
    assert roles["tuft cell"] == "rejected_subtype"


def test_the_line_says_which_of_the_two_checks_dropped_the_refinement(
    monkeypatch, pool
):
    """Both checks can drop a refinement, and the record has to say which one did.

    Here the free gate is satisfied, so the only thing left to drop `tuft cell` is the
    judge -- and the warning has to carry its verdict rather than a marker count that
    was never the reason.
    """
    import pandas as pd

    monkeypatch.setattr(agents, "_subtype_demotion", lambda *args: ("", ""))
    record = agents._result("1", pool, {}, _demotion_outcome(), pd.DataFrame())
    assert record["subtype"] == ""
    assert record["annotation_qc"] == agents.QC_DEMOTED
    assert any(
        line.startswith("rejected_subtype tuft cell: quality_control not_separable")
        and "PROX1 is curated for the parent as well." in line
        for line in record["claim_warnings"]
    )


def test_the_column_separates_a_check_that_failed_from_one_that_never_ran():
    """`not_checked` cannot be spelled the same as `passed`, or a replay with no judge
    would report every cluster as having survived a check nobody ran."""
    unchecked = {"checked": False}
    assert agents._qc_status(unchecked, False) == agents.QC_UNCHECKED
    assert agents._qc_status({}, False) == agents.QC_UNCHECKED
    carried = {
        "checked": True,
        "selected_round": 1,
        "judged_label": "enterocyte",
        "annotation_verdict": {"verdict": "supported"},
    }
    assert agents._qc_status(carried, False) == agents.QC_PASSED
    assert agents._qc_status({**carried, "selected_round": 2}, False) == agents.QC_REVISED
    assert agents._qc_status(carried, True) == agents.QC_DEMOTED
    refused = {**carried, "annotation_verdict": {"verdict": "contradicted"}}
    assert agents._qc_status(refused, False) == agents.QC_FAILED


def test_a_thin_panel_is_not_asked_about_twice(monkeypatch, pool):
    """`insufficient_evidence` buys nothing on a second pass.

    The panel that was too thin to decide is the same panel the next answer reads, so
    another round is a turn spent to be told the same thing. The verdict is still
    recorded -- it just does not trigger a retry.
    """
    judges = _Judges(annotation="insufficient_evidence")
    outcome, script = _run(monkeypatch, pool, [_final_text()], judges=judges)
    assert len(script.prompts) == 1
    assert _qc_of(outcome)["rounds"] == 1
    assert _qc_of(outcome)["retry"] is False
    assert _qc_of(outcome)["annotation_verdict"]["verdict"] == "insufficient_evidence"


def test_a_refinement_dropped_to_a_supported_parent_costs_no_extra_turn(
    monkeypatch, pool
):
    """Falling back to the parent IS the answer, so nothing is re-opened.

    The refinement was refused and the parent it fell back to is carried by the evidence.
    Spending another annotator turn there would re-litigate a call the evidence already
    supports, which is the difference between a check and a haggle.
    """
    judges = _Judges(subtype="not_separable", annotation="supported")
    outcome, script = _run(
        monkeypatch, pool, [_final_with_subtype()], judges=judges, cluster="2"
    )
    assert len(script.prompts) == 1
    review = _qc_of(outcome)
    assert review["demoted_subtype"] == "tuft cell"
    assert review["judged_label"] == "enteroendocrine cell"
    assert review["passed"] is True
    assert review["retry"] is False
    # The column still says a refinement was dropped, even though nothing was retried.
    assert agents._qc_status({**review, "selected_round": 1}, True) == agents.QC_DEMOTED
