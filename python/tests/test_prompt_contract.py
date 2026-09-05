from __future__ import annotations

import json
from pathlib import Path

from scmarkeragent import cluster_annotation as agents
from scmarkeragent import annotator_pool as pool_api

ROOT = Path(__file__).resolve().parents[1]
PROMPTS = ROOT / "src/scmarkeragent/prompts"


def _final(**overrides):
    """A well-formed answer for cluster 1 of the fixture dataset."""
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
        "reason": "FABP2 is detected in 92.4% of these cells.",
    }
    value.update(overrides)
    return value


def test_only_one_prompt_ever_sees_the_candidates_together():
    """One conversation decides one cluster; the others only check what it decided.

    The rule this protects was never "one file". It is that exactly one reading gets to
    compare the candidates, because a reading that happens BEFORE that comparison can
    only narrow what the comparison is allowed to conclude. The quality-control templates
    run after the answer exists and are handed one label with its own panel -- no
    `candidates` list, no retrieval order, and nothing the annotator said about why. They
    can send a question back; they cannot answer it.
    """
    assert sorted(path.name for path in PROMPTS.glob("*.txt")) == [
        "annotation_judge.txt",
        "cluster_annotator.txt",
        "cross_species_note.txt",
        "panel_reading.txt",
        "subtype_judge.txt",
    ]
    decides = agents._prompt("cluster_annotator")
    assert "`candidates` are the curated cell types retrieved" in decides
    for name in ("subtype_judge", "annotation_judge"):
        # The rendered template, so the shared block is checked where it is actually read.
        text = agents._prompt(name)
        assert "label_under_test" in text
        # Neither the pool nor the ranking, so a judge cannot become a second annotator.
        assert "`candidates`" not in text
        assert "retrieval" not in text
        # The judgment is bound to the supplied evidence, and there is no author label in
        # this package for it to reach for even if it were not.
        assert "Nothing else is supplied, and nothing else may be used." in text


def test_the_packet_is_the_last_thing_in_the_template():
    """Everything before the packet is identical across every call.

    That fixed head, plus the opening packet, is the only part a provider can serve from
    a cached prefix -- and it is why a tool result is APPENDED rather than folded back
    into the packet.
    """
    text = (PROMPTS / "cluster_annotator.txt").read_text(encoding="utf-8")
    head, placeholder, tail = text.partition("{{EVIDENCE_PACKET_JSON}}")
    assert placeholder
    assert not tail.strip()
    assert len(head) > 1000


def test_the_template_names_exactly_the_tools_the_dispatcher_implements():
    text = (PROMPTS / "cluster_annotator.txt").read_text(encoding="utf-8")
    for tool in pool_api.TOOLS:
        assert tool in text
    for retired in ("gene_stats", "panel_detail", "veto", "independent review"):
        assert retired not in text


def test_each_judge_can_ask_and_is_offered_only_the_tools_it_may_call():
    """A judge that cannot ask leaves code as the last word on what evidence exists.

    Which is how a fourteen-row panel cut came to hide 10% of the genes the answers rested
    on with nothing in the record to show it. The template has to offer the query turn, and
    it has to offer exactly the tools the judge dispatcher will answer -- naming one it
    cannot call spends a turn on an error.
    """
    for name in ("annotation_judge", "subtype_judge"):
        text = agents._prompt(name)
        flowed = " ".join(text.split())
        assert '"action":"query"' in text
        assert '"verdict"' in text
        for tool in pool_api.JUDGE_TOOLS:
            assert tool in text
        for withheld in set(pool_api.TOOLS) - set(pool_api.JUDGE_TOOLS):
            assert withheld not in text
        # And the panel it is handed is whole, so a sentence it has not read is the only
        # thing left to ask for.
        assert "WHOLE curated definition" in flowed
        assert "does NOT cover the whole panel" in flowed


def test_the_annotator_can_report_a_menu_it_could_not_answer_from():
    """Retrieval hands over 15 of a median 67 admitted cell types.

    Before `pool_search` and this field, a menu that missed the answer was indistinguishable
    from one that did not: `selected` was always a candidate, `Unknown` was returned 0/109,
    and no tool could see past the cut. The template has to offer both halves -- the way to
    look, and the field to say what was found -- and has to say it is not an answer.
    """
    text = agents._prompt("cluster_annotator")
    assert "pool_search" in text
    assert "unlisted_identity" in text
    assert "not a candidate" in text and "cannot be `selected`" in text
    # The schema line the agent copies has to carry the field, or it will never fill it.
    assert '"unlisted_identity":' in text


def test_no_template_invites_a_reference_label_into_the_output():
    """The zero-leak contract, enforced on the one field that could carry a breach.

    `state` is free text, which makes it the only place a reference label could enter a
    reference-free package -- and the annotator template used to say in as many words that
    it could be "a study's own cluster label". Nothing had leaked when this was found (546
    populated states across the 809 published clusters, none of them a cluster id), which
    is exactly why it needed a test: the invitation was there, and only the models' habits
    were keeping it shut.

    The check is on the templates rather than on output, because output only shows what a
    given model happened to do on a given day.
    """
    invitations = (
        "study's own cluster label",
        "the study's own label",
        "author annotation",
        "author label",
        "reference label",
        "ground truth",
        "original cluster label",
    )
    offences = []
    for path in sorted(PROMPTS.glob("*.txt")):
        text = path.read_text(encoding="utf-8")
        lowered = text.lower()
        for phrase in invitations:
            index = lowered.find(phrase.lower())
            while index != -1:
                line_no = text.count("\n", 0, index) + 1
                line = text.splitlines()[line_no - 1]
                # Saying one is NOT supplied, or must not be used, is the contract being
                # stated -- the opposite of an invitation. Only a permissive mention counts.
                permissive = not any(
                    marker in line.lower()
                    for marker in ("no ", "not ", "never", "must not", "cannot", "none")
                )
                if permissive:
                    offences.append(f"{path.name}:{line_no}: {line.strip()}")
                index = lowered.find(phrase.lower(), index + 1)
    assert not offences, "a template invites a reference label: " + "; ".join(offences)


def test_the_template_states_both_readings_of_an_unraised_marker():
    """Without this the agent cannot use the panels it is given.

    A high `pct_in` beside a high `pct_out` is a dominant population, not an absent
    identity, and the two are the same `avg_log2FC <= 0` on the page.
    """
    text = (PROMPTS / "cluster_annotator.txt").read_text(encoding="utf-8")
    assert "avg_log2FC > 0" in text
    assert "near zero is evidence AGAINST" in text
    assert "much of\nthe dataset drives its own defining markers to a flat fold change" in text


def test_every_template_that_reads_a_panel_reads_it_the_same_way():
    """One reading of a marker row, or the templates argue past each other.

    The judges were given the annotator's FIELDS without the annotator's READING: they
    used "raised" as the pivot of every rule and never defined it, and they carried one
    reading of an unraised marker where the annotator has three. A label's own
    best-published marker could therefore count against it in one template and be argued
    away in another, and the difference showed up as a verdict that moved when an
    unrelated rule was reworded. The block is one file; this asserts the annotator holds
    it verbatim and that both judges actually receive it.
    """
    shared = (PROMPTS / "panel_reading.txt").read_text(encoding="utf-8").strip()
    assert shared.startswith("[R1]")
    # Held verbatim by the annotator, so including it elsewhere cannot fork the reading.
    assert shared in (PROMPTS / "cluster_annotator.txt").read_text(encoding="utf-8")
    for name in ("cluster_annotator", "subtype_judge", "annotation_judge"):
        rendered = agents._prompt(name)
        assert agents.PANEL_READING_PLACEHOLDER not in rendered
        assert rendered.count(shared) == 1


def test_the_annotation_judge_offers_no_verdict_it_cannot_reach():
    """`too_broad` asks this judge to name an identity it was never shown.

    It is handed one label and told to use nothing else, so a verdict meaning "a narrower
    identity fits better" has no evidence behind it -- and in 241 judgments it was never
    returned. The inverted refinement it was meant to catch is caught by the subtype
    judge, which is given the parent as well and can compare the two.
    """
    text = agents._prompt("annotation_judge")
    assert "too_broad" not in text
    assert "too_broad" not in agents.RETRYABLE_VERDICTS
    assert "too_broad" not in agents._VERDICT_RANK
    for verdict in agents.RETRYABLE_VERDICTS:
        assert verdict in text


def test_both_subtype_gates_require_the_same_number_of_exclusive_markers():
    """A deterministic gate at two and a judge at one is two gates with two thresholds."""
    assert pool_api.SUBTYPE_EXCLUSIVE_REQUIRED == 2
    assert "AT LEAST TWO markers" in agents._prompt("subtype_judge")
    assert "at least two markers" in agents._prompt("cluster_annotator")


def test_rendered_prompt_carries_the_packet_and_no_placeholder(pool):
    packet = pool_api.cluster_packet(pool, "1")
    rendered = agents._render(agents._prompt("cluster_annotator"), packet)
    assert "{{EVIDENCE_PACKET_JSON}}" not in rendered
    assert '"enterocyte"' in rendered
    assert (
        json.dumps(packet, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
        in rendered
    )


def test_an_observation_appends_only_itself(pool):
    """Re-sending the packet each turn would destroy the cacheable prefix and the budget."""
    observation = pool_api.run_tool(pool, "1", "gene_across_clusters", {"gene": "FABP2"})
    block = agents._observation_block(1, observation)
    assert block.startswith("\n\n# Observation 1\n")
    assert "candidates" not in block
    assert '"gene":"FABP2"' in block


def test_a_query_turn_and_a_final_turn_are_told_apart(pool):
    assert agents.valid_query({"action": "query", "tool": "sources", "args": {}})
    assert not agents.valid_query({"action": "final"})
    assert not agents.valid_query({"tool": "sources"})
    assert agents.valid_final(_final(), pool, "1")
    assert not agents.valid_final({**_final(), "action": "query"}, pool, "1")


def test_every_claimed_identity_must_quote_a_gene_from_its_own_panel(pool):
    """A claim supported by a gene the resource never tied to that cell type is not
    supported by the curated evidence, whatever the measurement says."""
    borrowed = _final(
        claim_evidence=[
            {
                "identity": "enterocyte",
                "decisive_gene": "CHGA",
                "pct_in": 2.3,
                "pct_out": 4.1,
            }
        ]
    )
    assert not agents.valid_final(borrowed, pool, "1")

    invented = _final(
        claim_evidence=[
            {
                "identity": "enterocyte",
                "decisive_gene": "MADEUP",
                "pct_in": 1.0,
                "pct_out": 0.0,
            }
        ]
    )
    assert not agents.valid_final(invented, pool, "1")


def test_a_quoted_measurement_must_be_the_measured_one(pool):
    """The one thing the code checks about a number is that the agent did not type it.

    Whether 92.4 against 68.5 supports the call is the agent's judgement. Whether the
    cluster really shows 92.4 is not a judgement at all.
    """
    assert agents.valid_final(
        _final(
            claim_evidence=[
                {
                    "identity": "enterocyte",
                    "decisive_gene": "FABP2",
                    "pct_in": 92.44,
                    "pct_out": 68.5,
                }
            ]
        ),
        pool,
        "1",
    )
    assert not agents.valid_final(
        _final(
            claim_evidence=[
                {
                    "identity": "enterocyte",
                    "decisive_gene": "FABP2",
                    "pct_in": 99.0,
                    "pct_out": 68.5,
                }
            ]
        ),
        pool,
        "1",
    )


def test_a_co_occurring_identity_needs_its_own_entry(pool):
    """Selected and co-occurring are the same kind of assertion, so the same rule holds."""
    missing = _final(
        co_occurring_identities=["enteroendocrine cell"],
        identity_groups=[
            {"identity": "absorptive", "candidates": ["enterocyte"]},
            {"identity": "endocrine", "candidates": ["enteroendocrine cell"]},
            {"identity": "tuft", "candidates": ["tuft cell"]},
        ],
    )
    assert not agents.valid_final(missing, pool, "1")

    supplied = {
        **missing,
        "claim_evidence": missing["claim_evidence"]
        + [
            {
                "identity": "enteroendocrine cell",
                "decisive_gene": "PROX1",
                "pct_in": 77.7,
                "pct_out": 21.7,
            }
        ],
    }
    assert agents.valid_final(supplied, pool, "1")


def test_a_weak_claim_is_flagged_and_never_removed(pool):
    """The audit does not correct the agent, so a supportable-looking claim survives it.

    An enteroendocrine claim quoting PROX1 passes every structural check: PROX1 is on
    that identity's panel, it is raised, and the number is real. What the audit adds is
    the line saying CHGA, the identity's 160-publication marker, sits at 2.3 against
    4.1 -- visible to a reader, and not acted on by the code.
    """
    value = {
        **_final(
            co_occurring_identities=["enteroendocrine cell"],
            identity_groups=[
                {"identity": "absorptive", "candidates": ["enterocyte"]},
                {"identity": "endocrine", "candidates": ["enteroendocrine cell"]},
                {"identity": "tuft", "candidates": ["tuft cell"]},
            ],
        )
    }
    value["claim_evidence"] = value["claim_evidence"] + [
        {
            "identity": "enteroendocrine cell",
            "decisive_gene": "PROX1",
            "pct_in": 77.7,
            "pct_out": 21.7,
        }
    ]
    assert agents.valid_final(value, pool, "1")
    warned = pool_api.claim_warnings(pool, "1", agents.claimed_identities(value))
    assert any(line.startswith("cooc enteroendocrine cell") for _, line in warned)
    assert any("CHGA(2.3/4.1)" in line for _, line in warned)


def test_grouping_must_cover_every_supplied_candidate_exactly_once(pool):
    assert not agents.valid_final(
        _final(identity_groups=[{"identity": "absorptive", "candidates": ["enterocyte"]}]),
        pool,
        "1",
    )
    assert not agents.valid_final(
        _final(
            identity_groups=[
                {"identity": "a", "candidates": ["enterocyte", "tuft cell"]},
                {
                    "identity": "b",
                    "candidates": ["enteroendocrine cell", "tuft cell"],
                },
            ]
        ),
        pool,
        "1",
    )


def test_selected_is_a_supplied_candidate_or_unknown(pool):
    unknown = _final(
        selected=agents.UNKNOWN,
        claim_evidence=[],
        support_markers=[],
    )
    assert agents.valid_final(unknown, pool, "1")
    for retired in ("Uncertain", "Mixed: enterocyte + tuft cell", "plasma cell", ""):
        assert not agents.valid_final(_final(selected=retired), pool, "1")


def test_a_subtype_is_a_finer_name_for_the_same_cell(pool):
    """A name from another group is a different cell, not a finer reading of this one."""
    value = _final(
        selected="enteroendocrine cell",
        subtype="tuft cell",
        claim_evidence=[
            {
                "identity": "enteroendocrine cell",
                "decisive_gene": "PROX1",
                "pct_in": 77.7,
                "pct_out": 21.7,
            }
        ],
        support_markers=[],
        identity_groups=[
            {"identity": "absorptive", "candidates": ["enterocyte"]},
            {
                "identity": "secretory",
                "candidates": ["enteroendocrine cell", "tuft cell"],
            },
        ],
    )
    # tuft cell shares a group here, so the shape is allowed but its own entry is required.
    assert not agents.valid_final(value, pool, "1")
    value["claim_evidence"].append(
        {
            "identity": "tuft cell",
            "decisive_gene": "POU2F3",
            "pct_in": 1.0,
            "pct_out": 40.0,
        }
    )
    assert agents.valid_final(value, pool, "1")


def test_support_markers_must_be_listed_for_the_selected_identity(pool):
    assert not agents.valid_final(_final(support_markers=["CHGA"]), pool, "1")
    assert not agents.valid_final(_final(support_markers=["MADEUP"]), pool, "1")


def test_the_annotator_names_one_model_and_the_key_follows_it():
    from scmarkeragent import llm_client

    assert agents.ANNOTATOR_MODEL
    assert llm_client.resolve_model(agents.ANNOTATOR_MODEL) == agents.ANNOTATOR_MODEL
    assert llm_client._vkey("p", "high", agents.ANNOTATOR_MODEL) != llm_client._vkey(
        "p", "high", "other-model"
    )


def test_the_trajectory_identifies_the_audit_record_and_not_the_cache_key():
    """A replay has to hit on the prompt alone, or the second run asks everything again."""
    from scmarkeragent import llm_client

    assert llm_client._vkey("p", "high", agents.ANNOTATOR_MODEL) == llm_client._vkey(
        "p", "high", agents.ANNOTATOR_MODEL
    )
    source = (
        ROOT / "src/scmarkeragent/llm_client.py"
    ).read_text(encoding="utf-8")
    key_body = source.split("def _vkey(")[1].split("def ")[0]
    assert "trace_id" not in key_body
    assert '"trace_id": trace_id' in source
