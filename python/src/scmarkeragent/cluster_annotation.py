#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Annotate each cluster with one agent that can ask the data questions.

One agent decides one cluster. It opens with every retrieved candidate's COMPLETE
curated panel as measured here, ordered by publication support, and it may then call
read-only tools for the curated source sentences, for a gene's behaviour in the other
clusters, and for which candidates claim a gene. It answers when it is ready.

Nothing is removed before it answers. Reading one candidate on its own cannot foreclose
the comparison that is the whole judgement, because there is no such reading: the agent
sees the candidates together, and every judgement about which identity the cells carry is
made once, in one place, by the only participant that can see the alternatives.

The code decides nothing about biology. It supplies measurements, executes tool calls,
and afterwards checks the numbers the agent quoted against the DE table and flags claims
whose curated markers are not raised. That check neither deletes nor renames a claim: an
error is made visible on the record instead of being silently corrected, because a rule
strong enough to correct one would be strong enough to remove a correct call.

When no model is reachable the stage falls back to the top of the retrieval order and
says so in every field that records how the label was produced.
"""

from __future__ import annotations

import json
import os
import pickle
import sys
import uuid
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import annotator_pool as pool_api
from . import llm_client as llm
from .annotator_pool import build_pool, cluster_packet, is_raised, run_tool
from .config import (
    CLUSTER_ANNOTATION_ENABLED,
    CACHE_DIR as CACHE,
    CROSS_SPECIES,
    DEFAULTS,
    LLM_SETTINGS,
    NOT_AVAILABLE,
    PROMPT_DIR,
)
from .candidate_scoring import UNSUPPORTED
from .marker_sources import SourceDB, SourceServer

_A = DEFAULTS["cluster_annotation"]
ANNOTATOR_MODEL = str(_A["annotator_model"])
ANNOTATOR_EFFORT = str(_A["annotator_reasoning_effort"])
ANNOTATOR_SCHEMA = str(_A["annotator_schema_version"])
SCHEMA_RETRIES = int(_A["schema_retries"])
MAX_TURNS = int(_A["max_turns"])
SOURCES_PER_MARKER = int(_A["sources_per_marker"])
SOURCE_BATCHES_PER_MARKER = int(_A["source_batches_per_marker"])
UNKNOWN = str(_A["unknown_token"])
MAX_RATIONALE = int(_A["max_rationale_chars"])
JUDGE_MODEL = str(_A["judge_model"])
JUDGE_EFFORT = str(_A["judge_reasoning_effort"])
JUDGE_MAX_ROUNDS = int(_A["judge_max_rounds"])
# How many times one judgment may ask for evidence the packet did not carry.
JUDGE_QUERY_TURNS = int(_A["judge_query_turns"])

CONFIDENCE_VALUES = ("high", "medium", "low")

# A quoted detection fraction has to be the measured one. The packet shows a percent to
# one decimal, so a faithful copy differs from it by less than a rounding step; anything
# further apart is a number the agent produced rather than read.
QUOTE_TOLERANCE = 0.15

# The agent may keep asking as long as it is learning something. What it may not do is
# spin: a repeated call with identical arguments returns the same observation, and after
# a few of those in a row it is told to answer. That is the brake, rather than a turn
# budget that would cut off a cluster still making progress.
MAX_REPEATED_QUERIES = 3

RESOLVED = "resolved"
MIXED = "mixed"
UNRESOLVED = "unresolved"

# The panel-reading rules every template shares, and the marker that pulls them in.
PANEL_READING = "panel_reading"
PANEL_READING_PLACEHOLDER = "{{PANEL_READING}}"
# The evidence-provenance note spliced into every template when the cross-species switch
# is ON, and only then: the OFF prompt stays byte-identical to what it always was.
CROSS_SPECIES_NOTE = "cross_species_note"


# ------------------------------------------------------------------ helpers ----
def _prompt(name: str) -> str:
    """One template, with the shared panel-reading rules spliced in.

    How a marker row is read -- what counts as raised, what an unraised one means, how
    publication counts weigh -- is one question, and every template that reads a panel
    has to answer it the same way. Held in one file and included, so the annotator and
    the judges cannot drift into two readings of the same number: that drift is what let
    a label's own best-published marker count against it in one template and be argued
    away in another.

    Under the cross-species switch the pooled evidence deliberately crosses species, and
    every template that reads it is told so in one shared note: the sentence a marker
    carries may name another species, and that is not grounds against a label. Injected
    for the annotator and both quality-control judges alike, so the QC cannot reject on
    species grounds an answer the annotator was told not to reject on species grounds.
    """
    text = Path(PROMPT_DIR, f"{name}.txt").read_text(encoding="utf-8")
    if PANEL_READING_PLACEHOLDER in text:
        shared = Path(PROMPT_DIR, f"{PANEL_READING}.txt").read_text(encoding="utf-8")
        text = text.replace(PANEL_READING_PLACEHOLDER, shared.strip())
    if CROSS_SPECIES:
        note = (
            Path(PROMPT_DIR, f"{CROSS_SPECIES_NOTE}.txt")
            .read_text(encoding="utf-8")
            .strip()
        )
        if "\n# Input\n" in text:
            text = text.replace("\n# Input\n", f"\n{note}\n\n# Input\n", 1)
        elif "EVIDENCE_PACKET_JSON=" in text:
            text = text.replace(
                "EVIDENCE_PACKET_JSON=", f"{note}\n\nEVIDENCE_PACKET_JSON=", 1
            )
    return text


def _round(value, digits=4):
    if value is None:
        return None
    number = float(value)
    return None if not np.isfinite(number) else round(number, digits)


def _clip(text: str, limit: int) -> str:
    value = " ".join(str(text or "").split())
    return value if len(value) <= limit else value[: limit - 1] + "\u2026"


def _render(template: str, packet: dict) -> str:
    return template.replace(
        "{{EVIDENCE_PACKET_JSON}}",
        json.dumps(packet, ensure_ascii=False, sort_keys=True, separators=(",", ":")),
    )


def _observation_block(index: int, observation: dict) -> str:
    """One tool result appended to the running prompt.

    Only this observation is appended: never the opening packet again, and never a
    result the agent did not ask for. The fixed head plus the packet therefore stay
    byte-identical from turn to turn, which is the part a provider can serve from a
    cached prefix.
    """
    payload = json.dumps(
        observation, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    )
    return f"\n\n# Observation {index}\n{payload}\n"


# ---------------------------------------------------------------- answer shape ----
def valid_query(value: Any) -> bool:
    return (
        isinstance(value, dict)
        and value.get("action") == "query"
        and isinstance(value.get("tool"), str)
        and isinstance(value.get("args", {}), dict)
    )


def _groups_cover(value: Any, names: list[str]) -> list[list[str]] | None:
    """The agent's grouping, if it places every supplied candidate exactly once.

    The grouping is what separates "several names for one cell" from "several different
    cells", so an answer that leaves a candidate ungrouped has not done the step the rest
    of the answer rests on.
    """
    groups = value.get("identity_groups")
    if not isinstance(groups, list) or not groups:
        return None
    members: list[list[str]] = []
    grouped: list[str] = []
    for group in groups:
        if not isinstance(group, dict) or not str(group.get("identity", "")).strip():
            return None
        entries = group.get("candidates")
        if not isinstance(entries, list) or not entries:
            return None
        members.append([str(entry) for entry in entries])
        grouped.extend(str(entry) for entry in entries)
    return members if sorted(grouped) == sorted(names) else None


def claimed_identities(value: dict) -> list[tuple[str, str]]:
    """(role, identity) for every identity the answer asserts is present in the cluster.

    Selected, subtype and each co-occurring name are one kind of thing -- a claim that
    these cells carry that identity -- so the evidence requirement and the warning are
    applied to all of them identically.
    """
    claims: list[tuple[str, str]] = []
    selected = str(value.get("selected") or "")
    if selected and selected != UNKNOWN:
        claims.append(("selected", selected))
    subtype = str(value.get("subtype") or "")
    if subtype:
        claims.append(("subtype", subtype))
    for name in value.get("co_occurring_identities") or []:
        claims.append(("cooc", str(name)))
    return claims


def valid_final(value: Any, pool: dict, cluster: str) -> bool:
    """Whether one answer is answerable from the evidence this cluster was given.

    Three kinds of check, and none of them is a threshold. The answer has to be
    structurally complete; every name in it has to be a candidate that was supplied; and
    every claim has to quote a gene from its own identity's curated panel with the
    measurement the packet showed. Whether that measurement is good enough for the claim
    is the agent's judgement and stays the agent's judgement.
    """
    if not isinstance(value, dict) or value.get("action") != "final":
        return False
    if value.get("schema_version") != ANNOTATOR_SCHEMA:
        return False
    if value.get("confidence") not in CONFIDENCE_VALUES:
        return False
    if not isinstance(value.get("reason", ""), str):
        return False
    if not isinstance(value.get("state", ""), str):
        return False

    names = pool_api.candidate_names(pool, cluster)
    groups = _groups_cover(value, names)
    if groups is None:
        return False

    selected = value.get("selected")
    if not isinstance(selected, str) or not selected.strip():
        return False
    if selected != UNKNOWN and selected not in names:
        return False

    subtype = value.get("subtype", "")
    if not isinstance(subtype, str):
        return False
    if subtype:
        if selected == UNKNOWN or subtype == selected or subtype not in names:
            return False
        home = next((group for group in groups if selected in group), [])
        if subtype not in home:
            return False

    others = value.get("co_occurring_identities", [])
    if not isinstance(others, list):
        return False
    if any(
        not isinstance(name, str) or name not in names or name == selected
        for name in others
    ):
        return False
    # A co-occurring identity has to be a DIFFERENT identity. Another name from the
    # selected candidate's own group is the same cell, and reporting it as also present
    # would turn a synonym into a second population.
    home = next((group for group in groups if selected in group), [])
    if any(name in home for name in others):
        return False
    if selected == UNKNOWN and others:
        return False

    # An identity retrieval never offered is reportable but not answerable. Reportable
    # because the cut is real -- 15 of a median 67 admitted cell types reach the agent, and
    # until now a menu that missed the answer looked exactly like a menu that did not.
    # Not answerable because a name with no supplied panel has no `claim_evidence` to
    # quote, so asserting it would put an unverifiable label in the same column as verified
    # ones. Naming a candidate here instead is a contradiction: that one could just be
    # answered with.
    unlisted = value.get("unlisted_identity", "")
    if not isinstance(unlisted, str) or unlisted in names:
        return False

    listed = {
        str(marker["gene"]).upper()
        for marker in (pool_api.find_candidate(pool, cluster, selected) or {}).get(
            "markers", []
        )
    }
    support = value.get("support_markers", [])
    if not isinstance(support, list) or any(
        not isinstance(gene, str) or gene.upper() not in listed for gene in support
    ):
        return False

    return _quotes_match(value, pool, cluster)


def _quotes_match(value: dict, pool: dict, cluster: str) -> bool:
    """Every claim names a gene from its own panel and copies its measurement.

    The gene must belong to the claimed identity's curated panel, not merely to the
    dataset: a claim supported by a gene the resource never tied to that cell type is
    not a claim the curated evidence can carry.
    """
    entries = value.get("claim_evidence")
    if not isinstance(entries, list):
        return False
    quoted: dict[str, dict] = {}
    for entry in entries:
        if not isinstance(entry, dict):
            return False
        identity = str(entry.get("identity") or "")
        gene = str(entry.get("decisive_gene") or "").upper()
        candidate = pool_api.find_candidate(pool, cluster, identity)
        if candidate is None or not gene:
            return False
        marker = next(
            (
                row
                for row in candidate["markers"]
                if str(row["gene"]).upper() == gene
            ),
            None,
        )
        if marker is None:
            return False
        for field, measured in (
            ("pct_in", marker["pct_in"]),
            ("pct_out", marker["pct_out"]),
        ):
            try:
                stated = float(entry.get(field))
            except (TypeError, ValueError):
                return False
            if measured is None or abs(stated - float(measured)) > QUOTE_TOLERANCE:
                return False
        quoted[identity] = entry
    return all(name in quoted for _role, name in claimed_identities(value))


# ------------------------------------------------------------------ the loop ----
def _turn(prompt: str, api_key: str, api_url: str, trace_id: str, turn: int):
    return llm.cached_call_llm(
        prompt,
        api_url,
        api_key,
        reasoning_effort=ANNOTATOR_EFFORT,
        model=ANNOTATOR_MODEL,
        trace_id=trace_id,
        turn_index=turn,
    )


# ------------------------------------------------------------- quality control ----
SUBTYPE_PASS = "separable"
ANNOTATION_PASS = "supported"
# The verdicts a second answer could plausibly change, and only those. A label pitched
# finer than the evidence reaches can be answered from candidates the agent already
# holds, and a contradiction names the gene to answer against. `insufficient_evidence` is
# not there on purpose: the panel that was too thin to decide is the same panel the next
# answer would read, so asking again buys the same reply at the price of a turn. Nor is a
# refinement that was merely demoted -- the parent it fell back to is a real answer, and
# re-opening a call the evidence already carries is how a check turns into a haggle.
#
# There is no `too_broad` here because this judge cannot reach it: saying a label is too
# broad means naming a narrower identity, and this judge is handed one label and told to
# use nothing else. Across 241 judgments it never once returned it. An inverted refinement
# -- a "finer" name wider than its parent -- is caught upstream by the subtype judge's
# direction rule, which is given both labels and can compare them.
RETRYABLE_VERDICTS = ("too_specific", "contradicted")
# How good an answer was, for choosing between the versions the rounds produced. A label
# only pitched too fine still identified the lineage; a contradicted label did not, and a
# panel too thin to decide settled nothing either way.
_VERDICT_RANK = {
    ANNOTATION_PASS: 3,
    "too_specific": 2,
    "insufficient_evidence": 1,
    "contradicted": 0,
}


def _judge(
    name: str,
    packet: dict,
    pool: dict,
    cluster: str,
    sources,
    api_key: str,
    api_url: str,
    trace_id: str,
    turn: int,
) -> dict | None:
    """One quality-control judgment, with the turns it asked for, or None if none answered.

    Its own model and its own prompt, so the cache key differs from the annotator's and
    a judgment is never served the annotator's cached reply.

    A judge used to get one call and whatever evidence code had decided to send. That made
    code the last word on what evidence exists, and code was wrong about it: a fourteen-row
    cut hid 10% of the genes the answers rested on. It now reads the whole panel and may
    ask for what the packet did not carry -- the same growing-prompt loop the annotator
    uses, so the head stays byte-identical and a provider can serve it from cache. The
    budget is small on purpose: a check that could ask forever would be negotiating.
    """
    if not packet:
        return None
    labels = tuple(
        str(packet[key])
        for key in ("label_under_test", "parent_label")
        if str(packet.get(key) or "")
    )
    prompt = _render(_prompt(name), packet)
    asked = 0
    unusable = 0
    seen: dict[str, dict] = {}
    while True:
        text, _error, _replayed = llm.cached_call_llm(
            prompt,
            api_url,
            api_key,
            reasoning_effort=JUDGE_EFFORT,
            model=JUDGE_MODEL,
            trace_id=trace_id,
            turn_index=turn,
        )
        parsed = llm.parse_json(text) if text else None
        if isinstance(parsed, dict) and parsed.get("verdict"):
            return parsed
        answerable = (
            isinstance(parsed, dict)
            and valid_query(parsed)
            and asked < JUDGE_QUERY_TURNS
        )
        if not answerable:
            if unusable >= SCHEMA_RETRIES:
                return None
            # Neither a verdict nor a question that can still be answered, so it is asked
            # again rather than filed as unchecked. `not_checked` has to keep meaning "no
            # judge was reachable": reporting a label as never checked because one reply
            # was badly formatted tells a reader the check did not run, which is a
            # different fact about the run.
            unusable += 1
            prompt += f"\n\n# Retry {unusable}: return the verdict as one JSON object.\n"
            continue
        call = json.dumps(
            {"tool": parsed["tool"], "args": parsed.get("args", {})},
            sort_keys=True,
            separators=(",", ":"),
        )
        # A repeated call returns the same answer, so it is answered once and marked. The
        # budget is what stops the loop; this only stops it buying the same row twice.
        observation = None if parsed["tool"] == "sources" else seen.get(call)
        if observation is None:
            observation = pool_api.run_judge_tool(
                pool, cluster, labels, parsed["tool"], parsed.get("args", {}), sources
            )
            seen[call] = observation
        else:
            observation = {**observation, "duplicate_query": True}
        asked += 1
        prompt += _observation_block(asked, observation)
        if asked >= JUDGE_QUERY_TURNS:
            prompt += "\n# No further queries. Return the verdict now.\n"


def _review(
    pool: dict,
    cluster: str,
    final: dict,
    sources,
    api_key: str,
    api_url: str,
    trace_id: str,
    turn: int,
) -> dict[str, Any]:
    """Both quality-control questions asked of one answer.

    Two judges rather than one, because they are different questions asked of different
    evidence: whether the finer name is separable from ITS PARENT, and whether the
    delivered label is carried by its own panel. A single verdict would have to mean both
    at once, and a refinement dropped for being unseparable would be indistinguishable
    from a lineage call that was simply wrong.

    A refinement is demoted whenever it is not affirmatively separable, so an unchecked
    refinement costs granularity rather than passing by default. Where no judge answered
    at all, nothing is demoted and nothing is retried: silently coarsening every label of
    an offline replay would be a worse failure than leaving the answer unchecked.
    """
    selected = str(final.get("selected") or "")
    subtype = str(final.get("subtype") or "")
    review: dict[str, Any] = {
        "subtype_verdict": None,
        "annotation_verdict": None,
        "demoted_subtype": "",
        "judged_label": "",
        "checked": False,
        "passed": False,
        "retry": False,
        "rank": 0,
    }
    if subtype and selected and selected != UNKNOWN:
        verdict = _judge(
            "subtype_judge",
            pool_api.judge_packet(pool, cluster, subtype, sources, parent=selected),
            pool,
            cluster,
            sources,
            api_key,
            api_url,
            trace_id,
            turn,
        )
        review["subtype_verdict"] = verdict
        if verdict is not None:
            review["checked"] = True
            if str(verdict.get("verdict")) != SUBTYPE_PASS:
                review["demoted_subtype"] = subtype
    label = subtype if (subtype and not review["demoted_subtype"]) else selected
    review["judged_label"] = label
    if not label or label == UNKNOWN:
        # `Unknown` claims no identity, so there is nothing for a judge to carry or to
        # refute, and the answer stands as given. That counts as checked: `not_checked` is
        # reserved for a judge that could not be reached, and an abstention reported the
        # same way would be indistinguishable from a run whose model was down.
        review["checked"] = True
        review["passed"] = True
        review["rank"] = _VERDICT_RANK[ANNOTATION_PASS]
        return review
    verdict = _judge(
        "annotation_judge",
        pool_api.judge_packet(pool, cluster, label, sources),
        pool,
        cluster,
        sources,
        api_key,
        api_url,
        trace_id,
        turn,
    )
    review["annotation_verdict"] = verdict
    if verdict is not None:
        answer = str(verdict.get("verdict"))
        review["checked"] = True
        review["rank"] = _VERDICT_RANK.get(answer, 0)
        review["passed"] = answer == ANNOTATION_PASS
        review["retry"] = answer in RETRYABLE_VERDICTS
    return review


def _qc_feedback(review: dict) -> dict:
    """What the annotator is told, and only that.

    The packet is already at the top of this conversation, so none of it is sent again:
    what comes back is the verdict, the genes it turned on, and the direction the level
    was wrong in. That is the whole point of asking -- an error found and not carried to
    the next step is an error nobody acted on -- but it is also the whole of it. Which
    candidate to name is the agent's decision, made against the candidates it already
    holds, and re-deciding it here would put two annotators on one cluster.
    """
    block: dict[str, Any] = {}
    verdict = review.get("subtype_verdict")
    if verdict:
        block["subtype"] = {
            "name": review.get("demoted_subtype") or "",
            "verdict": verdict.get("verdict"),
            "conflicting_markers": verdict.get("conflicting_markers") or [],
            "reason": verdict.get("reason"),
        }
    verdict = review.get("annotation_verdict")
    if verdict:
        block["annotation"] = {
            "label": review.get("judged_label"),
            "verdict": verdict.get("verdict"),
            "conflicting_markers": verdict.get("conflicting_markers") or [],
            "reason": verdict.get("reason"),
        }
    block["note"] = (
        "quality control on the answer above, judged from that identity's own panel and "
        "its curated sentences alone. A refinement that is not separable has already "
        "been dropped to the parent. Answer again against the candidates you already "
        "hold: keep the call where the evidence still carries it, or name the candidate "
        "the evidence does carry. Query only where the verdict turns on something not "
        "already on this page."
    )
    return block


def annotate_cluster(
    pool: dict,
    cluster: str,
    template: str,
    api_key: str,
    api_url: str,
    sources=None,
) -> dict[str, Any]:
    """Run one cluster's conversation to a validated answer.

    Returns the parsed answer with a transcript of what was asked, or an error. The
    conversation is a single growing prompt string rather than a message list, so every
    turn is one ordinary call that the prompt-hash cache can replay exactly.
    """
    trace_id = uuid.uuid4().hex[:16]
    prompt = _render(template, cluster_packet(pool, cluster))
    transcript: list[dict[str, Any]] = []
    seen: dict[str, dict] = {}
    turn = 0
    schema_failures = 0
    repeated = 0
    forced = False
    audited: list[dict] | None = None
    attempts: list[dict[str, Any]] = []

    while True:
        turn += 1
        text, error, _replayed = _turn(prompt, api_key, api_url, trace_id, turn)
        if text is None:
            if error == "cache_miss_no_credentials":
                return {"error": error, "turns": turn, "transcript": transcript}
            schema_failures += 1
            if schema_failures > SCHEMA_RETRIES:
                return {
                    "error": error or "no response",
                    "turns": turn,
                    "transcript": transcript,
                }
            prompt += f"\n\n# Retry {schema_failures}\n"
            continue

        parsed = llm.parse_json(text)
        if valid_final(parsed, pool, cluster):
            # The audit the run already performs, returned in time to be answered. It
            # asserts nothing about the identity: the agent reads the exclusion against
            # the sources it already holds and answers again, keeping or changing the
            # call. Asked once, so a cluster costs at most one extra turn.
            pending = (
                []
                if audited is not None or forced
                else pool_api.binding_exclusions(
                    pool,
                    cluster,
                    str(parsed.get("selected") or ""),
                    parsed.get("identity_groups"),
                    parsed.get("co_occurring_identities"),
                )
            )
            if pending:
                audited = pending
                prompt += _observation_block(
                    len(transcript) + 1,
                    {
                        "exclusion_audit": pending,
                        "note": "curated exclusions the cluster detects on the identity "
                                "you selected and on each you reported as co-occurring; "
                                "an empty list is a measured absence, not a gap",
                    },
                )
                continue
            audited = audited or []
            # Quality control, asked of the answer rather than of the conversation. A
            # failing verdict goes back in as one observation and the agent answers
            # again; the rounds are bounded because a judge that could ask forever would
            # be negotiating, not checking.
            review = _review(
                pool, cluster, parsed, sources, api_key, api_url, trace_id, turn
            )
            attempts.append({"final": parsed, "review": review})
            if review["retry"] and not forced and len(attempts) < JUDGE_MAX_ROUNDS:
                prompt += _observation_block(
                    len(transcript) + 1, {"quality_control": _qc_feedback(review)}
                )
                continue
            # Every round produced a complete answer, so the one delivered is the one
            # the judges rated highest rather than whichever came last.
            best = max(
                range(len(attempts)),
                key=lambda index: (attempts[index]["review"]["rank"], index),
            )
            chosen = attempts[best]
            return {
                "final": chosen["final"],
                "turns": turn,
                "transcript": transcript,
                "trace_id": trace_id,
                "turn_budget_exhausted": forced,
                "exclusion_audit": audited or [],
                "quality_control": {
                    "rounds": len(attempts),
                    "selected_round": best + 1,
                    **chosen["review"],
                },
            }
        if valid_query(parsed) and not forced:
            schema_failures = 0
            call = json.dumps(
                {"tool": parsed["tool"], "args": parsed.get("args", {})},
                sort_keys=True,
                separators=(",", ":"),
            )
            # `sources` is the one tool whose answer changes when it is asked again:
            # it serves the next sentences rather than the same ones. Short-circuiting it
            # as a repeat is what used to make "these three did not settle it" an
            # unanswerable request. Its own per-pair limit is the brake instead.
            streaming = parsed["tool"] == "sources"
            if call in seen and not streaming:
                repeated += 1
                observation = {**seen[call], "duplicate_query": True}
            else:
                repeated = 0
                observation = run_tool(
                    pool, cluster, parsed["tool"], parsed.get("args", {}), sources
                )
                seen[call] = observation
            transcript.append(
                {"turn": turn, "call": json.loads(call), "duplicate": repeated > 0}
            )
            prompt += _observation_block(len(transcript), observation)
            if repeated >= MAX_REPEATED_QUERIES or (
                MAX_TURNS and turn >= MAX_TURNS
            ):
                forced = True
                prompt += (
                    "\n# No further queries. Return the final answer now.\n"
                )
            continue

        schema_failures += 1
        if schema_failures > SCHEMA_RETRIES:
            return {
                "error": "schema violation",
                "turns": turn,
                "transcript": transcript,
                "trace_id": trace_id,
            }
        prompt += f"\n\n# Retry {schema_failures}: return one valid JSON object.\n"


# --------------------------------------------------------------- record building ----
def _measurements(entry: dict, percentile: dict, cluster: str, limit: int = 30):
    """The measured markers reported for one claimed identity.

    The raised positive ones, best-corroborated first, and every curated exclusion this
    cluster detects. On the fraction scale the published tables use: the scale changes
    exactly here, at the boundary between what the agent read and what a reader reads, and
    nowhere else.

    The exclusions are here because the delivered reason argues from them. Of 327 detected
    exclusions carried by a delivered identity, 262 were named in that cluster's own
    rationale -- "the detected CD24 exclusion (67.4% versus 17.0%)" -- while the evidence
    table held no CD24 row at all, so a reader could not check the number or reach its
    source. A detected exclusion does not by itself refute the call, which is why the
    audit turn asks the agent to say what it means rather than deciding for it; but a
    reason that argues from a row belongs beside that row.

    They keep their own budget. A detected exclusion is at most a handful per identity
    against thirty positives, and letting the two compete for one cap is how the strongest
    evidence gets dropped by whatever sorts last.
    """
    positives, negatives = [], []
    for marker in entry["markers"]:
        polarity = marker["polarity"]
        pct_in = marker["pct_in"] or 0.0
        if polarity == "positive":
            if not is_raised(marker):
                continue
            target = positives
        elif polarity == "negative":
            if pct_in < pool_api.NEGATIVE_SOURCE_MIN_PCT_IN:
                continue
            target = negatives
        else:
            continue
        gene = str(marker["gene"])
        target.append(
            {
                "gene": gene,
                "polarity": polarity,
                "detection_fraction_in": _round(pct_in / 100.0),
                "detection_fraction_out": _round((marker["pct_out"] or 0.0) / 100.0),
                "avg_log2FC": marker["avg_log2FC"],
                "cross_cluster_percentile": percentile.get((cluster, gene.upper())),
                "publication_support": marker["n_pub"],
                "evidence_tier": marker["tier"],
            }
        )
    return positives[:limit] + negatives


def _candidate_entry(
    entry: dict,
    cluster: str,
    percentile: dict,
    role: str,
    evidence: dict | None,
    warning: str,
) -> dict:
    """One candidate's full audit record: what it curates, what was measured, what was claimed."""
    return {
        "cell_type": entry["cell_type"],
        "retrieval_rank": entry["retrieval_rank"],
        "claim_role": role,
        "claim_evidence": evidence or {},
        "claim_warning": warning,
        "panel": entry["markers"],
        "unmeasured_curated_genes": entry["unmeasured_curated_genes"],
        "single_cell_program": {
            "in_cluster_median": entry["program"]["median_in"],
            "out_of_cluster_median": entry["program"]["median_out"],
        },
        "decisive_marker_measurements": _measurements(entry, percentile, cluster),
    }


def _fallback(cluster: str, frame: pd.DataFrame, reason: str, detail: str) -> dict:
    """Top of the retrieval order, labelled as exactly that.

    A fallback label is not a model judgement, so there is nothing for quality control to
    have carried or refuted; the column says `not_checked` rather than leaving a reader
    to infer from a blank whether a check passed or never ran.
    """
    if frame.empty:
        return {
            "annotation_qc": QC_UNCHECKED,
            "quality_control": {},
            "cluster_id": cluster,
            "annotation": UNKNOWN,
            "subtype": "",
            "state": "",
            "co_occurring_identities": [],
            "confidence": NOT_AVAILABLE,
            "rationale": "no candidate's curated positive markers are significantly "
            "up-regulated in this cluster",
            "resolution_status": UNRESOLVED,
            "resolution_detail": "empty_candidate_pool",
            "annotation_source": "no_candidate",
            "llm_status": reason,
            "support_markers": [],
            "claim_warnings": [],
            "unlisted_identity": "",
        }
    top = frame.iloc[0]
    return {
        "annotation_qc": QC_UNCHECKED,
        "quality_control": {},
        "cluster_id": cluster,
        "annotation": str(top["candidate"]),
        "subtype": "",
        "state": "",
        "co_occurring_identities": [],
        "confidence": NOT_AVAILABLE,
        "rationale": (
            "no model judgement was available; this candidate leads the joint retrieval "
            "order over marker-level, cluster-level and single-cell-level evidence with "
            f"{int(top['hits'])} of {int(top['panel_size'])} measured positive markers "
            "significantly up-regulated here. The retrieval order is not a confidence."
        ),
        "resolution_status": RESOLVED,
        "resolution_detail": detail,
        "annotation_source": "relative_score_fallback",
        "llm_status": reason,
        "support_markers": [],
        "claim_warnings": [],
        "unlisted_identity": "",
    }


def _subtype_demotion(
    pool: dict, cluster: str, selected: str, subtype: str
) -> tuple[str, str]:
    """The subtype to drop and the line saying why, or two empty strings.

    Shared lineage markers cannot establish a finer identity, so the refinement is kept
    only when the cluster raises enough markers that the parent's curated definition does
    not contain. The line records the count against the requirement and names whichever
    exclusive markers there were, so a demotion can be checked rather than trusted.
    """
    if not subtype or not selected or selected == UNKNOWN:
        return "", ""
    exclusive = pool_api.subtype_exclusive_markers(pool, cluster, selected, subtype)
    required = pool_api.SUBTYPE_EXCLUSIVE_REQUIRED
    if len(exclusive) >= required:
        return "", ""
    found = ",".join(exclusive) if exclusive else "none"
    return subtype, (
        f"rejected_subtype {subtype}: exclusive_defining_markers "
        f"{len(exclusive)}/{required} | raised panel is shared with {selected}"
        f" | exclusive: {found}"
    )


QC_PASSED = "passed"
QC_REVISED = "passed_after_revision"
QC_DEMOTED = "demoted_to_parent"
QC_FAILED = "failed"
QC_UNCHECKED = "not_checked"
QC_VALUES = (QC_PASSED, QC_REVISED, QC_DEMOTED, QC_FAILED, QC_UNCHECKED)


def _qc_status(review: dict, demoted: bool) -> str:
    """How the delivered label left quality control, in one word.

    `not_checked` is its own answer rather than a pass: a replay with no reachable judge
    would otherwise report every cluster as having survived a check nobody ran.
    """
    if not review or not review.get("checked"):
        return QC_UNCHECKED
    label = str(review.get("judged_label") or "")
    verdict = str((review.get("annotation_verdict") or {}).get("verdict") or "")
    carried = verdict == ANNOTATION_PASS or not label or label == UNKNOWN
    if not carried:
        return QC_FAILED
    if demoted:
        return QC_DEMOTED
    return QC_REVISED if int(review.get("selected_round") or 1) > 1 else QC_PASSED


def _result(
    cluster: str,
    pool: dict,
    percentile: dict,
    outcome: dict,
    frame: pd.DataFrame,
) -> dict[str, Any]:
    """One cluster's record, with the single post-hoc audit applied to the answer."""
    names = pool_api.candidate_names(pool, cluster)
    final = outcome["final"]
    selected = str(final["selected"])
    subtype = str(final.get("subtype", ""))
    others = [str(name) for name in final.get("co_occurring_identities", [])]

    # The finer level has to earn its own name. A subtype whose raised markers all belong
    # to the parent's curated definition is the parent program under a narrower label, so
    # the reported identity falls back to the parent rather than claiming a refinement the
    # measurements do not separate. The rejected name keeps its evidence and its role: why
    # the cluster is NOT the finer type is as auditable as why it is the parent.
    rejected_subtype, demotion = _subtype_demotion(pool, cluster, selected, subtype)
    # The free gate and the judge answer the same question from different evidence -- one
    # from whether the parent's panel already contains the genes, the other from what the
    # curated sentences make of them -- and a refinement has to survive both. Either
    # verdict alone is enough to drop it, and whichever spoke first keeps its own line so
    # the record says which check the refinement failed.
    review = dict(outcome.get("quality_control") or {})
    judged_out = str(review.get("demoted_subtype") or "")
    if judged_out and not rejected_subtype:
        verdict = review.get("subtype_verdict") or {}
        rejected_subtype = judged_out
        demotion = (
            f"rejected_subtype {judged_out}: quality_control "
            f"{verdict.get('verdict') or 'unavailable'} | {verdict.get('reason') or ''}"
        ).strip()
    if rejected_subtype:
        subtype = ""
    claims = claimed_identities({**final, "subtype": subtype})
    if rejected_subtype:
        claims.append(("rejected_subtype", rejected_subtype))

    # The one audit. It verifies nothing about biology: the quoted numbers were already
    # checked against the DE table by `valid_final`, and this adds the short line naming
    # the claimed identity's curated markers the cluster does not raise.
    warned = pool_api.claim_warnings(pool, cluster, claims)
    if demotion:
        warned.append((rejected_subtype, demotion))
    warning_of: dict[str, str] = {}
    for name, line in warned:
        warning_of[name] = (
            f"{warning_of[name]} || {line}" if name in warning_of else line
        )
    warnings = [line for _, line in warned]
    role_of = {name: role for role, name in claims}
    evidence_of = {
        str(item.get("identity")): item
        for item in final.get("claim_evidence", [])
        if isinstance(item, dict)
    }

    entries = [
        _candidate_entry(
            entry,
            cluster,
            percentile,
            role_of.get(str(entry["cell_type"]), ""),
            evidence_of.get(str(entry["cell_type"])),
            warning_of.get(str(entry["cell_type"]), ""),
        )
        for entry in pool["clusters"][cluster]["candidates"]
    ]

    if selected == UNKNOWN:
        status, detail = UNRESOLVED, "agent_established_no_identity"
    elif others:
        status, detail = MIXED, "agent_majority_of_several_identities"
    else:
        status, detail = RESOLVED, "agent_selected"
    return {
        "cluster_id": cluster,
        "annotation": selected,
        "subtype": subtype,
        "state": str(final.get("state", "")),
        "co_occurring_identities": others,
        "confidence": str(final["confidence"]),
        "rationale": _clip(final.get("reason", ""), MAX_RATIONALE),
        "resolution_status": status,
        "resolution_detail": detail,
        "annotation_source": "cluster_annotation",
        "llm_status": "annotated",
        "support_markers": [
            str(gene) for gene in final.get("support_markers", []) if str(gene)
        ],
        "claim_warnings": warnings,
        # What the agent says the retrieved menu was missing, when it says so. A column
        # rather than a sentence buried in the rationale, because "how often did retrieval
        # not offer the answer" is a rate somebody has to be able to count.
        "unlisted_identity": str(final.get("unlisted_identity", "")),
        "candidates": names,
        "candidate_entries": entries,
        "identity_groups": final.get("identity_groups", []),
        "turns": int(outcome.get("turns", 0)),
        "tool_calls": [item["call"] for item in outcome.get("transcript", [])],
        "turn_budget_exhausted": bool(outcome.get("turn_budget_exhausted")),
        "exclusion_audit": list(outcome.get("exclusion_audit") or []),
        "annotation_qc": _qc_status(review, bool(rejected_subtype)),
        "quality_control": review,
    }


def _percentiles(de: pd.DataFrame) -> dict[tuple[str, str], float]:
    """Where each cluster ranks among all clusters for one gene's fold change.

    Reported in the marker table only. It never reaches the agent: it saturates at 1.0
    for most of the markers shown, so as evidence it separates almost nothing while
    costing a column of prompt.
    """
    frame = de.copy()
    frame["gene_key"] = frame["feature"].astype(str).str.upper()
    ranked = frame.pivot_table(
        index="gene_key", columns="group", values="avg_log2FC", aggfunc="first"
    ).rank(axis=1, pct=True, method="average")
    return {
        (str(cluster), str(gene)): _round(value)
        for gene, row in ranked.iterrows()
        for cluster, value in row.items()
        if value == value
    }


# --------------------------------------------------------------------- stage ----
def annotate(tag: str, enabled: bool | None = None) -> dict[str, dict[str, Any]]:
    enabled = CLUSTER_ANNOTATION_ENABLED if enabled is None else bool(enabled)
    with open(os.path.join(CACHE, f"{tag}_candidate_scoring.pkl"), "rb") as handle:
        scoring = pickle.load(handle)
    with open(os.path.join(CACHE, f"{tag}_de_meta.pkl"), "rb") as handle:
        prep = pickle.load(handle)

    de = prep["de"].copy()
    de["group"] = de["group"].astype(str)
    percentile = _percentiles(de)

    context = scoring["context"]
    sources = SourceDB().context(
        context["species"], context["tissue"], context["disease"]
    )
    pool = build_pool(scoring, prep, sources)
    scored = scoring["scored"]
    clusters = sorted(scoring["clusters"], key=lambda value: (len(value), value))

    api_key, api_url = llm.resolve_api()
    if not enabled:
        llm_status = "disabled"
    elif not api_key or not api_url:
        llm_status = "skipped_no_credentials"
    else:
        llm_status = "enabled"
    print(
        f"==== annotate(py) {tag}: {len(clusters)} clusters, "
        f"top {scoring['top_candidates']} candidates each, llm={llm_status} ====\n"
        f"  annotator: {llm.resolve_model(ANNOTATOR_MODEL)} ({ANNOTATOR_EFFORT})"
    )

    results: dict[str, dict[str, Any]] = {}
    runnable: list[str] = []
    for cluster in clusters:
        frame = scored[scored["cluster"] == cluster].reset_index(drop=True)
        if scoring["clusters"][cluster]["status"] == UNSUPPORTED:
            results[cluster] = _fallback(
                cluster, frame.head(0), llm_status, "empty_candidate_pool"
            )
            continue
        if llm_status != "enabled":
            results[cluster] = _fallback(
                cluster, frame, llm_status, "relative_score_fallback_top1"
            )
            results[cluster]["candidates"] = pool_api.candidate_names(pool, cluster)
            results[cluster]["candidate_entries"] = [
                _candidate_entry(entry, cluster, percentile, "", None, "")
                for entry in pool["clusters"][cluster]["candidates"]
            ]
            continue
        runnable.append(cluster)

    template = _prompt("cluster_annotator")

    def run_one(cluster: str):
        # One server per cluster: "already shown" is a fact about this conversation, and
        # sharing it across clusters would let one cluster's reading exhaust another's.
        server = SourceServer(
            sources, batch=SOURCES_PER_MARKER, max_batches=SOURCE_BATCHES_PER_MARKER
        )
        pool_api.register_packet_sources(server, pool, cluster)
        return cluster, annotate_cluster(
            pool, cluster, template, api_key, api_url, server
        )

    outcomes: dict[str, dict] = {}
    if runnable:
        # Each cluster is one sequential conversation; the clusters run beside each other.
        with ThreadPoolExecutor(
            max_workers=min(LLM_SETTINGS.threads, len(runnable))
        ) as executor:
            for cluster, outcome in executor.map(run_one, runnable):
                outcomes[cluster] = outcome

    for cluster in runnable:
        outcome = outcomes[cluster]
        frame = scored[scored["cluster"] == cluster].reset_index(drop=True)
        if "final" not in outcome:
            error = outcome.get("error") or "no answer"
            record = _fallback(
                cluster, frame, f"failed:{error}", "annotation_failed"
            )
            record["candidates"] = pool_api.candidate_names(pool, cluster)
            record["candidate_entries"] = [
                _candidate_entry(entry, cluster, percentile, "", None, "")
                for entry in pool["clusters"][cluster]["candidates"]
            ]
            record["turns"] = int(outcome.get("turns", 0))
            results[cluster] = record
            continue
        results[cluster] = _result(cluster, pool, percentile, outcome, frame)

    payload = {
        "tag": tag,
        "context": context,
        "llm_status": llm_status,
        "models": {
            "annotator": llm.resolve_model(ANNOTATOR_MODEL),
            "annotator_reasoning_effort": ANNOTATOR_EFFORT,
        },
        "results": results,
        "transcripts": {
            cluster: outcome.get("transcript", []) for cluster, outcome in outcomes.items()
        },
    }
    with open(os.path.join(CACHE, f"{tag}_annotations.pkl"), "wb") as handle:
        pickle.dump(payload, handle, protocol=4)
    llm.flush_response_cache()

    counts = {status: 0 for status in (RESOLVED, MIXED, UNRESOLVED)}
    for value in results.values():
        counts[value["resolution_status"]] += 1
    turns = [value.get("turns", 0) for value in results.values() if value.get("turns")]
    warned = sum(1 for value in results.values() if value.get("claim_warnings"))
    print(
        f"[done] {tag}: {counts[RESOLVED]} resolved, {counts[MIXED]} mixed, "
        f"{counts[UNRESOLVED]} unresolved of {len(results)} clusters; "
        f"turns median {int(np.median(turns)) if turns else 0}, max {max(turns, default=0)}; "
        f"{warned} clusters carry a claim warning; "
        f"{sum(1 for v in results.values() if v['annotation_source'] != 'cluster_annotation')} "
        "not model-decided"
    )
    return results


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit("usage: python -m scmarkeragent.cluster_annotation <tag>")
    annotate(sys.argv[1])
