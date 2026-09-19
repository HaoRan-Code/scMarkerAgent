#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
from . import borrowed_context
from . import llm_client as llm
from . import stages
from .annotator_pool import build_pool, cluster_packet, is_raised, run_tool
from .config import (
    CLUSTER_ANNOTATION_ENABLED,
    CACHE_DIR as CACHE,
    CROSS_SPECIES,
    CROSS_SPECIES_BORROW,
    DEFAULTS,
    LLM_SETTINGS,
    NOT_AVAILABLE,
    PROMPT_DIR,
)
from .candidate_scoring import UNSUPPORTED
from .marker_sources import CompositeSources, SourceDB, SourceServer

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

REVIEW_MODEL = str(_A["review_model"])
REVIEW_EFFORT = str(_A["review_reasoning_effort"])
REVIEW_MAX_ROUNDS_PER_TIER = int(_A["review_max_rounds_per_tier"])
REVIEW_QUERY_TURNS = int(_A["review_query_turns"])
ARBITER_MODEL = str(_A["arbiter_model"])
ARBITER_EFFORT = str(_A["arbiter_reasoning_effort"])

BORROW_SCREEN_MODEL = str(_A["borrow_screen_model"])
BORROW_SCREEN_EFFORT = str(_A["borrow_screen_reasoning_effort"])
BORROW_SCREEN_PROMPT = "borrow_tissue_screen"

REVIEW_CITATIONS_REQUIRED = 3

CONFIDENCE_VALUES = ("high", "medium", "low")

QUOTE_TOLERANCE = 0.15

MAX_REPEATED_QUERIES = 3

RESOLVED = "resolved"
MIXED = "mixed"
UNRESOLVED = "unresolved"

IDENTITY_PASS = "accept"
IDENTITY_VALUES = ("accept", "reject", "insufficient_evidence")
SUBTYPE_PASS = "separable"
SUBTYPE_VALUES = ("separable", "not_separable", "")

PANEL_READING = "panel_reading"
PANEL_READING_PLACEHOLDER = "{{PANEL_READING}}"
CROSS_SPECIES_NOTE = "cross_species_note"


def _prompt(name: str) -> str:
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


def _context_header(pool: dict, cluster: str) -> str:
    ctx = pool.get("context") or {}
    state = pool["clusters"][str(cluster)]
    disease = ctx.get("disease")
    disease = ", ".join(str(d) for d in disease) if isinstance(disease, list) else str(disease)
    tissue = ctx.get("tissue")
    tissue = ", ".join(str(t) for t in tissue) if isinstance(tissue, (list, tuple)) else str(tissue)
    return (
        "# Dataset context\n"
        f"species: {ctx.get('species')} | tissue: {tissue} | disease: {disease} | "
        f"development stage: {ctx.get('development_stage')} | cluster {cluster} of "
        f"{ctx.get('clusters_in_dataset')} clusters | {int(state['n_cells']):,} cells in this cluster\n\n"
    )


def _observation_block(index: int, observation: dict) -> str:
    payload = json.dumps(
        observation, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    )
    return f"\n\n# Observation {index}\n{payload}\n"


def valid_query(value: Any) -> bool:
    return (
        isinstance(value, dict)
        and value.get("action") == "query"
        and isinstance(value.get("tool"), str)
        and isinstance(value.get("args", {}), dict)
    )


def claimed_identities(value: dict) -> list[tuple[str, str]]:
    claims: list[tuple[str, str]] = []
    selected = str(value.get("selected") or "")
    if selected and selected != UNKNOWN:
        claims.append(("selected", selected))
    subtype = str(value.get("subtype") or "")
    if subtype and subtype != selected:
        claims.append(("subtype", subtype))
    for name in value.get("co_occurring_identities") or []:
        claims.append(("cooc", str(name)))
    return claims


def _listed_genes(pool: dict, cluster: str, selected: str, subtype: str) -> set[str]:
    return {
        str(marker["gene"]).upper()
        for name in {selected, subtype} if name
        for marker in (pool_api.find_candidate(pool, cluster, name) or {}).get("markers", [])
    }


def sanitize_support_markers(value: Any, pool: dict, cluster: str) -> tuple[Any, list[str]]:
    if not isinstance(value, dict) or value.get("action") != "final":
        return value, []
    selected = str(value.get("selected") or "")
    subtype = str(value.get("subtype") or "")
    support = value.get("support_markers")
    if not isinstance(support, list) or not support:
        return value, []
    listed = _listed_genes(pool, cluster, selected, subtype)
    kept = [g for g in support if isinstance(g, str) and g.upper() in listed]
    dropped = [g for g in support if not (isinstance(g, str) and g.upper() in listed)]
    if not dropped:
        return value, []
    return {**value, "support_markers": kept}, dropped


def valid_final(value: Any, pool: dict, cluster: str, tier: int = 1,
                program: dict | None = None) -> bool:
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

    names = pool_api.candidate_names(pool, cluster, tier=tier)

    selected = value.get("selected")
    if not isinstance(selected, str) or not selected.strip():
        return False
    if selected != UNKNOWN and selected not in names:
        return False

    subtype = value.get("subtype", "")
    if not isinstance(subtype, str):
        return False
    if subtype and (selected == UNKNOWN or subtype == selected or subtype not in names):
        return False

    lineage = value.get("lineage", "")
    if not isinstance(lineage, str):
        return False
    if lineage and lineage != selected:
        if selected == UNKNOWN or lineage not in names:
            return False

    others = value.get("co_occurring_identities", [])
    if not isinstance(others, list):
        return False
    if any(
        not isinstance(name, str) or name not in names or name == selected
        for name in others
    ):
        return False
    if selected == UNKNOWN and others:
        return False

    refinements = value.get("possible_refinements", [])
    if not isinstance(refinements, list):
        return False
    if any(not isinstance(name, str) or name not in names for name in refinements):
        return False

    listed = _listed_genes(pool, cluster, selected, subtype)
    support = value.get("support_markers", [])
    if not isinstance(support, list) or any(
        not isinstance(gene, str) or gene.upper() not in listed for gene in support
    ):
        return False

    if _program_fit_problems(value, program):
        return False

    return _quotes_match(value, pool, cluster)


def _program_fit_problems(value: dict, program: dict | None) -> list[str]:
    programs = [str(b.get("function") or "") for b in (program or {}).get("programs") or []]
    entries = value.get("program_fit", [])
    if not isinstance(entries, list):
        return ["program_fit must be a list, one entry per identity-bearing program of "
                "program_reading"]
    if not programs or str(value.get("selected") or "") == UNKNOWN:
        return []
    allowed = {name for _role, name in claimed_identities(value)} | {"unexplained"}
    problems, seen = [], set()
    for entry in entries:
        if not isinstance(entry, dict):
            problems.append("each program_fit entry must be an object")
            continue
        name = str(entry.get("program") or "")
        if name not in programs:
            problems.append(f"program_fit names a program that program_reading did not: {name!r}")
            continue
        if str(entry.get("accounted_for_by") or "") not in allowed:
            problems.append(f"{name!r}: accounted_for_by must be one of "
                            f"{sorted(allowed)} exactly as written")
        seen.add(name)
    missing = [p for p in programs if p not in seen]
    if missing:
        problems.append(f"every identity-bearing program needs one program_fit entry; "
                        f"missing {missing[:3]}")
    return problems


def _quotes_match(value: dict, pool: dict, cluster: str) -> bool:
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
            (row for row in candidate["markers"] if str(row["gene"]).upper() == gene),
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


def final_problems(value: Any, pool: dict, cluster: str, tier: int = 1,
                   program: dict | None = None) -> list[str]:
    problems: list[str] = []
    if not isinstance(value, dict):
        return problems
    names = pool_api.candidate_names(pool, cluster, tier=tier)
    selected = str(value.get("selected") or "")

    entries = value.get("claim_evidence")
    if not isinstance(entries, list):
        problems.append("claim_evidence must be a list with one entry per claimed identity")
        entries = []
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        identity = str(entry.get("identity") or "")
        gene = str(entry.get("decisive_gene") or "")
        candidate = pool_api.find_candidate(pool, cluster, identity)
        if candidate is None:
            problems.append(f"'{identity}' is not a supplied candidate")
            continue
        marker = next(
            (row for row in candidate["markers"] if str(row["gene"]).upper() == gene.upper()),
            None,
        )
        if marker is None:
            marker = next(
                (row for row in (candidate.get("definers") or [])
                 if str(row.get("gene") or "").upper() == gene.upper() and "pct_in" in row),
                None,
            )
        if marker is None:
            problems.append(f"{gene} is not on the curated panel of '{identity}'")
            continue
        for field in ("pct_in", "pct_out"):
            try:
                stated = float(entry.get(field))
            except (TypeError, ValueError):
                problems.append(f"{field} for {gene} must be the number shown")
                continue
            measured = marker.get(field)
            if measured is None or abs(stated - float(measured)) > QUOTE_TOLERANCE:
                problems.append(
                    f"{field} for {gene} under '{identity}' must be copied exactly "
                    f"({measured} shown)"
                )
    claimed = {name for _role, name in claimed_identities(value)}
    quoted = {str(e.get("identity") or "") for e in entries if isinstance(e, dict)}
    for name in sorted(claimed - quoted):
        problems.append(f"'{name}' is claimed but has no claim_evidence entry")
    problems.extend(_program_fit_problems(value, program))

    if not selected:
        problems.append("selected must be a supplied candidate or Unknown")
    elif selected != UNKNOWN and selected not in names:
        problems.append(f"selected '{selected}' is not a supplied candidate")
    if value.get("schema_version") != ANNOTATOR_SCHEMA:
        problems.append(f"schema_version must be '{ANNOTATOR_SCHEMA}'")
    if value.get("confidence") not in CONFIDENCE_VALUES:
        problems.append(f"confidence must be one of {list(CONFIDENCE_VALUES)}")
    for field in ("reason", "state", "subtype", "lineage"):
        if not isinstance(value.get(field, ""), str):
            problems.append(f"{field} must be a string")
    if str(value.get("subtype") or "") and str(value.get("subtype")) == selected:
        problems.append(
            "subtype must differ from selected (leave it empty when there is no finer name)"
        )
    for field in ("subtype", "lineage"):
        name = str(value.get(field) or "")
        if name and name != selected:
            if selected == UNKNOWN:
                problems.append(f"{field} cannot be given with selected Unknown")
            elif name not in names:
                problems.append(f"{field} '{name}' is not a supplied candidate")
    others = value.get("co_occurring_identities")
    if others is not None and not isinstance(others, list):
        problems.append("co_occurring_identities must be a list")
    for name in others if isinstance(others, list) else []:
        if not isinstance(name, str) or name not in names:
            problems.append(f"co-occurring '{name}' is not a supplied candidate")
        elif name == selected:
            problems.append(f"co-occurring '{name}' is the selected identity itself")
        elif selected == UNKNOWN:
            problems.append("co_occurring_identities must be empty with selected Unknown")
    refinements = value.get("possible_refinements")
    if refinements is not None and not isinstance(refinements, list):
        problems.append("possible_refinements must be a list")
    for name in refinements if isinstance(refinements, list) else []:
        if not isinstance(name, str) or name not in names:
            problems.append(
                f"possible_refinements '{name}' is not a supplied candidate; only "
                "candidate names can go there"
            )
    subtype_name = str(value.get("subtype") or "")
    listed_genes = _listed_genes(pool, cluster, selected, subtype_name)
    bad_support = [
        g for g in (value.get("support_markers") or []) if str(g).upper() not in listed_genes
    ]
    if bad_support:
        where = f"'{selected}'" + (
            f" or '{subtype_name}'" if subtype_name and subtype_name != selected else ""
        )
        problems.append(f"support_markers not on {where} panel: {bad_support[:6]}")
    return problems


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


def _packet_citations(packet: dict) -> set[tuple[str, str, str]]:
    out: set[tuple[str, str, str]] = set()
    for key in ("claimed_identities", "contested_identities"):
        for block in packet.get(key) or []:
            identity = str(block.get("cell_type") or "")
            for gene, records in (block.get("sources") or {}).items():
                for record in records or []:
                    for id_key in ("pmcid", "pmid"):
                        value = str(record.get(id_key) or "").strip()
                        if value and value != NOT_AVAILABLE:
                            out.add((identity, str(gene).upper(), value))
    return out


def _packet_sentence_count(packet: dict) -> int:
    count = 0
    for key in ("claimed_identities", "contested_identities"):
        for block in packet.get(key) or []:
            for records in (block.get("sources") or {}).values():
                for record in records or []:
                    pmcid = str(record.get("pmcid") or "").strip()
                    pmid = str(record.get("pmid") or "").strip()
                    if (pmcid and pmcid != NOT_AVAILABLE) or (pmid and pmid != NOT_AVAILABLE):
                        count += 1
    return count


def _citations_required(packet: dict) -> int:
    return min(REVIEW_CITATIONS_REQUIRED, _packet_sentence_count(packet))


def _cited_ok(verdict: dict, packet: dict) -> tuple[bool, list[str]]:
    available = _packet_citations(packet)
    required = _citations_required(packet)
    entries = verdict.get("evidence_cited")
    if not isinstance(entries, list):
        return False, ["evidence_cited must be a list of {identity, gene, citation}"]
    seen: set[tuple[str, str, str]] = set()
    problems: list[str] = []
    for entry in entries:
        if not isinstance(entry, dict):
            problems.append("each evidence_cited entry must be an object")
            continue
        key = (
            str(entry.get("identity") or ""),
            str(entry.get("gene") or "").upper(),
            str(entry.get("citation") or "").strip(),
        )
        if key in available:
            seen.add(key)
        else:
            problems.append(
                f"evidence_cited {key[2] or '(empty)'} for {key[1] or '(no gene)'} under "
                f"'{key[0]}' is not a sentence this packet supplied"
            )
    if len(seen) < required:
        problems.append(
            f"cite at least {required} distinct sentences from `sources` in this packet "
            f"(identity, gene and the citation exactly as shown); {len(seen)} of them were real"
        )
    return (not problems), problems


def _absent_definers(packet: dict) -> list[str]:
    selected = str((packet.get("delivered") or {}).get("selected") or "")
    audit = (packet.get("definer_audit") or {}).get(selected) or {}
    return [str(g) for g in (audit.get("definers_absent") or []) if str(g)]


ABSENT_DISPOSITIONS = ("dropout_single_gene", "tissue_context_absent", "refutes")


def _disposition_ok(verdict: dict, packet: dict) -> list[str]:
    absent = {g.upper() for g in _absent_definers(packet)}
    if not absent:
        return []
    entries = verdict.get("absent_definers_disposition")
    if not isinstance(entries, list):
        return [f"absent_definers_disposition must be a list, one entry per gene the "
                f"definer audit reads as absent here: {sorted(absent)[:6]}"]
    problems: list[str] = []
    seen = set()
    for entry in entries:
        if not isinstance(entry, dict):
            problems.append("each absent_definers_disposition entry must be an object")
            continue
        gene = str(entry.get("gene") or "").upper()
        if gene not in absent:
            problems.append(f"absent_definers_disposition names '{gene}', which the audit "
                            f"does not read as absent here")
            continue
        seen.add(gene)
        if entry.get("disposition") not in ABSENT_DISPOSITIONS:
            problems.append(f"{gene}: disposition must be one of {list(ABSENT_DISPOSITIONS)}")
        if not str(entry.get("reason") or "").strip():
            problems.append(f"{gene}: reason must say what the rows and the sentences show")
    missing = sorted(absent - seen)
    if missing:
        problems.append(f"every gene the audit reads as absent needs a disposition; "
                        f"missing {missing[:6]}")
    if (str(verdict.get("identity_verdict")) == IDENTITY_PASS
            and any(isinstance(e, dict) and e.get("disposition") == "refutes" for e in entries)):
        problems.append(
            "a definer marked `refutes` and an `accept` are the verdict disagreeing with "
            "itself: return 'reject' (with the candidate the evidence carries in "
            "better_candidate, where this packet holds it), or say what the rows show that "
            "makes that gene a dropout or absent from this tissue's context rather than a "
            "refutation"
        )
    return problems


def _detected_exclusions(packet: dict) -> list[str]:
    selected = str((packet.get("delivered") or {}).get("selected") or "")
    for block in packet.get("claimed_identities") or []:
        if str(block.get("cell_type") or "") == selected and str(block.get("role") or "") == "selected":
            return sorted(pool_api.compact_row_genes(block.get("detected_exclusions") or []))
    return []


EXCLUSION_DISPOSITIONS = ("background", "minority", "refutes")


def _exclusion_disposition_ok(verdict: dict, packet: dict) -> list[str]:
    detected = set(_detected_exclusions(packet))
    if not detected:
        return []
    entries = verdict.get("detected_exclusions_disposition")
    if not isinstance(entries, list):
        return [f"detected_exclusions_disposition must be a list, one entry per curated "
                f"exclusion of the delivered identity this cluster detects: {sorted(detected)[:6]}"]
    problems: list[str] = []
    seen = set()
    for entry in entries:
        if not isinstance(entry, dict):
            problems.append("each detected_exclusions_disposition entry must be an object")
            continue
        gene = str(entry.get("gene") or "").upper()
        if gene not in detected:
            problems.append(f"detected_exclusions_disposition names '{gene}', which is not a "
                            f"detected exclusion of the delivered identity here")
            continue
        seen.add(gene)
        if entry.get("disposition") not in EXCLUSION_DISPOSITIONS:
            problems.append(f"{gene}: disposition must be one of {list(EXCLUSION_DISPOSITIONS)}")
        if not str(entry.get("reason") or "").strip():
            problems.append(f"{gene}: reason must say what the rows and the sentences show")
    missing = sorted(detected - seen)
    if missing:
        problems.append(f"every detected exclusion of the delivered identity needs a "
                        f"disposition; missing {missing[:6]}")
    if (str(verdict.get("identity_verdict")) == IDENTITY_PASS
            and any(isinstance(e, dict) and e.get("disposition") == "refutes" for e in entries)):
        problems.append(
            "a detected exclusion marked `refutes` and an `accept` are the verdict disagreeing "
            "with itself: return 'reject' (with the candidate the evidence carries in "
            "better_candidate, where this packet holds it), or say what the rows show that "
            "makes that gene the tissue's background or a minority here rather than a "
            "refutation"
        )
    return problems


def _separating_sentences(packet: dict) -> list[str]:
    """Citable genes of the delivered identity that no contested name also curates."""
    selected = str((packet.get("delivered") or {}).get("selected") or "")
    shared = {str(gene).upper() for gene in (packet.get("shared_with_contested") or {})}
    return sorted({
        gene for identity, gene, _citation in _packet_citations(packet)
        if identity == selected and gene not in shared
    })


def _distinguishing_citation(verdict: dict, packet: dict) -> list[str]:
    contested = [str(b.get("cell_type")) for b in (packet.get("contested_identities") or [])]
    if not contested:
        return []
    shared = {str(gene).upper() for gene in (packet.get("shared_with_contested") or {})}
    selected = str((packet.get("delivered") or {}).get("selected") or "")
    own = [
        str(entry.get("gene") or "").upper()
        for entry in (verdict.get("evidence_cited") or [])
        if isinstance(entry, dict) and str(entry.get("identity") or "") == selected
    ]
    if any(gene and gene not in shared for gene in own):
        return []
    # Two names whose curated definitions the resource shares gene for gene leave nothing
    # to cite, and the demand would then refuse every verdict. What the evidence reaches
    # for such a pair is their common parent or a `mixed` answer, which is the reviewer's
    # reading to give.
    if not _separating_sentences(packet):
        return []
    return [
        f"every sentence cited for '{selected}' is for a gene that "
        f"{', '.join(repr(n) for n in contested)} also curates, so the citations do not "
        f"reach the distinction this verdict settles: cite at least one sentence for a "
        f"gene of '{selected}' that the contested identity does not claim"
    ]


def _valid_review(verdict: Any, packet: dict, subtype: str) -> tuple[bool, list[str]]:
    problems: list[str] = []
    if not isinstance(verdict, dict):
        return False, ["return one JSON object with the verdict"]
    if verdict.get("identity_verdict") not in IDENTITY_VALUES:
        problems.append(f"identity_verdict must be one of {list(IDENTITY_VALUES)}")
    elif (
        verdict.get("identity_verdict") == IDENTITY_PASS
        and str(verdict.get("better_candidate") or "").strip()
    ):
        problems.append(
            f"identity_verdict is '{IDENTITY_PASS}' but better_candidate is not empty -- "
            "naming a better candidate is a rejection of the delivered identity: return "
            "'reject' with that candidate in better_candidate, or clear better_candidate "
            "if the delivered identity is in fact the one to accept"
        )
    value = verdict.get("subtype_verdict", "")
    if value not in SUBTYPE_VALUES:
        problems.append(f"subtype_verdict must be one of {list(SUBTYPE_VALUES)}")
    elif subtype and not value:
        problems.append(
            f"a finer name ('{subtype}') was delivered, so subtype_verdict must be "
            "'separable' or 'not_separable'"
        )
    if not isinstance(verdict.get("reason", ""), str) or not str(verdict.get("reason") or "").strip():
        problems.append("reason must name genes with their measured percentages")
    ok, cited = _cited_ok(verdict, packet)
    problems.extend(cited)
    if ok:
        problems.extend(_distinguishing_citation(verdict, packet))
    problems.extend(_disposition_ok(verdict, packet))
    problems.extend(_exclusion_disposition_ok(verdict, packet))
    problems.extend(_population_ok(verdict, packet))
    return (not problems), problems


def _population_ok(verdict: dict, packet: dict) -> list[str]:
    required = stages.population_names(packet)
    entries = verdict.get("population_verdicts")
    if entries is None:
        entries = []
        if required:
            return [f"population_verdicts is required: one entry for each of {required[:6]}"]
    if not isinstance(entries, list):
        return ["population_verdicts must be a list"]
    problems: list[str] = []
    known = stages.packet_genes(packet)
    seen: set[str] = set()
    for entry in entries:
        if not isinstance(entry, dict):
            problems.append("each population_verdicts entry must be an object")
            continue
        name = str(entry.get("identity") or "")
        if name not in required:
            problems.append(
                f"population_verdicts names '{name}', which is neither a co-occurring "
                f"identity of the answer nor a contested identity of this packet"
            )
            continue
        seen.add(name)
        if entry.get("verdict") not in stages.POPULATION_VERDICTS:
            problems.append(f"{name}: verdict must be one of {list(stages.POPULATION_VERDICTS)}")
        genes = entry.get("genes")
        if not isinstance(genes, list) or not genes:
            problems.append(f"{name}: genes must list the measured genes the verdict rests on")
        else:
            bad = [g for g in genes if str(g).upper() not in known]
            if bad:
                problems.append(f"{name}: genes names genes not in this packet: {bad[:4]}")
        if not str(entry.get("reason") or "").strip():
            problems.append(f"{name}: reason must say what the rows and the sentences show")
    missing = [n for n in required if n not in seen]
    if missing:
        problems.append(f"every name the packet asks about needs a population verdict; missing {missing[:4]}")
    if str(verdict.get("identity_verdict")) == IDENTITY_PASS:
        undivided = [
            str(e.get("identity") or "") for e in entries
            if isinstance(e, dict) and e.get("verdict") == stages.NEITHER_SEPARATES
        ]
        if undivided:
            problems.append(
                f"'{stages.NEITHER_SEPARATES}' for {undivided[:3]} and an `accept` are the "
                f"verdict disagreeing with itself: where no gene in this packet tells the "
                f"delivered identity from that name, neither is established alone -- return "
                f"'reject' with their common parent in better_candidate where this packet "
                f"holds it (empty otherwise), or read the pair as one of the other verdicts"
            )
    return problems


def _review_passed(verdict: dict | None, subtype: str) -> bool:
    if not verdict:
        return False
    if str(verdict.get("identity_verdict")) != IDENTITY_PASS:
        return False
    return not subtype or str(verdict.get("subtype_verdict")) == SUBTYPE_PASS


def _review(
    pool: dict,
    cluster: str,
    final: dict,
    sources,
    tier: int,
    api_key: str,
    api_url: str,
    trace_id: str,
    turn: int,
    program: dict | None = None,
) -> tuple[dict | None, dict, str]:
    server = sources
    if sources is not None and hasattr(sources, "context"):
        server = SourceServer(
            sources.context, batch=SOURCES_PER_MARKER, max_batches=SOURCE_BATCHES_PER_MARKER
        )
    packet = pool_api.review_packet(
        pool, cluster, final, server, tier=tier, unknown_token=UNKNOWN
    )
    if program:
        packet["program_reading"] = program
    if not packet.get("claimed_identities"):
        return None, packet, "unreachable"
    required = stages.population_names(packet)
    if required:
        packet["population_verdicts_required"] = required
    if packet.get("contested_identities"):
        separating = _separating_sentences(packet)
        if separating:
            packet["separating_genes"] = separating
    subtype = str(final.get("subtype") or "")
    labels = tuple(
        str(block.get("cell_type"))
        for key in ("claimed_identities", "contested_identities")
        for block in (packet.get(key) or [])
        if block.get("cell_type")
    )
    prompt = _context_header(pool, cluster) + _render(_prompt("evidence_review"), packet)
    asked = 0
    unusable = 0
    last_verdict: dict | None = None
    seen: dict[str, dict] = {}
    while True:
        text, _error, _replayed = llm.cached_call_llm(
            prompt,
            api_url,
            api_key,
            reasoning_effort=REVIEW_EFFORT,
            model=REVIEW_MODEL,
            trace_id=trace_id,
            turn_index=turn,
        )
        parsed = llm.parse_json(text) if text else None
        if isinstance(parsed, dict) and parsed.get("identity_verdict"):
            ok, problems = _valid_review(parsed, packet, subtype)
            if ok:
                parsed["citations_required"] = _citations_required(packet)
                if packet.get("contested_identities"):
                    parsed["separating_sentences"] = _separating_sentences(packet)
                return parsed, packet, "ok"
            if str(parsed.get("identity_verdict")) in IDENTITY_VALUES:
                last_verdict = parsed
            if unusable >= SCHEMA_RETRIES:
                if last_verdict is not None:
                    last_verdict["citations_required"] = _citations_required(packet)
                    if packet.get("contested_identities"):
                        last_verdict["separating_sentences"] = _separating_sentences(packet)
                    last_verdict["uncited"] = True
                    last_verdict["uncited_problems"] = problems[:4]
                    return last_verdict, packet, "uncited"
                return None, packet, "unreachable"
            unusable += 1
            prompt += (
                f"\n\n# Retry {unusable}: the verdict was rejected -- "
                f"{'; '.join(problems[:4])}. Return one valid verdict object.\n"
            )
            continue
        answerable = (
            isinstance(parsed, dict) and valid_query(parsed) and asked < REVIEW_QUERY_TURNS
        )
        if not answerable:
            if unusable >= SCHEMA_RETRIES:
                return None, packet, "unreachable"
            unusable += 1
            prompt += f"\n\n# Retry {unusable}: return the verdict as one JSON object.\n"
            continue
        call = json.dumps(
            {"tool": parsed["tool"], "args": parsed.get("args", {})},
            sort_keys=True,
            separators=(",", ":"),
        )
        observation = None if parsed["tool"] == "sources" else seen.get(call)
        if observation is None:
            observation = pool_api.run_review_tool(
                pool, cluster, labels, parsed["tool"], parsed.get("args", {}), server
            )
            seen[call] = observation
        else:
            observation = {**observation, "duplicate_query": True}
        asked += 1
        prompt += _observation_block(asked, observation)
        if asked >= REVIEW_QUERY_TURNS:
            prompt += "\n# No further queries. Return the verdict now.\n"


def _rounds_note(tier: int, rounds_left: int, tiers_left: int) -> str:
    grow = (f" Only where that decision is Unknown are the next {tiers_left * 15} "
            "candidates added." if tiers_left else "")
    return (
        f"Candidates {tier * 15} of the retrieval order are on this page. "
        + (
            f"{rounds_left} further attempt(s) at this tier, then the answer is decided "
            "from the rounds so far." + grow
            if rounds_left
            else "This is the last attempt at this tier; after it the answer is decided "
                 "from the rounds so far." + grow
        )
    )


def _pair_tables(pool: dict, cluster: str, selected: str, rivals: list[str]) -> list[dict]:
    sel = pool_api.find_candidate(pool, cluster, selected)
    if sel is None:
        return []
    out = []
    for name in rivals:
        rentry = pool_api.find_candidate(pool, cluster, str(name))
        if rentry is not None:
            out.append(pool_api.pair_partition(sel, rentry))
    return out


def _contested_feedback(
    verdict: dict, reasons: list[str], audits: dict[str, dict],
    tier: int, rounds_left: int, tiers_left: int,
    pair_tables: list[dict] | None = None,
    population: list[str] | None = None,
) -> dict:
    block = {
        "evidence_review": {
            "identity_verdict": verdict.get("identity_verdict"),
            "subtype_verdict": verdict.get("subtype_verdict") or "",
            "supporting_markers": verdict.get("supporting_markers") or [],
            "conflicting_markers": verdict.get("conflicting_markers") or [],
            "evidence_cited": verdict.get("evidence_cited") or [],
            "population_verdicts": verdict.get("population_verdicts") or [],
            "reason": verdict.get("reason"),
            "flags": verdict.get("flags") or [],
        },
        "contested": reasons,
        "population_conflicts": list(population or []),
        "definer_audits": {
            name: {k: v for k, v in audit.items() if k != "species_definers"}
            for name, audit in (audits or {}).items()
        },
        "note": (
            "the evidence review accepted your identity. What did not close is listed "
            "above: under `contested`, the definer audit of the name you delivered, made on "
            "its OWN rows, contradicting that delivery; under `population_conflicts`, the "
            "reviewer's reading of each co-occurring identity you reported and of each "
            "contested identity you did not (`population_verdicts`) where it disagrees with "
            "your list. Answer again from the candidates you already hold, by [R8] and "
            "[A5]; `pair_partitions` below lists, for each disputed pair, what each side "
            "curates that the other does not and what both do, as measured here. "
            + _rounds_note(tier, rounds_left, tiers_left)
        ),
    }
    if pair_tables:
        block["pair_partitions"] = pair_tables
    return block


def _review_feedback(
    verdict: dict, tier: int, rounds_left: int, tiers_left: int,
    audits: dict[str, dict] | None = None, better_audit: dict | None = None,
    pair_tables: list[dict] | None = None,
) -> dict:
    block = {
        "evidence_review": {
            "identity_verdict": verdict.get("identity_verdict"),
            "subtype_verdict": verdict.get("subtype_verdict") or "",
            "better_candidate": str(verdict.get("better_candidate") or ""),
            "supporting_markers": verdict.get("supporting_markers") or [],
            "conflicting_markers": verdict.get("conflicting_markers") or [],
            "evidence_cited": verdict.get("evidence_cited") or [],
            "population_verdicts": verdict.get("population_verdicts") or [],
            "reason": verdict.get("reason"),
            "flags": verdict.get("flags") or [],
        },
    }
    if audits:
        block["definer_audit_of_your_claims"] = {
            name: {k: v for k, v in a.items() if k != "species_definers"}
            for name, a in audits.items()
        }
    if pair_tables:
        block["pair_partitions"] = pair_tables
    if better_audit:
        block["definer_audit_of_better_candidate"] = {
            k: v for k, v in better_audit.items() if k != "species_definers"
        }
        block["better_candidate_note"] = (
            "the reviewer's better_candidate was put through the same definer audit. A name "
            "the reviewer prefers is a name to check, not to adopt: take it as `selected` "
            "only where that audit reads its own program as `dominant` here on its OWN "
            "lineage genes. `minority` means at most a co-occurring identity beside a "
            "`selected` whose own program does occupy these cells; `absent` means the page "
            "does not carry that name at all, and the reviewer's objection to your answer "
            "still has to be met from the rows -- keep your name where its own lineage genes "
            "carry it, move to the common parent where neither name separates, or answer "
            "Unknown"
        )
    has_better = bool(str(verdict.get("better_candidate") or "").strip())
    legal = (
        "A refusal is not an instruction to move to another leaf of the page. The "
        "answers the evidence can carry are: (1) the name you gave, where its own lineage "
        "genes stand in these cells and you can answer the reviewer's objection from the "
        "rows and the sentences -- say which rows; "
        + ("(2) the reviewer's better_candidate, only under the conditions in "
           "`better_candidate_note`; " if has_better else
           "(2) -- the reviewer named no better candidate, so none is on offer; ")
        + "(3) another candidate whose OWN lineage genes occupy these cells at cluster "
        "level -- it is audited on the same terms before any review, and a name whose "
        "own program reads minority or absent there comes back refused; "
        "(4) their common parent, where this page holds it and neither separates; "
        "(5) a finer name dropped for the name it was delivered under, where the finer "
        "one is not separable; or (6) `Unknown`, where nothing on this page is carried "
        "by its own evidence -- say which measurements nothing accounts for. Moving to "
        "a candidate whose own program is weaker here than the one refused is none of "
        "these, and a co-occurring identity is reported by [A5] on its own lineage genes "
        "in its own share of the cells, not as a way of answering a refusal. "
    )
    return {
        **block,
        "note": (
            "your answer was checked against the curated source sentences behind the whole "
            "panel of every identity you claimed -- positive rows and curated exclusions "
            "alike -- and it did not survive. The reviewer read the sentences it cites "
            "above; it did not see your candidate ranking or your confidence. Answer again "
            "against the candidates you already hold. "
            + legal
            + _rounds_note(tier, rounds_left, tiers_left)
        ),
    }


def _arbitrate(
    pool: dict,
    cluster: str,
    attempts: list[dict],
    packet: dict,
    tier: int,
    api_key: str,
    api_url: str,
    trace_id: str,
    turn: int,
    contested: list[str] | None = None,
    audits: dict[str, dict] | None = None,
    pair_tables: list[dict] | None = None,
) -> dict | None:
    if not attempts:
        return None
    claimed: dict[str, dict] = {}
    for block in packet.get("claimed_identities") or []:
        claimed[str(block.get("cell_type"))] = block
    for item in attempts:
        for _role, name in claimed_identities(item["final"]):
            if name in claimed:
                continue
            entry = pool_api.find_candidate(pool, cluster, name)
            if entry is None:
                continue
            claimed[name] = {
                "cell_type": name,
                "panel": pool_api.compact_panel(pool_api._review_panel(entry)),
                "facts": pool_api.packet_facts(pool_api.candidate_facts(pool, cluster, entry), entry),
                "detected_exclusions": pool_api.compact_panel(
                    pool_api._detected_exclusion_markers(entry)
                ),
            }
    rounds = [
        {
            "round": index + 1,
            "tier": int(item["tier"]),
            "answer": {
                key: item["final"].get(key)
                for key in (
                    "selected",
                    "subtype",
                    "lineage",
                    "state",
                    "co_occurring_identities",
                    "possible_refinements",
                    "support_markers",
                    "claim_evidence",
                    "confidence",
                    "reason",
                )
            },
            "review": (
                {
                    key: (item["review"] or {}).get(key)
                    for key in (
                        "identity_verdict",
                        "subtype_verdict",
                        "better_candidate",
                        "supporting_markers",
                        "conflicting_markers",
                        "evidence_cited",
                        "population_verdicts",
                        "reason",
                        "flags",
                        "uncited",
                    )
                }
                if item.get("review")
                else {"identity_verdict": "not_checked"}
            ),
            "definer_audit": {
                name: {k: v for k, v in a.items() if k != "species_definers"}
                for name, a in (item["final"].get("definer_audit") or {}).items()
            },
            "structural_note": item.get("structural_note") or "",
        }
        for index, item in enumerate(attempts)
    ]
    audits_with_rows = dict(audits or {})
    for item in attempts:
        for name, a in (item["final"].get("definer_audit") or {}).items():
            audits_with_rows.setdefault(name, a)
    for block in packet.get("contested_identities") or []:
        claimed.setdefault(str(block.get("cell_type")), block)
    arbiter_packet = {
        "query": packet.get("query") or {},
        "rounds": rounds,
        "row_format": pool_api.MARKER_ROW_FORMAT,
        "identities_claimed": list(claimed.values()),
        "definer_audits": audits_with_rows,
        "program_reading": packet.get("program_reading") or {},
        "candidates_not_claimed": packet.get("candidates_not_claimed") or [],
        "cluster_top_enriched": packet.get("cluster_top_enriched") or [],
        "cluster_top_detection_gap": packet.get("cluster_top_detection_gap") or [],
        "cluster_top_depleted": packet.get("cluster_top_depleted") or [],
        "candidates_available": pool_api.candidate_names(pool, cluster, tier=tier),
    }
    screened = pool_api.screened_out_claimants(pool, cluster)
    if screened:
        arbiter_packet["screened_out_claimants"] = screened
    if contested:
        arbiter_packet["contested"] = contested
    if pair_tables:
        arbiter_packet["pair_partitions"] = pair_tables
    prompt = _context_header(pool, cluster) + _render(
        _prompt("final_arbiter"), arbiter_packet
    )
    names = set(arbiter_packet["candidates_available"])
    verdicts_of = _population_verdicts_by_name(attempts)
    for attempt in range(SCHEMA_RETRIES + 1):
        text, _error, _replayed = llm.cached_call_llm(
            prompt,
            api_url,
            api_key,
            reasoning_effort=ARBITER_EFFORT,
            model=ARBITER_MODEL,
            trace_id=trace_id,
            turn_index=turn,
        )
        parsed = llm.parse_json(text) if text else None
        problems = _arbiter_problems(parsed, names, verdicts_of)
        if not problems:
            return parsed
        if attempt >= SCHEMA_RETRIES:
            return None
        prompt += (
            f"\n\n# Retry {attempt + 1}: the decision was rejected -- "
            f"{'; '.join(problems[:4])}. Return one valid decision object.\n"
        )
    return None


def _population_verdicts_by_name(attempts: list[dict]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for item in attempts:
        for entry in ((item.get("review") or {}).get("population_verdicts") or []):
            if isinstance(entry, dict) and entry.get("identity"):
                out.setdefault(str(entry["identity"]), []).append(str(entry.get("verdict") or ""))
    return out


def _arbiter_problems(value: Any, names: set[str],
                      verdicts_of: dict[str, list[str]] | None = None) -> list[str]:
    if not isinstance(value, dict):
        return ["return one JSON object with the decision"]
    problems: list[str] = []
    selected = value.get("selected")
    if not isinstance(selected, str) or not selected.strip():
        problems.append("selected must be a candidate name or Unknown")
    elif selected != UNKNOWN and selected not in names:
        problems.append(f"selected '{selected}' is not one of the candidates offered")
    for field in ("subtype", "lineage"):
        name = value.get(field, "")
        if not isinstance(name, str):
            problems.append(f"{field} must be a string")
        elif name and name != selected and name not in names:
            problems.append(f"{field} '{name}' is not one of the candidates offered")
    others = value.get("co_occurring_identities", [])
    if not isinstance(others, list):
        problems.append("co_occurring_identities must be a list")
    else:
        for name in others:
            if not isinstance(name, str) or (name not in names):
                problems.append(f"co-occurring '{name}' is not one of the candidates offered")
                continue
            read = (verdicts_of or {}).get(name) or []
            if read and stages.ESTABLISHED not in read:
                problems.append(
                    f"co-occurring '{name}' was read by the rounds' reviews as "
                    f"{sorted(set(read))} and by none as an established second population; "
                    f"a co-occurring identity is delivered only on a review that read it so "
                    f"([F2c]) -- leave it out"
                )
    if value.get("confidence") not in CONFIDENCE_VALUES:
        problems.append(f"confidence must be one of {list(CONFIDENCE_VALUES)}")
    if not isinstance(value.get("reason", ""), str) or not str(value.get("reason") or "").strip():
        problems.append("reason must say which evidence decided it")
    return problems


def _arbiter_final(base: dict, decision: dict) -> dict:
    selected = str(decision.get("selected") or "")
    subtype = str(decision.get("subtype") or "")
    keep = {name for name in (selected, subtype) if name and name != UNKNOWN} | {
        str(name) for name in (decision.get("co_occurring_identities") or [])
    }
    return {
        **base,
        "selected": selected,
        "subtype": subtype,
        "lineage": str(decision.get("lineage") or "") or selected,
        "co_occurring_identities": [
            str(name) for name in (decision.get("co_occurring_identities") or [])
        ],
        "state": str(decision.get("state") or base.get("state") or ""),
        "confidence": str(decision.get("confidence") or base.get("confidence") or "low"),
        "reason": str(decision.get("reason") or ""),
        "claim_evidence": [
            item
            for item in (base.get("claim_evidence") or [])
            if isinstance(item, dict) and str(item.get("identity") or "") in keep
        ],
        "support_markers": list(base.get("support_markers") or []),
    }


def annotate_cluster(
    pool: dict,
    cluster: str,
    template: str,
    api_key: str,
    api_url: str,
    sources=None,
    consistency_note: dict | None = None,
) -> dict[str, Any]:
    trace_id = uuid.uuid4().hex[:16]
    tiers = pool_api.tiers_available(pool, cluster)
    tier = 1
    header = _context_header(pool, cluster)
    program = stages.program_reading(
        pool, cluster, tier, _prompt("program_reading"), header, api_key, api_url, trace_id, 0
    )
    opening = cluster_packet(pool, cluster, tier=tier)
    opening["program_reading"] = program or {}
    if consistency_note:
        opening["dataset_consistency_note"] = consistency_note
    prompt = header + _render(template, opening)
    audit_template = _prompt("definer_audit")
    transcript: list[dict[str, Any]] = []
    attempts: list[dict[str, Any]] = []
    seen: dict[str, dict] = {}
    turn = 0
    schema_failures = 0
    repeated = 0
    rounds_in_tier = 0
    forced = False
    last_packet: dict = {}
    last_refused_answer: tuple = ()
    review_unreachable = False
    structural_notes: list[str] = []
    audit_cache: dict[tuple[int, str], dict] = {}

    def _answer_key(final: dict) -> tuple:
        return (
            str(final.get("selected") or ""),
            str(final.get("subtype") or ""),
            tuple(sorted(str(n) for n in (final.get("co_occurring_identities") or []))),
        )

    def _escalate(reason: str) -> bool:
        nonlocal tier, rounds_in_tier, prompt
        if tier >= tiers:
            return False
        tier += 1
        rounds_in_tier = 0
        block = pool_api.tier_packet(pool, cluster, tier)
        block["note"] = (
            f"{reason} Candidates ranked {(tier - 1) * pool_api.CANDIDATE_TIER_SIZE + 1} "
            f"and below in the retrieval order are added here, with their whole measured "
            f"panels and facts. Everything already on this page is still available: the "
            f"earlier candidates were not withdrawn, and one of them may still be the "
            f"answer. Answer again."
        )
        prompt += _observation_block(len(transcript) + len(attempts) + 1, block)
        return True

    def _deliver(final: dict, arbitration: dict | None,
                 contested: list[str] | None = None) -> dict[str, Any]:
        capped, why = None, ""
        if str(final.get("selected") or "") != UNKNOWN:
            cap, why = stages.confidence_cap(final, final.get("definer_audit") or {})
            capped = cap if stages.apply_cap(final, cap) else None
        last = attempts[-1] if attempts else None
        last_review = (last or {}).get("review")
        if (
            last_review
            and str(final.get("selected") or "") == str(last["final"].get("selected") or "")
        ):
            minority = [
                str(e.get("gene") or "") for e in (last_review.get("detected_exclusions_disposition") or [])
                if isinstance(e, dict) and e.get("disposition") == "minority"
            ]
            if minority and stages.apply_cap(final, "medium"):
                capped, why = "medium", (
                    f"the review reads {minority[:4]} as detected exclusions of "
                    f"'{final.get('selected')}' that are real in a minority of these cells"
                )
        if review_unreachable:
            stages.apply_cap(final, "low")
        if contested:
            stages.apply_cap(final, "low")
        passed = bool(
            last is not None
            and last_review is not None
            and not last_review.get("uncited")
            and _review_passed(last_review, str(last["final"].get("subtype") or ""))
            and arbitration is None
            and not last.get("structural_note")
        )
        return {
            "final": final,
            "agent_selected": str((attempts[-1]["final"] if attempts else final).get("selected") or ""),
            "arbitration": arbitration,
            "turns": turn,
            "transcript": transcript,
            "trace_id": trace_id,
            "turn_budget_exhausted": forced,
            "delivered_tier": tier,
            "tiers_available": tiers,
            "program_reading": program or {},
            "definer_audit": final.get("definer_audit") or {},
            "contested_rivals": list(final.get("contested_rivals") or []),
            "contested": list(contested or []),
            "confidence_cap": {"cap": capped, "why": why} if capped else None,
            "review_unreachable": review_unreachable,
            "structural_notes": list(structural_notes),
            "consistency_note": consistency_note,
            "review": {
                "rounds": [
                    {
                        "round": index + 1,
                        "tier": int(item["tier"]),
                        "selected": str(item["final"].get("selected") or ""),
                        "subtype": str(item["final"].get("subtype") or ""),
                        "verdict": item.get("review") or None,
                        "structural_note": item.get("structural_note") or "",
                    }
                    for index, item in enumerate(attempts)
                ],
                "checked": any(item.get("review") for item in attempts),
                "passed": passed,
            },
        }

    def _audit_delivery(final: dict, offset: int) -> list[str]:
        audits = dict(final.get("definer_audit") or {})
        missing = [n for _r, n in claimed_identities(final) if n not in audits]
        if missing:
            audits.update(stages.audit_all(pool, cluster, final, tier, sources, program,
                                           audit_template, header, api_key, api_url,
                                           trace_id, turn + offset, cache=audit_cache))
        final["definer_audit"] = audits
        final["contested_rivals"] = stages.contested_rival_names(final, audits)
        return stages.delivery_contested(final, audits)

    def _from_decision(parsed: dict, decision: dict) -> dict:
        if str(decision.get("selected") or "") == UNKNOWN:
            return {**parsed, "selected": UNKNOWN, "subtype": "", "lineage": "",
                    "co_occurring_identities": [], "claim_evidence": [], "support_markers": [],
                    "reason": str(decision.get("reason") or ""),
                    "confidence": str(decision.get("confidence") or "low")}
        base = next(
            (item["final"] for item in reversed(attempts)
             if str(item["final"].get("selected") or "") == str(decision.get("selected") or "")),
            parsed,
        )
        return _arbiter_final(base, decision)

    def _arbitrate_now(parsed: dict) -> dict[str, Any] | None:
        nonlocal review_unreachable, last_refused_answer
        decision = _arbitrate(
            pool, cluster, attempts, last_packet, tier, api_key, api_url, trace_id, turn
        )
        if decision is None:
            review_unreachable = True
            structural_notes.append("arbiter unreachable; delivered unresolved at low confidence")
            return _deliver(parsed, None)
        if str(decision.get("selected") or "") == UNKNOWN and not forced and tier < tiers:
            why = str(decision.get("reason") or "").strip()
            note = (f"arbitration at tier {tier} answered Unknown: no candidate on this page "
                    "is carried by its own evidence; the next tier of candidates is added")
            structural_notes.append(note)
            attempts.append({"tier": tier, "review": None, "structural_note": note,
                             "final": {"selected": UNKNOWN, "subtype": "", "lineage": "",
                                       "co_occurring_identities": [], "claim_evidence": [],
                                       "support_markers": [], "confidence": "low",
                                       "reason": why, "arbitration_unknown": True}})
            _escalate(
                "The rounds at this tier were decided by arbitration, and the arbitration "
                "found no candidate on this page carried by its own evidence"
                + (f" ({why[:400]})" if why else "") + "."
            )
            last_refused_answer = ()
            return None
        final = _from_decision(parsed, decision)
        contested = _audit_delivery(final, 50)
        if contested:
            note = ("the arbitration was read again because " + "; ".join(contested))
            structural_notes.append(note)
            again = _arbitrate(
                pool, cluster, attempts, last_packet, tier, api_key, api_url, trace_id,
                turn + 60, contested=contested, audits=final.get("definer_audit"),
                pair_tables=_pair_tables(pool, cluster, str(final.get("selected") or ""),
                                         final.get("contested_rivals") or []),
            )
            if again is not None:
                decision, final = again, _from_decision(parsed, again)
                contested = _audit_delivery(final, 70)
        if contested:
            structural_notes.append(
                "delivered unresolved: " + "; ".join(contested)
            )
        return _deliver(final, decision, contested=contested)

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
        parsed, dropped_support = sanitize_support_markers(parsed, pool, cluster)
        if dropped_support:
            print(
                f"  [warn] cluster {cluster}: dropped support_markers not on "
                f"{parsed.get('selected')!r} panel: {dropped_support}"
            )
        if valid_final(parsed, pool, cluster, tier, program):
            schema_failures = 0
            if str(parsed.get("selected") or "") == UNKNOWN:
                if _escalate("You reported that no candidate on the page is carried."):
                    continue
                return _deliver(parsed, None)

            selected_now = str(parsed.get("selected") or "")
            parsed["definer_audit"] = stages.audit_all(
                pool, cluster, parsed, tier, sources, program, audit_template, header,
                api_key, api_url, trace_id, turn, cache=audit_cache,
            )
            parsed["contested_rivals"] = stages.contested_rival_names(
                parsed, parsed["definer_audit"]
            )

            if attempts and _answer_key(parsed) == last_refused_answer:
                note = (f"'{selected_now}' was re-delivered unchanged after a refusal; "
                        "rounds ended, decided by arbitration")
                structural_notes.append(note)
                attempts.append({"tier": tier, "final": parsed, "review": None,
                                 "structural_note": note})
                decided = _arbitrate_now(parsed)
                if decided is None:
                    continue
                return decided

            review, last_packet, status = _review(
                pool, cluster, parsed, sources, tier, api_key, api_url, trace_id, turn,
                program=program,
            )
            attempt = {"tier": tier, "final": parsed, "review": review}
            attempts.append(attempt)
            subtype = str(parsed.get("subtype") or "")
            if status == "unreachable" or review is None:
                review_unreachable = True
                structural_notes.append("evidence review unreachable; delivered unresolved at low confidence")
                return _deliver(parsed, None)
            if status == "uncited":
                note = (f"review verdict '{review.get('identity_verdict')}' could not be "
                        "validated (citations short after retries); decided by arbitration")
                structural_notes.append(note)
                attempt["structural_note"] = note
                decided = _arbitrate_now(parsed)
                if decided is None:
                    continue
                return decided
            if _review_passed(review, subtype):
                contested = stages.delivery_contested(parsed, parsed["definer_audit"])
                population = stages.population_conflicts(parsed, review)
                if not contested and not population:
                    return _deliver(parsed, None)
                note = ("review accepted, but " + "; ".join(contested + population)
                        + "; sent back to the annotator")
                structural_notes.append(note)
                attempt["structural_note"] = note
                rounds_in_tier += 1
                rounds_left = max(0, REVIEW_MAX_ROUNDS_PER_TIER - rounds_in_tier)
                if rounds_left and not forced:
                    last_refused_answer = _answer_key(parsed)
                    prompt += _observation_block(
                        len(transcript) + len(attempts) + 1,
                        _contested_feedback(review, contested, parsed["definer_audit"],
                                            tier, rounds_left, tiers - tier,
                                            _pair_tables(pool, cluster, selected_now,
                                                         parsed.get("contested_rivals") or []),
                                            population),
                    )
                    continue
                decided = _arbitrate_now(parsed)
                if decided is None:
                    continue
                return decided
            last_refused_answer = _answer_key(parsed)
            better = str(review.get("better_candidate") or "")
            better_audit = None
            if better and better != UNKNOWN and better not in parsed["definer_audit"]:
                better_audit = stages.definer_audit(
                    pool, cluster, better, "better_candidate",
                    [n for _r, n in claimed_identities(parsed)] + [better], tier, sources,
                    program, audit_template, header, api_key, api_url, trace_id, turn + 70,
                )
            rounds_in_tier += 1
            rounds_left = max(0, REVIEW_MAX_ROUNDS_PER_TIER - rounds_in_tier)
            versus = [n for n in ([better] if better and better != UNKNOWN else [])
                      + list(parsed.get("contested_rivals") or []) if n != selected_now]
            tables = _pair_tables(pool, cluster, selected_now, versus)
            if rounds_left and not forced:
                prompt += _observation_block(
                    len(transcript) + len(attempts) + 1,
                    _review_feedback(review, tier, rounds_left, tiers - tier,
                                     parsed["definer_audit"], better_audit, tables),
                )
                continue
            decided = _arbitrate_now(parsed)
            if decided is None:
                continue
            return decided

        if valid_query(parsed) and not forced:
            schema_failures = 0
            call = json.dumps(
                {"tool": parsed["tool"], "args": parsed.get("args", {})},
                sort_keys=True,
                separators=(",", ":"),
            )
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
            prompt += _observation_block(len(transcript) + len(attempts), observation)
            if repeated >= MAX_REPEATED_QUERIES or (MAX_TURNS and turn >= MAX_TURNS):
                forced = True
                prompt += "\n# No further queries. Return the final answer now.\n"
            continue

        schema_failures += 1
        if schema_failures > SCHEMA_RETRIES:
            return {
                "error": "schema violation",
                "turns": turn,
                "transcript": transcript,
                "trace_id": trace_id,
            }
        problems = (
            final_problems(parsed, pool, cluster, tier, program)
            if isinstance(parsed, dict) and parsed.get("action") == "final"
            else []
        )
        if problems:
            listed = "; ".join(problems[:6])
            prompt += (
                f"\n\n# Retry {schema_failures}: the final answer was rejected -- {listed}. "
                "Return one valid final JSON object.\n"
            )
        else:
            prompt += f"\n\n# Retry {schema_failures}: return one valid JSON object.\n"


def _measurements(entry: dict, percentile: dict, cluster: str, limit: int = 30):
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
                "auc": marker.get("auc"),
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
    return {
        "cell_type": entry["cell_type"],
        "retrieval_rank": entry["retrieval_rank"],
        "tier": pool_api.tier_of(entry),
        "borrowed_context": entry.get("borrowed_context") or None,
        "tissue_context": list(entry.get("tissue_context") or []),
        "retrieval_context": str(entry.get("retrieval_context") or ""),
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
    if frame.empty:
        return {
            "annotation_qc": QC_UNCHECKED,
            "review": {},
            "cluster_id": cluster,
            "annotation": UNKNOWN,
            "subtype": "",
            "lineage": "",
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
        }
    top = frame.iloc[0]
    return {
        "annotation_qc": QC_UNCHECKED,
        "review": {},
        "cluster_id": cluster,
        "annotation": str(top["candidate"]),
        "subtype": "",
        "lineage": str(top["candidate"]),
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
    }


QC_PASSED = "passed"
QC_REVISED = "passed_after_revision"
QC_ARBITRATED = "arbitrated"
QC_FAILED = "failed"
QC_UNCHECKED = "not_checked"
QC_VALUES = (QC_PASSED, QC_REVISED, QC_ARBITRATED, QC_FAILED, QC_UNCHECKED)


def _qc_status(outcome: dict) -> str:
    review = outcome.get("review") or {}
    rounds = review.get("rounds") or []
    if outcome.get("review_unreachable") or not review.get("checked"):
        return QC_UNCHECKED
    if outcome.get("arbitration"):
        return QC_ARBITRATED
    if not review.get("passed"):
        return QC_FAILED
    return QC_REVISED if len(rounds) > 1 else QC_PASSED


def _delivered_facts(pool: dict, cluster: str, names: set[str]) -> dict[str, dict]:
    precomputed = (pool.get("facts") or {}).get(cluster, [])
    out: dict[str, dict] = {}
    for name in names:
        if not name or name == UNKNOWN:
            continue
        row = next((r for r in precomputed if r["cell_type"] == name), None)
        if row is None:
            entry = pool_api.find_candidate(pool, cluster, name)
            row = pool_api.candidate_facts(pool, cluster, entry) if entry else None
        if row is not None:
            out[name] = row
    return out


def _result(
    cluster: str,
    pool: dict,
    percentile: dict,
    outcome: dict,
    frame: pd.DataFrame,
) -> dict[str, Any]:
    names = pool_api.candidate_names(pool, cluster)
    final = outcome["final"]
    selected = str(final["selected"])
    subtype = str(final.get("subtype", ""))
    others = [str(name) for name in final.get("co_occurring_identities", [])]
    lineage = str(final.get("lineage") or "") or selected

    claims = claimed_identities(final)
    warned = pool_api.claim_warnings(pool, cluster, claims)
    arbitration = outcome.get("arbitration") or None
    if arbitration:
        agent = str(outcome.get("agent_selected") or "")
        if agent and agent != selected:
            warned.append((
                agent,
                f"arbitration {agent} -> {selected}: the evidence review refused the "
                f"annotator's answer and the arbitration turn delivered this one | "
                f"{str(arbitration.get('reason') or '')[:600]}",
            ))
    warning_of: dict[str, str] = {}
    for name, line in warned:
        warning_of[name] = f"{warning_of[name]} || {line}" if name in warning_of else line
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

    contested = list(outcome.get("contested") or [])
    if selected == UNKNOWN:
        status, detail = UNRESOLVED, "no_candidate_carried_by_its_evidence"
    elif outcome.get("review_unreachable"):
        status, detail = UNRESOLVED, "evidence_review_unreachable"
    elif contested:
        status, detail = UNRESOLVED, "contested_after_rounds"
    elif others:
        status, detail = MIXED, "agent_majority_of_several_identities"
    else:
        status, detail = RESOLVED, "agent_selected"
    cap = outcome.get("confidence_cap") or None
    if cap:
        warnings.append(f"confidence capped at {cap['cap']}: {cap['why']}")
    for note in outcome.get("structural_notes") or []:
        warnings.append(f"structural: {note}")
    return {
        "program_reading": dict(outcome.get("program_reading") or {}),
        "definer_audit": dict(outcome.get("definer_audit") or {}),
        "contested": contested,
        "contested_rivals": list(outcome.get("contested_rivals") or []),
        "confidence_cap": cap,
        "structural_notes": list(outcome.get("structural_notes") or []),
        "consistency_note": outcome.get("consistency_note") or None,
        "cluster_id": cluster,
        "annotation": selected,
        "subtype": subtype,
        "lineage": lineage,
        "agent_selected": str(outcome.get("agent_selected") or final["selected"]),
        "arbitration": arbitration,
        "possible_refinements": sorted(
            {str(x) for x in (final.get("possible_refinements") or []) if str(x)}
        ),
        "candidate_facts": _delivered_facts(pool, cluster, {selected, subtype}),
        "borrowed_context": (pool.get("borrowed") or {}).get(cluster) or {},
        "delivered_borrowed": bool(
            (pool_api.find_candidate(pool, cluster, subtype or selected) or {}).get(
                "borrowed_context"
            )
        )
        if selected != UNKNOWN
        else False,
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
        "candidates": names,
        "candidate_entries": entries,
        "delivered_tier": int(outcome.get("delivered_tier") or 1),
        "tiers_available": int(outcome.get("tiers_available") or 1),
        "turns": int(outcome.get("turns", 0)),
        "tool_calls": [item["call"] for item in outcome.get("transcript", [])],
        "turn_budget_exhausted": bool(outcome.get("turn_budget_exhausted")),
        "annotation_qc": _qc_status(outcome),
        "review": dict(outcome.get("review") or {}),
        "review_flags": [
            {
                "round": item.get("round"),
                "label": item.get("selected"),
                "text": str(text)[:300],
            }
            for item in (outcome.get("review") or {}).get("rounds") or []
            for text in ((item.get("verdict") or {}).get("flags") or [])
            if str(text or "").strip()
        ],
    }


def _percentiles(de: pd.DataFrame) -> dict[tuple[str, str], float]:
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


def _screen_problems(value: Any, names: list[str]) -> list[str]:
    if not isinstance(value, dict) or not isinstance(value.get("screen"), list):
        return ["return one JSON object with a `screen` list"]
    problems: list[str] = []
    seen: set[str] = set()
    for entry in value["screen"]:
        if not isinstance(entry, dict):
            problems.append("each screen entry must be an object")
            continue
        name = str(entry.get("name") or "")
        if name not in names:
            problems.append(f"'{name}' is not one of the names given")
            continue
        seen.add(name)
        if entry.get("verdict") not in borrowed_context.SCREEN_VERDICTS:
            problems.append(f"{name}: verdict must be one of {list(borrowed_context.SCREEN_VERDICTS)}")
    missing = [n for n in names if n not in seen]
    if missing:
        problems.append(f"every name needs a verdict; missing {missing[:6]}")
    return problems


def borrow_tissue_screen(context: dict, api_key: str, api_url: str, record: dict):
    template = Path(PROMPT_DIR, f"{BORROW_SCREEN_PROMPT}.txt").read_text(encoding="utf-8")
    disease = context.get("disease")
    disease = ", ".join(str(d) for d in disease) if isinstance(disease, list) else str(disease)
    tissue = context.get("tissue")
    tissue = ", ".join(str(t) for t in tissue) if isinstance(tissue, (list, tuple)) else str(tissue)
    dataset_context = {
        "species": str(context.get("species")),
        "tissue": tissue,
        "disease": disease,
        "development_stage": str(context.get("development_stage") or ""),
    }
    record.update({
        "model": llm.resolve_model(BORROW_SCREEN_MODEL),
        "reasoning_effort": BORROW_SCREEN_EFFORT,
        "status": "not_run",
        "names": [],
        "foreign": [],
        "verdicts": {},
    })

    def screen(names: list[str]) -> dict[str, dict]:
        names = [str(n) for n in names]
        record["names"] = list(names)
        if not names:
            record["status"] = "nothing_to_screen"
            return {}
        prompt = _render(template, {"dataset_context": dataset_context, "names": names})
        trace_id = uuid.uuid4().hex[:16]
        parsed = None
        for attempt in range(SCHEMA_RETRIES + 1):
            text, error, _replayed = llm.cached_call_llm(
                prompt, api_url, api_key,
                reasoning_effort=BORROW_SCREEN_EFFORT, model=BORROW_SCREEN_MODEL,
                trace_id=trace_id, turn_index=attempt,
            )
            candidate = llm.parse_json(text) if text else None
            problems = _screen_problems(candidate, names)
            if not problems:
                parsed = candidate
                break
            if attempt >= SCHEMA_RETRIES:
                record["status"] = f"unavailable: {error or '; '.join(problems[:3])}"
                break
            prompt += (
                f"\n\n# Retry {attempt + 1}: the answer was rejected -- "
                f"{'; '.join(problems[:4])}. Return one valid screen object.\n"
            )
        if parsed is None:
            return {}
        verdicts = {
            str(e["name"]): {"verdict": str(e["verdict"]), "reason": str(e.get("reason") or "")}
            for e in parsed["screen"]
        }
        record["status"] = "ok"
        record["verdicts"] = verdicts
        record["foreign"] = sorted(
            n for n, v in verdicts.items() if v["verdict"] == borrowed_context.SCREEN_FOREIGN
        )
        return verdicts

    return screen


def prepare_pool(scoring: dict, prep: dict, markers_all: pd.DataFrame | None, source_db: SourceDB | None = None,
                 borrow_screen=None):
    from .marker_database import MarkerDatabase

    context = scoring["context"]
    source_db = source_db or SourceDB()
    native_sources = source_db.context(context["species"], context["tissue"], context["disease"])
    pool = build_pool(scoring, prep, native_sources, markers_all=markers_all)
    if len(pool.get("tissue_contexts") or []) > 1:
        by_retrieval = {}
        for name in pool["candidate_retrieval_context"].values():
            by_retrieval[name] = by_retrieval.get(name, 0) + 1
        print(f"  tissue contexts (operator list): {pool['tissue_contexts']}; {pool.get('retrieval_semantics')}; "
              f"eligible candidates by the context that admits them: {by_retrieval}")
    if markers_all is not None:
        pool["genome_de"] = pool_api.genome_de_summary(markers_all)
    else:
        pool["genome_de"] = {}
        print("  [warn] markers_all not found; the evidence review runs without genome-wide DE and nothing is borrowed")
    dbf = MarkerDatabase().load()
    sources = native_sources
    pool["borrowed"] = {}
    donor_species_by_name: dict[str, set[str]] = {}
    cross_species_kwargs: dict = {}
    if markers_all is not None:
        measured = {str(g).upper() for g in prep["genes"]}
        resource = borrowed_context.resource_slice(dbf, context["species"], context["disease"])
        if CROSS_SPECIES_BORROW:
            from .ortho_map import OrthoMap, SPECIES as ORTHO_SPECIES

            other_species = [s for s in ORTHO_SPECIES if s != str(context["species"])]
            cross_species_kwargs = {
                "species": str(context["species"]),
                "cross_species_resource": {
                    other: borrowed_context.resource_slice(dbf, other, context["disease"])
                    for other in other_species
                },
                "ortho": OrthoMap(),
            }
        pool["borrowed"] = borrowed_context.borrow_candidates(
            pool, resource, measured, markers_all, sources=None, screen=borrow_screen,
            **cross_species_kwargs)
        names = borrowed_context.borrowed_names(pool["borrowed"])
        donor_species_by_name = borrowed_context.borrowed_donor_species(pool["borrowed"])
        screened_out = sorted({
            str(item["cell_type"])
            for log in pool["borrowed"].values()
            for item in (log.get("considered") or [])
            if str(item.get("rejected") or "").startswith(borrowed_context.SCREEN_REJECTED_PREFIX)
        })
        if screened_out:
            print(f"  borrowed-name tissue screen: {len(screened_out)} name(s) read foreign to "
                  f"{context['tissue']} and not admitted ({', '.join(screened_out)})")
        if names:
            extra = source_db.cell_types_across_tissues(
                context["species"], names, context["disease"],
                donor_species_by_name=donor_species_by_name, ortho=cross_species_kwargs.get("ortho"),
            )
            sources = CompositeSources(native_sources, extra, names)
            for state in pool["clusters"].values():
                for entry in state["candidates"]:
                    if entry.get("borrowed_context"):
                        entry["exclusion_sources"] = pool_api._exclusion_sources(entry["cell_type"], entry["markers"], sources)
        n_events = sum(len(v.get("borrowed") or []) for v in pool["borrowed"].values())
        print(f"  borrowed-context retrieval: {n_events} candidates borrowed in "
              f"{sum(1 for v in pool['borrowed'].values() if v.get('borrowed'))} clusters ({', '.join(names) or 'none'})")
    recf = dbf[(dbf["species"].astype(str) == str(context["species"])) & (dbf["marker_polarity"] == "positive")
               & (dbf["is_recommended_marker"].astype(str).str.upper() == "TRUE")]
    rec_pairs = set(zip(recf["cell_type"].astype(str), recf["gene_symbol"].astype(str).str.upper()))
    marked = pool_api.attach_recommended(pool, rec_pairs)
    species_rows = MarkerDatabase.species_positive_rows(dbf, str(context["species"]))
    all_names = sorted({str(e["cell_type"]) for s in pool["clusters"].values() for e in s["candidates"]})
    definers = MarkerDatabase.species_definers(species_rows, all_names, k=pool_api.DEFINERS_PER_CANDIDATE)
    with_definers = pool_api.attach_definers(pool, definers)
    pool["species_sources"] = CompositeSources(
        native_sources,
        source_db.cell_types_across_tissues(
            context["species"], all_names, context["disease"],
            donor_species_by_name=donor_species_by_name, ortho=cross_species_kwargs.get("ortho"),
        ),
        all_names,
    )
    pool["species_claimants"] = MarkerDatabase.species_gene_claimants(species_rows)
    shared_rows = pool_api.attach_page_sharing(pool, pool["species_claimants"])
    pool_api.attach_genome_claims(
        pool,
        pool["species_claimants"],
        MarkerDatabase.known_genes(dbf),
        markers_all,
    )
    print(f"  evidence: species-wide definers on {with_definers} candidate entries "
          f"({len(definers)} names); {shared_rows} rows marked as shared with another page candidate")
    pool_api.assign_tiers(pool, scoring)
    pool["facts"] = {c: pool_api.candidate_facts_table(pool, c) for c in pool["clusters"]}
    per_tier: dict[int, int] = {}
    for state in pool["clusters"].values():
        for entry in state["candidates"]:
            key = pool_api.tier_of(entry)
            per_tier[key] = per_tier.get(key, 0) + 1
    print(f"  facts layer: {len(rec_pairs)} recommended (cell type, gene) pairs for {context['species']}; "
          f"{marked} panel rows marked; candidates by delivery tier {dict(sorted(per_tier.items()))}")
    return pool, sources


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
    markers_all_path = os.path.join(CACHE, f"{tag}_markers_all.pkl")
    markers_all = None
    if os.path.exists(markers_all_path):
        with open(markers_all_path, "rb") as handle:
            markers_all = pickle.load(handle)

    api_key, api_url = llm.resolve_api()
    if not enabled:
        llm_status = "disabled"
    elif not api_key or not api_url:
        llm_status = "skipped_no_credentials"
    else:
        llm_status = "enabled"
    borrow_screen_record: dict[str, Any] = {"status": "not_run"}
    borrow_screen = (
        borrow_tissue_screen(context, api_key, api_url, borrow_screen_record)
        if llm_status == "enabled" else None
    )
    pool, sources = prepare_pool(scoring, prep, markers_all, borrow_screen=borrow_screen)
    pool["borrow_screen"] = borrow_screen_record
    scored = scoring["scored"]
    clusters = sorted(scoring["clusters"], key=lambda value: (len(value), value))

    print(
        f"==== annotate(py) {tag}: {len(clusters)} clusters, up to "
        f"{scoring['top_candidates']} candidates each in {pool_api.CANDIDATE_TIERS} tiers "
        f"of {pool_api.CANDIDATE_TIER_SIZE}, llm={llm_status} ====\n"
        f"  annotator: {llm.resolve_model(ANNOTATOR_MODEL)} ({ANNOTATOR_EFFORT}); "
        f"review: {llm.resolve_model(REVIEW_MODEL)} ({REVIEW_EFFORT}), "
        f"{REVIEW_MAX_ROUNDS_PER_TIER} rounds per tier; "
        f"arbiter: {llm.resolve_model(ARBITER_MODEL)} ({ARBITER_EFFORT}); "
        f"borrow screen: {llm.resolve_model(BORROW_SCREEN_MODEL)} ({BORROW_SCREEN_EFFORT}), "
        f"{borrow_screen_record.get('status')}"
    )

    subset_raw = os.environ.get("SCMA_CLUSTER_SUBSET", "").strip()
    subset = {s.strip() for s in subset_raw.split(",") if s.strip()} if subset_raw else None

    results: dict[str, dict[str, Any]] = {}
    runnable: list[str] = []
    for cluster in clusters:
        frame = scored[scored["cluster"] == cluster].reset_index(drop=True)
        if scoring["clusters"][cluster]["status"] == UNSUPPORTED:
            results[cluster] = _fallback(
                cluster, frame.head(0), llm_status, "empty_candidate_pool"
            )
            continue
        if subset is not None and str(cluster) not in subset:
            results[cluster] = _fallback(
                cluster, frame, "not_run_in_subset", "subset_skipped"
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
        server = SourceServer(
            sources, batch=SOURCES_PER_MARKER, max_batches=SOURCE_BATCHES_PER_MARKER
        )
        pool_api.register_packet_sources(server, pool, cluster)
        return cluster, annotate_cluster(
            pool, cluster, template, api_key, api_url, server
        )

    outcomes: dict[str, dict] = {}
    if runnable:
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
            record = _fallback(cluster, frame, f"failed:{error}", "annotation_failed")
            record["candidates"] = pool_api.candidate_names(pool, cluster)
            record["candidate_entries"] = [
                _candidate_entry(entry, cluster, percentile, "", None, "")
                for entry in pool["clusters"][cluster]["candidates"]
            ]
            record["turns"] = int(outcome.get("turns", 0))
            results[cluster] = record
            continue
        results[cluster] = _result(cluster, pool, percentile, outcome, frame)

    consistency: list[dict] = []
    reopened: dict[str, dict] = {}
    if runnable and llm_status == "enabled" and subset is not None:
        print(f"  dataset consistency skipped: SCMA_CLUSTER_SUBSET annotates "
              f"{len(runnable)} of {len(results)} clusters")
    if runnable and llm_status == "enabled" and subset is None:
        consistency = stages.dataset_consistency(
            pool, results, _prompt("dataset_consistency"),
            f"# Dataset context\nspecies: {context['species']} | tissue: {context['tissue']} | "
            f"disease: {context['disease']}\n\n",
            api_key, api_url, uuid.uuid4().hex[:16],
        )
        to_reopen: dict[str, dict] = {}
        for item in consistency:
            reading = item.get("reading") or {}
            for member in reading.get("members") or []:
                if member.get("status") in ("weaker_member", "conflict"):
                    cid = str(member.get("cluster_id"))
                    if cid in results and cid in runnable:
                        to_reopen[cid] = {
                            "status": member["status"],
                            "why_grouped": item["group"]["why_grouped"],
                            "note": str(member.get("note_for_review") or ""),
                            "group_reason": str(reading.get("reason") or ""),
                        }
        if to_reopen:
            print(f"  dataset consistency: {len(consistency)} group reading(s); "
                  f"re-opening {sorted(to_reopen)}")

            def reopen_one(cluster: str):
                server = SourceServer(
                    sources, batch=SOURCES_PER_MARKER, max_batches=SOURCE_BATCHES_PER_MARKER
                )
                pool_api.register_packet_sources(server, pool, cluster)
                return cluster, annotate_cluster(
                    pool, cluster, template, api_key, api_url, server,
                    consistency_note=to_reopen[cluster],
                )

            with ThreadPoolExecutor(max_workers=min(LLM_SETTINGS.threads, len(to_reopen))) as executor:
                for cluster, outcome in executor.map(reopen_one, sorted(to_reopen)):
                    if "final" not in outcome:
                        continue
                    frame = scored[scored["cluster"] == cluster].reset_index(drop=True)
                    before = results[cluster]
                    after = _result(cluster, pool, percentile, outcome, frame)
                    after["reopened_by_consistency"] = {
                        **to_reopen[cluster],
                        "before": {"annotation": before.get("annotation"),
                                   "confidence": before.get("confidence"),
                                   "co_occurring_identities": before.get("co_occurring_identities")},
                    }
                    results[cluster] = after
                    reopened[cluster] = after["reopened_by_consistency"]
                    outcomes[cluster] = outcome

    payload = {
        "dataset_consistency": consistency,
        "reopened_by_consistency": reopened,
        "tag": tag,
        "context": context,
        "llm_status": llm_status,
        "models": {
            "annotator": llm.resolve_model(ANNOTATOR_MODEL),
            "annotator_reasoning_effort": ANNOTATOR_EFFORT,
            "review": llm.resolve_model(REVIEW_MODEL),
            "review_reasoning_effort": REVIEW_EFFORT,
            "arbiter": llm.resolve_model(ARBITER_MODEL),
            "arbiter_reasoning_effort": ARBITER_EFFORT,
            "borrow_screen": llm.resolve_model(BORROW_SCREEN_MODEL),
            "borrow_screen_reasoning_effort": BORROW_SCREEN_EFFORT,
        },
        "borrow_screen": borrow_screen_record,
        "results": results,
        "transcripts": {
            cluster: outcome.get("transcript", []) for cluster, outcome in outcomes.items()
        },
        "tissue_contexts": list(pool.get("tissue_contexts") or []),
        "candidate_tissue_contexts": pool.get("candidate_tissue_contexts") or {},
        "candidate_retrieval_context": pool.get("candidate_retrieval_context") or {},
        "retrieval_semantics": str(pool.get("retrieval_semantics") or ""),
        "candidate_tier_size": pool_api.CANDIDATE_TIER_SIZE,
        "candidate_tiers": pool_api.CANDIDATE_TIERS,
    }
    with open(os.path.join(CACHE, f"{tag}_annotations.pkl"), "wb") as handle:
        pickle.dump(payload, handle, protocol=4)
    llm.flush_response_cache()

    counts = {status: 0 for status in (RESOLVED, MIXED, UNRESOLVED)}
    for value in results.values():
        counts[value["resolution_status"]] += 1
    turns = [value.get("turns", 0) for value in results.values() if value.get("turns")]
    warned = sum(1 for value in results.values() if value.get("claim_warnings"))
    qc = {status: 0 for status in QC_VALUES}
    for value in results.values():
        qc[str(value.get("annotation_qc") or QC_UNCHECKED)] += 1
    tiers = {}
    for value in results.values():
        key = int(value.get("delivered_tier") or 0)
        if key:
            tiers[key] = tiers.get(key, 0) + 1
    print(
        f"[done] {tag}: {counts[RESOLVED]} resolved, {counts[MIXED]} mixed, "
        f"{counts[UNRESOLVED]} unresolved of {len(results)} clusters; "
        f"turns median {int(np.median(turns)) if turns else 0}, max {max(turns, default=0)}; "
        f"evidence review {dict(sorted(qc.items()))}; delivered from tier {dict(sorted(tiers.items()))}; "
        f"{warned} clusters carry a claim warning; "
        f"{sum(1 for v in results.values() if v['annotation_source'] != 'cluster_annotation')} "
        "not model-decided"
    )
    return results


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit("usage: python -m scmarkeragent.cluster_annotation <tag>")
    annotate(sys.argv[1])
