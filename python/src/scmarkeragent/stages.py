#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import json
import re
from typing import Any

from . import annotator_pool as pool_api
from . import llm_client as llm
from .config import DEFAULTS

_A = DEFAULTS["cluster_annotation"]
STAGE_MODEL = str(_A.get("stage_model") or _A["review_model"])
STAGE_EFFORT = str(_A.get("stage_reasoning_effort") or _A["review_reasoning_effort"])
SCHEMA_RETRIES = int(_A["schema_retries"])
AUDIT_PANEL_ROWS = int(_A.get("audit_panel_rows", 15))
AUDIT_SENTENCES_PER_GENE = int(_A.get("audit_sentences_per_gene", 2))
AUDIT_RIVALS_MAX = int(_A.get("audit_rivals_max", 4))
AUDIT_RIVAL_PANEL_ROWS = int(_A.get("audit_rival_panel_rows", 10))
AUDIT_GENOME_ROWS = int(_A.get("audit_genome_rows", 15))
CONSISTENCY_GROUP_MAX = int(_A.get("consistency_group_max", 8))

OWN_PROGRAM_VALUES = ("dominant", "minority", "absent", "single_study_only")
CONSISTENCY_VALUES = ("consistent", "weaker_member", "conflict")
LINEAGE_SPECIFICITY_VALUES = ("one_lineage", "several_lineages", "uncertain")

RELATION_VALUES = ("coarser_parent", "finer_child", "sibling_or_pole", "different_lineage",
                   "state_label")
OCCUPANCY_VALUES = ("identity_only", "rival_only", "both_cluster_level", "both_minority",
                    "neither")
POLE_RELATIONS = ("sibling_or_pole", "different_lineage")
CONTESTED_OCCUPANCY = ("rival_only", "both_cluster_level")


def _render(template: str, packet: dict) -> str:
    return template.replace(
        "{{EVIDENCE_PACKET_JSON}}",
        json.dumps(packet, ensure_ascii=False, sort_keys=True, separators=(",", ":")),
    )


def _call(prompt: str, api_key: str, api_url: str, trace_id: str, turn: int,
          validate, retry_note: str) -> dict | None:
    for attempt in range(SCHEMA_RETRIES + 1):
        text, _error, _replayed = llm.cached_call_llm(
            prompt, api_url, api_key, reasoning_effort=STAGE_EFFORT, model=STAGE_MODEL,
            trace_id=trace_id, turn_index=turn,
        )
        parsed = llm.parse_json(text) if text else None
        problems = validate(parsed)
        if not problems:
            return parsed
        if attempt >= SCHEMA_RETRIES:
            return None
        prompt += f"\n\n# Retry {attempt + 1}: {retry_note} -- {'; '.join(problems[:4])}.\n"
    return None


def compact_definers(rows: list[dict]) -> list[str]:
    return pool_api.compact_definer_rows(rows or [])


def program_packet(pool: dict, cluster: str, tier: int) -> dict:
    state = pool["clusters"][str(cluster)]
    return {
        "query": {**pool["context"], "cluster_id": str(cluster),
                  "cells_in_cluster": int(state["n_cells"])},
        **pool_api.genome_pages(pool, cluster, tier),
    }


def _page_genes(packet: dict) -> set[str]:
    genes = set()
    for key in ("cluster_top_enriched", "cluster_top_detection_gap", "cluster_top_depleted"):
        genes.update(str(r["gene"]).upper() for r in packet.get(key) or [])
    return genes


def _mentions_candidate(text: str, names: list[str]) -> str | None:
    low = str(text or "").lower()
    for name in names:
        n = name.lower().strip()
        if len(n) >= 6 and n in low:
            return name
    return None


def validate_program(value: Any, page_genes: set[str], names: list[str]) -> list[str]:
    if not isinstance(value, dict):
        return ["return one JSON object"]
    problems = []
    programs = value.get("programs")
    if not isinstance(programs, list):
        problems.append("programs must be a list")
        programs = []
    for block in programs:
        if not isinstance(block, dict):
            problems.append("each program must be an object")
            continue
        for g in block.get("genes") or []:
            if str(g).upper() not in page_genes:
                problems.append(f"gene {g} is not on the genome-wide pages")
                break
        hit = _mentions_candidate(block.get("function", ""), names)
        if hit:
            problems.append(f"'function' names a candidate cell type ({hit}); describe the function")
        if block.get("lineage_specificity") not in LINEAGE_SPECIFICITY_VALUES:
            problems.append(
                f"each program needs lineage_specificity, one of {list(LINEAGE_SPECIFICITY_VALUES)}"
            )
        hit = _mentions_candidate(block.get("lineage_note", ""), names)
        if hit:
            problems.append(f"'lineage_note' names a candidate cell type ({hit})")
    hit = _mentions_candidate(value.get("mixture_note", ""), names)
    if hit:
        problems.append(f"mixture_note names a candidate cell type ({hit})")
    if not isinstance(value.get("lineage_programs_count"), int):
        problems.append("lineage_programs_count must be an integer")
    return problems


def program_reading(pool: dict, cluster: str, tier: int, template: str, header: str,
                    api_key: str, api_url: str, trace_id: str, turn: int) -> dict | None:
    packet = program_packet(pool, cluster, tier)
    page_genes = _page_genes(packet)
    names = pool_api.candidate_names(pool, cluster)
    prompt = header + _render(template, packet)
    return _call(prompt, api_key, api_url, trace_id, turn,
                 lambda v: validate_program(v, page_genes, names),
                 "the reading was rejected")


def _sentences_for(sources, identity: str, genes: list[str]) -> dict[str, list[dict]]:
    out: dict[str, list[dict]] = {}
    if sources is None:
        return out
    for gene in genes:
        key = str(gene).upper()
        if key in out:
            continue
        records = sources.records_for_marker(identity, key, k=AUDIT_SENTENCES_PER_GENE)
        if records:
            out[key] = [pool_api._sentence_record(r) for r in records]
    return out


def _rival_names(pool: dict, cluster: str, identity: str, claimed: list[str], tier: int) -> list[str]:
    rivals = [n for n in claimed if n != identity]
    genome = (pool.get("genome_de") or {}).get(str(cluster)) or {}
    on_page = set(pool_api.candidate_names(pool, cluster, tier))
    votes: dict[str, int] = {}
    for row in (genome.get("top_enriched") or [])[:AUDIT_GENOME_ROWS]:
        for c in row.get("claimed_by_species") or []:
            name = c["cell_type"]
            if name in on_page and name != identity and name not in rivals:
                votes[name] = votes.get(name, 0) + 1
    for name, _n in sorted(votes.items(), key=lambda kv: (-kv[1], kv[0])):
        if len(rivals) >= AUDIT_RIVALS_MAX:
            break
        rivals.append(name)
    return rivals[:AUDIT_RIVALS_MAX]


def audit_packet(pool: dict, cluster: str, identity: str, role: str, claimed: list[str],
                 tier: int, sources, program: dict | None) -> dict | None:
    entry = pool_api.find_candidate(pool, cluster, identity)
    if entry is None:
        return None
    state = pool["clusters"][str(cluster)]
    definers = list(entry.get("definers") or [])
    panel = pool_api._review_panel(entry)
    positives = [r for r in panel if r.get("polarity") == "positive"][:AUDIT_PANEL_ROWS]
    negatives = [r for r in panel if r.get("polarity") == "negative"]
    species_sources = pool.get("species_sources") or sources
    sentences = _sentences_for(species_sources, identity, [str(d["gene"]) for d in definers])
    rivals = []
    for name in _rival_names(pool, cluster, identity, claimed, tier):
        rentry = pool_api.find_candidate(pool, cluster, name)
        if rentry is None:
            continue
        rdefiners = list(rentry.get("definers") or [])
        rpanel = [r for r in pool_api._review_panel(rentry) if r.get("polarity") == "positive"]
        rival = {
            "cell_type": name,
            "claimed_here": name in claimed,
            "species_definers": compact_definers(rdefiners),
            "panel": pool_api.compact_panel(rpanel[:AUDIT_RIVAL_PANEL_ROWS]),
            "sentences": _sentences_for(species_sources, name,
                                        [str(d["gene"]) for d in rdefiners]),
            "versus_identity": pool_api.pair_partition(entry, rentry),
        }
        if rentry.get("borrowed_context"):
            rival["borrowed_context"] = rentry["borrowed_context"]
        rivals.append(rival)
    genome = pool_api.genome_pages(pool, cluster, tier)
    packet = {
        "query": {**pool["context"], "cluster_id": str(cluster),
                  "cells_in_cluster": int(state["n_cells"])},
        "identity": identity,
        "role": role,
        "row_format": pool_api.MARKER_ROW_FORMAT,
        "species_definers": definers,
        "own_vs_shared": pool_api.lineage_partition(entry),
        "panel": pool_api.compact_panel(positives + negatives),
        "sentences": sentences,
        "detected_exclusions": pool_api.compact_panel(pool_api._detected_exclusion_markers(entry)),
        "rivals": rivals,
        "program_reading": program or {},
        "program_genes_that_are_not_identity": program_non_identity_genes(program),
        "cluster_top_enriched": genome["cluster_top_enriched"][:AUDIT_GENOME_ROWS],
    }
    if entry.get("borrowed_context"):
        packet["borrowed_context"] = entry["borrowed_context"]
    return packet


def program_non_identity_genes(program: dict | None) -> dict:
    out = {"state": [], "several_lineages": []}
    if not isinstance(program, dict):
        return out
    seen: set[str] = set()
    for block in program.get("state_programs") or []:
        for g in (block or {}).get("genes") or []:
            if str(g).upper() not in seen:
                seen.add(str(g).upper())
                out["state"].append(str(g))
    for block in program.get("programs") or []:
        if str((block or {}).get("lineage_specificity") or "") != "several_lineages":
            continue
        for g in block.get("genes") or []:
            if str(g).upper() not in seen:
                seen.add(str(g).upper())
                out["several_lineages"].append(str(g))
    return out


_GENE_LIST_KEYS = {
    "genes", "state", "several_lineages", "background_or_shared", "unplaced_top_genes",
}
_COMPACT_FIRST = re.compile(r"^(\S+) n(?:\d+|-)")


def _genes_in(obj: Any, out: set[str], key: str | None = None) -> None:
    if isinstance(obj, dict):
        gene = obj.get("gene")
        if isinstance(gene, str):
            out.add(gene.upper())
        if key in ("sentences", "sources", "exclusion_sources"):
            out.update(str(k).upper() for k in obj.keys())
        for k, v in obj.items():
            _genes_in(v, out, str(k))
        return
    if isinstance(obj, list):
        for item in obj:
            if isinstance(item, str) and key in _GENE_LIST_KEYS and " " not in item:
                out.add(item.upper())
            else:
                _genes_in(item, out, key)
        return
    if isinstance(obj, str):
        m = _COMPACT_FIRST.match(obj)
        if m:
            out.add(m.group(1).upper())


def packet_genes(packet: dict) -> set[str]:
    out: set[str] = set()
    _genes_in(packet, out)
    return out


def identity_side_genes(packet: dict) -> set[str]:
    out: set[str] = set()
    identity = str(packet.get("identity") or "")
    _genes_in({k: v for k, v in packet.items() if k != "rivals"}, out)
    for rival in packet.get("rivals") or []:
        versus = rival.get("versus_identity") or {}
        _genes_in(versus.get(f"only_{identity}"), out, "genes")
        _genes_in(versus.get("both"), out, "genes")
    return out


def validate_audit(value: Any, packet: dict) -> list[str]:
    if not isinstance(value, dict):
        return ["return one JSON object"]
    problems = []
    if value.get("own_program") not in OWN_PROGRAM_VALUES:
        problems.append(f"own_program must be one of {list(OWN_PROGRAM_VALUES)}")
    known = identity_side_genes(packet)
    rival_known = packet_genes(packet)
    for key in ("definers_present", "definers_absent"):
        genes = value.get(key)
        if not isinstance(genes, list):
            problems.append(f"{key} must be a list")
            continue
        bad = [g for g in genes if str(g).upper() not in known]
        if bad:
            problems.append(f"{key} names genes not in this packet: {bad[:4]}")
    rival_names = [str(r["cell_type"]) for r in packet.get("rivals") or []]
    read = value.get("rivals_read")
    if not isinstance(read, list):
        problems.append("rivals_read must be a list with one entry per rival")
        read = []
    seen = set()
    for entry in read:
        if not isinstance(entry, dict):
            problems.append("each rivals_read entry must be an object")
            continue
        name = str(entry.get("rival") or "")
        if name not in rival_names:
            problems.append(f"rivals_read names '{name}', which is not a rival in this packet")
            continue
        seen.add(name)
        if entry.get("relation") not in RELATION_VALUES:
            problems.append(f"{name}: relation must be one of {list(RELATION_VALUES)}")
        if entry.get("whose_program_occupies") not in OCCUPANCY_VALUES:
            problems.append(f"{name}: whose_program_occupies must be one of {list(OCCUPANCY_VALUES)}")
        for key in ("separating_genes", "rival_separating_genes"):
            genes = entry.get(key)
            if not isinstance(genes, list):
                problems.append(f"{name}: {key} must be a list")
                continue
            bad = [g for g in genes if str(g).upper() not in rival_known]
            if bad:
                problems.append(f"{name}: {key} names genes not in this packet: {bad[:4]}")
        if not str(entry.get("note") or "").strip():
            problems.append(f"{name}: note must quote the measured percentages of both sides")
    missing = [n for n in rival_names if n not in seen]
    if missing:
        problems.append(f"every rival in this packet needs one rivals_read entry; missing {missing[:4]}")
    if not str(value.get("reason") or "").strip():
        problems.append("reason must quote the measured percentages")
    return problems


def definer_audit(pool: dict, cluster: str, identity: str, role: str, claimed: list[str],
                  tier: int, sources, program: dict | None, template: str, header: str,
                  api_key: str, api_url: str, trace_id: str, turn: int) -> dict | None:
    packet = audit_packet(pool, cluster, identity, role, claimed, tier, sources, program)
    if packet is None:
        return None
    prompt = header + _render(template, packet)
    out = _call(prompt, api_key, api_url, trace_id, turn,
                lambda v: validate_audit(v, packet), "the audit was rejected")
    if out is not None:
        out["identity"] = identity
        out["role"] = role
        out["species_definers"] = compact_definers(packet["species_definers"])
    return out


def cached_audit(cache: dict | None, tier: int, identity: str, need_rivals: list[str]) -> dict | None:
    if cache is None:
        return None
    audit = cache.get((tier, identity))
    if audit is None:
        return None
    read = {str(e.get("rival") or "") for e in rivals_read(audit)}
    if any(name not in read for name in need_rivals if name != identity):
        return None
    return audit


def audit_all(pool: dict, cluster: str, final: dict, tier: int, sources, program: dict | None,
              template: str, header: str, api_key: str, api_url: str, trace_id: str,
              turn: int, extra: dict[str, str] | None = None,
              cache: dict | None = None) -> dict[str, dict]:
    from .cluster_annotation import claimed_identities

    claims = list(claimed_identities(final))
    claimed_names = {n for _r, n in claims}
    for name, role in (extra or {}).items():
        if name and name not in claimed_names:
            claims.append((role, name))
            claimed_names.add(name)
    names = [n for _r, n in claims]
    out: dict[str, dict] = {}
    for index, (role, name) in enumerate(claims):
        hit = cached_audit(cache, tier, name, _rival_names(pool, cluster, name, names, tier))
        if hit is not None:
            out[name] = {**hit, "role": role}
            continue
        audit = definer_audit(pool, cluster, name, role, names, tier, sources, program,
                              template, header, api_key, api_url, trace_id, turn * 100 + index)
        if audit is not None:
            out[name] = audit
            if cache is not None:
                cache[(tier, name)] = audit
    return out


def rivals_read(audit: dict | None) -> list[dict]:
    return [e for e in ((audit or {}).get("rivals_read") or []) if isinstance(e, dict)]


def contested_rival_names(final: dict, audits: dict[str, dict]) -> list[str]:
    from .cluster_annotation import claimed_identities

    selected = str(final.get("selected") or "")
    claimed = {name for _role, name in claimed_identities(final)}
    out = []
    for entry in rivals_read(audits.get(selected)):
        name = str(entry.get("rival") or "")
        if not name or name in claimed or name == selected:
            continue
        if str(entry.get("relation")) in POLE_RELATIONS and name not in out:
            out.append(name)
    return out


def _pair_reading(audit: dict | None, other: str) -> dict | None:
    for entry in rivals_read(audit):
        if str(entry.get("rival") or "") == other:
            return entry
    return None


def delivery_contested(final: dict, audits: dict[str, dict]) -> list[str]:
    selected = str(final.get("selected") or "")
    audit = audits.get(selected)
    if not audit:
        return []
    reasons: list[str] = []
    own = str(audit.get("own_program") or "")
    if own and own != "dominant":
        reasons.append(
            f"the definer audit of '{selected}' reads its own defining program as {own} here"
            + (f" ({', '.join(str(g) for g in (audit.get('definers_absent') or []))[:160]})"
               if audit.get("definers_absent") else "")
        )
    for name in contested_rival_names(final, audits):
        entry = _pair_reading(audit, name)
        if entry is None:
            continue
        relation = str(entry.get("relation") or "")
        occupancy = str(entry.get("whose_program_occupies") or "")
        if occupancy == "rival_only":
            reasons.append(
                f"its own audit reads the pair '{selected}' / '{name}' ({relation}) as "
                f"{occupancy} here: {str(entry.get('note') or '')[:200]}"
            )
            continue
        if occupancy == "neither":
            reasons.append(
                f"its own audit reads the pair '{selected}' / '{name}' ({relation}) as "
                f"neither: no gene in the packet tells them apart, so the delivered pole is "
                f"not established against '{name}' and the level the rows reach is their "
                f"common parent where this page holds it: {str(entry.get('note') or '')[:200]}"
            )
            continue
    return reasons


POPULATION_VERDICTS = (
    "established_second_population",
    "shared_program_of_same_population",
    "neither_separates",
    "not_present",
)
ESTABLISHED = "established_second_population"
NEITHER_SEPARATES = "neither_separates"


def population_names(packet: dict) -> list[str]:
    delivered = packet.get("delivered") or {}
    out: list[str] = []
    for name in delivered.get("co_occurring_identities") or []:
        if name and name not in out:
            out.append(str(name))
    for block in packet.get("contested_identities") or []:
        name = str(block.get("cell_type") or "")
        if name and name not in out:
            out.append(name)
    return out


def population_conflicts(final: dict, review: dict | None) -> list[str]:
    if not review:
        return []
    claimed = {str(n) for n in (final.get("co_occurring_identities") or []) if n}
    out: list[str] = []
    for entry in review.get("population_verdicts") or []:
        if not isinstance(entry, dict):
            continue
        name = str(entry.get("identity") or "")
        verdict = str(entry.get("verdict") or "")
        why = str(entry.get("reason") or "")[:220]
        if name in claimed and verdict != ESTABLISHED:
            out.append(
                f"'{name}' is reported as a co-occurring identity, and the review reads it "
                f"as {verdict}: {why}"
            )
        elif name not in claimed and verdict == ESTABLISHED:
            out.append(
                f"'{name}' is not reported, and the review reads it as an established "
                f"second population here: {why}"
            )
    return out


def confidence_cap(final: dict, audits: dict[str, dict]) -> tuple[str | None, str]:
    selected = str(final.get("selected") or "")
    audit = audits.get(selected)
    if not audit:
        return None, ""
    if audit.get("own_program") == "single_study_only":
        return "low", f"'{selected}' is carried only by single-study genes in its definer audit"
    if audit.get("own_program") != "dominant":
        return "medium", f"the definer audit reads '{selected}' as {audit.get('own_program')}, not dominant"
    for entry in rivals_read(audit):
        if (str(entry.get("relation")) in ("coarser_parent", "finer_child")
                and str(entry.get("whose_program_occupies")) == "rival_only"):
            return "medium", (
                f"the audit reads the distinction '{selected}' draws against "
                f"'{entry.get('rival')}' as carried by '{entry.get('rival')}' here"
            )
    return None, ""


CONF_ORDER = {"high": 3, "medium": 2, "low": 1}


def apply_cap(final: dict, cap: str | None) -> bool:
    if cap is None:
        return False
    current = str(final.get("confidence") or "low")
    if CONF_ORDER.get(current, 1) > CONF_ORDER[cap]:
        final["confidence"] = cap
        return True
    return False


def _batched(clusters: list[str], size: int) -> list[list[str]]:
    if len(clusters) <= size:
        return [clusters]
    return [clusters[i:i + size] for i in range(0, len(clusters), size)]


def consistency_groups(results: dict[str, dict]) -> list[dict]:
    unknown = str(_A["unknown_token"])
    by_name: dict[str, list[str]] = {}
    for cluster, record in results.items():
        name = str(record.get("annotation") or "")
        if not name or name == unknown:
            continue
        by_name.setdefault(name, []).append(str(cluster))

    groups: list[dict] = []
    for name, clusters in sorted(by_name.items()):
        if len(clusters) < 2:
            continue
        batches = _batched(sorted(clusters, key=lambda c: (len(c), c)), CONSISTENCY_GROUP_MAX)
        for index, batch in enumerate(batches):
            suffix = f" (batch {index + 1} of {len(batches)})" if len(batches) > 1 else ""
            groups.append({"why_grouped": f"same delivered name: {name}{suffix}",
                           "clusters": sorted(batch, key=lambda c: (len(c), c))})

    seen_pairs: set[tuple[str, ...]] = set()
    for cluster, record in sorted(results.items(), key=lambda kv: (len(kv[0]), kv[0])):
        name = str(record.get("annotation") or "")
        audit = (record.get("definer_audit") or {}).get(name)
        for entry in rivals_read(audit):
            rival = str(entry.get("rival") or "")
            if (str(entry.get("relation")) not in POLE_RELATIONS
                    or str(entry.get("whose_program_occupies")) not in CONTESTED_OCCUPANCY):
                continue
            others = sorted(
                (c for c, r in results.items() if str(r.get("annotation") or "") == rival),
                key=lambda c: (len(c), c),
            )
            if not others:
                continue
            members = sorted({str(cluster), *others}, key=lambda c: (len(c), c))
            key = tuple(members)
            if len(members) < 2 or key in seen_pairs:
                continue
            seen_pairs.add(key)
            for batch in _batched(members, CONSISTENCY_GROUP_MAX):
                if len(batch) >= 2:
                    groups.append({
                        "why_grouped": (f"contested pair: c{cluster} was delivered {name!r} while "
                                        f"its own audit read {rival!r} as "
                                        f"{entry.get('whose_program_occupies')} here"),
                        "clusters": batch,
                    })
    return groups


def consistency_packet(pool: dict, group: dict, results: dict[str, dict]) -> dict:
    names = []
    for cluster in group["clusters"]:
        name = str((results.get(cluster) or {}).get("annotation") or "")
        if name and name not in names:
            names.append(name)
    members = []
    for cluster in group["clusters"]:
        record = results[cluster]
        name = str(record.get("annotation") or "")
        genome = (pool.get("genome_de") or {}).get(str(cluster)) or {}
        audits = record.get("definer_audit") or {}
        definers_here = {}
        for other in names:
            entry = pool_api.find_candidate(pool, cluster, other)
            if entry is not None:
                definers_here[other] = compact_definers(entry.get("definers") or [])
        members.append({
            "cluster_id": str(cluster),
            "n_cells": int(pool["clusters"][str(cluster)]["n_cells"]),
            "delivered": {
                "selected": name,
                "co_occurring_identities": list(record.get("co_occurring_identities") or []),
                "confidence": str(record.get("confidence") or ""),
            },
            "species_definers_of_each_delivered_name": definers_here,
            "definer_audit": {k: v for k, v in (audits.get(name) or {}).items()
                              if k in ("own_program", "definers_present", "definers_absent",
                                       "rivals_read")},
            "top_enriched": [f"{r['gene']} {r['pct_in']}/{r['pct_out']}"
                             for r in (genome.get("top_enriched") or [])[:12]],
        })
    return {"query": dict(pool["context"]),
            "group": {"why_grouped": group["why_grouped"], "delivered_names": names},
            "members": members}


def validate_consistency(value: Any, clusters: list[str]) -> list[str]:
    if not isinstance(value, dict):
        return ["return one JSON object"]
    problems = []
    members = value.get("members")
    if not isinstance(members, list):
        return ["members must be a list"]
    seen = set()
    for m in members:
        if not isinstance(m, dict):
            problems.append("each member must be an object")
            continue
        cid = str(m.get("cluster_id") or "")
        if cid not in clusters:
            problems.append(f"cluster_id '{cid}' is not in this group")
        seen.add(cid)
        if m.get("status") not in CONSISTENCY_VALUES:
            problems.append(f"status must be one of {list(CONSISTENCY_VALUES)}")
        if m.get("status") != "consistent" and not str(m.get("note_for_review") or "").strip():
            problems.append(f"cluster {cid}: a marked member needs note_for_review")
    missing = [c for c in clusters if c not in seen]
    if missing:
        problems.append(f"every member needs a status; missing {missing[:4]}")
    return problems


def dataset_consistency(pool: dict, results: dict[str, dict], template: str, header: str,
                        api_key: str, api_url: str, trace_id: str) -> list[dict]:
    readings = []
    groups = consistency_groups(results)
    for index, group in enumerate(groups):
        packet = consistency_packet(pool, group, results)
        prompt = header + _render(template, packet)
        out = _call(prompt, api_key, api_url, trace_id, 900 + index,
                    lambda v, cl=group["clusters"]: validate_consistency(v, cl),
                    "the reading was rejected")
        readings.append({"group": group, "reading": out})
    return readings
