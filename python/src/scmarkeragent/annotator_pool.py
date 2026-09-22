#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import os
import re
from typing import Any

import numpy as np
import pandas as pd

from . import evidence_gate
from .config import DEFAULTS, NOT_AVAILABLE

_A = DEFAULTS["cluster_annotation"]
SOURCES_PER_GENE = int(_A["sources_per_marker"])
SOURCE_BATCHES_PER_MARKER = int(_A["source_batches_per_marker"])
SOURCE_GENES_PER_QUERY = int(_A["source_genes_per_query"])
SOURCE_MAX_CHARS = int(_A["source_sentence_max_chars"])
UNMEASURED_SHOWN = int(_A["unmeasured_genes_shown"])
WARNING_GENES_SHOWN = int(_A["warning_genes_shown"])
CANDIDATE_TIER_SIZE = int(_A["candidate_tier_size"])
CANDIDATE_TIERS = int(_A["candidate_tiers"])
DECISIVE_MIN_PUB = int(_A["decisive_min_publications"])
SINGLE_STUDY_PUB = int(_A["single_study_publications"])
SPARSE_LABEL_MAX_PAPERS = int(_A["sparse_label_max_papers"])
EXCLUSION_MIN_PUB = int(_A["exclusion_min_publications"])
DM_GATE_TOP = int(_A["defining_marker_gate_top"])
DM_GATE_FLOOR_PCT = float(_A["defining_marker_gate_floor_pct"])
DM_GATE_DOMINANT_PCT = float(_A["defining_marker_gate_dominant_pct"])
GENOME_ENRICHED_ROWS = int(_A["evidence_judge_enriched_rows"])
GENOME_DEPLETED_ROWS = int(_A["evidence_judge_depleted_rows"])
REVIEW_OTHER_CANDIDATE_ROWS = int(_A["review_other_candidate_rows"])
REVIEW_SOURCES_PER_MARKER = int(_A["review_sources_per_marker"])

_LABEL_PAPERS: dict[tuple[str, str], int] | None = None


def label_papers(species: str, label: str) -> int | None:
    global _LABEL_PAPERS
    if _LABEL_PAPERS is None:
        import csv as _csv
        from pathlib import Path as _Path

        path = os.environ.get("SCMA_LABEL_PAPERS", "") or str(
            _Path(__file__).resolve().parent / "resources" / "label_papers_v4.csv"
        )
        table: dict[tuple[str, str], int] = {}
        if os.path.exists(path):
            with open(path, newline="", encoding="utf-8") as handle:
                for row in _csv.DictReader(handle):
                    table[(str(row["species"]), str(row["cell_type"]))] = int(row["n_papers"])
        _LABEL_PAPERS = table
    return _LABEL_PAPERS.get((str(species), str(label)))
NEGATIVE_SOURCE_MIN_PCT_IN = float(_A["negative_source_min_pct_in"])
CLUSTER_MARKER_ROWS = int(_A["cluster_marker_rows"])

TOOLS = ("sources", "gene_across_clusters", "candidates_with_gene")

REVIEW_TOOLS = ("sources", "gene_across_clusters", "candidates_with_gene")

PCT_DECIMALS = 1
LOGFC_DECIMALS = 3
SPECIFICITY_DECIMALS = 3
AUC_DECIMALS = 3

PACKET_MARKER_FIELDS = (
    "gene",
    "polarity",
    "n_pub",
    "recommended",
    "specificity",
    "pct_in",
    "pct_out",
    "avg_log2FC",
    "auc",
    "significant",
    "shared_on_page",
)
SPARSE_LABEL_PAPERS_MAX = SPARSE_LABEL_MAX_PAPERS


def _num(value, decimals):
    if value is None:
        return None
    number = float(value)
    return None if not np.isfinite(number) else round(number, decimals)


def _clip(text: str, limit: int) -> str:
    value = " ".join(str(text or "").split())
    return value if len(value) <= limit else value[: limit - 1] + "\u2026"


def _identifier(value) -> str:
    text = str(value or "").strip()
    return text if text and text.lower() not in {"nan", "none", "na"} else NOT_AVAILABLE


def _sentence_record(record: dict) -> dict:
    n_pub = record.get("n_pub")
    n_pub = int(n_pub) if n_pub is not None else None
    out = {
        "pmid": _identifier(record.get("pmid")),
        "pmcid": _identifier(record.get("pmcid")),
        "sentence": _clip(record.get("sentence"), SOURCE_MAX_CHARS),
    }
    if n_pub is not None:
        out["n_pub"] = n_pub
        out["single_study"] = n_pub <= SINGLE_STUDY_PUB
    return out


def is_raised(marker: dict) -> bool:
    value = marker.get("avg_log2FC")
    if value is None or float(value) <= 0.0:
        return False
    auc = marker.get("auc")
    return auc is not None and float(auc) > evidence_gate.SIG_AUC


def _significant(measured: dict) -> bool:
    pct_in = measured.get("pct_in")
    auc = measured.get("auc")
    padj = measured.get("adjusted_p_value")
    if (
        measured.get("avg_log2FC") is None
        or pct_in is None
        or auc is None
        or padj is None
    ):
        return False
    return bool(
        evidence_gate.sig_pass(
            float(measured["avg_log2FC"]),
            float(pct_in) / 100.0,
            float(padj),
            float(auc),
        )
    )


def build_pool(
    scoring: dict, prep: pd.DataFrame | dict, sources=None, markers_all: pd.DataFrame | None = None
) -> dict[str, Any]:
    de = (prep["de"] if isinstance(prep, dict) else prep).copy()
    de["group"] = de["group"].astype(str)
    de["gene_key"] = de["feature"].astype(str).str.upper()

    stats: dict[tuple[str, str], dict[str, float]] = {}
    genome_stats: dict[tuple[str, str], dict[str, float]] = {}
    by_gene: dict[str, list[tuple[str, float, float]]] = {}
    native: dict[str, str] = {}
    for row in de.itertuples(index=False):
        cluster, gene = str(row.group), str(row.gene_key)
        pct_in = _num(row.pct_in, PCT_DECIMALS)
        pct_out = _num(row.pct_out, PCT_DECIMALS)
        stats[(cluster, gene)] = {
            "pct_in": pct_in,
            "pct_out": pct_out,
            "avg_log2FC": _num(row.avg_log2FC, LOGFC_DECIMALS),
            "auc": _num(row.auc, AUC_DECIMALS),
            "adjusted_p_value": _num(row.padj, 6),
        }
        by_gene.setdefault(gene, []).append((cluster, pct_in, pct_out))
        native.setdefault(gene, str(row.feature))
    if markers_all is not None and len(markers_all):
        ga = markers_all
        for row in ga.itertuples(index=False):
            cluster, gene = str(row.group), str(row.feature).upper()
            if (cluster, gene) in stats:
                continue
            pct_in = _num(row.pct_in, PCT_DECIMALS)
            pct_out = _num(row.pct_out, PCT_DECIMALS)
            genome_stats[(cluster, gene)] = {
                "pct_in": pct_in,
                "pct_out": pct_out,
                "avg_log2FC": _num(row.avg_log2FC, LOGFC_DECIMALS),
                "auc": _num(row.auc, AUC_DECIMALS),
                "adjusted_p_value": _num(row.padj, 6),
            }
            by_gene.setdefault(gene, []).append((cluster, pct_in, pct_out))
            native.setdefault(gene, str(row.feature))
    measured_genes = (
        {str(g).upper() for g in prep.get("genes", [])} if isinstance(prep, dict) else set()
    )

    panel_records = scoring["panel_records"]
    meta: dict[tuple[str, str], dict[str, Any]] = {}
    carriers: dict[str, list[tuple[str, int, str, str]]] = {}
    curated: dict[str, set[str]] = {}
    for row in panel_records.itertuples(index=False):
        candidate, gene = str(row.candidate), str(row.gene_key).upper()
        polarity = str(row.marker_polarity)
        n_pub = int(row.n_pub) if pd.notna(row.n_pub) else 0
        tier = str(row.tier) if pd.notna(row.tier) else NOT_AVAILABLE
        meta[(candidate, gene)] = {
            "n_pub": n_pub,
            "tier": tier,
            "polarity": polarity,
        }
        carriers.setdefault(gene, []).append((candidate, n_pub, tier, polarity))
        if polarity == "positive":
            curated.setdefault(candidate, set()).add(gene)

    context = scoring["context"]
    scored = scoring["scored"]
    positive_panels = scoring["measured_panels"]
    negative_panels = scoring["negative_panels"]
    specificity = scoring.get("marker_specificity") or {}
    tissue_contexts = list(scoring.get("tissue_contexts") or [])
    candidate_contexts = scoring.get("candidate_tissue_contexts") or {}
    retrieval_context = scoring.get("candidate_retrieval_context") or {}

    clusters: dict[str, dict[str, Any]] = {}
    for cluster, state in scoring["clusters"].items():
        cluster = str(cluster)
        frame = scored[scored["cluster"] == cluster]
        candidates = []
        for name in state.get("candidates") or []:
            name = str(name)
            match = frame[frame["candidate"] == name]
            if match.empty:
                continue
            row = match.iloc[0]
            markers = _panel_rows(
                name,
                positive_panels.get(name, []),
                negative_panels.get(name, []),
                cluster,
                stats,
                meta,
                native,
                specificity,
            )
            measured = {marker["gene"].upper() for marker in markers}
            unmeasured = sorted(curated.get(name, set()) - measured)
            candidates.append(
                {
                    "cell_type": name,
                    "retrieval_rank": int(row["retrieval_rank"]),
                    "tissue_context": list(candidate_contexts.get(name) or []),
                    "retrieval_context": str(retrieval_context.get(name) or ""),
                    "markers": markers,
                    "exclusion_sources": _exclusion_sources(name, markers, sources),
                    "unmeasured_curated_genes": {
                        "count": len(unmeasured),
                        "genes": [native.get(gene, gene) for gene in unmeasured][
                            :UNMEASURED_SHOWN
                        ],
                    },
                    "program": {
                        "median_in": _num(row["program_in_median"], 4),
                        "median_out": _num(row["program_out_median"], 4),
                    },
                }
            )
        clusters[cluster] = {
            "cluster_id": cluster,
            "status": state.get("status"),
            "n_cells": int(state.get("n_cells", 0)),
            "candidates": candidates,
        }

    return {
        "context": {
            "species": context["species"],
            "tissue": context["tissue"],
            "disease": context["disease"],
            "development_stage": context["development_stage"] or NOT_AVAILABLE,
            "clusters_in_dataset": int(context["n_clusters"]),
        },
        "clusters": clusters,
        "tissue_contexts": tissue_contexts,
        "candidate_tissue_contexts": {str(k): list(v) for k, v in candidate_contexts.items()},
        "candidate_retrieval_context": {str(k): str(v) for k, v in retrieval_context.items()},
        "retrieval_semantics": str(scoring.get("retrieval_semantics") or ""),
        "gene_clusters": by_gene,
        "gene_carriers": carriers,
        "native_gene": native,
        "native_menu": sorted(curated),
        "native_menu_genes": sorted({g for genes in curated.values() for g in genes}),
        "de_stats": stats,
        "genome_stats": genome_stats,
        "measured_genes": measured_genes,
        "marker_specificity": specificity,
    }


def _opening_records(sources, label: str, gene: str, k: int) -> list[dict]:
    if sources is None:
        return []
    if hasattr(sources, "opening"):
        return sources.opening(label, gene, k=k)
    return sources.records_for_marker(label, gene, k=k)


def register_packet_sources(server, pool: dict, cluster: str) -> None:
    if server is None or not hasattr(server, "note_delivered"):
        return
    for entry in pool["clusters"][str(cluster)]["candidates"]:
        candidate = str(entry["cell_type"])
        for gene, records in (entry.get("exclusion_sources") or {}).items():
            server.note_delivered(candidate, gene, records)


def _exclusion_sources(candidate: str, markers: list[dict], sources) -> dict:
    out = {}
    for marker in markers:
        if marker["polarity"] != "negative":
            continue
        pct_in = marker["pct_in"]
        if pct_in is None or pct_in < NEGATIVE_SOURCE_MIN_PCT_IN:
            continue
        gene = str(marker["gene"]).upper()
        records = _opening_records(sources, candidate, gene, SOURCES_PER_GENE)
        out[gene] = [_sentence_record(record) for record in records]
    return out


def _panel_rows(
    candidate: str,
    positive: list[str],
    negative: list[str],
    cluster: str,
    stats: dict,
    meta: dict,
    native: dict,
    specificity: dict,
) -> list[dict]:
    rows = []
    seen: set[str] = set()
    for gene, polarity in [(gene, "positive") for gene in positive] + [
        (gene, "negative") for gene in negative
    ]:
        gene = str(gene).upper()
        if gene in seen:
            continue
        measured = stats.get((cluster, gene))
        if measured is None:
            continue
        curated = meta.get((candidate, gene), {})
        if (
            str(curated.get("polarity") or polarity) == "negative"
            and int(curated.get("n_pub", 0)) < EXCLUSION_MIN_PUB
        ):
            continue
        seen.add(gene)
        rows.append(
            {
                "gene": native.get(gene, gene),
                "polarity": str(curated.get("polarity") or polarity),
                "n_pub": int(curated.get("n_pub", 0)),
                "tier": str(curated.get("tier") or NOT_AVAILABLE),
                "pct_in": measured["pct_in"],
                "pct_out": measured["pct_out"],
                "avg_log2FC": measured["avg_log2FC"],
                "auc": measured.get("auc"),
                "significant": _significant(measured),
                "specificity": _num(specificity.get(gene), SPECIFICITY_DECIMALS),
                "recommended": False,
            }
        )
    rows.sort(key=lambda row: (-(row["auc"] if row["auc"] is not None else -1.0), row["gene"]))
    return rows


def assign_tiers(pool: dict, scoring: dict) -> None:
    for cluster, state in pool["clusters"].items():
        scored_state = (scoring.get("clusters") or {}).get(cluster) or {}
        extras = {
            str(name)
            for block in (scored_state.get("extra_contexts") or {}).values()
            for name in (block.get("candidates") or [])
        }
        seen_primary = 0
        for entry in state["candidates"]:
            if str(entry["cell_type"]) in extras or entry.get("borrowed_context"):
                entry["tier"] = 1
                continue
            entry["tier"] = seen_primary // CANDIDATE_TIER_SIZE + 1
            seen_primary += 1


def tier_of(entry: dict) -> int:
    return int(entry.get("tier") or 1)


def tiers_available(pool: dict, cluster: str) -> int:
    entries = pool["clusters"][str(cluster)]["candidates"]
    highest = max((tier_of(e) for e in entries), default=1)
    if highest > CANDIDATE_TIERS:
        raise ValueError(
            f"cluster {cluster}: candidates stamped tier {highest} above the "
            f"{CANDIDATE_TIERS}-tier ceiling; assign_tiers and retrieval.top_candidates "
            "disagree"
        )
    return min(CANDIDATE_TIERS, highest)


def candidate_names(pool: dict, cluster: str, tier: int | None = None) -> list[str]:
    return [
        str(entry["cell_type"])
        for entry in pool["clusters"][str(cluster)]["candidates"]
        if tier is None or tier_of(entry) <= int(tier)
    ]


def find_candidate(pool: dict, cluster: str, name: str) -> dict | None:
    for entry in pool["clusters"][str(cluster)]["candidates"]:
        if str(entry["cell_type"]) == str(name):
            return entry
    return None


def _packet_marker(row: dict) -> dict:
    return {field: row[field] for field in PACKET_MARKER_FIELDS if field in row}


MARKER_ROW_FORMAT = (
    "GENE n<publications>[R] <pct_in>/<pct_out>[!] fc<avg_log2FC> auc<auc> "
    "sp<specificity>[ sig][ neg][ ~Name1; Name2]: R the resource's recommended flag; "
    "! not raised here (avg_log2FC at or below 0, or auc at or below 0.5); sig the row "
    "clears the pipeline's significance gate (avg_log2FC above 0.25, pct_in above 10, "
    "auc above 0.5, adjusted p below 0.05); neg the resource curates the gene as a "
    "NEGATIVE marker of this identity, one its cells do not express; ~ the other "
    "candidates of this page the species' literature also ties the gene to, full names "
    "separated by '; '. A '-' where a number belongs means not measured."
)
_SHARED_SEP = "; "


def _fmt_num(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, (bool, np.bool_)):
        return "1" if value else "0"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if number != number:
        return "-"
    return repr(number)


def compact_marker_row(row: dict) -> str:
    n_pub = row.get("n_pub")
    tag = f"{row['gene']} n{'-' if n_pub is None else int(n_pub)}"
    tag += "R" if row.get("recommended") else ""
    tag += f" {_fmt_num(row.get('pct_in'))}/{_fmt_num(row.get('pct_out'))}"
    tag += "" if is_raised(row) else "!"
    tag += f" fc{_fmt_num(row.get('avg_log2FC'))} auc{_fmt_num(row.get('auc'))}"
    tag += f" sp{_fmt_num(row.get('specificity'))}"
    if row.get("significant"):
        tag += " sig"
    if str(row.get("polarity") or "positive") == "negative":
        tag += " neg"
    shared = [str(s) for s in (row.get("shared_on_page") or []) if str(s)]
    if shared:
        tag += " ~" + _SHARED_SEP.join(shared)
    return tag


_NUM = r"(?:-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?|-)"
_COMPACT_ROW = re.compile(
    r"^(?P<gene>\S+) n(?P<pub>\d+|-)(?P<rec>R?) (?P<pin>" + _NUM + r")/(?P<pout>" + _NUM
    + r")(?P<unraised>!?) fc(?P<fc>" + _NUM + r") auc(?P<auc>" + _NUM + r") sp(?P<sp>"
    + _NUM + r")(?P<flags>(?: sig| neg)*)(?: ~(?P<shared>.*))?$"
)


def _parse_num(text: str) -> float | None:
    return None if text == "-" else float(text)


def parse_compact_marker_row(text: str) -> dict:
    m = _COMPACT_ROW.match(str(text))
    if not m:
        raise ValueError(f"not a compact marker row: {text!r}")
    flags = m.group("flags") or ""
    shared = m.group("shared")
    return {
        "gene": m.group("gene"),
        "polarity": "negative" if " neg" in flags else "positive",
        "n_pub": None if m.group("pub") == "-" else int(m.group("pub")),
        "recommended": bool(m.group("rec")),
        "specificity": _parse_num(m.group("sp")),
        "pct_in": _parse_num(m.group("pin")),
        "pct_out": _parse_num(m.group("pout")),
        "avg_log2FC": _parse_num(m.group("fc")),
        "auc": _parse_num(m.group("auc")),
        "significant": " sig" in flags,
        "shared_on_page": [s for s in shared.split(_SHARED_SEP) if s] if shared else [],
        "raised": not m.group("unraised"),
    }


def compact_panel(rows: list[dict]) -> list[str]:
    return [compact_marker_row(row) for row in rows]


def compact_row_genes(rows: list) -> set[str]:
    out: set[str] = set()
    for row in rows or []:
        if isinstance(row, dict):
            gene = row.get("gene")
        else:
            m = re.match(r"^(\S+)", str(row))
            gene = m.group(1) if m else None
        if gene:
            out.add(str(gene).upper())
    return out


def _detected_exclusion_markers(entry: dict) -> list[dict]:
    rows = [
        marker for marker in entry["markers"]
        if marker["polarity"] == "negative" and marker["pct_in"] is not None
        and float(marker["pct_in"]) >= NEGATIVE_SOURCE_MIN_PCT_IN
    ]
    return sorted(rows, key=lambda row: -float(row["pct_in"] or 0.0))


def packet_facts(facts: dict, entry: dict) -> dict:
    out = dict(facts)
    positives = [m for m in entry["markers"] if m["polarity"] == "positive"]
    by_pub = sorted(positives, key=lambda m: (-int(m.get("n_pub") or 0), str(m["gene"])))
    out["recommended_markers"] = [str(m["gene"]) for m in positives if m.get("recommended")]
    out["best_published_markers"] = [str(m["gene"]) for m in by_pub[:2]]
    out["detected_exclusions"] = compact_panel(_detected_exclusion_markers(entry))
    return out


def cluster_marker_table(pool: dict, cluster: str, tier: int | None = None) -> list[dict]:
    rows: dict[str, dict] = {}
    for entry in pool["clusters"][str(cluster)]["candidates"]:
        if tier is not None and tier_of(entry) > int(tier):
            continue
        name = str(entry["cell_type"])
        for marker in entry["markers"]:
            if marker["polarity"] != "positive" or not marker["significant"]:
                continue
            gene = str(marker["gene"])
            row = rows.setdefault(
                gene,
                {
                    "gene": gene,
                    "pct_in": marker["pct_in"],
                    "pct_out": marker["pct_out"],
                    "specificity": marker["specificity"],
                    "claimed_by": [],
                },
            )
            row["claimed_by"].append({"cell_type": name, "n_pub": marker["n_pub"]})
    here = set(candidate_names(pool, cluster, tier))
    claimants = pool.get("species_claimants") or {}
    for gene, row in rows.items():
        row["claimed_by"].sort(key=lambda claim: (-claim["n_pub"], claim["cell_type"]))
        named = {c["cell_type"] for c in row["claimed_by"]}
        also = [name for name, _n in (claimants.get(str(gene).upper()) or [])
                if name in here and name not in named]
        if also:
            row["also_curated_species_wide"] = also[:6]
    ordered = sorted(
        rows.values(),
        key=lambda row: (-((row["pct_in"] or 0.0) - (row["pct_out"] or 0.0)), row["gene"]),
    )
    return ordered[:CLUSTER_MARKER_ROWS]


def candidate_block(pool: dict, entry: dict) -> dict:
    species = str(pool["context"].get("species") or "")
    block = {
        "cell_type": entry["cell_type"],
        "markers": compact_panel(entry["markers"]),
        "exclusion_sources": entry["exclusion_sources"],
        "unmeasured_curated_genes": entry["unmeasured_curated_genes"],
        "program": entry["program"],
    }
    if entry.get("borrowed_context"):
        block["borrowed_context"] = entry["borrowed_context"]
        block["flat_or_undetected_curated_genes"] = entry.get(
            "flat_or_undetected_curated_genes"
        )
    papers = label_papers(species, str(entry["cell_type"]))
    if papers is not None:
        block["n_papers"] = int(papers)
    if entry.get("definers") is not None:
        block["species_definers"] = compact_definer_rows(entry.get("definers") or [])
    return block


def compact_definer_rows(rows: list[dict]) -> list[str]:
    out = []
    for r in rows:
        tag = f"{r['gene']} n{r['n_pub']}{'R' if r.get('recommended') else ''}"
        tag += "" if r.get("in_this_tissue_panel") else "*"
        if "pct_in" in r:
            tag += f" {r['pct_in']}/{r['pct_out']}" + ("" if r.get("raised") else "!")
            peaks = r.get("peak_clusters") or []
            if peaks:
                tag += f" peak={peaks[0]}"
        else:
            tag += f" {r.get('measurement', 'unmeasured')}"
        shared = r.get("shared_with_page") or []
        if shared:
            tag += " ~" + ",".join(shared[:3]) + (f",+{len(shared) - 3}" if len(shared) > 3 else "")
        out.append(tag)
    return out


def cluster_packet(pool: dict, cluster: str, tier: int = 1) -> dict:
    state = pool["clusters"][str(cluster)]
    entries = [e for e in state["candidates"] if tier_of(e) <= int(tier)]
    packet = {
        "query": {
            **pool["context"],
            "cluster_id": str(cluster),
            "cells_in_cluster": int(state["n_cells"]),
        },
        **genome_pages(pool, cluster, tier),
        "cluster_markers": cluster_marker_table(pool, cluster, tier=tier),
        "row_format": MARKER_ROW_FORMAT,
        "candidate_facts": [
            packet_facts(candidate_facts(pool, cluster, entry), entry) for entry in entries
        ],
        "candidates": [
            candidate_block(pool, entry)
            for entry in sorted(entries, key=lambda entry: str(entry["cell_type"]))
        ],
    }
    return packet


def tier_packet(pool: dict, cluster: str, tier: int) -> dict:
    entries = [
        e for e in pool["clusters"][str(cluster)]["candidates"] if tier_of(e) == int(tier)
    ]
    return {
        "additional_candidates_tier": int(tier),
        "cluster_markers": cluster_marker_table(pool, cluster, tier=tier),
        "candidate_facts": [
            packet_facts(candidate_facts(pool, cluster, entry), entry) for entry in entries
        ],
        "candidates": [
            candidate_block(pool, entry)
            for entry in sorted(entries, key=lambda entry: str(entry["cell_type"]))
        ],
    }


def _review_panel(entry: dict) -> list[dict]:
    positives = [row for row in entry["markers"] if row["polarity"] == "positive"]
    negatives = [row for row in entry["markers"] if row["polarity"] == "negative"]
    return [_packet_marker(row) for row in positives + negatives]


def _review_sources(label: str, panel: list[dict], sources) -> dict:
    if sources is None:
        return {}
    found: dict[str, list[dict]] = {}
    for row in panel:
        gene = str(row["gene"]).upper()
        if gene in found:
            continue
        records = _opening_records(sources, label, gene, REVIEW_SOURCES_PER_MARKER)
        if records:
            found[gene] = [_sentence_record(record) for record in records]
    return found


def review_packet(
    pool: dict,
    cluster: str,
    final: dict,
    sources=None,
    tier: int = 1,
    unknown_token: str = "Unknown",
) -> dict:
    state = pool["clusters"][str(cluster)]
    selected = str(final.get("selected") or "")
    subtype = str(final.get("subtype") or "")
    others = [str(name) for name in (final.get("co_occurring_identities") or [])]

    claims: list[tuple[str, str]] = []
    if selected and selected != unknown_token:
        claims.append(("selected", selected))
    if subtype and subtype != selected:
        claims.append(("finer_name", subtype))
    claimed_names = {name for _r, name in claims}
    for name in others:
        if name and name not in claimed_names:
            claims.append(("co_occurring", name))
            claimed_names.add(name)

    claimed: list[dict] = []
    for role, name in claims:
        block = _claimed_block(pool, cluster, name, role, sources)
        if block is not None:
            claimed.append(block)

    contested = [name for name in (final.get("contested_rivals") or [])
                 if name and name not in claimed_names]
    contested_blocks = [
        block for block in (
            _claimed_block(pool, cluster, name, "contested_rival", sources, final)
            for name in contested
        ) if block is not None
    ]
    covered = {name for _role, name in claims} | {b["cell_type"] for b in contested_blocks}

    packet = {
        "query": {
            **pool["context"],
            "cluster_id": str(cluster),
            "cells_in_cluster": int(state["n_cells"]),
        },
        "delivered": {
            "selected": selected,
            "subtype": subtype,
            "co_occurring_identities": others,
            "state": str(final.get("state") or ""),
            "reason": str(final.get("reason") or ""),
        },
        "row_format": MARKER_ROW_FORMAT,
        "claimed_identities": claimed,
        **genome_pages(pool, cluster, tier),
        "candidates_not_claimed": _other_candidate_blocks(pool, cluster, covered, tier),
    }
    if contested_blocks:
        packet["contested_identities"] = contested_blocks
        packet["shared_with_contested"] = _shared_with(
            pool, cluster, selected, [b["cell_type"] for b in contested_blocks]
        )
        sel_entry = find_candidate(pool, cluster, selected)
        if sel_entry is not None:
            packet["pair_partitions"] = [
                pair_partition(sel_entry, rentry)
                for rentry in (find_candidate(pool, cluster, b["cell_type"]) for b in contested_blocks)
                if rentry is not None
            ]
    if final.get("definer_audit"):
        packet["definer_audit"] = final["definer_audit"]
    screened = screened_out_claimants(pool, cluster)
    if screened:
        packet["screened_out_claimants"] = screened
    return packet


SCREENED_OUT_SHOWN = 5
SCREENED_OUT_ANCHORS_SHOWN = 12


def screened_out_claimants(pool: dict, cluster: str) -> list[dict]:
    log = (pool.get("borrowed") or {}).get(str(cluster)) or {}
    native_sym = dict(pool.get("native_gene") or {})
    prefix = "tissue screen: "
    out = []
    for item in log.get("considered") or []:
        rejected = str(item.get("rejected") or "")
        if not rejected.startswith(prefix):
            continue
        anchors = [str(g) for g in (item.get("anchor_genes") or [])]
        out.append({
            "cell_type": str(item.get("cell_type") or ""),
            "borrow_rank": int(item.get("borrow_rank") or 0),
            "anchors_covered": int(item.get("anchors_covered") or len(anchors)),
            "anchor_genes": [native_sym.get(g, g) for g in anchors][:SCREENED_OUT_ANCHORS_SHOWN],
            "recommended_anchors": [native_sym.get(str(g), str(g))
                                    for g in (item.get("recommended_anchors") or [])],
            "screen_reason": rejected[len(prefix):],
        })
    out.sort(key=lambda r: (r["borrow_rank"] or 10**6, -r["anchors_covered"], r["cell_type"]))
    return out[:SCREENED_OUT_SHOWN]


def _shared_with(pool: dict, cluster: str, identity: str, others: list[str]) -> dict[str, list[str]]:
    entry = find_candidate(pool, cluster, identity)
    if entry is None or not others:
        return {}
    genes = {str(m["gene"]).upper() for m in entry["markers"] if m["polarity"] == "positive"}
    genes |= {str(d["gene"]).upper() for d in entry.get("definers") or []}
    claimants = pool.get("species_claimants") or {}
    others_set = set(others)
    out = {}
    for gene in sorted(genes):
        also = [name for name, _n_pub in (claimants.get(gene) or []) if name in others_set]
        if also:
            out[gene] = also
    return out


def _claimed_block(pool: dict, cluster: str, name: str, role: str, sources,
                   final: dict | None = None) -> dict | None:
    entry = find_candidate(pool, cluster, name)
    if entry is None:
        return None
    panel = _review_panel(entry)
    block = {
        "cell_type": name,
        "role": role,
        "panel": compact_panel(panel),
        "sources": _review_sources(name, panel, sources),
        "detected_exclusions": compact_panel(_detected_exclusion_markers(entry)),
        "unmeasured_curated_genes": entry["unmeasured_curated_genes"],
        "facts": packet_facts(candidate_facts(pool, cluster, entry), entry),
    }
    if entry.get("definers") is not None:
        block["species_definers"] = definer_block(entry)
    if entry.get("borrowed_context"):
        block["borrowed_context"] = entry["borrowed_context"]
    audit = ((final or {}).get("definer_audit") or {}).get(name)
    if audit:
        block["definer_audit"] = audit
    return block


def _other_candidate_blocks(
    pool: dict, cluster: str, claimed: set[str], tier: int
) -> list[dict]:
    blocks = []
    for entry in pool["clusters"][str(cluster)]["candidates"]:
        name = str(entry["cell_type"])
        if name in claimed or tier_of(entry) > int(tier):
            continue
        rows = _review_panel(entry)
        positives = [r for r in rows if r.get("polarity") == "positive"][
            :REVIEW_OTHER_CANDIDATE_ROWS
        ]
        negatives = [r for r in rows if r.get("polarity") == "negative"]
        block = {"cell_type": name, "panel": compact_panel(positives + negatives)}
        if entry.get("borrowed_context"):
            block["borrowed_context"] = entry["borrowed_context"]
        papers = label_papers(str(pool["context"].get("species") or ""), name)
        if papers is not None:
            block["n_papers"] = int(papers)
        blocks.append(block)
    return blocks


def _carried(marker: dict) -> bool:
    pct_in = float(marker.get("pct_in") or 0.0)
    raised = is_raised(marker) and bool(marker.get("significant")) and pct_in >= DM_GATE_FLOOR_PCT
    return raised or pct_in >= DM_GATE_DOMINANT_PCT


def definers_absent(rows: list[dict], k: int) -> bool:
    top = sorted(rows, key=lambda m: (-int(m.get("n_pub") or 0), str(m["gene"])))[:k]
    return len(top) == k and all(float(m.get("pct_in") or 0.0) < DM_GATE_FLOOR_PCT for m in top)


def genome_de_summary(markers_all: pd.DataFrame) -> dict[str, dict]:
    out: dict[str, dict] = {}
    if markers_all is None:
        return out
    df = markers_all
    for cluster, sub in df.groupby(df["group"].astype(str)):
        up = sub[(sub["padj"] < 0.05) & (sub["avg_log2FC"] > 0) & (sub["pct_in"] >= 10.0)]
        up = up.sort_values("avg_log2FC", ascending=False).head(GENOME_ENRICHED_ROWS)
        down = sub[(sub["avg_log2FC"] < 0) & (sub["pct_out"] >= 20.0)]
        down = down.sort_values("avg_log2FC").head(GENOME_DEPLETED_ROWS)
        anchors = sub[(sub["padj"] < 0.05) & (sub["avg_log2FC"] > 0) & (sub["pct_in"] >= DM_GATE_DOMINANT_PCT)]
        anchors = anchors.sort_values("avg_log2FC", ascending=False).head(GENOME_ENRICHED_ROWS)
        out[str(cluster)] = {
            "top_enriched": [
                {"gene": str(r.feature), "log2FC": _num(r.avg_log2FC, 2), "pct_in": _num(r.pct_in, 1),
                 "pct_out": _num(r.pct_out, 1), "auc": _num(r.auc, 3)}
                for r in up.itertuples()
            ],
            "top_depleted": [
                {"gene": str(r.feature), "log2FC": _num(r.avg_log2FC, 2), "pct_in": _num(r.pct_in, 1),
                 "pct_out": _num(r.pct_out, 1)}
                for r in down.itertuples()
            ],
            "anchor_genes": [
                {"gene": str(r.feature), "log2FC": _num(r.avg_log2FC, 2), "pct_in": _num(r.pct_in, 1),
                 "pct_out": _num(r.pct_out, 1)}
                for r in anchors.itertuples()
            ],
        }
    return out


DEFINERS_PER_CANDIDATE = int(_A.get("definers_per_candidate", 8))
DEFINER_TOP_CLUSTERS = int(_A.get("definer_top_clusters", 3))
GENOME_DETECTION_GAP_ROWS = int(_A.get("genome_detection_gap_rows", 20))
GENOME_CLAIMANTS_SHOWN = int(_A.get("genome_claimants_shown", 3))


def _measurement_for(pool: dict, cluster: str, gene: str) -> tuple[dict | None, str]:
    key = (str(cluster), str(gene).upper())
    row = pool.get("de_stats", {}).get(key)
    if row is not None:
        return row, "de_table"
    row = pool.get("genome_stats", {}).get(key)
    if row is not None:
        return row, "genome_table"
    if key[1] in (pool.get("measured_genes") or set()):
        return None, "below_reporting_screen"
    return None, "not_measured"


def _top_clusters_for(pool: dict, gene: str, k: int) -> list[dict]:
    rows = pool.get("gene_clusters", {}).get(str(gene).upper()) or []
    ranked = sorted(rows, key=lambda item: -(item[1] or 0.0))[:k]
    return [{"cluster_id": item[0], "pct_in": item[1], "pct_out": item[2]} for item in ranked]


def attach_definers(pool: dict, definers_by_name: dict[str, list[dict]]) -> int:
    n = 0
    for cluster, state in pool["clusters"].items():
        tissue_claimants: dict[str, set[str]] = {}
        for e in state["candidates"]:
            e_name = str(e["cell_type"])
            for m in e["markers"]:
                if m["polarity"] == "positive":
                    tissue_claimants.setdefault(str(m["gene"]).upper(), set()).add(e_name)
        for entry in state["candidates"]:
            name = str(entry["cell_type"])
            panel_genes = {
                str(m["gene"]).upper() for m in entry["markers"] if m["polarity"] == "positive"
            }
            rows = []
            for d in definers_by_name.get(name) or []:
                gene_u = str(d["gene"]).upper()
                measured, how = _measurement_for(pool, cluster, gene_u)
                in_panel = gene_u in panel_genes
                row = {
                    "gene": pool.get("native_gene", {}).get(gene_u, d["gene"]),
                    "n_pub": int(d["n_pub"]),
                    "recommended": bool(d["recommended"]),
                    "in_this_tissue_panel": in_panel,
                }
                if not in_panel:
                    others = sorted(tissue_claimants.get(gene_u, set()) - {name})
                    if others:
                        row["shared_with_page"] = others
                if measured is not None:
                    row.update(
                        {
                            "pct_in": measured["pct_in"],
                            "pct_out": measured["pct_out"],
                            "auc": measured.get("auc"),
                            "raised": is_raised(measured),
                            "peak_clusters": [
                                f"{c['cluster_id']}:{c['pct_in']}"
                                for c in _top_clusters_for(pool, gene_u, DEFINER_TOP_CLUSTERS)
                            ],
                        }
                    )
                else:
                    row["measurement"] = how
                rows.append(row)
            entry["definers"] = rows
            n += int(bool(rows))
    return n


PARTITION_OWN_SHOWN = int(_A.get("partition_own_shown", 15))
PARTITION_SHARED_SHOWN = int(_A.get("partition_shared_shown", 8))


def attach_page_sharing(pool: dict, species_claimants: dict[str, list[tuple[str, int]]]) -> int:
    n = 0
    for cluster, state in pool["clusters"].items():
        here = {str(e["cell_type"]) for e in state["candidates"]}
        cache: dict[str, list[str]] = {}

        def others_for(gene_u: str) -> list[str]:
            if gene_u not in cache:
                cache[gene_u] = [name for name, _n in (species_claimants.get(gene_u) or [])
                                 if name in here]
            return cache[gene_u]

        for entry in state["candidates"]:
            name = str(entry["cell_type"])
            for marker in entry["markers"]:
                if marker.get("polarity") != "positive":
                    marker.pop("shared_on_page", None)
                    continue
                shared = [o for o in others_for(str(marker["gene"]).upper()) if o != name]
                if shared:
                    marker["shared_on_page"] = shared
                    n += 1
                else:
                    marker.pop("shared_on_page", None)
            for row in entry.get("definers") or []:
                shared = [o for o in others_for(str(row["gene"]).upper()) if o != name]
                prior = [o for o in (row.get("shared_with_page") or []) if o != name]
                merged = sorted(set(shared) | set(prior))
                if merged:
                    row["shared_with_page"] = merged
                    n += 1
                else:
                    row.pop("shared_with_page", None)
    return n


def _partition_string(row: dict, with_sharers: bool = False) -> str:
    tag = f"{row['gene']} n{int(row.get('n_pub') or 0)}"
    if row.get("pct_in") is not None:
        tag += f" {row['pct_in']}/{row['pct_out']}" + ("" if row.get("_raised") else "!")
    else:
        tag += f" {row.get('measurement') or 'unmeasured'}"
    shared = row.get("shared") or []
    if with_sharers and shared:
        tag += " ~" + ",".join(shared[:3]) + (f",+{len(shared) - 3}" if len(shared) > 3 else "")
    return tag


def _gene_rows(entry: dict) -> list[dict]:
    by_gene: dict[str, dict] = {}
    for marker in entry.get("markers") or []:
        if marker.get("polarity") != "positive":
            continue
        gene_u = str(marker["gene"]).upper()
        by_gene[gene_u] = {
            "gene": marker["gene"], "n_pub": marker.get("n_pub"),
            "pct_in": marker.get("pct_in"), "pct_out": marker.get("pct_out"),
            "auc": marker.get("auc"), "_raised": is_raised(marker),
            "shared": list(marker.get("shared_on_page") or []),
        }
    for row in entry.get("definers") or []:
        gene_u = str(row["gene"]).upper()
        if gene_u in by_gene:
            if int(row.get("n_pub") or 0) > int(by_gene[gene_u].get("n_pub") or 0):
                by_gene[gene_u]["n_pub"] = row.get("n_pub")
            continue
        by_gene[gene_u] = {
            "gene": row["gene"], "n_pub": row.get("n_pub"),
            "pct_in": row.get("pct_in"), "pct_out": row.get("pct_out"),
            "auc": row.get("auc"), "_raised": bool(row.get("raised")),
            "measurement": row.get("measurement"),
            "shared": list(row.get("shared_with_page") or []),
        }

    def order(r: dict):
        auc = r.get("auc")
        return (0 if auc is not None else 1, -(float(auc) if auc is not None else 0.0), str(r["gene"]))

    return sorted(by_gene.values(), key=order)


def lineage_partition(entry: dict) -> dict:
    rows = _gene_rows(entry)
    own = [r for r in rows if not r["shared"]]
    shared_rows = [r for r in rows if r["shared"]]
    return {
        "own": [_partition_string(r) for r in own[:PARTITION_OWN_SHOWN]],
        "own_total": len(own),
        "shared": [_partition_string(r, with_sharers=True)
                   for r in shared_rows[:PARTITION_OWN_SHOWN]],
        "shared_total": len(shared_rows),
    }


def pair_partition(entry_a: dict, entry_b: dict) -> dict:
    name_a = str(entry_a["cell_type"])
    name_b = str(entry_b["cell_type"])
    rows_a = {str(r["gene"]).upper(): r for r in _gene_rows(entry_a)}
    rows_b = {str(r["gene"]).upper(): r for r in _gene_rows(entry_b)}
    both_keys = set()
    for g, r in rows_a.items():
        if g in rows_b or name_b in (r.get("shared") or []):
            both_keys.add(g)
    for g, r in rows_b.items():
        if g in rows_a or name_a in (r.get("shared") or []):
            both_keys.add(g)
    only_a = [r for g, r in rows_a.items() if g not in both_keys]
    only_b = [r for g, r in rows_b.items() if g not in both_keys]
    both = [rows_a.get(g) or rows_b.get(g) for g in both_keys]

    def order(r: dict):
        auc = r.get("auc")
        return (0 if auc is not None else 1, -(float(auc) if auc is not None else 0.0), str(r["gene"]))

    return {
        "pair": [name_a, name_b],
        f"only_{name_a}": [_partition_string(r) for r in sorted(only_a, key=order)[:PARTITION_OWN_SHOWN]],
        f"only_{name_b}": [_partition_string(r) for r in sorted(only_b, key=order)[:PARTITION_OWN_SHOWN]],
        "both": [_partition_string(r) for r in sorted(both, key=order)[:PARTITION_OWN_SHOWN]],
        "totals": {name_a: len(only_a), name_b: len(only_b), "both": len(both)},
    }


def attach_genome_claims(
    pool: dict,
    claimants: dict[str, list[tuple[str, int]]],
    known_genes: set[str],
    markers_all: pd.DataFrame | None,
) -> None:
    genome = pool.get("genome_de") or {}
    gaps: dict[str, list[dict]] = {}
    if markers_all is not None and len(markers_all):
        df = markers_all
        for cluster, sub in df.groupby(df["group"].astype(str)):
            dom = sub[sub["pct_in"] >= DM_GATE_DOMINANT_PCT].copy()
            dom["_gap"] = dom["pct_in"] - dom["pct_out"]
            dom = dom.sort_values("_gap", ascending=False).head(GENOME_DETECTION_GAP_ROWS * 2)
            gaps[str(cluster)] = [
                {"gene": str(r.feature), "pct_in": _num(r.pct_in, 1), "pct_out": _num(r.pct_out, 1),
                 "log2FC": _num(r.avg_log2FC, 2), "auc": _num(r.auc, 3)}
                for r in dom.itertuples()
            ]

    def _claims(gene: str) -> list[dict]:
        return [
            {"cell_type": name, "n_pub": n_pub}
            for name, n_pub in (claimants.get(str(gene).upper()) or [])[:GENOME_CLAIMANTS_SHOWN]
        ]

    def _filter(rows: list[dict], cap: int) -> tuple[list[dict], int]:
        kept, dropped = [], 0
        for row in rows:
            if str(row["gene"]).upper() not in known_genes:
                dropped += 1
                continue
            kept.append({**row, "claimed_by_species": _claims(row["gene"])})
            if len(kept) >= cap:
                break
        return kept, dropped

    for cluster in pool["clusters"]:
        block = genome.setdefault(str(cluster), {})
        enriched, d1 = _filter(block.get("top_enriched", []), GENOME_ENRICHED_ROWS)
        depleted, d2 = _filter(block.get("top_depleted", []), GENOME_DEPLETED_ROWS)
        gap, d3 = _filter(gaps.get(str(cluster), []), GENOME_DETECTION_GAP_ROWS)
        block["top_enriched"] = enriched
        block["top_depleted"] = depleted
        block["top_detection_gap"] = gap
        block["rows_not_curated_anywhere"] = {"enriched": d1, "depleted": d2, "detection_gap": d3}
    pool["genome_de"] = genome


def _mark_on_page(rows: list[dict], on_page: set[str]) -> list[dict]:
    out = []
    for row in rows:
        claims = [
            {**c, "on_page": c["cell_type"] in on_page} for c in row.get("claimed_by_species", [])
        ]
        out.append({**row, "claimed_by_species": claims})
    return out


def genome_pages(pool: dict, cluster: str, tier: int | None) -> dict:
    block = (pool.get("genome_de") or {}).get(str(cluster)) or {}
    on_page = set(candidate_names(pool, cluster, tier))
    return {
        "cluster_top_enriched": _mark_on_page(block.get("top_enriched", []), on_page),
        "cluster_top_detection_gap": _mark_on_page(block.get("top_detection_gap", []), on_page),
        "cluster_top_depleted": _mark_on_page(block.get("top_depleted", []), on_page),
        "genome_rows_not_curated_anywhere": block.get("rows_not_curated_anywhere", {}),
    }


def definer_block(entry: dict) -> list[dict]:
    return list(entry.get("definers") or [])


def attach_recommended(pool: dict, recommended_pairs: set[tuple[str, str]]) -> int:
    n = 0
    for state in pool["clusters"].values():
        for entry in state["candidates"]:
            name = str(entry["cell_type"])
            for marker in entry["markers"]:
                flag = (name, str(marker["gene"]).upper()) in recommended_pairs
                marker["recommended"] = bool(flag)
                n += int(flag)
    return n


def _measured_row(marker: dict) -> dict:
    return {
        "gene": marker["gene"],
        "n_pub": int(marker.get("n_pub") or 0),
        "recommended": bool(marker.get("recommended")),
        "pct_in": marker["pct_in"],
        "pct_out": marker["pct_out"],
        "avg_log2FC": marker["avg_log2FC"],
        "auc": marker.get("auc"),
        "significant": bool(marker.get("significant")),
    }


def detected_exclusions(entry: dict) -> list[dict]:
    return [_measured_row(marker) for marker in _detected_exclusion_markers(entry)]


def candidate_facts(pool: dict, cluster: str, entry: dict) -> dict:
    species = str(pool["context"].get("species") or "")
    name = str(entry["cell_type"])
    positives = [m for m in entry["markers"] if m["polarity"] == "positive"]
    by_pub = sorted(positives, key=lambda m: (-int(m.get("n_pub") or 0), str(m["gene"])))
    raised = [m for m in positives if is_raised(m) and m.get("significant")]
    papers = label_papers(species, name)
    return {
        "cell_type": name,
        "retrieval_rank": int(entry.get("retrieval_rank") or 0),
        "tier": tier_of(entry),
        "n_papers": None if papers is None else int(papers),
        "recommended_markers": [_measured_row(m) for m in positives if m.get("recommended")],
        "best_published_markers": [_measured_row(m) for m in by_pub[:2]],
        "positive_markers_measured": len(positives),
        "positive_markers_raised_and_significant": len(raised),
        "detected_exclusions": detected_exclusions(entry),
        "own_vs_shared": lineage_partition(entry),
    }


def candidate_facts_table(pool: dict, cluster: str, tier: int | None = None) -> list[dict]:
    return [
        candidate_facts(pool, cluster, entry)
        for entry in pool["clusters"][str(cluster)]["candidates"]
        if tier is None or tier_of(entry) <= int(tier)
    ]


def borrow_admission_row(pool: dict, cluster: str, entry: dict) -> dict:
    name = str(entry["cell_type"])
    positives = [m for m in entry["markers"] if m["polarity"] == "positive"]
    rec = [m for m in positives if m.get("recommended")]
    rec_carried = [m for m in rec if _carried(m)]
    top2 = sorted(positives, key=lambda m: (-int(m.get("n_pub") or 0), str(m["gene"])))[:2]
    top2_carried = [m for m in top2 if _carried(m)]
    if rec:
        gate, gate_basis = len(rec_carried) >= 1, "recommended"
    elif len(top2) >= 2:
        gate, gate_basis = len(top2_carried) >= 1, "top2_n_pub"
    else:
        gate, gate_basis = True, "too_thin_to_gate"
    if gate and definers_absent(positives, DM_GATE_TOP):
        gate, gate_basis = False, "top2_absent"
    return {
        "cell_type": name,
        "recommended_measured": [m["gene"] for m in rec],
        "recommended_carried": [m["gene"] for m in rec_carried],
        "top2_n_pub": [m["gene"] for m in top2],
        "top2_carried": [m["gene"] for m in top2_carried],
        "identity_gate": "pass" if gate else "fail",
        "identity_gate_basis": gate_basis,
    }


def run_tool(pool: dict, cluster: str, tool: str, args: Any, sources=None) -> dict:
    args = args if isinstance(args, dict) else {}
    if tool == "sources":
        return _tool_sources(pool, cluster, args, sources)
    if tool == "gene_across_clusters":
        return _tool_gene_across_clusters(pool, args)
    if tool == "candidates_with_gene":
        return _tool_candidates_with_gene(pool, cluster, args)
    return {"tool": str(tool), "error": f"unknown tool; available: {', '.join(TOOLS)}"}


def run_review_tool(
    pool: dict, cluster: str, labels: tuple[str, ...], tool: str, args: Any, sources=None
) -> dict:
    args = args if isinstance(args, dict) else {}
    if tool not in REVIEW_TOOLS:
        return {
            "tool": str(tool),
            "error": f"unknown tool; available: {', '.join(REVIEW_TOOLS)}",
        }
    if tool == "gene_across_clusters":
        return _tool_gene_across_clusters(pool, args)
    if tool == "candidates_with_gene":
        return _tool_candidates_with_gene(pool, cluster, args)
    asked = str(args.get("label") or args.get("candidate") or "")
    if asked and asked not in labels:
        return {
            "tool": "sources",
            "error": (
                f"'{asked}' is not under test here; "
                f"available: {', '.join(labels)}"
            ),
        }
    label = asked or (labels[0] if labels else "")
    return _tool_sources(pool, cluster, {**args, "candidate": label}, sources)


def _tool_sources(pool: dict, cluster: str, args: dict, sources) -> dict:
    candidate = str(args.get("candidate") or "")
    entry = find_candidate(pool, cluster, candidate)
    if entry is None:
        return {
            "tool": "sources",
            "error": f"'{candidate}' is not a candidate of this cluster",
        }
    requested = args.get("genes")
    requested = requested if isinstance(requested, list) else [requested]
    genes = [str(gene).strip() for gene in requested if str(gene or "").strip()]
    truncated = len(genes) > SOURCE_GENES_PER_QUERY
    result = {}
    left = {}
    for gene in genes[:SOURCE_GENES_PER_QUERY]:
        key = gene.upper()
        if sources is None:
            result[key] = []
            continue
        if hasattr(sources, "take"):
            answer = sources.take(candidate, key)
            records = answer["sources"]
            if answer["remaining"] or answer["limit_reached"]:
                left[key] = (
                    "limit reached"
                    if answer["limit_reached"]
                    else f"{answer['remaining']} more"
                )
        else:
            records = sources.records_for_marker(candidate, key, k=SOURCES_PER_GENE)
        result[key] = [_sentence_record(record) for record in records]
    answer = {
        "tool": "sources",
        "candidate": candidate,
        "sources": result,
        "truncated": truncated,
        "also_curated_on_this_page": _page_claimants(pool, cluster, list(result)),
    }
    if left:
        answer["not_yet_shown"] = left
    return answer


def _page_claimants(pool: dict, cluster: str, genes: list[str]) -> dict[str, list[dict]]:
    here = set(candidate_names(pool, cluster))
    out: dict[str, list[dict]] = {}
    for gene in genes:
        rows = [
            {"cell_type": name, "n_pub": int(n_pub)}
            for name, n_pub in (pool.get("species_claimants", {}).get(str(gene).upper()) or [])
            if name in here
        ]
        if len(rows) > 1:
            out[str(gene).upper()] = rows[:6]
    return out


def _tool_gene_across_clusters(pool: dict, args: dict) -> dict:
    gene = str(args.get("gene") or "").strip().upper()
    rows = pool["gene_clusters"].get(gene)
    if not rows:
        if gene in (pool.get("measured_genes") or set()):
            return {
                "tool": "gene_across_clusters",
                "gene": gene,
                "measured": True,
                "note": "measured in this dataset, but below the reporting screen in every "
                        "cluster: no cluster separates it (|log2FC| under the floor, or "
                        "detected in under 1% of cells inside and outside)",
                "clusters": [],
            }
        return {
            "tool": "gene_across_clusters",
            "gene": gene,
            "measured": False,
            "error": "not measured in this dataset",
        }
    return {
        "tool": "gene_across_clusters",
        "gene": pool["native_gene"].get(gene, gene),
        "measured": True,
        "clusters": [
            {"cluster_id": item[0], "pct_in": item[1], "pct_out": item[2]}
            for item in sorted(rows, key=lambda item: -(item[1] or 0.0))
        ],
    }


def _tool_candidates_with_gene(pool: dict, cluster: str, args: dict) -> dict:
    gene = str(args.get("gene") or "").strip().upper()
    here = set(candidate_names(pool, cluster))
    rows = [
        {"cell_type": name, "n_pub": n_pub, "tier": tier, "polarity": polarity}
        for name, n_pub, tier, polarity in pool["gene_carriers"].get(gene, [])
        if name in here
    ]
    rows.sort(key=lambda row: (-row["n_pub"], row["cell_type"]))
    return {
        "tool": "candidates_with_gene",
        "gene": pool["native_gene"].get(gene, gene),
        "candidates": rows,
        "candidates_species_wide": _page_claimants(pool, cluster, [gene]).get(gene, []),
    }


def claim_warnings(
    pool: dict, cluster: str, claims: list[tuple[str, str]]
) -> list[tuple[str, str]]:
    lines = []
    for role, name in claims:
        entry = find_candidate(pool, cluster, name)
        if entry is None:
            continue
        positive = [
            marker for marker in entry["markers"] if marker["polarity"] == "positive"
        ]
        unraised = [marker for marker in positive if not is_raised(marker)]
        if not positive or not unraised:
            continue
        shown = unraised[:WARNING_GENES_SHOWN]
        listed = ",".join(
            f"{marker['gene']}({marker['pct_in']}/{marker['pct_out']})"
            for marker in shown
        )
        lines.append((
            name,
            f"{role} {name}: not_raised {len(unraised)}/{len(positive)}"
            f" | top_n_pub: {listed} | +{len(unraised) - len(shown)} more",
        ))
    return lines
