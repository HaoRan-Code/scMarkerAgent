#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import collections
import csv
import hashlib
import math
import os
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import annotator_pool as ap
from .candidate_scoring import (
    EPS,
    _percentile_rank,
    cross_cluster_percentile,
    marker_specificity,
)
from .config import DEFAULTS
from .evidence_gate import significant_genes_by_cluster
from .marker_database import normalize_disease_query

_B = DEFAULTS["borrowed_context"]
BORROW_MIN_ANCHORS = int(_B["min_anchor_genes"])
BORROW_MIN_PUB = ap.DECISIVE_MIN_PUB
BORROW_MAX_PER_CLUSTER = int(_B["max_per_cluster"])
BORROW_ANCHOR_ROWS = ap.GENOME_ENRICHED_ROWS
BORROW_ANCHOR_MIN_PCT = ap.DM_GATE_DOMINANT_PCT
BORROW_QUALIFIER_PROMISCUITY_MIN = int(_B["qualifier_promiscuity_min"])
GENERIC_HEAD_NOUNS = frozenset(str(value) for value in _B["generic_head_nouns"])
BORROW_OWN_SIGNIFICANT_MIN = int(_B["own_significant_markers_min"])


def norm_name(value) -> str:
    text = str(value).strip().lower().replace(",", " ")
    return re.sub(r"\s+", " ", text).strip()


def _label_papers_normalized() -> dict[tuple[str, str], int]:
    path = os.environ.get("SCMA_LABEL_PAPERS", "") or str(
        Path(__file__).resolve().parent / "resources" / "label_papers_v4.csv"
    )
    table: dict[tuple[str, str], int] = {}
    if os.path.exists(path):
        with open(path, newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                table[(str(row["species"]), norm_name(row["cell_type"]))] = int(row["n_papers"])
    return table


def _lineage_classes(parents: list[str]) -> int:
    items = sorted(set(parents))
    root = list(range(len(items)))

    def find(a: int) -> int:
        while root[a] != a:
            root[a] = root[root[a]]
            a = root[a]
        return a

    toks = [set(p.split(" ")) - GENERIC_HEAD_NOUNS for p in items]
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            if toks[i] & toks[j]:
                ri, rj = find(i), find(j)
                if ri != rj:
                    root[ri] = rj
    return len({find(i) for i in range(len(items))})


def qualifier_lexicon(db: pd.DataFrame, cutoff: int = BORROW_QUALIFIER_PROMISCUITY_MIN) -> dict[str, Any]:
    papers = _label_papers_normalized()
    frame = db[["species", "cell_type"]].dropna().astype(str).drop_duplicates()
    by_species: dict[str, set[str]] = collections.defaultdict(set)
    for species, name in zip(frame["species"], frame["cell_type"]):
        by_species[species].add(norm_name(name))
    raw: dict[str, set[str]] = collections.defaultdict(set)
    for species, names in by_species.items():
        for name in names:
            tokens = name.split(" ")
            if len(tokens) < 2:
                continue
            parent = " ".join(tokens[1:])
            if parent in names and parent not in GENERIC_HEAD_NOUNS and papers.get((species, parent), 0) > ap.SPARSE_LABEL_PAPERS_MAX:
                raw[tokens[0]].add(parent)
    tokens = sorted(w for w, ps in raw.items() if _lineage_classes(sorted(ps)) >= cutoff)
    digest = hashlib.sha256("\n".join(tokens).encode("utf-8")).hexdigest()
    return {"tokens": tokens, "size": len(tokens), "sha256": digest, "cutoff": int(cutoff),
            "vocabulary_pairs": int(len(frame)), "modifier_tokens": len(raw)}


def qualified_form_of(claimant: str, eligible: list[str], lexicon: set[str]) -> tuple[str, str] | None:
    tokens = norm_name(claimant).split(" ")
    local = {norm_name(x): x for x in eligible if norm_name(x) not in GENERIC_HEAD_NOUNS}
    for k in range(1, len(tokens)):
        extra, rest = tokens[:k], " ".join(tokens[k:])
        if rest in local and all(t in lexicon for t in extra):
            return local[rest], " ".join(extra)
    return None


def resource_slice(db: pd.DataFrame, species: str, disease) -> pd.DataFrame:
    terms = normalize_disease_query(disease)
    frame = db[(db["species"].astype(str) == str(species)) & db["marker_polarity"].isin(["positive", "negative"])]
    qc_ok = (frame["gene_qc_pass"] == "TRUE") | frame["gene_qc_pass"].isna()
    mask = qc_ok & (~frame["is_in_vitro"]) & frame["gene_symbol"].notna() & frame["cell_type"].notna()
    if terms is not None:
        mask = mask & frame["_disease_norm"].isin(terms)
    out = frame.loc[mask, ["tissue_type", "gene_symbol", "cell_type", "marker_polarity", "n_pub_support",
                           "confidence_tier", "is_recommended_marker"]].copy()
    out["gu"] = out["gene_symbol"].astype(str).str.upper()
    out["cell_type"] = out["cell_type"].astype(str)
    out["recommended"] = out["is_recommended_marker"].astype(str).str.upper() == "TRUE"
    return out


def translate_resource(resource: pd.DataFrame, gene_map: dict[str, str]) -> pd.DataFrame:
    out = resource.copy()
    mapped = out["gu"].map(gene_map)
    out = out[mapped.notna()].copy()
    out["donor_gene"] = out["gene_symbol"]
    out["gene_symbol"] = mapped
    out["gu"] = out["gene_symbol"].astype(str).str.upper()
    return out


def _type_tables(positive: pd.DataFrame) -> tuple[dict, dict, dict]:
    by_type_genes = positive.groupby("cell_type")["gu"].apply(set).to_dict()
    rec_by_type = positive[positive["recommended"]].groupby("cell_type")["gu"].apply(set).to_dict()
    npub_by_pair = positive.groupby(["cell_type", "gu"])["n_pub_support"].max().to_dict()
    return by_type_genes, rec_by_type, npub_by_pair


def _panel_meta(rows: pd.DataFrame) -> tuple[dict, list[str], list[str], dict]:
    tier_rank = {"high": 0, "medium": 1, "low": 2}
    rows = rows.assign(_tr=rows["confidence_tier"].astype(str).str.lower().map(tier_rank).fillna(3).astype(int))
    rows = rows.sort_values(["n_pub_support", "_tr"], ascending=[False, True], kind="stable")
    meta: dict = {}
    positive: list[str] = []
    negative: list[str] = []
    tissues: dict[str, set] = {}
    name = str(rows["cell_type"].iloc[0])
    for r in rows.itertuples(index=False):
        gene, pol = str(r.gu), str(r.marker_polarity)
        tissues.setdefault(gene, set()).add(str(r.tissue_type))
        if (name, gene) in meta:
            continue
        meta[(name, gene)] = {"n_pub": int(r.n_pub_support), "tier": str(r.confidence_tier), "polarity": pol}
        (positive if pol == "positive" else negative).append(gene)
    return meta, positive, negative, tissues


SCREEN_FOREIGN = "foreign"
SCREEN_PLAUSIBLE = "plausible"
SCREEN_VERDICTS = (SCREEN_PLAUSIBLE, SCREEN_FOREIGN)
SCREEN_REJECTED_PREFIX = "tissue screen: "


def borrow_candidates(pool: dict, resource: pd.DataFrame, measured_genes: set[str],
                      markers_all: pd.DataFrame, sources=None,
                      qualifier_tokens: set[str] | None = None,
                      species: str | None = None,
                      cross_species_resource: dict[str, pd.DataFrame] | None = None,
                      ortho: Any = None,
                      screen: Any = None) -> dict[str, list[dict]]:
    native_genes = set(pool.get("native_menu_genes") or [])
    all_positive = resource[resource["marker_polarity"] == "positive"]
    positive = all_positive[all_positive["n_pub_support"] >= BORROW_MIN_PUB]
    by_type_genes, rec_by_type, npub_by_pair = _type_tables(positive)
    panel_genes_by_type = all_positive.groupby("cell_type")["gu"].apply(set).to_dict()
    tissues_by_type = resource.groupby("cell_type")["tissue_type"].nunique().to_dict()

    cross_tables: dict[str, dict] = {}
    if species and cross_species_resource and ortho is not None:
        for other, other_resource in cross_species_resource.items():
            gene_map = ortho.one_to_one(other, species)
            donor_positive = other_resource[
                (other_resource["marker_polarity"] == "positive") & (other_resource["n_pub_support"] >= BORROW_MIN_PUB)
            ]
            translated = translate_resource(donor_positive, gene_map)
            translated = translated[~translated["cell_type"].isin(by_type_genes)]
            by_type_genes_x, rec_by_type_x, npub_by_pair_x = _type_tables(translated)
            cross_tables[other] = {
                "translated": translated,
                "by_type_genes": by_type_genes_x,
                "rec_by_type": rec_by_type_x,
                "npub_by_pair": npub_by_pair_x,
                "tissues_by_type": other_resource.groupby("cell_type")["tissue_type"].nunique().to_dict(),
                "ortholog_pairs": len(gene_map),
            }

    ma = markers_all.copy()
    ma["group"] = ma["group"].astype(str)
    ma["gu"] = ma["feature"].astype(str).str.upper()

    ma["gene_key"] = ma["gu"]
    significant = significant_genes_by_cluster(markers_all)
    relative = cross_cluster_percentile(ma)
    relative.columns = relative.columns.astype(str)
    universe = pd.concat(
        [all_positive[["cell_type", "gu"]]]
        + [tab["translated"][["cell_type", "gu"]] for tab in cross_tables.values()],
        ignore_index=True,
    ).rename(columns={"cell_type": "cell", "gu": "gene_key"})
    claimant_specificity = marker_specificity(universe)
    specificity = pool.get("marker_specificity") or {}

    events: dict[str, list[dict]] = {}
    pending: list[dict] = []
    for cluster, state in pool["clusters"].items():
        genome = (pool.get("genome_de") or {}).get(cluster) or {}
        anchors = [str(a["gene"]).upper() for a in genome.get("anchor_genes") or []]
        anchor_set = set(anchors)
        if not anchor_set:
            continue
        unclaimed = [g for g in anchors if g not in native_genes]
        native_best = 0
        native_best_name = ""
        for entry in state["candidates"]:
            covered = sum(1 for m in entry["markers"] if m["polarity"] == "positive" and str(m["gene"]).upper() in anchor_set)
            if covered > native_best:
                native_best, native_best_name = covered, str(entry["cell_type"])
        admitted = {str(e["cell_type"]) for e in state["candidates"]}
        log = {
            "anchor_genes": anchors,
            "unclaimed_anchor_genes": unclaimed,
            "native_best": {"cell_type": native_best_name, "anchors_covered": native_best},
            "considered": [],
            "borrowed": [],
            "name_gated": [],
            "qualifier_gate": "removed",
        }
        events[cluster] = log
        if len(unclaimed) < BORROW_MIN_ANCHORS:
            continue
        unclaimed_set = set(unclaimed)
        hits_here = significant.get(cluster, frozenset())
        relative_col = relative[cluster] if cluster in relative.columns else None

        def _scored_claimant(name, genes, rec_set, npub_pairs, panel_set, donor):
            if name in admitted:
                return None
            hit = sorted(genes & unclaimed_set)
            if len(hit) < BORROW_MIN_ANCHORS:
                return None
            panel = sorted(g for g in panel_set if g in measured_genes)
            sig_hits = [g for g in panel if g in hits_here]
            if not panel or not sig_hits:
                return None
            w_hit = np.array([claimant_specificity.get(g, 0.0) for g in sig_hits], dtype=np.float64)
            w_panel = np.array([claimant_specificity.get(g, 0.0) for g in panel], dtype=np.float64)
            if relative_col is None:
                perc = np.zeros(len(panel), dtype=np.float64)
            else:
                perc = np.array(
                    [
                        (0.0 if pd.isna(v) else float(v))
                        for v in (
                            relative_col.at[g] if g in relative_col.index else 0.0
                            for g in panel
                        )
                    ],
                    dtype=np.float64,
                )
            marker_level = float(w_hit.sum() / math.sqrt(len(panel)))
            cluster_level = (
                float((w_panel * perc).sum() / max(float(w_panel.sum()), EPS))
                if w_panel.sum() > 0
                else float(perc.mean())
            )
            return {
                "cell_type": name,
                "anchors_covered": len(hit),
                "anchor_genes": hit,
                "recommended_anchors": sorted(set(hit) & rec_set.get(name, set())),
                "anchor_n_pub_sum": int(sum(npub_pairs.get((name, g), 0) for g in hit)),
                "donor_species": donor,
                "marker_level": marker_level,
                "cluster_level": cluster_level,
                "panel_measured": len(panel),
                "significant_hits": len(sig_hits),
            }

        scored = []
        for name, genes in by_type_genes.items():
            row = _scored_claimant(name, genes, rec_by_type, npub_by_pair,
                                   panel_genes_by_type.get(name, set()), None)
            if row is not None:
                scored.append(row)
        for other, tab in cross_tables.items():
            for name, genes in tab["by_type_genes"].items():
                row = _scored_claimant(name, genes, tab["rec_by_type"], tab["npub_by_pair"],
                                       genes, other)
                if row is not None:
                    scored.append(row)
        if not scored:
            continue
        rank_m = _percentile_rank(np.array([r["marker_level"] for r in scored], dtype=np.float64))
        rank_c = _percentile_rank(np.array([r["cluster_level"] for r in scored], dtype=np.float64))
        for row, rm, rc in zip(scored, rank_m, rank_c):
            row["borrow_score"] = float(math.sqrt(float(rm) * float(rc)))
        scored.sort(key=lambda r: (-r["borrow_score"], -r["cluster_level"], -r["marker_level"], r["cell_type"]))
        for position, row in enumerate(scored, start=1):
            row["borrow_rank"] = position
        pending.append({"cluster": cluster, "state": state, "log": log, "scored": scored,
                        "native_best": native_best, "native_best_name": native_best_name})

    screened: dict[str, dict] = {}
    if screen is not None and pending:
        names = sorted({str(item["cell_type"]) for p in pending for item in p["scored"]})
        screened = dict(screen(names) or {})

    for p in pending:
        cluster, state, log, scored = p["cluster"], p["state"], p["log"], p["scored"]
        native_best, native_best_name = p["native_best"], p["native_best_name"]
        for item in scored:
            name = item["cell_type"]
            hit = item["anchor_genes"]
            donor_species = item["donor_species"]
            source_frame = resource if donor_species is None else cross_tables[donor_species]["translated"]
            source_tissues_by_type = tissues_by_type if donor_species is None else cross_tables[donor_species]["tissues_by_type"]
            item["tissue_contexts"] = int(source_tissues_by_type.get(name, 0))
            verdict = screened.get(name) or {}
            if str(verdict.get("verdict") or "") == SCREEN_FOREIGN:
                item["rejected"] = SCREEN_REJECTED_PREFIX + str(verdict.get("reason") or "").strip()
                log["considered"].append(item)
                continue
            if len(log["borrowed"]) >= BORROW_MAX_PER_CLUSTER:
                item["rejected"] = f"cap of {BORROW_MAX_PER_CLUSTER} borrowed candidates reached"
                log["considered"].append(item)
                continue
            rows = source_frame[source_frame["cell_type"] == name]
            meta, pos_genes, neg_genes, tissues = _panel_meta(rows)
            here = ma[ma["group"] == cluster]
            stats = {}
            native_sym = dict(pool["native_gene"])
            for r in here.itertuples(index=False):
                stats[(cluster, str(r.gu))] = {
                    "pct_in": ap._num(r.pct_in, ap.PCT_DECIMALS), "pct_out": ap._num(r.pct_out, ap.PCT_DECIMALS),
                    "avg_log2FC": ap._num(r.avg_log2FC, ap.LOGFC_DECIMALS),
                    "auc": ap._num(r.auc, ap.AUC_DECIMALS), "adjusted_p_value": ap._num(r.padj, 6),
                }
                native_sym.setdefault(str(r.gu), str(r.feature))
            stats.update({k: v for k, v in pool["de_stats"].items() if k[0] == cluster})
            markers = ap._panel_rows(name, [g for g in pos_genes if g in measured_genes],
                                     [g for g in neg_genes if g in measured_genes],
                                     cluster, stats, meta, native_sym, specificity)
            in_table = {str(m["gene"]).upper() for m in markers}
            unmeasured = sorted(g for g in pos_genes if g not in measured_genes)
            flat = sorted((g for g in pos_genes if g in measured_genes and g not in in_table),
                          key=lambda g: (-int(meta[(name, g)]["n_pub"]), g))
            entry = {
                "cell_type": name,
                "retrieval_rank": len(state["candidates"]) + 1,
                "markers": markers,
                "exclusion_sources": ap._exclusion_sources(name, markers, sources),
                "unmeasured_curated_genes": {"count": len(unmeasured),
                                             "genes": [native_sym.get(g, g) for g in unmeasured][:ap.UNMEASURED_SHOWN]},
                "flat_or_undetected_curated_genes": {"count": len(flat),
                                                     "genes": [native_sym.get(g, g) for g in flat][:ap.UNMEASURED_SHOWN]},
                "program": {"median_in": None, "median_out": None},
                "borrowed_context": {
                    "source_tissues": sorted({t for g in hit for t in tissues.get(g, set())}),
                    "tissue_contexts_documenting_type": int(source_tissues_by_type.get(name, 0)),
                    "donor_species": donor_species,
                    "anchor_genes": hit,
                    "recommended_anchors": item["recommended_anchors"],
                    "native_best": {"cell_type": native_best_name, "anchors_covered": native_best},
                },
            }
            rec_set = rec_by_type if donor_species is None else cross_tables[donor_species]["rec_by_type"]
            for m in markers:
                m["recommended"] = bool(meta.get((name, str(m["gene"]).upper()), {}).get("n_pub")) and \
                    str(m["gene"]).upper() in rec_set.get(name, set())
            own_sig = [m for m in markers if m["polarity"] == "positive" and ap._carried(m)]
            if len(own_sig) < BORROW_OWN_SIGNIFICANT_MIN:
                item["rejected"] = (
                    f"fewer than {BORROW_OWN_SIGNIFICANT_MIN} of its own positive markers "
                    f"carried here ({len(own_sig)} of "
                    f"{len([m for m in markers if m['polarity'] == 'positive'])} measured)"
                )
                log["considered"].append(item)
                continue
            item["own_markers_significant"] = [m["gene"] for m in own_sig]
            item["panel_measured"] = len(markers)
            log["borrowed"].append(item)
            state["candidates"].append(entry)
            for m in markers:
                gene = str(m["gene"]).upper()
                carriers = pool["gene_carriers"].setdefault(gene, [])
                if not any(c[0] == name for c in carriers):
                    carriers.append((name, int(m["n_pub"]), str(m["tier"]), m["polarity"]))
                if gene not in pool["gene_clusters"]:
                    rows_g = ma[ma["gu"] == gene]
                    pool["gene_clusters"][gene] = [(str(r.group), ap._num(r.pct_in, 1), ap._num(r.pct_out, 1))
                                                   for r in rows_g.itertuples(index=False)]
                pool["native_gene"].setdefault(gene, native_sym.get(gene, gene))
    return events


def borrowed_names(events: dict[str, Any]) -> list[str]:
    names = set()
    for log in events.values():
        for item in log.get("borrowed") or []:
            names.add(str(item["cell_type"]))
    return sorted(names)


def borrowed_donor_species(events: dict[str, Any]) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for log in events.values():
        for item in log.get("borrowed") or []:
            donor = item.get("donor_species")
            if donor:
                out.setdefault(str(item["cell_type"]), set()).add(str(donor))
    return out
