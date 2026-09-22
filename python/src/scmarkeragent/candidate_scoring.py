#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import os
import pickle
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import rankdata

from . import evidence_gate
from .config import (
    CACHE_DIR as CACHE,
    CORROBORATION_ONLY,
    CROSS_SPECIES,
    DEFAULTS,
)
from .marker_database import TISSUE_CONTEXT_SEP, MarkerDatabase, tissue_list

_R = DEFAULTS["retrieval"]
TOP_CANDIDATES = int(_R["top_candidates"])
MIN_HITS = int(_R["min_significant_hits"])
MIN_POOL_FLOOR = int(_R["min_pool_floor"])
EPS = float(_R["epsilon"])
RANK_DECIMALS = int(_R["rank_quantization_decimals"])
if MIN_HITS < 1 or MIN_POOL_FLOOR < 1 or TOP_CANDIDATES < 1:
    raise ValueError("retrieval admission controls must be positive integers")

UNSUPPORTED = "unsupported_empty_candidate_pool"


def marker_specificity(positive_panel: pd.DataFrame) -> dict[str, float]:
    total = positive_panel["cell"].nunique()
    holders = positive_panel.groupby("gene_key")["cell"].nunique()
    scale = np.log(total + 1.0)
    return {
        str(gene): float(np.log((total + 1.0) / (count + 1.0)) / scale)
        for gene, count in holders.items()
    }


def cross_cluster_percentile(de: pd.DataFrame) -> pd.DataFrame:
    wide = de.pivot_table(
        index="gene_key", columns="group", values="avg_log2FC", aggfunc="first"
    )
    return wide.rank(axis=1, pct=True, method="average")


def _expression_percentile(norm: sp.csc_matrix) -> tuple[sp.csc_matrix, np.ndarray]:
    n_cells = norm.shape[0]
    baseline = np.zeros(norm.shape[1], dtype=np.float64)
    delta = norm.copy().astype(np.float64)
    for column in range(norm.shape[1]):
        start, end = delta.indptr[column], delta.indptr[column + 1]
        detected = end - start
        zeros = n_cells - detected
        zero_percentile = ((zeros + 1.0) / 2.0) / n_cells
        baseline[column] = zero_percentile
        if detected:
            ranks = rankdata(delta.data[start:end], method="average") + zeros
            delta.data[start:end] = ranks / n_cells - zero_percentile
    return delta, baseline


def _program_scores(
    delta: sp.csc_matrix,
    baseline: np.ndarray,
    column_of: dict[str, int],
    panels: dict[str, list[str]],
    weights: dict[str, float],
) -> tuple[np.ndarray, list[str]]:
    names = sorted(panels)
    scores = np.zeros((delta.shape[0], len(names)), dtype=np.float64)
    for position, name in enumerate(names):
        genes = panels[name]
        columns = [column_of[gene] for gene in genes]
        weight = np.array([weights.get(gene, 0.0) for gene in genes], dtype=np.float64)
        if weight.sum() <= 0:
            weight = np.ones(len(genes), dtype=np.float64)
        selector = sp.csc_matrix(
            (weight, (columns, np.zeros(len(columns), dtype=np.int64))),
            shape=(delta.shape[1], 1),
        )
        contributed = np.asarray((delta @ selector).todense()).ravel()
        constant = float(baseline[columns] @ weight)
        scores[:, position] = (constant + contributed) / weight.sum()
    return scores, names


def _percentile_rank(values: np.ndarray) -> np.ndarray:
    if len(values) == 1:
        return np.ones(1, dtype=np.float64)
    quantized = np.round(values, RANK_DECIMALS)
    return rankdata(quantized, method="average") / len(quantized)


def admit_by_hits(frame: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    threshold = MIN_HITS
    while threshold > 1 and int((frame["hits"] >= threshold).sum()) < MIN_POOL_FLOOR:
        threshold -= 1
    saturated = (frame["hits"] == frame["panel_size"]) & (frame["panel_size"] < MIN_HITS)
    admitted = frame[(frame["hits"] >= threshold) | saturated].reset_index(drop=True)
    return admitted, threshold


EXTRA_CONTEXT_TOP = int(_R.get("extra_context_top_candidates", 5))
if EXTRA_CONTEXT_TOP < 0:
    raise ValueError("retrieval.extra_context_top_candidates must be >= 0")
RETRIEVAL_SEMANTICS = "primary_preserved_v1"


def _score_context(
    database: MarkerDatabase,
    species: str,
    tissue_name: str,
    disease,
    cross_species: bool,
    measured: set,
    column_of: dict,
    delta: sp.csc_matrix,
    baseline: np.ndarray,
    significant: dict,
    relative: pd.DataFrame,
    meta: pd.DataFrame,
    clusters: list,
) -> dict:
    context = database.context_subset(
        species, tissue_name, disease, cross_species=cross_species
    )
    panel = database.panel(context, tissue_contexts=[tissue_name])
    positive, eligible = database.eligible_candidates(
        panel, measured, corrob_only=CORROBORATION_ONLY
    )
    positive = positive.copy()
    positive["gene_key"] = positive["g"].astype(str).str.upper()
    specificity = marker_specificity(positive)
    negatives = panel[
        (panel["marker_polarity"] == "negative") & (panel["g"].astype(str).isin(measured))
    ].copy()
    negatives["gene_key"] = negatives["g"].astype(str).str.upper()
    negatives = negatives[negatives["cell"].isin(set(eligible))]
    print(f"  [{tissue_name}] eligible free-text candidates: {len(eligible)}")

    measured_panels = {
        str(name): sorted({gene for gene in group["gene_key"] if gene in column_of})
        for name, group in positive.groupby("cell")
    }
    measured_panels = {
        name: genes for name, genes in measured_panels.items() if genes
    }
    negative_panels = {
        str(name): sorted({gene for gene in group["gene_key"]})
        for name, group in negatives.groupby("cell")
    }
    panel_records = (
        pd.concat(
            [
                positive.assign(marker_polarity="positive"),
                negatives.assign(marker_polarity="negative"),
            ],
            ignore_index=True,
        )[["cell", "g", "gene_key", "marker_polarity", "n_pub", "tier", "tissue_context"]]
        .rename(columns={"cell": "candidate", "g": "gene"})
        .drop_duplicates(["candidate", "gene_key", "marker_polarity"])
        .reset_index(drop=True)
    )

    program, program_names = _program_scores(
        delta, baseline, column_of, measured_panels, specificity
    )
    program_index = {name: position for position, name in enumerate(program_names)}
    program_rank = np.empty(program.shape, dtype=np.float64)
    for position in range(program.shape[1]):
        program_rank[:, position] = rankdata(
            np.round(program[:, position], RANK_DECIMALS), method="average"
        )

    cluster_of = meta["cluster"].to_numpy()
    frames: dict[str, pd.DataFrame | None] = {}
    before_gate: dict[str, int] = {}
    thresholds: dict[str, int | None] = {}
    for cluster in clusters:
        inside = cluster_of == cluster
        n_in = int(inside.sum())
        n_out = int(len(meta) - n_in)
        hits_here = significant.get(cluster)
        if hits_here is None:
            raise ValueError(f"cluster {cluster} is absent from the DE table")
        records = []
        for name, genes in measured_panels.items():
            intersection = [gene for gene in genes if gene in hits_here]
            if not intersection:
                continue
            weight_hit = np.array(
                [specificity.get(gene, 0.0) for gene in intersection]
            )
            weight_panel = np.array([specificity.get(gene, 0.0) for gene in genes])
            percentile = np.array(
                [
                    float(relative.at[gene, cluster]) if gene in relative.index else 0.0
                    for gene in genes
                ]
            )
            position = program_index[name]
            rank_in = program_rank[inside, position]
            marker_level = float(weight_hit.sum() / np.sqrt(len(genes)))
            cluster_level = float(
                (weight_panel * percentile).sum() / max(weight_panel.sum(), EPS)
                if weight_panel.sum() > 0
                else percentile.mean()
            )
            single_cell_level = (
                float((rank_in.sum() - n_in * (n_in + 1.0) / 2.0) / (n_in * n_out))
                if n_out > 0
                else float("nan")
            )
            values = program[inside, position]
            outside_values = program[~inside, position]
            records.append(
                dict(
                    cluster=cluster,
                    candidate=name,
                    hits=len(intersection),
                    panel_size=len(genes),
                    marker_level=marker_level,
                    cluster_level=cluster_level,
                    single_cell_level=single_cell_level,
                    program_in_q25=float(np.quantile(values, 0.25)),
                    program_in_median=float(np.median(values)),
                    program_in_q75=float(np.quantile(values, 0.75)),
                    program_out_q25=(
                        float(np.quantile(outside_values, 0.25)) if n_out else float("nan")
                    ),
                    program_out_median=(
                        float(np.median(outside_values)) if n_out else float("nan")
                    ),
                    program_out_q75=(
                        float(np.quantile(outside_values, 0.75)) if n_out else float("nan")
                    ),
                )
            )
        before_gate[cluster] = len(records)
        if not records:
            frames[cluster] = None
            thresholds[cluster] = None
            continue
        frame, hits_threshold = admit_by_hits(pd.DataFrame(records))
        single = frame["single_cell_level"].to_numpy()
        if not np.isfinite(single).all():
            single = np.zeros(len(frame))
        frame["rank_marker_level"] = _percentile_rank(
            frame["marker_level"].to_numpy()
        )
        frame["rank_cluster_level"] = _percentile_rank(
            frame["cluster_level"].to_numpy()
        )
        frame["rank_single_cell_level"] = _percentile_rank(single)
        frame["retrieval_score"] = np.cbrt(
            frame["rank_marker_level"]
            * frame["rank_cluster_level"]
            * frame["rank_single_cell_level"]
        )
        frame = frame.sort_values(
            by=[
                "retrieval_score",
                "cluster_level",
                "single_cell_level",
                "marker_level",
                "candidate",
            ],
            ascending=[False, False, False, False, True],
            kind="stable",
        ).reset_index(drop=True)
        frame["retrieval_rank"] = np.arange(1, len(frame) + 1)
        frames[cluster] = frame
        thresholds[cluster] = int(hits_threshold)
    return dict(
        tissue=tissue_name,
        context=context,
        eligible=set(str(name) for name in eligible),
        specificity=specificity,
        measured_panels=measured_panels,
        negative_panels=negative_panels,
        panel_records=panel_records,
        candidate_cl_id=database.candidate_cl_id(context),
        frames=frames,
        before_gate=before_gate,
        thresholds=thresholds,
    )


def compute_candidate_scores(
    tag: str,
    top_candidates: int | None = None,
    cross_species: bool | None = None,
) -> pd.DataFrame:
    top_candidates = int(
        TOP_CANDIDATES if top_candidates is None else top_candidates
    )
    if cross_species is None:
        cross_species = CROSS_SPECIES
    with open(os.path.join(CACHE, f"{tag}_de_meta.pkl"), "rb") as handle:
        prep = pickle.load(handle)
    de = prep["de"].copy()
    meta = prep["meta"].copy()
    meta["cluster"] = meta["cluster"].astype(str)
    de["feature"] = de["feature"].astype(str)
    de["group"] = de["group"].astype(str)
    de["gene_key"] = de["feature"].str.upper()
    species, tissue, disease = prep["species"], prep["tissue"], prep["disease"]
    tissue_contexts = tissue_list(tissue)
    if not tissue_contexts:
        raise ValueError(f"{tag}: the preprocessing payload names no tissue context")
    primary = tissue_contexts[0]
    extras = tissue_contexts[1:]
    measured = set(prep["genes"])
    menu_genes = list(prep["menu_genes"])
    column_of = {str(gene).upper(): index for index, gene in enumerate(menu_genes)}
    norm = sp.csc_matrix(sp.load_npz(os.path.join(CACHE, f"{tag}_norm.npz")))
    clusters = sorted(meta["cluster"].unique())
    print(
        f"==== retrieve(py) {tag}: {len(meta)} cells, {len(clusters)} clusters "
        f"({species}/{tissue}/{disease}) cross_species={cross_species} ===="
    )
    if extras:
        print(f"  tissue contexts (operator list): primary {primary!r}; extra {extras} "
              f"(<= {EXTRA_CONTEXT_TOP} attributed candidates each; {RETRIEVAL_SEMANTICS})")

    significant = evidence_gate.significant_genes_by_cluster(de)
    relative = cross_cluster_percentile(de)
    print("  per-cell marker-program scores (unsmoothed expression percentiles) ...")
    delta, baseline = _expression_percentile(norm)

    database = MarkerDatabase()
    scored_by = {}
    for name in tissue_contexts:
        scored_by[name] = _score_context(
            database, species, name, disease, cross_species, measured, column_of,
            delta, baseline, significant, relative, meta, clusters,
        )
    P = scored_by[primary]

    attributed: dict[str, str] = {}
    for name in tissue_contexts:
        for candidate in sorted(scored_by[name]["eligible"]):
            attributed.setdefault(candidate, name)
    eligible_in = {
        candidate: [name for name in tissue_contexts if candidate in scored_by[name]["eligible"]]
        for candidate in attributed
    }

    def _label(candidate: str) -> str:
        return TISSUE_CONTEXT_SEP.join(eligible_in[str(candidate)])

    rows = []
    per_cluster = {}
    for cluster in clusters:
        n_in = int((meta["cluster"].to_numpy() == cluster).sum())
        frame = P["frames"][cluster]
        parts = []
        offset = 0
        if frame is None:
            print(f"  cluster {cluster}: empty candidate pool -> {UNSUPPORTED}")
            state = dict(n_cells=n_in, status=UNSUPPORTED, candidates=[])
        else:
            frame = frame.copy()
            frame["tissue_contexts"] = [_label(c) for c in frame["candidate"]]
            frame["retrieval_context"] = primary
            frame["context_rank"] = frame["retrieval_rank"]
            parts.append(frame)
            offset = len(frame)
            selected = frame.head(top_candidates)["candidate"].tolist()
            state = dict(
                n_cells=n_in,
                status="pool",
                candidates=selected,
                pool_size=int(len(frame)),
                pool_size_before_hits_gate=int(P["before_gate"][cluster]),
                hits_threshold_applied=int(P["thresholds"][cluster]),
            )
            head = ", ".join(f"{name[:28]}" for name in selected[:3])
            relaxed = "" if P["thresholds"][cluster] == MIN_HITS else f" (gate relaxed to {P['thresholds'][cluster]})"
            print(
                f"  cluster {cluster:<3} n={n_in:<6} pool={len(frame):<4} "
                f"of {P['before_gate'][cluster]:<4}{relaxed} top: {head}"
            )
        if extras:
            extra_state = {}
            merged = list(state["candidates"])
            for name in extras:
                E = scored_by[name]
                fe = E["frames"][cluster]
                chosen: list[str] = []
                n_attributed = 0
                if fe is not None:
                    own = fe[[attributed.get(str(c)) == name for c in fe["candidate"]]].copy()
                    n_attributed = len(own)
                    if len(own):
                        own["tissue_contexts"] = [_label(c) for c in own["candidate"]]
                        own["retrieval_context"] = name
                        own["context_rank"] = own["retrieval_rank"]
                        own["retrieval_rank"] = np.arange(offset + 1, offset + 1 + len(own))
                        parts.append(own)
                        offset += len(own)
                        chosen = own.head(EXTRA_CONTEXT_TOP)["candidate"].tolist()
                extra_state[name] = dict(
                    candidates=chosen,
                    attributed_pool_size=int(n_attributed),
                    context_pool_size=int(len(fe)) if fe is not None else 0,
                    context_pool_size_before_hits_gate=int(E["before_gate"][cluster]),
                    context_hits_threshold_applied=E["thresholds"][cluster],
                )
                merged.extend(chosen)
                if chosen:
                    print(f"  cluster {cluster:<3} + {name!r}: {len(chosen)} of {n_attributed} attributed "
                          f"(context pool {len(fe)} of {E['before_gate'][cluster]}): {', '.join(c[:28] for c in chosen)}")
            state["primary_status"] = state["status"]
            state["extra_contexts"] = extra_state
            state["candidates"] = merged
            if merged and state["status"] == UNSUPPORTED:
                state["status"] = "pool"
        per_cluster[cluster] = state
        if parts:
            rows.append(pd.concat(parts, ignore_index=True))

    scored = (
        pd.concat(rows, ignore_index=True)
        if rows
        else pd.DataFrame(
            columns=[
                "cluster",
                "candidate",
                "hits",
                "panel_size",
                "marker_level",
                "cluster_level",
                "single_cell_level",
                "retrieval_score",
                "retrieval_rank",
                "tissue_contexts",
                "retrieval_context",
                "context_rank",
            ]
        )
    )

    measured_panels = dict(P["measured_panels"])
    negative_panels = dict(P["negative_panels"])
    candidate_cl_id = dict(P["candidate_cl_id"])
    specificity = dict(P["specificity"])
    record_parts = [P["panel_records"]]
    for name in extras:
        E = scored_by[name]
        mine = {c for c, ctx in attributed.items() if ctx == name}
        measured_panels.update({c: g for c, g in E["measured_panels"].items() if c in mine})
        negative_panels.update({c: g for c, g in E["negative_panels"].items() if c in mine})
        candidate_cl_id.update({c: v for c, v in E["candidate_cl_id"].items() if c in mine and c not in candidate_cl_id})
        for gene, value in E["specificity"].items():
            specificity.setdefault(gene, value)
        pr = E["panel_records"]
        record_parts.append(pr[pr["candidate"].isin(mine)])
    panel_records = pd.concat(record_parts, ignore_index=True)
    keys_by_context = {
        name: set(zip(scored_by[name]["panel_records"]["candidate"], scored_by[name]["panel_records"]["gene_key"],
                      scored_by[name]["panel_records"]["marker_polarity"]))
        for name in tissue_contexts
    }
    panel_records["tissue_context"] = [
        TISSUE_CONTEXT_SEP.join(name for name in tissue_contexts if key in keys_by_context[name])
        for key in zip(panel_records["candidate"], panel_records["gene_key"], panel_records["marker_polarity"])
    ]
    payload = dict(
        context=dict(
            species=species,
            tissue=tissue,
            disease=list(disease) if not isinstance(disease, str) else [disease],
            development_stage=str(prep.get("development_stage") or ""),
            n_cells=int(len(meta)),
            n_clusters=len(clusters),
        ),
        scored=scored,
        clusters=per_cluster,
        marker_specificity=specificity,
        measured_panels=measured_panels,
        negative_panels=negative_panels,
        panel_records=panel_records,
        top_candidates=top_candidates,
        min_significant_hits=MIN_HITS,
        min_pool_floor=MIN_POOL_FLOOR,
        candidate_cl_id=candidate_cl_id,
        tissue_contexts=tissue_contexts,
        candidate_tissue_contexts={c: eligible_in[c] for c in sorted(attributed)},
        candidate_retrieval_context={c: attributed[c] for c in sorted(attributed)},
        eligible_by_context={name: sorted(scored_by[name]["eligible"]) for name in tissue_contexts},
        retrieval_semantics=RETRIEVAL_SEMANTICS,
        extra_context_top_candidates=EXTRA_CONTEXT_TOP,
    )
    with open(os.path.join(CACHE, f"{tag}_candidate_scoring.pkl"), "wb") as handle:
        pickle.dump(payload, handle, protocol=4)
    print(
        f"\n[done] {tag}: {len(scored)} scored candidate-cluster pairs across "
        f"{len(clusters)} clusters -> {tag}_candidate_scoring.pkl"
    )
    return scored


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit("usage: python -m scmarkeragent.candidate_scoring <tag>")
    compute_candidate_scores(
        sys.argv[1],
        top_candidates=int(sys.argv[2]) if len(sys.argv) > 2 else None,
    )
