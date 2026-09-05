#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Retrieve the free-text cell-type candidates a cluster could plausibly be.

Identity is the curated free-text cell-type string for a species. No ontology is read,
resolved, expanded or compared here: two curated names are two candidates, and whether
they denote the same biology is a question for the annotating agent, not for a
merge rule.

A candidate enters the pool for a cluster when enough of its measured positive markers
are significantly up-regulated there -- see `admit_by_hits` for what "enough" means and
why. The pool is then ordered by three scores, each answering a different "relative to
what" question about the candidate's whole measured panel:

  M  marker level       how much specific curated marker evidence for this identity is
                        actually present, length-normalized so a large panel neither wins
                        by size nor a two-gene panel by scarcity
  C  cluster level      does this panel single out THIS cluster among all clusters
  S  single-cell level  do the cells of this cluster carry the panel's program more than
                        the cells outside it

The three are percentile-ranked within the cluster and combined by their geometric mean
into one retrieval order. That order decides which candidates are shown to the agents; it
is not a confidence, not a probability, and not a biological conclusion.

Negative markers never enter retrieval or any of the three scores. Their measured
behaviour is collected here only so it can be shown as evidence.
"""

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
from .marker_database import MarkerDatabase

_R = DEFAULTS["retrieval"]
TOP_CANDIDATES = int(_R["top_candidates"])
MIN_HITS = int(_R["min_significant_hits"])
MIN_POOL_FLOOR = int(_R["min_pool_floor"])
EPS = float(_R["epsilon"])
RANK_DECIMALS = int(_R["rank_quantization_decimals"])
if MIN_HITS < 1 or MIN_POOL_FLOOR < 1 or TOP_CANDIDATES < 1:
    raise ValueError("retrieval admission controls must be positive integers")

# A cluster with no candidate whose panel intersects its significant DE has nothing to
# retrieve. Its cells are kept so coverage stays honest, and the annotation stage emits
# an explicit Unknown for it.
UNSUPPORTED = "unsupported_empty_candidate_pool"


def marker_specificity(positive_panel: pd.DataFrame) -> dict[str, float]:
    """m_g in [0, 1]: how few of this context's cell types claim gene g as a marker.

    Computed over the free-text cell-type universe of the context, so it measures the
    curated database's own discriminative power for the gene and nothing about the data.
    Publication counts and evidence tiers deliberately stay out of it; they are shown to
    the agents as separate evidence rather than folded into one number.
    """
    total = positive_panel["cell"].nunique()
    holders = positive_panel.groupby("gene_key")["cell"].nunique()
    scale = np.log(total + 1.0)
    return {
        str(gene): float(np.log((total + 1.0) / (count + 1.0)) / scale)
        for gene, count in holders.items()
    }


def cross_cluster_percentile(de: pd.DataFrame) -> pd.DataFrame:
    """Per gene, the average-rank percentile of its avg_log2FC across the clusters.

    1.0 means no cluster separates this gene more strongly than this one. Computed over
    every cluster the DE table tested, so it is a genuine between-cluster contrast rather
    than a restatement of the gene's own effect size.
    """
    wide = de.pivot_table(
        index="gene_key", columns="group", values="avg_log2FC", aggfunc="first"
    )
    return wide.rank(axis=1, pct=True, method="average")


def _expression_percentile(norm: sp.csc_matrix) -> tuple[sp.csc_matrix, np.ndarray]:
    """Per-gene expression percentile over ALL cells, in sparse form.

    Returns `(delta, baseline)` where `baseline[j]` is the percentile every non-detecting
    cell gets for gene j and `delta` holds `percentile - baseline` on the detected
    entries only. A program score is then `baseline . w + delta @ w`, which is exactly the
    dense average-rank ECDF without ever building a dense cells-by-genes matrix.

    The expression is the log-normalized value as measured. There is no neighbourhood
    smoothing: a smoothed matrix would make the single-cell score a function of the kNN
    graph, and the whole point of the third axis is that it reads individual cells.
    """
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
    """Cells-by-candidates marker-program score, one column per candidate.

    Kept in double precision: the R arm accumulates in double, and a single-precision
    sum here left the two arms differing by ~1e-5 on this axis -- small, but the axis is
    percentile-ranked without further rounding, so a near-tie could resolve differently
    in the two arms for no reason a reader could act on.
    """
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
    """Average-rank percentile within one cluster's candidate pool."""
    if len(values) == 1:
        return np.ones(1, dtype=np.float64)
    quantized = np.round(values, RANK_DECIMALS)
    return rankdata(quantized, method="average") / len(quantized)


def admit_by_hits(frame: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Which candidates the retrieval order is built from, and the gate that let them in.

    A candidate needs `min_significant_hits` of its OWN positive markers significantly
    up-regulated here. One or two hits is overwhelmingly a retrieval false positive:
    measured over the benchmark, 20-30% of one-hit candidates survive their own
    independent review against 79-93% of candidates carrying five or more. Those
    candidates also crowd out the real ones, because the shortlist the review tier can
    afford to read is short.

    The requirement is relaxed one step at a time where enforcing it would leave fewer
    than `min_pool_floor` candidates to choose between. A thin curated menu is a fact
    about the database's coverage of that organ -- three of the benchmark's contexts hold
    under 110 eligible cell types -- not evidence that the cluster has no identity, and a
    gate that empties a pool would turn a coverage gap into a confident Unknown.

    The threshold actually applied travels with the cluster, so a reader never has to
    infer whether a pool was gated at three hits or relaxed to one.
    """
    threshold = MIN_HITS
    while threshold > 1 and int((frame["hits"] >= threshold).sum()) < MIN_POOL_FLOOR:
        threshold -= 1
    admitted = frame[frame["hits"] >= threshold].reset_index(drop=True)
    return admitted, threshold


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
    de["feature"] = de["feature"].astype(str)  # matrix-native symbol, no case mutation
    de["group"] = de["group"].astype(str)
    de["gene_key"] = de["feature"].str.upper()
    species, tissue, disease = prep["species"], prep["tissue"], prep["disease"]
    measured = set(prep["genes"])
    menu_genes = list(prep["menu_genes"])
    column_of = {str(gene).upper(): index for index, gene in enumerate(menu_genes)}
    norm = sp.csc_matrix(sp.load_npz(os.path.join(CACHE, f"{tag}_norm.npz")))
    clusters = sorted(meta["cluster"].unique())
    print(
        f"==== retrieve(py) {tag}: {len(meta)} cells, {len(clusters)} clusters "
        f"({species}/{tissue}/{disease}) cross_species={cross_species} ===="
    )

    # ---- curated candidate universe for this context ------------------------
    database = MarkerDatabase()
    context = database.context_subset(
        species, tissue, disease, cross_species=cross_species
    )
    panel = database.panel(context)
    positive, eligible = database.eligible_candidates(
        panel, measured, corrob_only=CORROBORATION_ONLY
    )
    positive = positive.copy()
    # `gene_key` stays uppercase: by the time it is built, `g` IS a measured symbol
    # (`eligible_candidates` matches exactly), and measured symbols are unique
    # case-insensitively because preprocessing summed case duplicates. So the key is an
    # injective relabelling within one namespace, not a second crossing.
    positive["gene_key"] = positive["g"].astype(str).str.upper()
    specificity = marker_specificity(positive)
    # Negatives cross the same boundary and are gated the same way. They were only ever
    # displayed for genes the DE table carries, so restricting them here changes nothing
    # except that a curated symbol differing from the measured one by case no longer
    # reaches a panel the positive of the same name was refused from.
    negatives = panel[
        (panel["marker_polarity"] == "negative") & (panel["g"].astype(str).isin(measured))
    ].copy()
    negatives["gene_key"] = negatives["g"].astype(str).str.upper()
    negatives = negatives[negatives["cell"].isin(set(eligible))]
    print(f"  eligible free-text candidates: {len(eligible)}")

    # A candidate's panel is its measured positive markers. Every score reads the whole
    # panel; the significant subset only decides admission and what is shown as a hit.
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
    # Curated reliability metadata for every marker of every eligible candidate, carried
    # forward so the annotation stage never has to re-read the marker resource.
    panel_records = (
        pd.concat(
            [
                positive.assign(marker_polarity="positive"),
                negatives.assign(marker_polarity="negative"),
            ],
            ignore_index=True,
        )[["cell", "g", "gene_key", "marker_polarity", "n_pub", "tier"]]
        .rename(columns={"cell": "candidate", "g": "gene"})
        .drop_duplicates(["candidate", "gene_key", "marker_polarity"])
        .reset_index(drop=True)
    )

    significant = evidence_gate.significant_genes_by_cluster(de)
    relative = cross_cluster_percentile(de)

    print("  per-cell marker-program scores (unsmoothed expression percentiles) ...")
    delta, baseline = _expression_percentile(norm)
    program, program_names = _program_scores(
        delta, baseline, column_of, measured_panels, specificity
    )
    program_index = {name: position for position, name in enumerate(program_names)}
    program_rank = np.empty(program.shape, dtype=np.float64)
    for position in range(program.shape[1]):
        # Quantize before ranking. Two cells with the same expression across the panel
        # should tie, and whether their sums land bit-identically depends on accumulation
        # order; without this the two arms differ by exactly one tie's worth of rank.
        program_rank[:, position] = rankdata(
            np.round(program[:, position], RANK_DECIMALS), method="average"
        )

    cluster_of = meta["cluster"].to_numpy()
    rows = []
    per_cluster = {}
    for cluster in clusters:
        inside = cluster_of == cluster
        n_in = int(inside.sum())
        n_out = int(len(meta) - n_in)
        hits_here = significant.get(cluster)
        if hits_here is None:
            raise ValueError(f"{tag}: cluster {cluster} is absent from the DE table")
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
        if not records:
            print(f"  cluster {cluster}: empty candidate pool -> {UNSUPPORTED}")
            per_cluster[cluster] = dict(
                n_cells=n_in, status=UNSUPPORTED, candidates=[]
            )
            continue
        # The admission gate runs before the three axes are percentile-ranked, because a
        # percentile is a position within the pool: ranking first and filtering after
        # would leave every candidate's score describing a pool it is no longer in.
        pool_before_gate = len(records)
        frame, hits_threshold = admit_by_hits(pd.DataFrame(records))
        # A single-cluster dataset has no out-group, so the single-cell axis is undefined
        # rather than perfect; it then contributes nothing and the other two decide.
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
        rows.append(frame)
        selected = frame.head(top_candidates)["candidate"].tolist()
        per_cluster[cluster] = dict(
            n_cells=n_in,
            status="pool",
            candidates=selected,
            pool_size=int(len(frame)),
            pool_size_before_hits_gate=int(pool_before_gate),
            hits_threshold_applied=int(hits_threshold),
        )
        head = ", ".join(f"{name[:28]}" for name in selected[:3])
        relaxed = "" if hits_threshold == MIN_HITS else f" (gate relaxed to {hits_threshold})"
        print(
            f"  cluster {cluster:<3} n={n_in:<6} pool={len(frame):<4} "
            f"of {pool_before_gate:<4}{relaxed} top: {head}"
        )

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
            ]
        )
    )
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
        candidate_cl_id=database.candidate_cl_id(context),
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
