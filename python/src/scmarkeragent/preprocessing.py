#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepare raw single-cell counts for production annotation.

The stage performs quality control, normalization, native Leiden partitioning and
differential expression. Clustering is always recomputed here: an annotation run
never inherits a partition from the input object, so the Python and R arms enter
candidate scoring from the same standard workflow.
"""

import os
import argparse
import gc
import pickle
from pathlib import Path
from typing import Any
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import rankdata, norm as _norm

from .config import CACHE_DIR as CACHE
from .config import (
    QC_MIN_GENES,
    QC_MIN_CELLS,
    QC_MAX_MT,
    COMPUTE_UMAP,
    LEIDEN_RESOLUTION,
    DEFAULTS,
)
from .marker_database import MarkerDatabase

# The counts matrix is never guessed. A run states where its raw counts live, and the
# declaration is echoed into the run configuration, so a normalized matrix can never be
# silently annotated as if it were counts.
COUNTS_SOURCE_HELP = (
    "explicit location of the raw counts inside the .h5ad: "
    "'X', 'raw.X', or 'layers/<name>'"
)


# ============================ presto::wilcoxauc port ==========================
def wilcoxauc(X, genes, labels, gene_chunk=200):
    """One-vs-rest Wilcoxon AUC DE, port of presto::wilcoxauc.

    X: cells x genes LOG-NORMALIZED data, dense or scipy sparse.
    genes: list of gene symbols (columns of X).
    labels: array-like cluster labels (n_cells).
    Returns long DataFrame: feature, group, avgExpr, logFC, avg_log2FC, auc, pct_in,
    pct_out, scaled_mean, pval, padj.  auc = U1/(n1*n2) with average-tie ranks over all
    cells per gene; pval is the two-sided Mann-Whitney normal approximation with exact
    tie-corrected rank variance (== presto/scipy); padj is Benjamini-Hochberg applied
    WITHIN each cluster (Seurat FindMarkers per-cluster adjustment).

    `logFC` is the difference of mean log-normalized expression (the presto column) and
    is kept for audit only. `avg_log2FC` is the Seurat-style statistic the significance
    gate reads: de-log1p first, compare arithmetic means, add a pseudocount of one, take
    log2. Computing it here, on the same cells and in the same pass as the rest of the
    DE table, is what keeps the gate, the report and both arms on one number.

    Every statistic except the BH step is per gene, so genes are processed in chunks and
    only the chunk is ever dense. That bounds peak memory by `n_cells x gene_chunk` rather
    than by the matrix, which is what lets the same function serve both the small
    menu-restricted DE and the genome-wide FindAllMarkers-equivalent pass. BH is applied
    once at the end, over whichever gene universe was passed in.
    """
    sparse_input = sp.issparse(X)
    Xc = X.tocsc() if sparse_input else np.asarray(X, dtype=np.float32)
    n_cells, G = Xc.shape
    labels = np.asarray(labels).astype(str)
    groups = sorted(pd.unique(labels))
    K = len(groups)
    gi = {g: i for i, g in enumerate(groups)}
    lab_i = np.fromiter((gi[label] for label in labels), dtype=np.int64, count=n_cells)
    onehot = sp.csr_matrix(
        (np.ones(n_cells, dtype=np.float64), (np.arange(n_cells), lab_i)),
        shape=(n_cells, K),
    )
    onehotT = onehot.T.tocsr()
    n1 = np.asarray(onehot.sum(axis=0)).ravel()  # cells per group
    n2 = n_cells - n1

    mean_in = np.empty((K, G), dtype=np.float64)
    logFC = np.empty((K, G), dtype=np.float64)
    avg_log2FC = np.empty((K, G), dtype=np.float64)
    pct_in = np.empty((K, G), dtype=np.float64)
    pct_out = np.empty((K, G), dtype=np.float64)
    auc = np.empty((K, G), dtype=np.float64)
    pval = np.empty((K, G), dtype=np.float64)
    half = (n1 * (n1 + 1.0) / 2.0)[:, None]
    denom = (n1 * n2)[:, None]
    meanU = (n1 * n2 / 2.0)[:, None]
    for s in range(0, G, gene_chunk):
        e = min(s + gene_chunk, G)
        block = Xc[:, s:e]
        B = np.asarray(block.todense(), dtype=np.float32) if sparse_input else block

        sum_in = onehotT.dot(B).astype(np.float64)
        grand_sum = B.sum(axis=0, dtype=np.float64)
        m_in = sum_in / n1[:, None]
        mean_in[:, s:e] = m_in
        logFC[:, s:e] = m_in - (grand_sum[None, :] - sum_in) / n2[:, None]

        linear = np.expm1(B.astype(np.float64))
        lin_in = onehotT.dot(linear)
        lin_grand = linear.sum(axis=0)
        lin_mean_in = lin_in / n1[:, None]
        lin_mean_out = (lin_grand[None, :] - lin_in) / n2[:, None]
        avg_log2FC[:, s:e] = np.log2((lin_mean_in + 1.0) / (lin_mean_out + 1.0))
        del linear

        posmat = (B > 0).astype(np.float64)
        sum_in_pos = onehotT.dot(posmat)
        grand_pos = posmat.sum(axis=0)
        pct_in[:, s:e] = 100.0 * sum_in_pos / n1[:, None]
        pct_out[:, s:e] = 100.0 * (grand_pos[None, :] - sum_in_pos) / n2[:, None]
        del posmat

        R = rankdata(B, axis=0)  # cells x chunk, avg ties
        sum_rank_in = onehotT.dot(R)  # K x chunk
        U1 = sum_rank_in - half
        auc[:, s:e] = U1 / denom
        # two-sided Mann-Whitney normal approx with exact tie-corrected rank variance
        # (mean rank == (n+1)/2 exactly; sigma_U^2 = n1*n2/(n-1) * popVar(rank)).
        var_rank = R.var(axis=0)  # population var per gene (tie-aware)
        sigmaU = np.sqrt(
            (n1[:, None] * n2[:, None] / max(n_cells - 1.0, 1.0)) * var_rank[None, :]
        )
        # presto::compute_pval applies the standard 0.5 continuity correction:
        # z_num = U - mean(U) - sign(U-mean(U))*0.5. Keep it explicitly so
        # Python and R pval/padj columns are numerically identical.
        z_num = U1 - meanU
        z_num = z_num - np.sign(z_num) * 0.5
        with np.errstate(divide="ignore", invalid="ignore"):
            z = np.where(
                sigmaU > 0,
                z_num / np.where(sigmaU > 0, sigmaU, 1.0),
                0.0,
            )
        pval[:, s:e] = 2.0 * _norm.sf(np.abs(z))
        del R, B

    # scaled_mean per feature = min-max of avgExpr across groups (preprocessing.R)
    mmin = mean_in.min(axis=0)
    mmax = mean_in.max(axis=0)
    rng = mmax - mmin
    scaled = np.where(
        rng > 0, (mean_in - mmin[None, :]) / np.where(rng > 0, rng, 1.0), 0.0
    )

    pval = np.where(np.isnan(pval), 1.0, np.clip(pval, 0.0, 1.0))
    ga = np.array(genes)
    feat = np.repeat(ga, K)  # G*K, feature-major
    grp = np.tile(np.array(groups), G)
    de = pd.DataFrame(
        {
            "feature": feat,
            "group": grp,
            "avgExpr": mean_in.T.ravel(),
            "logFC": logFC.T.ravel(),
            "avg_log2FC": avg_log2FC.T.ravel(),
            "auc": auc.T.ravel(),
            "pct_in": pct_in.T.ravel(),
            "pct_out": pct_out.T.ravel(),
            "scaled_mean": scaled.T.ravel(),
            "pval": pval.T.ravel(),
        }
    )
    de["padj"] = de.groupby("group")["pval"].transform(_bh_adjust)
    return de


def _bh_adjust(p):
    """Benjamini-Hochberg FDR within one group (== p.adjust(method='BH'))."""
    p = np.asarray(p, dtype=float)
    m = len(p)
    if m == 0:
        return p
    order = np.argsort(p)
    ranked = np.empty(m, dtype=int)
    ranked[order] = np.arange(m)
    q = p[order] * m / (np.arange(m) + 1.0)
    q = np.minimum.accumulate(q[::-1])[::-1]
    return np.clip(q[ranked], 0.0, 1.0)


def _summarize_obs(obs, column, cap=8):
    """Summarize one sample-metadata column, or "" when it is absent or empty.

    Kept alongside the reference-free cell metadata because sample context (which
    developmental stages the tissue was taken from) constrains which identities
    are admissible, while carrying no author cell-type label.
    """
    if column not in obs.columns:
        return ""
    values = sorted(
        {
            str(value).strip()
            for value in obs[column].dropna().astype(str).unique()
            if str(value).strip()
            and str(value).strip().lower() not in {"nan", "none", "na", "unknown"}
        }
    )
    if not values:
        return ""
    if len(values) <= cap:
        return " | ".join(values)
    return " | ".join(values[:cap]) + f" | plus {len(values) - cap} more"


def screen_all_markers(table):
    """Apply the Seurat FindAllMarkers reporting screen to a genome-wide DE table.

    Seurat screens genes before testing (`min.pct`, `logfc.threshold`) and then returns
    only what it tested. Here every gene is tested, so the per-cluster BH family is the
    whole transcriptome and does not depend on where the screen was set, and the screen
    only decides which rows are reported. The R arm applies the identical screen to its
    presto table, so both arms deliver the same rows.

    `return.thresh` is deliberately NOT applied: this table is the "all markers" half of
    a pair, and its companion is the subset passing the project's significance gate, so
    pre-filtering it by p-value here would collapse the two into one.
    """
    screen = DEFAULTS["preprocessing"]["all_markers"]
    keep = (
        np.maximum(table["pct_in"].values, table["pct_out"].values)
        >= 100.0 * float(screen["min_pct"])
    ) & (np.abs(table["avg_log2FC"].values) >= float(screen["logfc_threshold"]))
    reported = table.loc[keep].reset_index(drop=True)
    print(
        f"  all-gene DE: {len(reported)} reported rows "
        f"(min_pct={screen['min_pct']}, logfc_threshold={screen['logfc_threshold']})"
    )
    return reported


def _sample_symbols(symbols, cap=4):
    """A few symbols, to show what a namespace actually looks like in a failure."""
    shown = sorted(str(symbol) for symbol in symbols)[:cap]
    return ", ".join(shown) if shown else "nothing"


def menu_measured_genes(measured, menu_positive, menu_negative):
    """The curated menu intersected with what this object measured, by EXACT symbol.

    This is the one place a curated gene symbol meets an uploaded one, and the match is
    case-sensitive because a gene symbol IS its case: HGNC PECAM1 and MGI Pecam1 are two
    species' names for two different measurements. Folding case here matched a mouse
    matrix against the human menu -- 1,029 of 14,877 mouse symbols "hit" a Human/blood
    menu that way, against 4 real human/mouse homographs -- and the run then annotated
    confidently against markers it had never actually measured.

    Because the match is exact, the surviving symbol is simultaneously the curated symbol
    and the matrix's own, so the DE and scoring labels stay in the matrix's native
    namespace with nothing to translate. Returns (positive, negative-only), both sorted.
    """
    measured = {str(gene) for gene in measured}
    positive = sorted(measured.intersection(str(gene) for gene in menu_positive))
    negative = sorted(
        measured.intersection(str(gene) for gene in menu_negative) - set(positive)
    )
    return positive, negative


def menu_de_table(table, menu_genes):
    """Rows of the genome-wide DE table for the curated menu genes, unfiltered.

    The menu table is a slice of the genome-wide pass rather than a second Wilcoxon run,
    so every statistic in it is the one the transcriptome-wide pass produced and `padj`
    carries the genome-wide per-cluster BH family that the retrieval gate reads. Rows are
    ordered feature-major over `menu_genes`, so both arms lay the table out identically.
    """
    order = {str(gene): position for position, gene in enumerate(menu_genes)}
    subset = table[table["feature"].astype(str).isin(order)].copy()
    subset["_feature_order"] = subset["feature"].astype(str).map(order)
    subset = (
        subset.sort_values(["_feature_order", "group"], kind="stable")
        .drop(columns="_feature_order")
        .reset_index(drop=True)
    )
    return subset


# ============================ h5ad loading ===================================
def _select_counts(adata, counts_source, path):
    """Resolve the declared counts location to (matrix, var) or fail loudly."""
    source = str(counts_source).strip()
    if source == "X":
        return adata.X, adata.var
    if source in {"raw.X", "raw/X"}:
        if adata.raw is None:
            raise ValueError(f"{path}: counts source 'raw.X' requested but .raw is absent")
        return adata.raw.X, adata.raw.var
    if source.startswith("layers/"):
        name = source.split("/", 1)[1]
        if name not in adata.layers:
            available = ", ".join(sorted(adata.layers.keys())) or "<none>"
            raise ValueError(
                f"{path}: counts source {source!r} requested but layer {name!r} is "
                f"absent; available layers: {available}"
            )
        return adata.layers[name], adata.var
    raise ValueError(
        f"{path}: unsupported counts source {source!r}; use 'X', 'raw.X' or "
        "'layers/<name>'"
    )


def load_counts_h5ad(path, counts_source):
    """Return (counts CSC cells x genes, gene_symbols unique-summed, obs).

    `counts_source` names the matrix explicitly. Nothing is inferred: an input whose
    declared source does not hold nonnegative integer counts stops the run instead of
    being coerced, because silently annotating a normalized matrix as counts is the
    failure this contract exists to prevent.
    """
    import scanpy as sc

    a = sc.read_h5ad(path)
    Xc, var = _select_counts(a, counts_source, path)
    Xc = sp.csr_matrix(Xc)
    if Xc.data.size:
        values = np.asarray(Xc.data, dtype=float)
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(
                f"{path}: counts source {counts_source!r} contains "
                "non-finite/negative values"
            )
        if np.max(np.abs(values - np.rint(values))) > 1e-6:
            raise ValueError(
                f"{path}: counts source {counts_source!r} is not integer-like; "
                "point --counts-source at the raw counts matrix, or round the input "
                "once at its origin"
            )
    if "feature_name" in var.columns:
        sym = (
            var["feature_name"].astype(str).to_numpy()
        )  # NATIVE species nomenclature (NO case mutation)
    else:
        sym = np.array([str(s) for s in var.index])
    # sum duplicate gene symbols (Seurat rownames are unique). Dedupe on a
    # CASE-INSENSITIVE key so a native/upper duplicate of the SAME gene (e.g. a
    # source carrying both Pisd and PISD) is summed, but KEEP the native symbol for
    # storage -- gene casing is never mutated (species symbol nomenclature preserved).
    key = np.array([s.upper() for s in sym])
    ukey, inv = np.unique(key, return_inverse=True)
    if len(ukey) != len(sym):
        M = sp.csr_matrix(
            (np.ones(len(sym)), (np.arange(len(sym)), inv)), shape=(len(sym), len(ukey))
        )
        Xc = Xc.dot(M)  # cells x uniq (sum dups)
        rep: dict[str, str] = {}
        for s, k in zip(sym, key):
            rep.setdefault(k, s)  # native representative per key (first seen)
        sym = np.array([rep[k] for k in ukey])
    Xc = sp.csc_matrix(Xc)
    return Xc, sym, a.obs.copy()


# ============================ main prep ======================================
def preprocess(
    tag: str,
    h5ad: str | Path,
    species: str,
    tissue: str,
    disease: str | list[str],
    counts_source: str,
    res: float | None = None,
    min_genes: int = QC_MIN_GENES,
    min_cells: int = QC_MIN_CELLS,
    max_mt: float = QC_MAX_MT,
    random_state: int | None = None,
    cross_species: bool = False,
    compute_umap: bool = COMPUTE_UMAP,
) -> dict[str, Any]:
    Path(CACHE).mkdir(parents=True, exist_ok=True)
    import scanpy as sc

    prep_cfg = DEFAULTS["preprocessing"]
    random_state = int(
        prep_cfg["random_state"] if random_state is None else random_state
    )
    res = float(LEIDEN_RESOLUTION if res is None else res)
    print(
        f"==== prep(py): {tag}  ({species}/{tissue}/{disease}) "
        f"counts_source={counts_source} cross_species={cross_species} "
        f"random_state={random_state} resolution={res} umap={compute_umap} ===="
    )
    counts, genes, obs = load_counts_h5ad(h5ad, counts_source)
    n0 = counts.shape[0]
    # Preserve the source AnnData observation names as the auditable cell IDs.
    # Older caches used reset_index(drop=True), leaving only positional surrogates;
    # the validated source snapshot resolves those positions back through the source h5ad.
    obs = obs.copy()
    obs.index = obs.index.astype(str)
    A = sc.AnnData(
        X=counts, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="symbol"))
    )
    A.var_names_make_unique()
    # ---- QC -----------------------------------------------------------------
    sc.pp.filter_cells(A, min_genes=min_genes)
    sc.pp.filter_genes(A, min_cells=min_cells)
    A.var["mt"] = A.var_names.str.upper().str.startswith(
        "MT-"
    )  # CI: catches human MT- + mouse mt- + rat Mt- (== preprocessing.R toupper)
    sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], inplace=True, percent_top=None)
    if A.var["mt"].sum() > 0:
        A = A[A.obs["pct_counts_mt"] < max_mt].copy()
    print(
        f"  QC: {n0} -> {A.n_obs} cells, {A.n_vars} genes "
        f"(min_genes={min_genes}, max_mt={max_mt})"
    )
    # ---- normalize (CP10k + log1p) -----------------------------------------
    A.layers["counts"] = A.X.copy()
    sc.pp.normalize_total(A, target_sum=float(prep_cfg["normalization_scale_factor"]))
    sc.pp.log1p(A)  # A.X = log-normalized
    # ---- HVG -> scale -> PCA -> neighbors -> Leiden -> UMAP ------------------
    # seurat_v3 is the variance-stabilizing selection run on RAW counts, matching the R
    # arm's FindVariableFeatures(selection.method="vst"). The older "seurat" flavor is a
    # different (dispersion-binned, log-scale) selection, which gave the two arms
    # different HVGs and therefore different PCA, clustering and kNN neighbourhoods.
    # Selection only writes .var columns, so it runs on A directly. Copying the whole
    # AnnData first would duplicate both the log-normalized matrix and the raw counts
    # layer, which on a 700k-cell input is tens of GB for nothing.
    sc.pp.highly_variable_genes(
        A,
        n_top_genes=int(prep_cfg["hvg_n"]),
        flavor="seurat_v3",
        layer="counts",
    )
    Ac = A[:, A.var["highly_variable"]].copy()
    del A.layers["counts"]  # selection is done; the raw copy is not needed again
    gc.collect()
    sc.pp.scale(Ac, max_value=float(prep_cfg["scale_max_value"]))
    sc.tl.pca(Ac, n_comps=int(prep_cfg["pca_n"]), svd_solver="arpack")
    sc.pp.neighbors(
        Ac,
        n_neighbors=int(prep_cfg["neighbors_k"]),
        n_pcs=int(prep_cfg["pca_n"]),
    )
    # Clustering is always computed here. A partition carried by the input object was
    # produced by an unknown upstream workflow, so reusing it would make two datasets --
    # and the two arms -- incomparable; the annotation contract is one standard workflow.
    ccol = f"leiden_res{res:g}"  # %g, matching preprocessing.R's label exactly
    sc.tl.leiden(
        Ac,
        resolution=res,
        random_state=random_state,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    clusters = Ac.obs["leiden"].astype(str).to_numpy()
    print(f"  clusters: Leiden res={res} -> {len(set(clusters))} clusters")
    if len(set(clusters)) < 2:
        # Every later stage reads one-vs-rest differential expression, which has no
        # "rest" for a single cluster; continuing would hand the retrieval gate NaN
        # statistics and report an unresolved cluster as if the evidence had been weighed.
        raise ValueError(
            f"Leiden clustering at resolution {res:g} put every cell in one cluster, and "
            "one-vs-rest differential expression is undefined for a single cluster. "
            "Raise --clustering-resolution, or check that the input holds more than one "
            "population."
        )
    score_pca_n = int(prep_cfg["score_pca_n"])
    emb = Ac.obsm["X_pca"][:, :score_pca_n]
    if compute_umap:
        sc.tl.umap(Ac, random_state=random_state)
        umap = np.asarray(Ac.obsm["X_umap"], dtype=np.float32)
    else:
        umap = None  # plots skipped for this run
    del Ac

    # ---- menu gene universe -------------------------------------------------
    # The menu holds both polarities. Positive markers are what candidate retrieval and
    # every score read; negative markers are measured only so that their actual detection
    # fraction in a cluster can be shown as evidence, never scored or used as a veto.
    db = MarkerDatabase()
    menu_positive = set(
        db.menu_genes(
            species, tissue, disease, cross_species=cross_species, polarity="positive"
        )
    )
    menu_negative = set(
        db.menu_genes(
            species, tissue, disease, cross_species=cross_species, polarity="negative"
        )
    )
    measured = {str(gene) for gene in A.var_names}
    positive_genes, negative_genes = menu_measured_genes(
        measured, menu_positive, menu_negative
    )
    de_genes = sorted(set(positive_genes) | set(negative_genes))
    if not positive_genes:
        # The symbols from both sides, so the reader can see WHY nothing matched instead
        # of being told to re-check three dropdowns. Nothing here infers what species the
        # object is: it reports what was measured, what is curated, and what to do.
        raise ValueError(
            f"{tag}: not one of the {len(menu_positive)} curated positive markers for "
            f"{species}/{tissue}/{disease} is among this object's {len(measured)} "
            "measured genes. Gene symbols are matched exactly, case included -- this "
            f"object measures {_sample_symbols(measured)}; the menu curates "
            f"{_sample_symbols(menu_positive)}. If those are the same genes written "
            "differently, this object is not written in the nomenclature the selected "
            "species is curated under: select the species these symbols belong to, or "
            "rewrite the symbols in the selected species' nomenclature. Otherwise the "
            "tissue or disease does not describe this object."
        )
    print(
        f"  menu genes: positive={len(menu_positive)} negative={len(menu_negative)}, "
        f"measured={len(measured)}, menu-measured={len(de_genes)} "
        f"({len(positive_genes)} positive + {len(negative_genes)} negative-only)"
    )

    # ---- per-cell expression matrix (ALL cells x measured menu genes) -------
    norm = sp.csc_matrix(A[:, de_genes].X)

    # Production metadata is reference-free by contract. Benchmark labels are
    # joined later, outside the annotation package, using stable cell IDs.
    meta = pd.DataFrame(
        {
            "cell": A.obs_names.to_numpy(),
            "cluster": clusters,
        }
    )

    development_stage = _summarize_obs(A.obs, "development_stage")

    # ---- one genome-wide one-vs-rest DE pass --------------------------------
    # Every gene is tested in a single pass, so the per-cluster BH family is the whole
    # transcriptome and the menu table is a slice of that same pass rather than a second
    # Wilcoxon run with its own, smaller family. Every cell is tested: subsampling made
    # the DE table -- and therefore the retrieval gate -- depend on a random draw.
    # Everything else has been read out of the AnnData by now, so it is released before
    # the pass rather than held alongside it.
    all_gene_names = list(A.var_names)
    all_gene_matrix = A.X.tocsc() if sp.issparse(A.X) else A.X
    del A
    gc.collect()
    print(
        f"  genome-wide DE: {len(all_gene_names)} genes x {len(set(clusters))} clusters ..."
    )
    full_de = wilcoxauc(all_gene_matrix, all_gene_names, clusters)
    del all_gene_matrix
    gc.collect()
    de = menu_de_table(full_de, de_genes)
    all_markers = screen_all_markers(full_de)
    del full_de
    gc.collect()

    payload = dict(
        de=de,
        meta=meta,
        genes=sorted(measured),
        menu_genes=de_genes,
        menu_positive_genes=positive_genes,
        menu_negative_genes=negative_genes,
        cells=meta["cell"].to_numpy(),
        cluster_col=ccol,
        umap=umap,
        leiden_resolution=float(res),
        counts_source=str(counts_source),
        source_path=os.path.abspath(h5ad),
        species=species,
        tissue=tissue,
        cross_species=bool(cross_species),
        disease=[disease] if isinstance(disease, str) else list(disease),
        development_stage=development_stage,
    )
    with open(os.path.join(CACHE, f"{tag}_de_meta.pkl"), "wb") as fh:
        pickle.dump(payload, fh, protocol=4)
    sp.save_npz(os.path.join(CACHE, f"{tag}_norm.npz"), norm)
    np.save(os.path.join(CACHE, f"{tag}_emb.npy"), np.asarray(emb, dtype=np.float32))
    with open(os.path.join(CACHE, f"{tag}_markers_all.pkl"), "wb") as fh:
        pickle.dump(all_markers, fh, protocol=4)
    print(
        f"  [done] -> {tag}_de_meta.pkl ({len(de)} DE rows) + {tag}_norm.npz "
        f"({norm.shape}) + {tag}_emb.npy ({emb.shape}) + {tag}_markers_all.pkl "
        f"({len(all_markers)} rows)"
    )
    return payload


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tag")
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--tissue", required=True)
    ap.add_argument("--disease", default="Normal")
    ap.add_argument("--counts-source", required=True, help=COUNTS_SOURCE_HELP)
    ap.add_argument(
        "--clustering-resolution",
        type=float,
        default=LEIDEN_RESOLUTION,
        help="Leiden resolution for the always-recomputed partition",
    )
    ap.add_argument("--cross-species", action="store_true")
    a = ap.parse_args()
    disease = a.disease.split("|") if "|" in a.disease else a.disease
    preprocess(
        a.tag,
        a.h5ad,
        a.species,
        a.tissue,
        disease,
        a.counts_source,
        res=a.clustering_resolution,
        cross_species=a.cross_species,
    )


if __name__ == "__main__":
    main()
