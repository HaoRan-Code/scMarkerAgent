#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import gc
import pickle
import threading
from concurrent.futures import ThreadPoolExecutor
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
from .marker_database import MarkerDatabase, parse_tissue_arg

COUNTS_SOURCE_HELP = (
    "explicit location of the raw counts inside the .h5ad: "
    "'X', 'raw.X', or 'layers/<name>'"
)


def resolve_wilcox_threads(n_jobs=None) -> int:
    if n_jobs is not None:
        return max(1, int(n_jobs))
    raw = str(os.environ.get("SCMA_WILCOX_THREADS") or "").strip()
    if raw:
        return max(1, int(raw))
    return 1


def wilcoxauc(X, genes, labels, gene_chunk=200, n_jobs=None):
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
    n1 = np.asarray(onehot.sum(axis=0)).ravel()
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
    workers = resolve_wilcox_threads(n_jobs)
    slices = [(s, min(s + gene_chunk, G)) for s in range(0, G, gene_chunk)]
    n_slices = len(slices)
    workers = max(1, min(workers, n_slices or 1))

    def _one_chunk(se):
        s, e = se
        block = Xc[:, s:e]
        B = np.asarray(block.todense(), dtype=np.float32) if sparse_input else np.asarray(block, dtype=np.float32)

        sum_in = onehotT.dot(B).astype(np.float64)
        grand_sum = B.sum(axis=0, dtype=np.float64)
        m_in = sum_in / n1[:, None]
        logfc = m_in - (grand_sum[None, :] - sum_in) / n2[:, None]

        linear = np.expm1(B.astype(np.float64))
        lin_in = onehotT.dot(linear)
        lin_grand = linear.sum(axis=0)
        lin_mean_in = lin_in / n1[:, None]
        lin_mean_out = (lin_grand[None, :] - lin_in) / n2[:, None]
        avg = np.log2((lin_mean_in + 1.0) / (lin_mean_out + 1.0))
        del linear

        posmat = (B > 0).astype(np.float64)
        sum_in_pos = onehotT.dot(posmat)
        grand_pos = posmat.sum(axis=0)
        pin = 100.0 * sum_in_pos / n1[:, None]
        pout = 100.0 * (grand_pos[None, :] - sum_in_pos) / n2[:, None]
        del posmat

        R = rankdata(B, axis=0)
        sum_rank_in = onehotT.dot(R)
        U1 = sum_rank_in - half
        auc_s = U1 / denom
        var_rank = R.var(axis=0)
        sigmaU = np.sqrt(
            (n1[:, None] * n2[:, None] / max(n_cells - 1.0, 1.0)) * var_rank[None, :]
        )
        z_num = U1 - meanU
        z_num = z_num - np.sign(z_num) * 0.5
        with np.errstate(divide="ignore", invalid="ignore"):
            z = np.where(
                sigmaU > 0,
                z_num / np.where(sigmaU > 0, sigmaU, 1.0),
                0.0,
            )
        pval_s = 2.0 * _norm.sf(np.abs(z))
        del R, B
        return s, e, m_in, logfc, avg, pin, pout, auc_s, pval_s

    done = 0
    done_lock = threading.Lock()
    mark_every = 1 if n_slices <= 10 else max(1, n_slices // 10)

    def _one_chunk_logged(se):
        nonlocal done
        out = _one_chunk(se)
        with done_lock:
            done += 1
            if done == 1 or done == n_slices or done % mark_every == 0:
                print(f"  wilcoxauc chunks {done}/{n_slices} n_jobs={workers}", flush=True)
        return out

    if workers == 1:
        chunk_rows = [_one_chunk_logged(se) for se in slices]
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            chunk_rows = list(pool.map(_one_chunk_logged, slices))

    for s, e, m_in, logfc, avg, pin, pout, auc_s, pval_s in chunk_rows:
        mean_in[:, s:e] = m_in
        logFC[:, s:e] = logfc
        avg_log2FC[:, s:e] = avg
        pct_in[:, s:e] = pin
        pct_out[:, s:e] = pout
        auc[:, s:e] = auc_s
        pval[:, s:e] = pval_s

    mmin = mean_in.min(axis=0)
    mmax = mean_in.max(axis=0)
    rng = mmax - mmin
    scaled = np.where(
        rng > 0, (mean_in - mmin[None, :]) / np.where(rng > 0, rng, 1.0), 0.0
    )

    pval = np.where(np.isnan(pval), 1.0, np.clip(pval, 0.0, 1.0))
    ga = np.array(genes)
    feat = np.repeat(ga, K)
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
    shown = sorted(str(symbol) for symbol in symbols)[:cap]
    return ", ".join(shown) if shown else "nothing"


def menu_measured_genes(measured, menu_positive, menu_negative):
    measured = {str(gene) for gene in measured}
    positive = sorted(measured.intersection(str(gene) for gene in menu_positive))
    negative = sorted(
        measured.intersection(str(gene) for gene in menu_negative) - set(positive)
    )
    return positive, negative


def menu_de_table(table, menu_genes):
    order = {str(gene): position for position, gene in enumerate(menu_genes)}
    subset = table[table["feature"].astype(str).isin(order)].copy()
    subset["_feature_order"] = subset["feature"].astype(str).map(order)
    subset = (
        subset.sort_values(["_feature_order", "group"], kind="stable")
        .drop(columns="_feature_order")
        .reset_index(drop=True)
    )
    return subset


def _select_counts(adata, counts_source, path):
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
        )
    else:
        sym = np.array([str(s) for s in var.index])
    key = np.array([s.upper() for s in sym])
    ukey, inv = np.unique(key, return_inverse=True)
    if len(ukey) != len(sym):
        M = sp.csr_matrix(
            (np.ones(len(sym)), (np.arange(len(sym)), inv)), shape=(len(sym), len(ukey))
        )
        Xc = Xc.dot(M)
        rep: dict[str, str] = {}
        for s, k in zip(sym, key):
            rep.setdefault(k, s)
        sym = np.array([rep[k] for k in ukey])
    Xc = sp.csc_matrix(Xc)
    return Xc, sym, a.obs.copy()


def preprocess(
    tag: str,
    h5ad: str | Path,
    species: str,
    tissue: str | list[str],
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
    obs = obs.copy()
    obs.index = obs.index.astype(str)
    A = sc.AnnData(
        X=counts, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="symbol"))
    )
    A.var_names_make_unique()
    sc.pp.filter_cells(A, min_genes=min_genes)
    sc.pp.filter_genes(A, min_cells=min_cells)
    A.var["mt"] = A.var_names.str.upper().str.startswith(
        "MT-"
    )
    sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], inplace=True, percent_top=None)
    if A.var["mt"].sum() > 0:
        A = A[A.obs["pct_counts_mt"] < max_mt].copy()
    print(
        f"  QC: {n0} -> {A.n_obs} cells, {A.n_vars} genes "
        f"(min_genes={min_genes}, max_mt={max_mt})"
    )
    A.layers["counts"] = A.X.copy()
    sc.pp.normalize_total(A, target_sum=float(prep_cfg["normalization_scale_factor"]))
    sc.pp.log1p(A)
    sc.pp.highly_variable_genes(
        A,
        n_top_genes=int(prep_cfg["hvg_n"]),
        flavor="seurat_v3",
        layer="counts",
    )
    Ac = A[:, A.var["highly_variable"]].copy()
    del A.layers["counts"]
    gc.collect()
    sc.pp.scale(Ac, max_value=float(prep_cfg["scale_max_value"]))
    sc.tl.pca(Ac, n_comps=int(prep_cfg["pca_n"]), svd_solver="arpack")
    sc.pp.neighbors(
        Ac,
        n_neighbors=int(prep_cfg["neighbors_k"]),
        n_pcs=int(prep_cfg["pca_n"]),
    )
    ccol = f"leiden_res{res:g}"
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
        umap = None
    del Ac

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

    norm = sp.csc_matrix(A[:, de_genes].X)

    meta = pd.DataFrame(
        {
            "cell": A.obs_names.to_numpy(),
            "cluster": clusters,
        }
    )

    development_stage = _summarize_obs(A.obs, "development_stage")

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
    ap.add_argument(
        "--tissue",
        required=True,
        help="tissue context; several names joined by '|' are one operator-fixed LIST "
        "(union of the listed tissues' closures), a single name is the default single context",
    )
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
    tissue = parse_tissue_arg(a.tissue)
    preprocess(
        a.tag,
        a.h5ad,
        a.species,
        tissue,
        disease,
        a.counts_source,
        res=a.clustering_resolution,
        cross_species=a.cross_species,
    )


if __name__ == "__main__":
    main()
