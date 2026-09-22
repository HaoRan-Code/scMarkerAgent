#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from .config import DEFAULTS

_G = DEFAULTS["evidence_gate"]
SIG_LOG2FC = float(_G["avg_log2fc_min_exclusive"])
SIG_PCT = float(_G["pct_in_min_exclusive"])
SIG_AUC = float(_G["auc_min_exclusive"])
SIG_PADJ = float(_G["padj_max_exclusive"])


def sig_pass(avg_log2fc, pct_in_frac, padj, auc):
    return (
        avg_log2fc > SIG_LOG2FC
        and pct_in_frac > SIG_PCT
        and auc > SIG_AUC
        and padj < SIG_PADJ
    )


def significant_genes_by_cluster(de):
    if de is None or not len(de):
        return {}
    columns = {"feature", "group", "avg_log2FC", "pct_in", "auc", "padj"}
    missing = columns - set(de.columns)
    if missing:
        raise ValueError(f"DE table lacks {sorted(missing)}; rebuild prep")
    hits: dict[str, set[str]] = {str(group): set() for group in de["group"].unique()}
    for row in de.itertuples(index=False):
        if sig_pass(
            float(row.avg_log2FC),
            float(row.pct_in) / 100.0,
            float(row.padj),
            float(row.auc),
        ):
            hits[str(row.group)].add(str(row.feature).upper())
    return {cluster: frozenset(genes) for cluster, genes in hits.items()}
