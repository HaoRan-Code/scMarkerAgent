#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The single marker-significance predicate.

Shared verbatim by every caller that needs it, so no two of them can drift apart on what
"significant" means.

It is NOT, however, a filter on the marker evidence handed to the agent or written to
`marker_evidence.csv`. That table carries the identity's whole measured panel, gate or no
gate, and that is deliberate: a curated marker sitting UNRAISED is the evidence
`claim_warnings` exists to report, and filtering it out would hide exactly the rows a
reader needs to weigh the call. On the 2026-08-07 freeze, 1,704 of 12,997 positive rows
(13.1%) are below this gate -- 1,542 detected in 10% of the cluster's cells or fewer, 441
with an adjusted p-value at or above 0.05. Retrieval is what the gate decides; the panel is
what the reader is shown.

A menu marker is significant in a cluster iff its one-versus-rest DE satisfies all of:
  avg_log2FC   >  SIG_LOG2FC (Seurat-style one-vs-rest log2 fold change)
  pct_in_frac  >  SIG_PCT    (detected in more than this FRACTION of the cluster's cells)
  padj         <  SIG_PADJ   (genome-wide per-cluster BH-adjusted significance)

The gate decides which candidates are retrieved, not which call is defensible; the
biological judgement is made downstream from the marker evidence itself. It is therefore
deliberately permissive, and carries no AUC term: for a gene with a positive one-vs-rest
log2 fold change the AUC condition was satisfied in all but pathological ties, so it only
removed recall without removing a decision.

`avg_log2FC` is the Seurat-style statistic written by preprocessing into the DE table:
de-log1p first, compare arithmetic means in versus out, add a pseudocount of one, take
log2. It is deliberately NOT the DE table's `logFC` column, which is a difference of mean
log-normalized expression on a natural-log scale and is retained for audit only.

`padj` carries the genome-wide per-cluster BH family, because the DE table is a slice of
one transcriptome-wide pass rather than a separate test over the menu alone.

pct_in_frac is a FRACTION in [0, 1]; a caller that stores pct as a PERCENT (0-100) must
divide by 100 before calling. The thresholds come from the configuration file, so a
deliberate sweep needs no code edit.
"""

from .config import DEFAULTS

_G = DEFAULTS["evidence_gate"]
SIG_LOG2FC = float(_G["avg_log2fc_min_exclusive"])
SIG_PCT = float(_G["pct_in_min_exclusive"])
SIG_PADJ = float(_G["padj_max_exclusive"])


def sig_pass(avg_log2fc, pct_in_frac, padj):
    """True iff a marker is significantly up-regulated in the cluster. pct_in_frac is a
    FRACTION (0-1). Shared verbatim by retrieval and by the displayed marker evidence.
    """
    return avg_log2fc > SIG_LOG2FC and pct_in_frac > SIG_PCT and padj < SIG_PADJ


def significant_genes_by_cluster(de):
    """{cluster: frozenset(GENE_UPPER)} of the genes passing `sig_pass` in that cluster.

    The cheap read of the same predicate: a caller that only needs "is this gene
    significantly up in this cluster" gets it without rebuilding the full DE lookup. The
    DE table stores pct as a PERCENT, so it is divided by 100 here, once.

    EVERY cluster the DE table tested gets a key, the empty set included, so an ABSENT
    key means the cluster was never tested. Callers rely on that distinction: a genuinely
    marker-silent cluster must abstain, but an untested one is a cache inconsistency and
    has to fail loudly instead of being reported as marker silence.
    """
    if de is None or not len(de):
        return {}
    columns = {"feature", "group", "avg_log2FC", "pct_in", "padj"}
    missing = columns - set(de.columns)
    if missing:
        raise ValueError(f"DE table lacks {sorted(missing)}; rebuild prep")
    hits: dict[str, set[str]] = {str(group): set() for group in de["group"].unique()}
    for row in de.itertuples(index=False):
        if sig_pass(
            float(row.avg_log2FC),
            float(row.pct_in) / 100.0,
            float(row.padj),
        ):
            hits[str(row.group)].add(str(row.feature).upper())
    return {cluster: frozenset(genes) for cluster, genes in hits.items()}
