#!/usr/bin/env Rscript
# =============================================================================
# evidence_gate.R -- the ONE marker-significance predicate. Single source of truth
# shared by candidate retrieval and by the marker evidence handed to the annotating
# agent, so the two can never drift.
#
# A menu marker is significant in a cluster iff its one-vs-rest DE is:
#   avg_log2FC   >  SIG_LOG2FC (Seurat-style one-vs-rest log2 fold change)
#   pct_in_frac  >  SIG_PCT    (detected in more than this FRACTION of the cluster's cells)
#   padj         <  SIG_PADJ   (genome-wide per-cluster BH-adjusted significance)
#
# The gate decides which candidates are retrieved, not which call is defensible, so it is
# deliberately permissive and carries no AUC term: for a gene with a positive one-vs-rest
# log2 fold change the AUC condition held in all but pathological ties.
#
# avg_log2FC is the Seurat-style statistic preprocessing.R writes into the DE table
# (de-log1p, compare arithmetic means in versus out, add a pseudocount of one, take log2).
# It is deliberately NOT presto's `logFC` column, which is a difference of mean
# log-normalized expression on a natural-log scale and is retained for audit only.
#
# pct_in_frac is a FRACTION in [0, 1]; a caller that stores pct as a PERCENT
# (0-100) must divide by 100 before calling.
# =============================================================================
SIG_LOG2FC <- as.numeric(CFG$evidence_gate$avg_log2fc_min_exclusive)
SIG_PCT <- as.numeric(CFG$evidence_gate$pct_in_min_exclusive)
SIG_PADJ <- as.numeric(CFG$evidence_gate$padj_max_exclusive)

# TRUE iff a marker is significantly up-regulated in the cluster. `pct_in_frac` is a
# FRACTION (0-1). Vector-safe and NA-safe: NA inputs (e.g. a gene absent from the DE
# table) yield FALSE, which is how every caller guards a missing row.
sig_pass <- function(avg_log2fc, pct_in_frac, padj) {
  ok <- avg_log2fc > SIG_LOG2FC & pct_in_frac > SIG_PCT & padj < SIG_PADJ
  ok[is.na(ok)] <- FALSE
  ok
}
