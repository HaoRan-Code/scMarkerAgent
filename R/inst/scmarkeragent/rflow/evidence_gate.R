#!/usr/bin/env Rscript
SIG_LOG2FC <- as.numeric(CFG$evidence_gate$avg_log2fc_min_exclusive)
SIG_PCT <- as.numeric(CFG$evidence_gate$pct_in_min_exclusive)
SIG_AUC <- as.numeric(CFG$evidence_gate$auc_min_exclusive)
SIG_PADJ <- as.numeric(CFG$evidence_gate$padj_max_exclusive)
sig_pass <- function(avg_log2fc, pct_in_frac, padj, auc) {
  ok <- avg_log2fc > SIG_LOG2FC & pct_in_frac > SIG_PCT & auc > SIG_AUC & padj < SIG_PADJ
  ok[is.na(ok)] <- FALSE
  ok
}
