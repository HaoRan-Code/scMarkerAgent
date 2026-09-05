#!/usr/bin/env Rscript
# =============================================================================
# ortho_map.R -- bidirectional ortholog SYMBOL mapping across Human/Mouse/Rat for
# the cross-species marker-sharing switch, built from the two human-hub CSVs
# (data/ortholog/ortho_Human_to_{Mouse,Rat}.csv produced by ortholog_export.R).
#
# Built EXACTLY as every cross-species arm builds it, so all of them use identical
# mappings:
#   Human<->Mouse, Human<->Rat : direct table (reversed by column swap)
#   Mouse<->Rat                : COMPOSED via Human (Mouse->Human->Rat)
# 1:many orthologs are ALL kept. Returns data.table(src_sym, tgt_sym) where src_sym is
# an UPPERCASE key and tgt_sym the target species' NATIVE spelling: the target symbol is
# written into the query species' symbol space and matched to that object's measured
# genes EXACTLY, so an all-upper PECAM1 handed to a mouse matrix would match nothing.
# The source side is only ever merged against symbols the caller uppercases first.
#   ortho_pairs(src, tgt)  ->  data.table mapping src-species symbols to tgt-species
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
})

.ortho_hub <- function(other, ortho_dir) {
  f <- file.path(ortho_dir, sprintf("ortho_Human_to_%s.csv", other))
  if (!file.exists(f)) stop("ortholog table missing: ", f, " (run rflow/ortholog_export.R)")
  d <- data.table::fread(f, colClasses = "character")
  d <- d[, 1:2]
  data.table::setnames(d, c("hum", "tgt"))
  d <- d[!is.na(hum) & nzchar(hum) & !is.na(tgt) & nzchar(tgt)]
  # Both spellings are kept as the table stores them, with an uppercase key beside each:
  # either column can be the TARGET of a direction, so neither may be folded in place.
  d[, hum_key := toupper(hum)][, tgt_key := toupper(tgt)]
  unique(d)
}

ortho_pairs <- function(src, tgt, ortho_dir) {
  sp <- c("Human", "Mouse", "Rat")
  if (src == tgt) {
    return(NULL)
  }
  stopifnot(src %in% sp, tgt %in% sp)
  if (src == "Human" && tgt %in% c("Mouse", "Rat")) {
    d <- .ortho_hub(tgt, ortho_dir)[, .(src_sym = hum_key, tgt_sym = tgt)]
  } else if (tgt == "Human" && src %in% c("Mouse", "Rat")) {
    d <- .ortho_hub(src, ortho_dir)[, .(src_sym = tgt_key, tgt_sym = hum)]
  } else { # Mouse <-> Rat composed via Human
    # Crossed on the uppercase key: the first leg now delivers Human in its native
    # spelling and the second keys on the folded one, so the 23 HGNC symbols carrying
    # lower case (C1orf43) would otherwise drop out of the composition.
    a <- ortho_pairs(src, "Human", ortho_dir)[, .(s = src_sym, h = toupper(tgt_sym))]
    b <- ortho_pairs("Human", tgt, ortho_dir)
    data.table::setnames(b, c("h", "t"))
    d <- merge(a, b, by = "h", allow.cartesian = TRUE)[, .(src_sym = s, tgt_sym = t)]
  }
  unique(d[!is.na(src_sym) & nzchar(src_sym) & !is.na(tgt_sym) & nzchar(tgt_sym)])
}
