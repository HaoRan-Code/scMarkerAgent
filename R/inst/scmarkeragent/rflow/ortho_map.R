suppressPackageStartupMessages({
  library(data.table)
})
ORTHO_SPECIES <- c("Human", "Mouse", "Rat")
.ortho_hub <- function(other, ortho_dir) {
  f <- file.path(ortho_dir, sprintf("ortho_Human_to_%s.csv", other))
  if (!file.exists(f)) stop("ortholog table missing: ", f, " (run rflow/ortholog_export.R)")
  d <- data.table::fread(f, colClasses = "character")
  d <- d[, 1:2]
  data.table::setnames(d, c("hum", "tgt"))
  d <- d[!is.na(hum) & nzchar(hum) & !is.na(tgt) & nzchar(tgt)]
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
  } else {
    a <- ortho_pairs(src, "Human", ortho_dir)[, .(s = src_sym, h = toupper(tgt_sym))]
    b <- ortho_pairs("Human", tgt, ortho_dir)
    data.table::setnames(b, c("h", "t"))
    d <- merge(a, b, by = "h", allow.cartesian = TRUE)[, .(src_sym = s, tgt_sym = t)]
  }
  unique(d[!is.na(src_sym) & nzchar(src_sym) & !is.na(tgt_sym) & nzchar(tgt_sym)])
}
ortho_one_to_one <- function(src, tgt, ortho_dir) {
  if (identical(src, tgt)) {
    return(setNames(character(0), character(0)))
  }
  d <- copy(ortho_pairs(src, tgt, ortho_dir))
  if (is.null(d) || !nrow(d)) {
    return(setNames(character(0), character(0)))
  }
  d[, tgt_key := toupper(tgt_sym)]
  fwd_n <- d[, .(n = uniqueN(tgt_key)), by = src_sym]
  rev_n <- d[, .(n = uniqueN(src_sym)), by = tgt_key]
  fwd_ok <- fwd_n[n == 1L, src_sym]
  rev_ok <- rev_n[n == 1L, tgt_key]
  d <- d[src_sym %in% fwd_ok & tgt_key %in% rev_ok]
  setNames(as.character(d$tgt_sym), as.character(d$src_sym))
}
