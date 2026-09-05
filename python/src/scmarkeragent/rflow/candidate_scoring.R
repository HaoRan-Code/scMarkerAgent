#!/usr/bin/env Rscript
# =============================================================================
# candidate_scoring.R -- retrieve the free-text cell-type candidates a cluster could be.
# Every arm applies the same formulas in the same order.
#
# Identity is the curated free-text cell-type string for a species. No ontology is read,
# resolved, expanded or compared: two curated names are two candidates, and whether they
# denote the same biology is a question for the annotating agent, not a merge rule.
#
# A candidate enters a cluster's pool when enough of its measured positive markers are
# significantly up-regulated there -- see .admit_by_hits for what "enough" means and why.
# The pool is then ordered by three scores, each answering a different "relative to what"
# question about the candidate's WHOLE measured panel:
#   M  marker level      length-normalized mass of specific curated marker evidence present
#   C  cluster level     does this panel single out THIS cluster among all clusters
#   S  single-cell level do this cluster's cells carry the program more than other cells
# The three are percentile-ranked within the cluster and combined by geometric mean into
# one retrieval order. That order decides what the agents see; it is not a confidence.
#
# Negative markers never enter retrieval or any score. Their measured behaviour is
# collected only so it can be shown as evidence.
#
# Output: cache/<tag>_candidate_scoring.rds
#
#   Usage: Rscript candidate_scoring.R <tag> [top_candidates]
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})
.script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(.script_argument) != 1L) {
  stop("candidate_scoring.R must be executed with Rscript")
}
.sd <- dirname(normalizePath(sub("^--file=", "", .script_argument)))
source(file.path(.sd, "config.R"))
source(file.path(.sd, "marker_database.R"))
sys.source(file.path(.sd, "evidence_gate.R"), envir = environment())

## A cluster with no candidate whose panel intersects its significant DE has nothing to
## retrieve. Its cells are kept so coverage stays honest, and annotation emits Unknown.
UNSUPPORTED <- "unsupported_empty_candidate_pool"

args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args) >= 1 && nzchar(args[1])) args[1] else INPUT_TAG
top_candidates <- if (length(args) >= 2) as.integer(args[2]) else TOP_CANDIDATES
eps <- as.numeric(CFG$retrieval$epsilon)
rank_decimals <- as.integer(CFG$retrieval$rank_quantization_decimals)
set.seed(as.integer(CFG$preprocessing$random_state))

## ---- frozen preprocessing ---------------------------------------------------
cc <- readRDS(file.path(CACHE, paste0(tag, "_de_meta.rds")))
meta <- as.data.table(cc$meta)
meta[, cluster := as.character(cluster)]
de <- as.data.table(cc$de)
de[, feature := as.character(feature)]
de[, group := as.character(group)]
de[, gene_key := toupper(feature)]
SP <- cc$species
TI <- cc$tissue
DI <- cc$disease
DEV <- if (length(cc$development_stage)) as.character(cc$development_stage)[1] else ""
if (length(SP) == 0 || is.na(SP[1]) || !nzchar(SP[1])) {
  stop("frozen prep metadata lacks species/tissue/disease context")
}
if (!("avg_log2FC" %in% names(de))) {
  stop(sprintf("%s: de cache has no 'avg_log2FC' column; rebuild prep", tag))
}
dat <- readRDS(file.path(CACHE, sprintf("%s_norm.rds", tag))) # genes x cells, log-norm
dat <- dat[, meta$cell, drop = FALSE]
menu_genes <- rownames(dat)
clusters <- sort(unique(meta$cluster))
cat(sprintf(
  "==== retrieve(r) %s: %d cells, %d clusters (%s/%s/%s)%s ====\n",
  tag, nrow(meta), length(clusters), SP, TI, paste(DI, collapse = "+"),
  if (CROSS_SPECIES) " [cross_species ON]" else ""
))

## ---- curated candidate universe for this context ----------------------------
db <- scma_load_marker_db(DB_CSV)
setDT(db)
db_audit <- attr(db, "scma_db_audit")
source(file.path(SELF, "uberon_ontology.R"))
UB <- uberon_load(OBO_UBERON)
.tissue_closure <- function(tissue, tissue_root = TISSUE_ROOT) {
  scma_tissue_closure(UB, tissue, tissue_root)
}
db[, .ttn := .ub_norm(tissue_type)]
.TC <- .tissue_closure(TI)
tt_match <- (db$tissue_type_uberon_id %in% .TC$ids) | (db[[".ttn"]] %in% .TC$names)
iv_ok <- if (EXCLUDE_IN_VITRO) !db$is_in_vitro else rep(TRUE, nrow(db))
cat(sprintf(
  "  [tissue_root=%s] menu pools %d uberon ids / %d names for '%s'\n",
  TISSUE_ROOT, length(.TC$ids), length(.TC$names), TI
))
## Disease policy: normalize both sides by trim+lower, then require exact equality.
## Deliberately no substring or descendant expansion: "glioblastoma" must not absorb
## "IDH-wildtype glioblastoma".
.dmatch <- .disease_exact_match(db$disease_normalized, DI)
sub_db <- db[species == SP & tt_match & iv_ok & .dmatch &
  (toupper(gene_qc_pass) %in% c("", "TRUE") | is.na(gene_qc_pass)) &
  !is.na(gene_symbol) & !(tolower(gene_symbol) %in% c("", "unknown")) &
  !is.na(cell_type) & cell_type != "" &
  marker_polarity %in% c("positive", "negative")]
sub_db[, g := as.character(gene_symbol)]
sub_db[, cell := cell_type]
if (CROSS_SPECIES) {
  source(file.path(SELF, "ortho_map.R"))
  add <- list()
  for (osp in setdiff(c("Human", "Mouse", "Rat"), SP)) {
    osub <- db[species == osp & tt_match & iv_ok & .dmatch &
      (toupper(gene_qc_pass) %in% c("", "TRUE") | is.na(gene_qc_pass)) &
      !is.na(gene_symbol) & !(tolower(gene_symbol) %in% c("", "unknown")) &
      !is.na(cell_type) & cell_type != "" &
      marker_polarity %in% c("positive", "negative")]
    if (!nrow(osub)) next
    osub[, g := toupper(gene_symbol)]
    osub[, cell := cell_type]
    pr <- ortho_pairs(osp, SP, ORTHO_DIR)
    if (is.null(pr) || !nrow(pr)) next
    osub <- merge(osub, pr, by.x = "g", by.y = "src_sym", allow.cartesian = TRUE)
    osub[, g := tgt_sym][, tgt_sym := NULL]
    add[[osp]] <- osub
  }
  if (length(add)) {
    sub_db <- rbindlist(c(list(sub_db), add), use.names = TRUE, fill = TRUE)
    sub_db <- unique(sub_db, by = c("cell", "g", "marker_polarity", "n_pub_support", "confidence_tier"))
  }
}
has_tier <- "confidence_tier" %in% names(sub_db)
## The curated CL id is carried forward for ONE display column and nothing else: it is
## never resolved against an ontology, expanded, or compared between two candidates.
sub_db[, clid := fifelse(
  !is.na(cell_type_canonical_cl_id) & cell_type_canonical_cl_id != "",
  cell_type_canonical_cl_id, cell_type_cl_id
)]
candidate_cl_id <- sub_db[!is.na(clid) & clid != "", .N,
  by = .(cell, clid)
][order(cell, -N)][, .(cl_id = clid[1]), by = cell]

## Exactly one row per (cell, g, marker_polarity), keeping the max n_pub and the best
## tier. Deduplicating on n_pub or tier as well would let the same curated marker survive
## as several rows and be double-counted downstream.
panel <- sub_db[, .(cell, g, marker_polarity,
  n_pub = suppressWarnings(as.integer(n_pub_support)),
  tier = if (has_tier) confidence_tier else NA_character_
)]
panel[is.na(n_pub), n_pub := 1L]
.tier_rank <- c(high = 0L, medium = 1L, low = 2L)
panel[, .tr := {
  r <- unname(.tier_rank[tolower(as.character(tier))])
  r[is.na(r)] <- 3L
  as.integer(r)
}]
setorder(panel, -n_pub, .tr)
panel <- panel[!duplicated(panel[, .(cell, g, marker_polarity)])]
panel[, .tr := NULL]
rm(db, sub_db)
gc(verbose = FALSE)

## EXACT symbol membership (== marker_database.eligible_candidates). This is where a
## curated symbol meets an uploaded one for the second time, and case-folding it let a
## mouse panel qualify on a human matrix and vice versa. `gene_key` stays uppercase: by
## the time it is built, `g` IS a measured symbol, and measured symbols are unique
## case-insensitively (duplicates were summed in preprocessing), so the key is an
## injective relabelling within one namespace rather than a second crossing.
measured_genes <- as.character(cc$genes)
pos <- panel[marker_polarity == "positive" & g %in% measured_genes]
pos[, gene_key := toupper(g)]
elig <- pos[, .(np = uniqueN(gene_key)), by = cell][np >= RU_MIN_ELIG_GENES, cell]
if (isTRUE(CORROBORATION_ONLY)) {
  reli <- pos[, .(
    mx = suppressWarnings(max(n_pub, na.rm = TRUE)),
    himed = sum(tier %in% CORROBORATING_TIERS, na.rm = TRUE)
  ), by = cell]
  reli[!is.finite(mx), mx := 1L]
  elig <- intersect(elig, reli[mx >= MIN_CORROBORATING_PUBLICATIONS | himed >= 1, cell])
}
pos <- pos[cell %in% elig]
## Negatives cross the same boundary and are gated the same way. They were only ever
## displayed for genes the DE table carries, so restricting them here changes nothing
## except that a curated symbol differing from the measured one by case no longer
## reaches a panel the positive of the same name was refused from.
negatives <- panel[marker_polarity == "negative" & cell %in% elig & g %in% measured_genes]
negatives[, gene_key := toupper(g)]
cat(sprintf("  eligible free-text candidates: %d\n", uniqueN(pos$cell)))

## ---- m_g: how few curated cell types claim the gene -------------------------
## Publication counts and evidence tiers deliberately stay out of it; they reach the
## agents as separate evidence rather than folded into one number.
total_candidates <- uniqueN(pos$cell)
dfp <- pos[, .(holders = uniqueN(cell)), by = gene_key]
dfp[, m := log((total_candidates + 1) / (holders + 1)) / log(total_candidates + 1)]
specificity <- setNames(dfp$m, dfp$gene_key)

## Panels are restricted to the measured menu; every score reads the whole panel, and the
## significant subset only decides admission and what is shown as a hit.
menu_upper <- toupper(menu_genes)
column_of <- setNames(seq_along(menu_genes), menu_upper)
measured_panels <- lapply(
  split(pos$gene_key, pos$cell),
  function(genes) sort(unique(genes[genes %in% menu_upper]))
)
measured_panels <- measured_panels[lengths(measured_panels) > 0]
negative_panels <- lapply(split(negatives$gene_key, negatives$cell), function(g) sort(unique(g)))

## ---- cross-cluster percentile of avg_log2FC ---------------------------------
## 1.0 means no cluster separates this gene more strongly than this one.
wide <- dcast(de, gene_key ~ group, value.var = "avg_log2FC")
rel_genes <- wide$gene_key
rel <- as.matrix(wide[, ..clusters])
rel_pct <- matrix(
  apply(rel, 1, function(row) frank(row, ties.method = "average") / length(row)),
  nrow = length(clusters), ncol = length(rel_genes)
)
rel_pct <- t(rel_pct)
dimnames(rel_pct) <- list(rel_genes, clusters)

## ---- significant genes per cluster ------------------------------------------
de[, significant := sig_pass(avg_log2FC, pct_in / 100, padj)]
significant_by_cluster <- lapply(
  clusters, function(cl) unique(de[group == cl & significant == TRUE, gene_key])
)
names(significant_by_cluster) <- clusters

## ---- per-cell marker-program scores -----------------------------------------
## Per-gene expression percentile over ALL cells, computed analytically on the sparse
## matrix: every non-detecting cell shares one average rank, so only the detected values
## need sorting and no dense cells-by-genes matrix is ever built. The expression is the
## log-normalized value as measured -- no neighbourhood smoothing, because the whole point
## of the third axis is that it reads individual cells rather than the kNN graph.
cat("  per-cell marker-program scores (unsmoothed expression percentiles) ...\n")
X <- as(t(dat), "CsparseMatrix") # cells x genes
n_cells <- nrow(X)
baseline <- numeric(ncol(X))
delta_x <- X@x
for (j in seq_len(ncol(X))) {
  from <- X@p[j] + 1L
  to <- X@p[j + 1L]
  zeros <- n_cells - (to - X@p[j])
  zero_pct <- ((zeros + 1) / 2) / n_cells
  baseline[j] <- zero_pct
  if (to >= from) {
    ranks <- frank(X@x[from:to], ties.method = "average") + zeros
    delta_x[from:to] <- ranks / n_cells - zero_pct
  }
}
Delta <- X
Delta@x <- delta_x
rm(X)
gc(verbose = FALSE)

candidate_names <- sort(names(measured_panels))
weight_rows <- rbindlist(lapply(seq_along(candidate_names), function(index) {
  genes <- measured_panels[[candidate_names[index]]]
  weights <- unname(specificity[genes])
  weights[is.na(weights)] <- 0
  if (sum(weights) <= 0) weights <- rep(1, length(genes))
  data.table(i = unname(column_of[genes]), j = index, x = weights)
}))
## A marker resource carrying nothing for this species leaves no candidates at all, and
## rbindlist() of an empty list drops the columns rather than keeping them empty, so the
## three slots arrive as NULL and sparseMatrix() cannot derive a type from them. Coercing
## keeps that case a zero-column W, which every step below already treats as a no-op, so
## an empty pool is reported instead of a failure.
W <- sparseMatrix(
  i = as.integer(weight_rows$i), j = as.integer(weight_rows$j),
  x = as.numeric(weight_rows$x),
  dims = c(ncol(Delta), length(candidate_names))
)
weight_total <- Matrix::colSums(W)
constant <- as.numeric(Matrix::crossprod(W, baseline))
program <- as.matrix(Delta %*% W)
program <- sweep(program, 2, constant, "+")
program <- sweep(program, 2, weight_total, "/")
colnames(program) <- candidate_names
## Quantize before ranking. Two cells with the same expression across the panel should
## tie, and whether their sums land bit-identically depends on accumulation order;
## without this the two arms differ by exactly one tie's worth of rank.
program_rank <- apply(program, 2, function(column) {
  frank(round(column, rank_decimals), ties.method = "average")
})
rm(Delta, W)
gc(verbose = FALSE)

## ---- per cluster: pool, three scores, retrieval order -----------------------
.percentile_rank <- function(values) {
  if (length(values) == 1L) {
    return(1)
  }
  frank(round(values, rank_decimals), ties.method = "average") / length(values)
}

## Which candidates the retrieval order is built from, and the gate that let them in.
## A candidate needs RU_MIN_HITS of its OWN positive markers significantly up-regulated
## here; one or two hits is overwhelmingly a retrieval false positive and crowds out the
## real candidates in a shortlist the review tier can afford to read. The requirement is
## relaxed one step at a time where enforcing it would leave fewer than RU_MIN_POOL_FLOOR
## candidates: a thin curated menu is a fact about the database's coverage of that organ,
## not evidence that the cluster has no identity. The threshold actually applied travels
## with the cluster so a reader never has to infer it.
.admit_by_hits <- function(frame) {
  threshold <- RU_MIN_HITS
  while (threshold > 1L && sum(frame$hits >= threshold) < RU_MIN_POOL_FLOOR) {
    threshold <- threshold - 1L
  }
  list(frame = frame[hits >= threshold], threshold = threshold)
}

scored_list <- list()
cluster_state <- list()
cluster_of <- meta$cluster
for (cl in clusters) {
  inside <- cluster_of == cl
  n_in <- sum(inside)
  n_out <- n_cells - n_in
  hits_here <- significant_by_cluster[[cl]]
  rows <- list()
  for (candidate in candidate_names) {
    genes <- measured_panels[[candidate]]
    intersection <- genes[genes %in% hits_here]
    if (!length(intersection)) next
    weight_hit <- unname(specificity[intersection])
    weight_hit[is.na(weight_hit)] <- 0
    weight_panel <- unname(specificity[genes])
    weight_panel[is.na(weight_panel)] <- 0
    percentile <- unname(rel_pct[genes, cl])
    percentile[is.na(percentile)] <- 0
    position <- match(candidate, candidate_names)
    rank_in <- program_rank[inside, position]
    values <- program[inside, position]
    outside <- program[!inside, position]
    rows[[length(rows) + 1L]] <- data.table(
      cluster = cl,
      candidate = candidate,
      hits = length(intersection),
      panel_size = length(genes),
      marker_level = sum(weight_hit) / sqrt(length(genes)),
      cluster_level = if (sum(weight_panel) > 0) {
        sum(weight_panel * percentile) / max(sum(weight_panel), eps)
      } else {
        mean(percentile)
      },
      single_cell_level = if (n_out > 0) {
        ## `sum()` over a logical yields an INTEGER, so both products below overflow
        ## .Machine$integer.max and return NA -- with a warning, not an error -- once
        ## a cluster is large enough: n_in * n_out passes 2^31-1 at roughly 93k cells
        ## evenly split, and n_in * (n_in + 1) at n_in = 46,341. The Python arm reads
        ## the same formula in float64, so leaving these as integers is also what made
        ## the two arms disagree above that size.
        (sum(rank_in) - as.numeric(n_in) * (as.numeric(n_in) + 1) / 2) /
          (as.numeric(n_in) * as.numeric(n_out))
      } else {
        NA_real_
      },
      program_in_q25 = unname(quantile(values, 0.25)),
      program_in_median = median(values),
      program_in_q75 = unname(quantile(values, 0.75)),
      program_out_q25 = if (n_out > 0) unname(quantile(outside, 0.25)) else NA_real_,
      program_out_median = if (n_out > 0) median(outside) else NA_real_,
      program_out_q75 = if (n_out > 0) unname(quantile(outside, 0.75)) else NA_real_
    )
  }
  if (!length(rows)) {
    cat(sprintf("  cluster %s: empty candidate pool -> %s\n", cl, UNSUPPORTED))
    cluster_state[[cl]] <- list(
      n_cells = as.integer(n_in), status = UNSUPPORTED, candidates = character(0)
    )
    next
  }
  ## The admission gate runs before the three axes are percentile-ranked, because a
  ## percentile is a position within the pool: ranking first and filtering after would
  ## leave every candidate's score describing a pool it is no longer in.
  pool_before_gate <- length(rows)
  .admitted <- .admit_by_hits(rbindlist(rows))
  frame <- .admitted$frame
  hits_threshold <- .admitted$threshold
  ## A single-cluster dataset has no out-group, so the single-cell axis is undefined
  ## rather than perfect; it then contributes nothing and the other two decide. Any
  ## OTHER NA means the axis failed to compute, and zeroing it drops a third of the
  ## ranking evidence for the whole cluster while still publishing a retrieval_score
  ## that looks comparable. That is how an integer overflow went unnoticed across two
  ## released example packages, so it now stops the run instead.
  single <- frame$single_cell_level
  if (anyNA(single)) {
    if (n_out > 0L) {
      stop(sprintf(
        paste0("cluster %s: single_cell_level is NA for %d of %d candidates while ",
               "the out-group holds %d cells -- the third retrieval axis failed to ",
               "compute; refusing to rank on two axes"),
        cl, sum(is.na(single)), length(single), n_out
      ))
    }
    single <- rep(0, nrow(frame))
  }
  frame[, rank_marker_level := .percentile_rank(marker_level)]
  frame[, rank_cluster_level := .percentile_rank(cluster_level)]
  frame[, rank_single_cell_level := .percentile_rank(single)]
  frame[, retrieval_score := (rank_marker_level * rank_cluster_level *
    rank_single_cell_level)^(1 / 3)]
  setorder(
    frame, -retrieval_score, -cluster_level, -single_cell_level, -marker_level, candidate
  )
  frame[, retrieval_rank := seq_len(.N)]
  scored_list[[cl]] <- frame
  selected <- frame$candidate[seq_len(min(top_candidates, nrow(frame)))]
  cluster_state[[cl]] <- list(
    n_cells = as.integer(n_in), status = "pool",
    candidates = selected, pool_size = as.integer(nrow(frame)),
    pool_size_before_hits_gate = as.integer(pool_before_gate),
    hits_threshold_applied = as.integer(hits_threshold)
  )
  cat(sprintf(
    "  cluster %-3s n=%-6d pool=%-4d of %-4d%s top: %s\n", cl, n_in, nrow(frame),
    pool_before_gate,
    if (hits_threshold == RU_MIN_HITS) "" else sprintf(" (gate relaxed to %d)", hits_threshold),
    paste(substr(selected[seq_len(min(3, length(selected)))], 1, 28), collapse = ", ")
  ))
}
## Every cluster can come back with an empty pool -- that is what a marker resource
## holding nothing for this species looks like, not a failure. The reporting layer still
## asks this table for a cluster's pool, so the empty case has to keep its columns; these
## are the ones every arm emits for the same case, so the layer below is handed the same
## frame whether or not anything scored.
scored <- if (length(scored_list)) {
  rbindlist(scored_list)
} else {
  data.table(
    cluster = character(0), candidate = character(0), hits = integer(0),
    panel_size = integer(0), marker_level = numeric(0), cluster_level = numeric(0),
    single_cell_level = numeric(0), retrieval_score = numeric(0),
    retrieval_rank = integer(0)
  )
}

panel_records <- rbindlist(list(
  pos[, .(candidate = cell, gene = g, gene_key, marker_polarity = "positive", n_pub, tier)],
  negatives[, .(candidate = cell, gene = g, gene_key, marker_polarity = "negative", n_pub, tier)]
))
panel_records <- panel_records[!duplicated(panel_records[, .(candidate, gene_key, marker_polarity)])]

saveRDS(
  list(
    context = list(
      species = SP, tissue = TI, disease = DI, development_stage = DEV,
      n_cells = nrow(meta), n_clusters = length(clusters)
    ),
    scored = scored,
    clusters = cluster_state,
    marker_specificity = specificity,
    measured_panels = measured_panels,
    negative_panels = negative_panels,
    panel_records = panel_records,
    top_candidates = top_candidates,
    min_significant_hits = RU_MIN_HITS,
    min_pool_floor = RU_MIN_POOL_FLOOR,
    candidate_cl_id = setNames(candidate_cl_id$cl_id, candidate_cl_id$cell),
    tissue_members = .TC$names,
    db_audit = db_audit
  ),
  file.path(CACHE, sprintf("%s_candidate_scoring.rds", tag))
)
cat(sprintf(
  "\n[done] %s: %d scored candidate-cluster pairs across %d clusters -> %s_candidate_scoring.rds\n",
  tag, nrow(scored), length(clusters), tag
))
