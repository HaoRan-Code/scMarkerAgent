#!/usr/bin/env Rscript
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
UNSUPPORTED <- "unsupported_empty_candidate_pool"
args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args) >= 1 && nzchar(args[1])) args[1] else INPUT_TAG
top_candidates <- if (length(args) >= 2) as.integer(args[2]) else TOP_CANDIDATES
eps <- as.numeric(CFG$retrieval$epsilon)
rank_decimals <- as.integer(CFG$retrieval$rank_quantization_decimals)
set.seed(as.integer(CFG$preprocessing$random_state))
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
  stop("preprocessing metadata lacks species/tissue/disease context")
}
if (!("avg_log2FC" %in% names(de))) {
  stop(sprintf("%s: de cache has no 'avg_log2FC' column; rebuild prep", tag))
}
dat <- readRDS(file.path(CACHE, sprintf("%s_norm.rds", tag)))
dat <- dat[, meta$cell, drop = FALSE]
menu_genes <- rownames(dat)
clusters <- sort(unique(meta$cluster))
cat(sprintf(
  "==== retrieve(r) %s: %d cells, %d clusters (%s/%s/%s)%s ====\n",
  tag, nrow(meta), length(clusters), SP, paste(TI, collapse = "|"), paste(DI, collapse = "+"),
  if (CROSS_SPECIES) " [cross_species ON]" else ""
))
source(file.path(SELF, "uberon_ontology.R"))
TI <- scma_tissue_list(TI)
if (!length(TI)) stop("preprocessing metadata names no tissue context")
PRIMARY <- TI[1]
EXTRAS <- TI[-1]
if (length(EXTRAS)) {
  cat(sprintf(
    "  tissue contexts (operator list): primary '%s'; extra [%s] (<= %d attributed candidates each; %s)\n",
    PRIMARY, paste(EXTRAS, collapse = ", "), EXTRA_CONTEXT_TOP, RETRIEVAL_SEMANTICS
  ))
}
db <- scma_load_marker_db(DB_CSV)
setDT(db)
db_audit <- attr(db, "scma_db_audit")
UB <- uberon_load(OBO_UBERON)
db[, .ttn := .ub_norm(tissue_type)]
.CLOSURES <- scma_tissue_closures(UB, TI, TISSUE_ROOT)
for (tname in TI) {
  cat(sprintf(
    "  [tissue_root=%s] menu pools %d uberon ids / %d names for '%s'\n",
    TISSUE_ROOT, length(.CLOSURES[[tname]]$ids), length(.CLOSURES[[tname]]$names), tname
  ))
}
iv_ok <- if (EXCLUDE_IN_VITRO) !db$is_in_vitro else rep(TRUE, nrow(db))
.dmatch <- .disease_exact_match(db$disease_normalized, DI)
.row_ok <- iv_ok & .dmatch &
  (toupper(db$gene_qc_pass) %in% c("", "TRUE") | is.na(db$gene_qc_pass)) &
  !is.na(db$gene_symbol) & !(tolower(db$gene_symbol) %in% c("", "unknown")) &
  !is.na(db$cell_type) & db$cell_type != "" &
  db$marker_polarity %in% c("positive", "negative")
measured_genes <- as.character(cc$genes)
if (CROSS_SPECIES) source(file.path(SELF, "ortho_map.R"))
wide <- dcast(de, gene_key ~ group, value.var = "avg_log2FC")
rel_genes <- wide$gene_key
rel <- as.matrix(wide[, ..clusters])
rel_pct <- matrix(
  apply(rel, 1, function(row) frank(row, ties.method = "average") / length(row)),
  nrow = length(clusters), ncol = length(rel_genes)
)
rel_pct <- t(rel_pct)
dimnames(rel_pct) <- list(rel_genes, clusters)
de[, significant := sig_pass(avg_log2FC, pct_in / 100, padj, auc)]
significant_by_cluster <- lapply(
  clusters, function(cl) unique(de[group == cl & significant == TRUE, gene_key])
)
names(significant_by_cluster) <- clusters
cat("  per-cell marker-program scores (unsmoothed expression percentiles) ...\n")
menu_upper <- toupper(menu_genes)
column_of <- setNames(seq_along(menu_genes), menu_upper)
X <- as(t(dat), "CsparseMatrix")
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
.percentile_rank <- function(values) {
  if (length(values) == 1L) {
    return(1)
  }
  frank(round(values, rank_decimals), ties.method = "average") / length(values)
}
.admit_by_hits <- function(frame) {
  threshold <- RU_MIN_HITS
  while (threshold > 1L && sum(frame$hits >= threshold) < RU_MIN_POOL_FLOOR) {
    threshold <- threshold - 1L
  }
  saturated <- (frame$hits == frame$panel_size) & (frame$panel_size < RU_MIN_HITS)
  list(frame = frame[(hits >= threshold) | saturated], threshold = threshold)
}
cluster_of <- meta$cluster
.score_context <- function(tname) {
  .TC <- .CLOSURES[[tname]]
  tt_match <- (db$tissue_type_uberon_id %in% .TC$ids) | (db[[".ttn"]] %in% .TC$names)
  sub_db <- db[species == SP & tt_match & .row_ok]
  sub_db[, g := as.character(gene_symbol)]
  sub_db[, cell := cell_type]
  if (CROSS_SPECIES) {
    add <- list()
    for (osp in setdiff(c("Human", "Mouse", "Rat"), SP)) {
      osub <- db[species == osp & tt_match & .row_ok]
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
  sub_db[, clid := fifelse(
    !is.na(cell_type_canonical_cl_id) & cell_type_canonical_cl_id != "",
    cell_type_canonical_cl_id, cell_type_cl_id
  )]
  candidate_cl_id <- sub_db[!is.na(clid) & clid != "", .N,
    by = .(cell, clid)
  ][order(cell, -N)][, .(cl_id = clid[1]), by = cell]
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
  negatives <- panel[marker_polarity == "negative" & cell %in% elig & g %in% measured_genes]
  negatives[, gene_key := toupper(g)]
  cat(sprintf("  [%s] eligible free-text candidates: %d\n", tname, uniqueN(pos$cell)))
  total_candidates <- uniqueN(pos$cell)
  dfp <- pos[, .(holders = uniqueN(cell)), by = gene_key]
  dfp[, m := log((total_candidates + 1) / (holders + 1)) / log(total_candidates + 1)]
  specificity <- setNames(dfp$m, dfp$gene_key)
  measured_panels <- lapply(
    split(pos$gene_key, pos$cell),
    function(genes) sort(unique(genes[genes %in% menu_upper]))
  )
  measured_panels <- measured_panels[lengths(measured_panels) > 0]
  negative_panels <- lapply(split(negatives$gene_key, negatives$cell), function(g) sort(unique(g)))
  candidate_names <- sort(names(measured_panels))
  weight_rows <- rbindlist(lapply(seq_along(candidate_names), function(index) {
    genes <- measured_panels[[candidate_names[index]]]
    weights <- unname(specificity[genes])
    weights[is.na(weights)] <- 0
    if (sum(weights) <= 0) weights <- rep(1, length(genes))
    data.table(i = unname(column_of[genes]), j = index, x = weights)
  }))
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
  program_rank <- apply(program, 2, function(column) {
    frank(round(column, rank_decimals), ties.method = "average")
  })
  rm(W)
  frames <- list()
  before_gate <- list()
  thresholds <- list()
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
    before_gate[[cl]] <- length(rows)
    if (!length(rows)) {
      frames[cl] <- list(NULL)
      thresholds[cl] <- list(NULL)
      next
    }
    .admitted <- .admit_by_hits(rbindlist(rows))
    frame <- .admitted$frame
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
    frames[[cl]] <- frame
    thresholds[[cl]] <- .admitted$threshold
  }
  panel_records <- rbindlist(list(
    pos[, .(candidate = cell, gene = g, gene_key, marker_polarity = "positive", n_pub, tier)],
    negatives[, .(candidate = cell, gene = g, gene_key, marker_polarity = "negative", n_pub, tier)]
  ))
  panel_records <- panel_records[!duplicated(panel_records[, .(candidate, gene_key, marker_polarity)])]
  panel_records[, tissue_context := tname]
  list(
    tissue = tname,
    eligible = sort(unique(as.character(pos$cell))),
    specificity = specificity,
    measured_panels = measured_panels,
    negative_panels = negative_panels,
    panel_records = panel_records,
    candidate_cl_id = setNames(candidate_cl_id$cl_id, candidate_cl_id$cell),
    frames = frames,
    before_gate = before_gate,
    thresholds = thresholds
  )
}
scored_by <- list()
for (tname in TI) scored_by[[tname]] <- .score_context(tname)
P <- scored_by[[PRIMARY]]
rm(Delta)
gc(verbose = FALSE)
attributed <- character(0)
for (tname in TI) {
  for (candidate in scored_by[[tname]]$eligible) {
    if (is.na(attributed[candidate])) attributed[candidate] <- tname
  }
}
attributed <- attributed[!is.na(attributed)]
eligible_in <- lapply(names(attributed), function(candidate) {
  TI[vapply(TI, function(tname) candidate %in% scored_by[[tname]]$eligible, TRUE)]
})
names(eligible_in) <- names(attributed)
.label <- function(candidates) {
  vapply(as.character(candidates), function(c) paste(eligible_in[[c]], collapse = ";"), "")
}
scored_list <- list()
cluster_state <- list()
for (cl in clusters) {
  n_in <- sum(cluster_of == cl)
  frame <- P$frames[[cl]]
  parts <- list()
  offset <- 0L
  if (is.null(frame)) {
    cat(sprintf("  cluster %s: empty candidate pool -> %s\n", cl, UNSUPPORTED))
    state <- list(
      n_cells = as.integer(n_in), status = UNSUPPORTED, candidates = character(0)
    )
  } else {
    frame <- copy(frame)
    frame[, tissue_contexts := .label(candidate)]
    frame[, retrieval_context := PRIMARY]
    frame[, context_rank := retrieval_rank]
    parts[[length(parts) + 1L]] <- frame
    offset <- nrow(frame)
    selected <- frame$candidate[seq_len(min(top_candidates, nrow(frame)))]
    state <- list(
      n_cells = as.integer(n_in), status = "pool",
      candidates = selected, pool_size = as.integer(nrow(frame)),
      pool_size_before_hits_gate = as.integer(P$before_gate[[cl]]),
      hits_threshold_applied = as.integer(P$thresholds[[cl]])
    )
    cat(sprintf(
      "  cluster %-3s n=%-6d pool=%-4d of %-4d%s top: %s\n", cl, n_in, nrow(frame),
      P$before_gate[[cl]],
      if (P$thresholds[[cl]] == RU_MIN_HITS) "" else sprintf(" (gate relaxed to %d)", P$thresholds[[cl]]),
      paste(substr(selected[seq_len(min(3, length(selected)))], 1, 28), collapse = ", ")
    ))
  }
  if (length(EXTRAS)) {
    extra_state <- list()
    merged <- as.character(state$candidates)
    for (tname in EXTRAS) {
      E <- scored_by[[tname]]
      fe <- E$frames[[cl]]
      chosen <- character(0)
      n_attributed <- 0L
      if (!is.null(fe)) {
        own <- fe[unname(attributed[as.character(fe$candidate)]) == tname]
        own <- own[!is.na(candidate)]
        n_attributed <- nrow(own)
        if (nrow(own)) {
          own <- copy(own)
          own[, tissue_contexts := .label(candidate)]
          own[, retrieval_context := tname]
          own[, context_rank := retrieval_rank]
          own[, retrieval_rank := offset + seq_len(.N)]
          parts[[length(parts) + 1L]] <- own
          offset <- offset + nrow(own)
          chosen <- own$candidate[seq_len(min(EXTRA_CONTEXT_TOP, nrow(own)))]
        }
      }
      extra_state[[tname]] <- list(
        candidates = chosen,
        attributed_pool_size = as.integer(n_attributed),
        context_pool_size = as.integer(if (is.null(fe)) 0L else nrow(fe)),
        context_pool_size_before_hits_gate = as.integer(E$before_gate[[cl]]),
        context_hits_threshold_applied = if (is.null(E$thresholds[[cl]])) NULL else as.integer(E$thresholds[[cl]])
      )
      merged <- c(merged, chosen)
      if (length(chosen)) {
        cat(sprintf(
          "  cluster %-3s + '%s': %d of %d attributed (context pool %d of %d): %s\n", cl, tname,
          length(chosen), n_attributed, if (is.null(fe)) 0L else nrow(fe), E$before_gate[[cl]],
          paste(substr(chosen, 1, 28), collapse = ", ")
        ))
      }
    }
    state$primary_status <- state$status
    state$extra_contexts <- extra_state
    state$candidates <- merged
    if (length(merged) && identical(state$status, UNSUPPORTED)) state$status <- "pool"
  }
  cluster_state[[cl]] <- state
  if (length(parts)) scored_list[[cl]] <- rbindlist(parts, use.names = TRUE)
}
scored <- if (length(scored_list)) {
  rbindlist(scored_list, use.names = TRUE)
} else {
  data.table(
    cluster = character(0), candidate = character(0), hits = integer(0),
    panel_size = integer(0), marker_level = numeric(0), cluster_level = numeric(0),
    single_cell_level = numeric(0), retrieval_score = numeric(0),
    retrieval_rank = integer(0), tissue_contexts = character(0),
    retrieval_context = character(0), context_rank = integer(0)
  )
}
measured_panels <- P$measured_panels
negative_panels <- P$negative_panels
candidate_cl_id <- P$candidate_cl_id
specificity <- P$specificity
record_parts <- list(P$panel_records)
for (tname in EXTRAS) {
  E <- scored_by[[tname]]
  mine <- names(attributed)[attributed == tname]
  for (c in intersect(names(E$measured_panels), mine)) measured_panels[[c]] <- E$measured_panels[[c]]
  for (c in intersect(names(E$negative_panels), mine)) negative_panels[[c]] <- E$negative_panels[[c]]
  for (c in setdiff(intersect(names(E$candidate_cl_id), mine), names(candidate_cl_id))) {
    candidate_cl_id[c] <- E$candidate_cl_id[[c]]
  }
  for (g in setdiff(names(E$specificity), names(specificity))) specificity[g] <- E$specificity[[g]]
  record_parts[[length(record_parts) + 1L]] <- E$panel_records[candidate %in% mine]
}
panel_records <- rbindlist(record_parts, use.names = TRUE)
.keys_of <- function(pr) paste(pr$candidate, pr$gene_key, pr$marker_polarity, sep = "\r")
keys_by_context <- lapply(TI, function(tname) .keys_of(scored_by[[tname]]$panel_records))
names(keys_by_context) <- TI
.row_keys <- .keys_of(panel_records)
panel_records[, tissue_context := vapply(.row_keys, function(key) {
  paste(TI[vapply(TI, function(tname) key %in% keys_by_context[[tname]], TRUE)], collapse = ";")
}, "")]
eligible_by_context <- lapply(TI, function(tname) scored_by[[tname]]$eligible)
names(eligible_by_context) <- TI
tissue_members <- scma_tissue_closure(UB, TI, TISSUE_ROOT)$names
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
    candidate_cl_id = candidate_cl_id,
    tissue_members = tissue_members,
    db_audit = db_audit,
    tissue_contexts = TI,
    candidate_tissue_contexts = eligible_in[sort(names(eligible_in))],
    candidate_retrieval_context = attributed[sort(names(attributed))],
    eligible_by_context = eligible_by_context,
    retrieval_semantics = RETRIEVAL_SEMANTICS,
    extra_context_top_candidates = EXTRA_CONTEXT_TOP
  ),
  file.path(CACHE, sprintf("%s_candidate_scoring.rds", tag))
)
cat(sprintf(
  "\n[done] %s: %d scored candidate-cluster pairs across %d clusters -> %s_candidate_scoring.rds\n",
  tag, nrow(scored), length(clusters), tag
))
