suppressPackageStartupMessages(library(data.table))
scma_load_marker_db <- function(path) {
  if (!file.exists(path)) {
    stop("curated marker resource not found: ", path)
  }
  db <- fread(path, colClasses = "character", showProgress = FALSE)
  required <- c(
    "species", "tissue_type", "tissue_type_uberon_id",
    "disease_normalized", "gene_symbol", "gene_qc_pass", "cell_type",
    "cell_type_cl_id", "cell_type_canonical_cl_id", "marker_polarity",
    "n_pub_support", "confidence_tier", "is_in_vitro"
  )
  missing <- setdiff(required, names(db))
  if (length(missing)) {
    stop(
      "curated marker resource missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
  db[, n_pub_support := {
    value <- suppressWarnings(as.integer(n_pub_support))
    value[is.na(value)] <- 1L
    value
  }]
  db[, is_in_vitro := tolower(trimws(as.character(is_in_vitro))) %in%
    c("true", "1", "yes", "t")]
  attr(db, "scma_db_audit") <- list(
    source = normalizePath(path),
    runtime_blocklist_required = FALSE,
    blocked_rows_remaining = 0L
  )
  db
}
scma_species_positive_rows <- function(db, species) {
  frame <- as.data.table(db)
  qc_ok <- frame$gene_qc_pass == "TRUE" | is.na(frame$gene_qc_pass)
  g_ok <- !is.na(frame$gene_symbol) & !(tolower(frame$gene_symbol) %in% c("", "unknown"))
  mask <- frame$species == as.character(species) & frame$marker_polarity == "positive" &
    qc_ok & g_ok & !is.na(frame$cell_type)
  if (EXCLUDE_IN_VITRO) mask <- mask & !frame$is_in_vitro
  frame[mask, .(
    cell_type = as.character(cell_type), gene_symbol = as.character(gene_symbol),
    n_pub_support = as.integer(n_pub_support), confidence_tier = as.character(confidence_tier),
    is_recommended_marker = as.character(is_recommended_marker),
    tissue_type = as.character(tissue_type)
  )]
}
scma_species_definers <- function(rows, names, k = 8L) {
  sub <- rows[cell_type %in% names]
  if (!nrow(sub)) return(list())
  tier_rank <- c(high = 0L, medium = 1L, low = 2L)
  sub <- copy(sub)
  sub[, rec_flag := as.integer(toupper(is_recommended_marker) == "TRUE")]
  sub[, tier_r := {
    r <- unname(tier_rank[tolower(confidence_tier)])
    r[is.na(r)] <- 3L
    as.integer(r)
  }]
  grp <- sub[, .(
    n_pub = max(n_pub_support), rec_flag = max(rec_flag), tier_rank = min(tier_r),
    tissue_contexts = uniqueN(tissue_type)
  ), by = .(cell_type, gene_symbol)]
  inv_tier <- c("0" = "high", "1" = "medium", "2" = "low", "3" = NOT_AVAILABLE)
  out <- list()
  for (name in unique(grp$cell_type)) {
    block <- grp[cell_type == name]
    ord <- order(-block$n_pub, -block$rec_flag, .cp(block$gene_symbol), method = "radix")
    block <- block[ord]
    head_n <- min(as.integer(k), nrow(block))
    head <- block[seq_len(head_n)]
    rest <- block[-seq_len(head_n)]
    extra <- rest[rec_flag == 1L][seq_len(min(4L, sum(rest$rec_flag == 1L)))]
    chosen <- rbind(head, extra)
    out[[as.character(name)]] <- lapply(seq_len(nrow(chosen)), function(i) {
      row <- chosen[i]
      list(
        gene = as.character(row$gene_symbol), n_pub = as.integer(row$n_pub),
        tier = unname(inv_tier[as.character(row$tier_rank)]),
        recommended = as.logical(row$rec_flag), tissue_contexts = as.integer(row$tissue_contexts)
      )
    })
  }
  out
}
scma_species_gene_claimants <- function(rows) {
  if (!nrow(rows)) return(list())
  grp <- rows[, .(n_pub = max(n_pub_support)), by = .(gene = toupper(gene_symbol), cell_type)]
  setorder(grp, gene, -n_pub, cell_type)
  by_gene <- split(grp[, .(cell_type, n_pub)], grp$gene)
  lapply(by_gene, function(block) {
    lapply(seq_len(nrow(block)), function(i) {
      list(as.character(block$cell_type[i]), as.integer(block$n_pub[i]))
    })
  })
}
scma_known_genes <- function(db) {
  g <- as.character(db$gene_symbol)
  g <- g[!is.na(g)]
  unique(toupper(g))
}
