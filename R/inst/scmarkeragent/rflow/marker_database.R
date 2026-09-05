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
