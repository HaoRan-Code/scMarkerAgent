suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})
#' The tie-break among sentences a marker's publication count cannot separate.
#' The same seed, separator, UTF-8 payload and SHA-256 on every arm, so all of them
#' serve the same sentences.
scma_source_order_key <- function(seed, pmcid, sentence) {
  payload <- enc2utf8(paste(
    as.character(seed), as.character(pmcid), as.character(sentence),
    sep = "\u001f"
  ))
  vapply(payload, function(one) {
    digest(one, algo = "sha256", serialize = FALSE)
  }, character(1), USE.NAMES = FALSE)
}
#' Species-scoped load of the curated source export.
scma_load_sources <- function(path, species) {
  if (!file.exists(path)) {
    stop(
      "curated source export not found: ", path,
      ". Point SCMA_RESOURCE_DIR at a resource bundle containing it.",
      call. = FALSE
    )
  }
  wanted <- c(
    "species", "tissue_type", "disease_normalized", "cell_type",
    "gene_symbol", "n_pub_support", "pmcid", "pmid", "source"
  )
  db <- fread(
    path,
    select = wanted, colClasses = "character", showProgress = FALSE
  )
  missing <- setdiff(wanted, names(db))
  if (length(missing)) {
    stop(
      "curated source export missing columns: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  wanted_species <- as.character(species)
  db <- db[db$species == wanted_species]
  db[, n_pub_support := {
    value <- suppressWarnings(as.integer(n_pub_support))
    value[is.na(value)] <- 1L
    value
  }]
  db[, tissue_norm := .ub_norm(tissue_type)]
  db[, disease_norm := .normalize_disease_name(disease_normalized)]
  db[, gene_upper := toupper(gene_symbol)]
  db[]
}
#' A (species, tissue, disease)-scoped slice with each pair's sentences in served order.
#'
#' Deduplicated on the sentence, then ordered by publication count descending and
#' `scma_source_order_key` within a count -- the WHOLE evidence for the pair, not a
#' shortlist, because a caller that has read the first few can come back for the rest.
#'
#' `cross_species` defaults to the run configuration (features.cross_species_markers),
#' so the sentence lookup tracks the SAME switch that pools the marker menu
#' (candidate_scoring.R). ON pools the OTHER species' source rows at the SAME tissue +
#' disease, their gene symbols ortholog-mapped into the query species' symbol space by
#' the same `ortho_pairs` the marker pooling uses -- a translated marker row therefore
#' carries its donor species' sentences instead of an empty lookup. OFF is byte-identical
#' to the within-species path.
scma_source_context <- function(ub, path, species, tissue, disease,
                                tissue_root = TISSUE_ROOT,
                                seed = SOURCE_ORDER_SEED,
                                cross_species = CROSS_SPECIES) {
  closure <- scma_tissue_closure(ub, tissue, tissue_root)
  target <- .normalize_disease_query(disease)
  .scope <- function(rows) {
    keep <- rows$tissue_norm %in% closure$names
    if (!is.null(target)) {
      keep <- keep & (rows$disease_norm %in% target)
    }
    rows[keep]
  }
  sub <- .scope(scma_load_sources(path, species))
  if (isTRUE(cross_species)) {
    if (!exists("ortho_pairs", mode = "function")) {
      source(file.path(PACKAGE_ROOT, "rflow", "ortho_map.R"))
    }
    add <- list()
    for (osp in setdiff(c("Human", "Mouse", "Rat"), as.character(species))) {
      donor <- .scope(scma_load_sources(path, osp))
      if (!nrow(donor)) next
      pr <- ortho_pairs(osp, as.character(species), ORTHO_DIR)
      if (is.null(pr) || !nrow(pr)) next
      donor <- merge(donor, pr,
        by.x = "gene_upper", by.y = "src_sym", allow.cartesian = TRUE
      )
      if (!nrow(donor)) next
      donor[, gene_symbol := tgt_sym]
      donor[, gene_upper := toupper(tgt_sym)]
      donor[, tgt_sym := NULL]
      add[[osp]] <- donor
    }
    if (length(add)) {
      sub <- rbindlist(c(list(sub), add), use.names = TRUE)
    }
  }
  if (nrow(sub)) {
    sub <- unique(sub, by = c("cell_type", "gene_upper", "source"))
    sub[, order_key := scma_source_order_key(seed, pmcid, source)]
    setorder(sub, cell_type, gene_upper, -n_pub_support, order_key)
    setkey(sub, cell_type, gene_upper)
  }
  list(rows = sub, seed = as.character(seed))
}
#' The sentences behind a few cell types in EVERY tissue context of one species,
#' disease-gated exactly as `scma_source_context` is.
#'
#' A borrowed candidate's panel was assembled from other tissues' rows, so its sentences
#' live there too; the tissue gate is dropped for these names only, the disease gate is
#' kept, and the served order is the same fixed order.
#'
#' `donor_species_by_name` (a named list, name -> character vector of donor species) and
#' `ortho_dir` cover the OTHER borrowed shape: a name `scma_borrow_candidates` placed on
#' the page from a different species entirely (`features.cross_species_borrow`), for
#' which no within-species row -- in this tissue or any other -- exists at all, by
#' construction. Its sentences only exist under the donor's own species and its own gene
#' spelling, so for each donor species this pools that donor's rows for its share of the
#' wanted names (disease-gated the same way) and translates `gene_symbol` through the
#' SAME one-to-one ortholog map `scma_borrow_candidates` used to build the candidate's
#' panel in the first place.
scma_source_types_across_tissues <- function(path, species, cell_types, disease,
                                             seed = SOURCE_ORDER_SEED,
                                             donor_species_by_name = NULL,
                                             ortho_dir = NULL) {
  wanted <- unique(as.character(cell_types))
  target <- .normalize_disease_query(disease)
  sub <- scma_load_sources(path, species)
  keep <- sub$cell_type %in% wanted
  if (!is.null(target)) {
    keep <- keep & (sub$disease_norm %in% target)
  }
  sub <- sub[keep]
  if (length(donor_species_by_name) && !is.null(ortho_dir)) {
    if (!exists("ortho_one_to_one", mode = "function")) {
      source(file.path(PACKAGE_ROOT, "rflow", "ortho_map.R"))
    }
    by_donor <- list()
    for (nm in names(donor_species_by_name)) {
      for (donor in donor_species_by_name[[nm]]) {
        by_donor[[donor]] <- union(by_donor[[donor]], nm)
      }
    }
    add <- list()
    for (donor in names(by_donor)) {
      gene_map <- ortho_one_to_one(donor, as.character(species), ortho_dir)
      if (!length(gene_map)) next
      donor_rows <- scma_load_sources(path, donor)
      donor_keep <- donor_rows$cell_type %in% by_donor[[donor]]
      if (!is.null(target)) {
        donor_keep <- donor_keep & (donor_rows$disease_norm %in% target)
      }
      donor_rows <- donor_rows[donor_keep]
      if (!nrow(donor_rows)) next
      map_dt <- data.table(src_sym = names(gene_map), tgt_sym = as.character(gene_map))
      donor_rows <- merge(donor_rows, map_dt,
        by.x = "gene_upper", by.y = "src_sym", allow.cartesian = TRUE
      )
      if (!nrow(donor_rows)) next
      donor_rows[, gene_symbol := tgt_sym]
      donor_rows[, gene_upper := toupper(tgt_sym)]
      donor_rows[, tgt_sym := NULL]
      add[[donor]] <- donor_rows
    }
    if (length(add)) {
      sub <- rbindlist(c(list(sub), add), use.names = TRUE)
    }
  }
  if (nrow(sub)) {
    sub <- unique(sub, by = c("cell_type", "gene_upper", "source"))
    sub[, order_key := scma_source_order_key(seed, pmcid, source)]
    setorder(sub, cell_type, gene_upper, -n_pub_support, order_key)
    setkey(sub, cell_type, gene_upper)
  }
  list(rows = sub, seed = as.character(seed))
}
#' The native context's sentences, plus other contexts' sentences for the borrowed
#' candidates only.
#'
#' Quacks like a context for `scma_source_records`, so the packet builders and the per-cluster
#' server need no branch: a name that is not borrowed is answered from the native slice alone.
scma_composite_sources <- function(native, borrowed, borrowed_names) {
  list(
    kind = "composite", native = native, borrowed = borrowed,
    borrowed_names = unique(as.character(borrowed_names)),
    rows = native$rows, seed = native$seed
  )
}
#' Every record backing (cell_type, gene) at this context, in the served order.
scma_source_records <- function(context, cell_type, gene) {
  if (identical(context$kind, "composite")) {
    records <- scma_source_records(context$native, cell_type, gene)
    if (as.character(cell_type) %in% context$borrowed_names) {
      seen <- vapply(records, function(r) as.character(r$sentence), "")
      extra <- Filter(
        function(r) !(as.character(r$sentence) %in% seen),
        scma_source_records(context$borrowed, cell_type, gene)
      )
      records <- c(records, extra)
    }
    return(records)
  }
  rows <- context$rows
  if (!nrow(rows)) {
    return(list())
  }
  key_label <- as.character(cell_type)
  key_gene <- toupper(as.character(gene))
  found <- rows[.(key_label, key_gene), nomatch = NULL]
  if (!nrow(found)) {
    return(list())
  }
  lapply(seq_len(nrow(found)), function(i) {
    list(
      pmid = as.character(found$pmid[i]),
      pmcid = as.character(found$pmcid[i]),
      sentence = as.character(found$source[i]),
      n_pub = as.integer(found$n_pub_support[i])
    )
  })
}
#' The first k records, for a packet that ships sentences rather than asking for them.
scma_source_head <- function(context, cell_type, gene, k) {
  every <- scma_source_records(context, cell_type, gene)
  if (!length(every)) {
    return(list())
  }
  every[seq_len(min(as.integer(k), length(every)))]
}
#' Serves one cluster's sentences without handing back what it already gave.
#'
#' State is per cluster because "already shown" is a fact about one conversation; sharing
#' it would let one cluster's reading exhaust another's.
scma_source_server <- function(context, batch = SOURCES_PER_MARKER,
                               max_batches = SOURCE_BATCHES_PER_MARKER) {
  self <- new.env(parent = emptyenv())
  self$context <- context
  self$batch <- as.integer(batch)
  self$max_batches <- as.integer(max_batches)
  self$served <- new.env(parent = emptyenv())
  self$draws <- new.env(parent = emptyenv())
  self$key <- function(cell_type, gene) {
    paste0(as.character(cell_type), "\u001f", toupper(as.character(gene)))
  }
  self$note_delivered <- function(cell_type, gene, records) {
    key <- self$key(cell_type, gene)
    prior <- if (exists(key, envir = self$served, inherits = FALSE)) {
      get(key, envir = self$served, inherits = FALSE)
    } else {
      character(0)
    }
    sentences <- vapply(records, function(r) as.character(r$sentence), character(1))
    assign(key, unique(c(prior, sentences)), envir = self$served)
    invisible(NULL)
  }
  self$opening <- function(cell_type, gene, k = self$batch) {
    records <- scma_source_head(self$context, cell_type, gene, k)
    self$note_delivered(cell_type, gene, records)
    records
  }
  self$take <- function(cell_type, gene) {
    key <- self$key(cell_type, gene)
    every <- scma_source_records(self$context, cell_type, gene)
    served <- if (exists(key, envir = self$served, inherits = FALSE)) {
      get(key, envir = self$served, inherits = FALSE)
    } else {
      character(0)
    }
    drawn <- if (exists(key, envir = self$draws, inherits = FALSE)) {
      get(key, envir = self$draws, inherits = FALSE)
    } else {
      0L
    }
    remaining <- Filter(function(r) !(as.character(r$sentence) %in% served), every)
    if (self$max_batches > 0L && drawn >= self$max_batches) {
      return(list(
        sources = list(),
        remaining = length(remaining),
        exhausted = !length(remaining),
        limit_reached = TRUE
      ))
    }
    take_n <- min(self$batch, length(remaining))
    batch <- if (take_n > 0L) remaining[seq_len(take_n)] else list()
    self$note_delivered(cell_type, gene, batch)
    assign(key, drawn + 1L, envir = self$draws)
    list(
      sources = batch,
      remaining = length(remaining) - take_n,
      exhausted = length(remaining) <= take_n,
      limit_reached = FALSE
    )
  }
  self
}
