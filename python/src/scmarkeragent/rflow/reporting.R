suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
  library(Matrix)
})
.rep_script_dir <- function() {
  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      value <- frame$ofile
      if (is.null(value)) "" else as.character(value)[1]
    },
    ""
  )
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files)) {
    return(dirname(normalizePath(frame_files[length(frame_files)])))
  }
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep("^--file=", args, value = TRUE)
  if (length(match)) {
    return(dirname(normalizePath(sub("^--file=", "", match[1]))))
  }
  normalizePath(getwd())
}
.REP_SELF <- .rep_script_dir()
.REP_PACKAGE <- dirname(.REP_SELF)
if (!exists("CFG")) {
  source(file.path(.REP_SELF, "config.R"))
}
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
if (!exists("uberon_load", mode = "function")) {
  source(file.path(.REP_SELF, "uberon_ontology.R"))
}
if (!exists("scma_source_context", mode = "function")) {
  source(file.path(.REP_SELF, "marker_sources.R"))
}
if (!exists("sig_pass", mode = "function")) {
  source(file.path(.REP_SELF, "evidence_gate.R"))
}
.REP_SCHEMA <- fromJSON(
  file.path(.REP_PACKAGE, "schemas", "output_schema.json"),
  simplifyVector = FALSE
)
.REP_DEFAULTS <- fromJSON(
  file.path(.REP_PACKAGE, "config", "defaults.json"),
  simplifyVector = FALSE
)
.REP_DECIMALS <- .REP_DEFAULTS$output$decimals
ND_CONF <- as.integer(.REP_DECIMALS$confidence)
ND_FRAC <- as.integer(.REP_DECIMALS$fraction)
ND_S <- as.integer(.REP_DECIMALS$specificity)
ND_LOGFC <- as.integer(.REP_DECIMALS$logfc)
TABLE_A_COLS <- as.character(unlist(.REP_SCHEMA$cluster_summary$columns))
TABLE_B_COLS <- as.character(unlist(.REP_SCHEMA$marker_evidence$columns))
CLUSTER_EVIDENCE_FILE <- as.character(.REP_SCHEMA$cluster_evidence$file)
SUMMARY_MARKER_CAP <- 12L
DOTPLOT_COLS <- c(
  "gene", "marker_slot", "gene_order", "gene_group", "gene_group_order",
  "cluster", "cluster_order", "cluster_celltype", "pct_exp", "avg_exp_scaled"
)
#' Numeric clusters sort as numbers, everything else as text after them.
#'
#' Returned as a sortable key rather than applied directly, because several builders order
#' by cluster and they must all order the same way.
.cluster_sort_key <- function(clusters) {
  text <- as.character(clusters)
  numeric_like <- grepl("^-?[0-9]+$", text)
  order(
    ifelse(numeric_like, 0L, 1L),
    ifelse(numeric_like, suppressWarnings(as.numeric(text)), 0),
    text
  )
}
#' The SHORTEST decimal that reads back as this double.
#'
#' p-values are written in the double's own notation rather than rounded to a fixed number
#' of places, which would read every strongly raised marker as exactly zero.
#'
#' The search asks R to parse each candidate back, and R's `as.numeric` is NOT correctly
#' rounded -- it lands 1 ULP off on about 0.014% of real values, so this can pick a string
#' one digit shorter than a correctly-rounded parser would. That is why the parity bar is
#' per field with numbers compared AS DOUBLES: the two strings denote the same double, and
#' the parity comparator passes them. Implementing correctly-rounded parsing here would
#' buy nothing a reader can see.
scma_shortest_double <- function(x) {
  if (is.null(x) || is.na(x) || !is.finite(x)) {
    return(NOT_AVAILABLE)
  }
  if (x == 0) {
    return("0.0")
  }
  significant <- 17L
  for (p in 1:17) {
    if (as.numeric(sprintf(paste0("%.", p, "g"), x)) == x) {
      significant <- p
      break
    }
  }
  canonical <- sprintf(paste0("%.", significant - 1L, "e"), x)
  pieces <- strsplit(canonical, "e", fixed = TRUE)[[1]]
  exponent <- as.integer(pieces[2])
  if (exponent >= -4L && exponent < 16L) {
    decimals <- max(0L, significant - 1L - exponent)
    out <- sprintf(paste0("%.", decimals, "f"), x)
    if (!grepl(".", out, fixed = TRUE)) {
      out <- paste0(out, ".0")
    }
    return(out)
  }
  mantissa <- sub("\\.$", "", sub("0+$", "", pieces[1]))
  paste0(
    mantissa, "e",
    if (exponent < 0L) "-" else "+",
    sprintf("%02d", abs(exponent))
  )
}
.fmt_frac <- function(x) sprintf(paste0("%.", ND_FRAC, "f"), as.numeric(x))
#' An adjusted p-value for a table cell, or the not-applicable token.
.fmt_pvalue <- function(value) {
  if (is.null(value) || !length(value)) {
    return(NOT_AVAILABLE)
  }
  number <- suppressWarnings(as.numeric(value[1]))
  if (is.na(number) || !is.finite(number)) {
    return(NOT_AVAILABLE)
  }
  scma_shortest_double(number)
}
#' A number formatted for a table cell, or the not-applicable token.
.fmt <- function(value, decimals) {
  if (is.null(value) || !length(value)) {
    return(NOT_AVAILABLE)
  }
  number <- suppressWarnings(as.numeric(value[1]))
  if (is.na(number) || !is.finite(number)) {
    return(NOT_AVAILABLE)
  }
  sprintf(paste0("%.", as.integer(decimals), "f"), number)
}
#' The other identities the cluster also carries, ';'-joined.
#'
#' What a reader needs from a summary row is which OTHER cells are in this cluster. The
#' full candidate set, each with its claim role and retrieval rank, is in
#' `cluster_evidence.jsonl`, where it can be queried instead of widening this column past
#' the point a reader can use it.
.fmt_alternative_candidates <- function(names) {
  parts <- trimws(as.character(unlist(names %||% list())))
  parts <- parts[nzchar(parts)]
  if (!length(parts)) NOT_AVAILABLE else paste(parts, collapse = "; ")
}
#' A named-vector lookup that yields NULL for an absent key.
#'
#' `x[[key]]` raises on a missing name and `x[key]` yields an NA element carrying an NA
#' name, either of which would turn "this label has no ontology id" into a crash or into a
#' cell reading `NA` instead of the shared token.
.lookup <- function(mapping, key) {
  if (is.null(mapping) || !length(mapping) || is.null(key) || !nzchar(key)) {
    return(NULL)
  }
  if (!(key %in% names(mapping))) {
    return(NULL)
  }
  mapping[[key]]
}
#' The audit's short note on claims whose curated markers this cluster does not raise.
#'
#' It is a flag for a reader, not a correction: the annotation is exactly what the agent
#' decided, and this column says which of its claims rest on an identity whose
#' best-corroborated markers sit unraised here.
.fmt_claim_warnings <- function(lines) {
  parts <- trimws(as.character(unlist(lines %||% list())))
  parts <- parts[nzchar(parts)]
  if (!length(parts)) NOT_AVAILABLE else paste(parts, collapse = " || ")
}
#' Load the saved run artifacts every builder reads. Nothing is re-decided here.
#'
#' The `.rds` counterpart of `_load_run_artifacts`. The sparse matrix is kept sparse: the
#' whole point of reading it here rather than through a bridge is that a genome-wide
#' matrix must never be materialised dense.
scma_load_run_artifacts <- function(tag, cache = CACHE) {
  read_one <- function(kind) {
    path <- file.path(cache, sprintf("%s_%s.rds", tag, kind))
    if (!file.exists(path)) {
      stop(sprintf("missing run artifact: %s", path), call. = FALSE)
    }
    readRDS(path)
  }
  annotations <- read_one("annotations")
  norm <- read_one("norm")
  if (!inherits(norm, "CsparseMatrix")) {
    norm <- as(as(norm, "CsparseMatrix"), "generalMatrix")
  }
  list(
    tag = tag,
    dm = read_one("de_meta"),
    scoring = read_one("candidate_scoring"),
    annotations = annotations,
    res = annotations$results,
    norm = norm
  )
}
#' Per-cluster annotation intermediate shared with the figure bundle.
#'
#' Both arms consume these labels rather than re-deriving them, so no plot can disagree
#' with the delivered table.
scma_build_clustermap <- function(fr) {
  res <- fr$res
  cl_id <- fr$scoring$candidate_cl_id %||% character(0)
  clusters <- names(res)
  clusters <- clusters[.cluster_sort_key(clusters)]
  rows <- lapply(clusters, function(cluster) {
    record <- res[[cluster]]
    label <- as.character(record$annotation %||% "")[1]
    if (!nzchar(label)) label <- NOT_AVAILABLE
    subtype <- as.character(record$subtype %||% "")[1]
    primary <- if (nzchar(subtype)) subtype else label
    data.table(
      cluster = as.character(cluster),
      n_cells = as.integer(fr$scoring$clusters[[as.character(cluster)]]$n_cells),
      cell_type_annotation = label,
      cell_subtype_annotation = na_display(record$subtype),
      primary_annotation = primary,
      cell_state = na_display(record$state),
      cell_ontology = na_display(.lookup(cl_id, label)),
      annotation_confidence = na_display(record$confidence),
      annotation_rationale = na_display(record$rationale),
      resolution_status = as.character(record$resolution_status %||% "")[1],
      resolution_detail = as.character(record$resolution_detail %||% "")[1],
      annotation_source = as.character(record$annotation_source %||% "")[1],
      llm_status = as.character(record$llm_status %||% "")[1],
      annotation_qc = as.character(record$annotation_qc %||% "")[1],
      cluster_celltype = sprintf("%s: %s", cluster, primary),
      alternative_candidates = .fmt_alternative_candidates(
        record$co_occurring_identities
      ),
      claim_warnings = .fmt_claim_warnings(record$claim_warnings),
      evidence_review_rounds = length((record$review %||% list())$rounds %||% list()),
      delivered_tier = as.integer(record$delivered_tier %||% 0L),
      arbitrated = isTRUE(length(record$arbitration) > 0L)
    )
  })
  clustermap <- rbindlist(rows)
  clustermap[, cluster_order := seq_len(.N) - 1L]
  setcolorder(clustermap, c("cluster", "cluster_order"))
  clustermap[]
}
#' What part this gene played in the claim the annotating agent made.
#'
#' `decisive` is the gene the agent quoted for that identity, `support` a gene it named
#' behind the selected label, and `panel` a curated marker of the identity that the
#' cluster raises without either having been cited.
.marker_role <- function(entry, result, gene) {
  upper <- toupper(as.character(gene))
  decisive <- toupper(as.character(entry$claim_evidence$decisive_gene %||% "")[1])
  if (identical(decisive, upper) && nzchar(upper)) {
    return("decisive")
  }
  support <- toupper(as.character(unlist(result$support_markers %||% list())))
  if (upper %in% support) {
    return("support")
  }
  "panel"
}
#' Table B -- the measured, source-backed markers behind every identity claimed.
#'
#' A row exists for each raised curated marker of every identity the agent asserted the
#' cluster carries: the assigned one, its subtype, and each co-occurring identity. The
#' candidates considered and not claimed keep their complete panels in
#' `cluster_evidence.jsonl`.
#'
#' One kind of unclaimed candidate does put rows here: a sibling in the assigned identity's
#' own group, and only its DETECTED EXCLUSIONS. The audit turn binds those to the assigned
#' identity, and the delivered reason then argues from them -- of the exclusions a reason
#' named, two thirds were curated under a group sibling rather than under the assigned
#' name, so leaving them out kept most of the exclusion argument unverifiable.
#'
#' A claimed identity's panel is not always native to this run's own tissue and species:
#' `scma_borrow_candidates` places a name on a cluster's page from another tissue of the
#' same species, or (`features.cross_species_borrow`) another species entirely, when the
#' cluster's own context does not cover it. Its sentences live under that other context,
#' not this one, so the default `sources` built below is widened for exactly those names
#' -- the same way `scma_prepare_pool` already widens `pool$species_sources` for the live
#' conversation -- to every tissue of this species, plus, for a name that crossed
#' species, the donor's own rows translated through the same one-to-one ortholog map
#' that built its panel. A name nothing here claims stays answered from the native,
#' tissue-scoped slice alone.
.borrowed_context_names <- function(fr) {
  names_found <- character(0)
  donor_species_by_name <- list()
  for (result in fr$res) {
    for (entry in result$candidate_entries %||% list()) {
      borrowed <- entry$borrowed_context
      if (is.null(borrowed) || !length(borrowed)) next
      nm <- as.character(entry$cell_type)
      names_found <- c(names_found, nm)
      donor <- borrowed$donor_species
      if (!is.null(donor) && nzchar(as.character(donor))) {
        donor_species_by_name[[nm]] <- union(donor_species_by_name[[nm]], as.character(donor))
      }
    }
  }
  list(
    names = sort(.cp(unique(names_found)), method = "radix"),
    donor_species_by_name = donor_species_by_name
  )
}
scma_build_marker_evidence <- function(tag, fr, sources = NULL) {
  context <- fr$scoring$context
  if (is.null(sources)) {
    native_sources <- scma_source_context(
      uberon_load(OBO_UBERON), DB_SOURCES,
      context$species, context$tissue, context$disease
    )
    borrowed <- .borrowed_context_names(fr)
    if (length(borrowed$names)) {
      extra <- scma_source_types_across_tissues(
        DB_SOURCES, context$species, borrowed$names, context$disease,
        donor_species_by_name = borrowed$donor_species_by_name, ortho_dir = ORTHO_DIR
      )
      sources <- scma_composite_sources(native_sources, extra, borrowed$names)
    } else {
      sources <- native_sources
    }
  }
  specificity <- fr$scoring$marker_specificity
  de <- as.data.table(fr$dm$de)
  padj_key <- paste(as.character(de$group), toupper(as.character(de$feature)), sep = "\x1f")
  padj_of <- setNames(as.numeric(de$padj), padj_key)
  clusters <- names(fr$res)
  clusters <- clusters[.cluster_sort_key(clusters)]
  rows <- list()
  for (cluster in clusters) {
    result <- fr$res[[cluster]]
    assigned <- unique(c(
      as.character(result$annotation %||% "")[1],
      as.character(result$subtype %||% "")[1]
    ))
    assigned <- assigned[nzchar(assigned)]
    reported <- as.character(unlist(result$possible_refinements %||% list()))
    for (entry in result$candidate_entries %||% list()) {
      candidate <- as.character(entry$cell_type)[1]
      is_selected <- candidate %in% assigned
      claimed <- is_selected || nzchar(as.character(entry$claim_role %||% "")[1])
      if (!claimed && !(candidate %in% reported)) next
      for (marker in entry$decisive_marker_measurements %||% list()) {
        polarity <- as.character(marker$polarity %||% "positive")[1]
        if (!claimed && !identical(polarity, "negative")) next
        gene <- as.character(marker$gene)[1]
        records <- scma_source_head(sources, candidate, gene, 1L)
        source_row <- if (length(records)) records[[1]] else list()
        rows[[length(rows) + 1L]] <- data.table(
          dataset = tag,
          cluster_id = as.character(cluster),
          candidate_annotation = candidate,
          is_selected_annotation = if (is_selected) "True" else "False",
          gene = gene,
          marker_polarity = polarity,
          marker_role = .marker_role(entry, result, gene),
          marker_provenance = "db_cited",
          cluster_detection_fraction = .fmt(marker$detection_fraction_in, ND_FRAC),
          out_of_cluster_detection_fraction = .fmt(
            marker$detection_fraction_out, ND_FRAC
          ),
          average_log2_fold_change = .fmt(marker$avg_log2FC, ND_LOGFC),
          auc = .fmt(marker$auc, ND_FRAC),
          adjusted_p_value = .fmt_pvalue(
            .lookup(
              padj_of,
              paste(as.character(cluster), toupper(gene), sep = "\x1f")
            )
          ),
          cross_cluster_percentile = .fmt(marker$cross_cluster_percentile, ND_FRAC),
          marker_specificity_weight = .fmt(
            .lookup(specificity, toupper(gene)), ND_S
          ),
          publication_support_count = {
            value <- marker$publication_support
            if (is.null(value) || !length(value) || is.na(value[1])) {
              NOT_AVAILABLE
            } else {
              format(value[1], scientific = FALSE, trim = TRUE)
            }
          },
          evidence_tier = na_display(marker$evidence_tier),
          pmid = na_display(source_row$pmid),
          pmcid = na_display(source_row$pmcid),
          source_sentence = na_display(source_row$sentence)
        )
      }
    }
  }
  if (!length(rows)) {
    empty <- data.table()
    for (column in TABLE_B_COLS) empty[, (column) := character(0)]
    return(empty[])
  }
  table_b <- rbindlist(rows)
  setcolorder(table_b, TABLE_B_COLS)
  cluster_rank <- match(
    as.character(table_b$cluster_id), unique(as.character(clusters))
  )
  table_b <- table_b[order(
    cluster_rank,
    -xtfrm(table_b$is_selected_annotation),
    .cp(table_b$candidate_annotation),
    .cp(table_b$gene),
    method = "radix"
  )]
  unique(table_b, by = c("cluster_id", "candidate_annotation", "gene"))[]
}
#' Rank a cluster's marker rows so a cap keeps the strongest, best-corroborated ones.
#'
#' Table B is written gene-ascending because it is a full record and a reader looks genes
#' up in it. Table A shows twelve, and taking the first twelve of that alphabetical order
#' dropped PECAM1 and VWF from an `endothelial cell` cluster in favour of a one-publication
#' gene, because P and V sort after F. So the summary is ordered the way every other panel
#' in this package is: publication support first, contrast to break ties.
.by_evidence <- function(frame) {
  ranked <- frame
  if ("marker_polarity" %in% names(ranked)) {
    ranked <- ranked[as.character(marker_polarity) != "negative"]
  }
  if (!nrow(ranked)) {
    return(ranked)
  }
  support <- suppressWarnings(as.numeric(ranked$publication_support_count))
  support[is.na(support)] <- 0
  detected <- suppressWarnings(as.numeric(ranked$cluster_detection_fraction))
  detected[is.na(detected)] <- 0
  outside <- suppressWarnings(as.numeric(ranked$out_of_cluster_detection_fraction))
  outside[is.na(outside)] <- 0
  ranked[order(-support, -(detected - outside), .cp(ranked$gene), method = "radix")]
}
.format_summary_markers <- function(rows, tag_identity) {
  if (!length(rows) || !nrow(rows)) {
    return(list(markers = "", pmcids = ""))
  }
  markers <- if (isTRUE(tag_identity)) {
    paste(
      sprintf(
        "%s|%s|%s", rows$gene, rows$candidate_annotation, rows$cluster_detection_fraction
      ),
      collapse = "; "
    )
  } else {
    paste(
      sprintf("%s(%s)", rows$gene, rows$cluster_detection_fraction),
      collapse = "; "
    )
  }
  pmcids <- as.character(rows$pmcid)
  pmcids <- pmcids[!(pmcids %in% c("", NOT_AVAILABLE))]
  list(markers = markers, pmcids = paste(unique(pmcids), collapse = "; "))
}
#' Shortlist rows for one summary column, fair-shared across the identities in `group`.
#'
#' A single global corroboration sort lets a high-n_pub parent -- or one co-occurring
#' identity -- fill the whole cap and hide the markers that actually define a subtype or a
#' second population. Each identity gets a share of the cap first; leftovers fill in the
#' same corroboration order.
.summary_marker_rows <- function(group, tag_identity) {
  if (is.null(group) || !nrow(group)) {
    return(list(markers = "", pmcids = ""))
  }
  identities <- unique(as.character(group$candidate_annotation))
  ranked_by_identity <- lapply(identities, function(identity) {
    .by_evidence(group[as.character(candidate_annotation) == identity])
  })
  names(ranked_by_identity) <- identities
  share <- max(1L, SUMMARY_MARKER_CAP %/% length(identities))
  chosen <- list()
  seen <- character(0)
  take <- function(rows, from, to) {
    if (!nrow(rows) || from > nrow(rows)) {
      return(FALSE)
    }
    for (i in seq(from, min(to, nrow(rows)))) {
      gene <- as.character(rows$gene[i])
      if (gene %in% seen) next
      seen <<- c(seen, gene)
      chosen[[length(chosen) + 1L]] <<- rows[i]
      if (length(chosen) >= SUMMARY_MARKER_CAP) {
        return(TRUE)
      }
    }
    FALSE
  }
  full <- FALSE
  for (identity in identities) {
    if (take(ranked_by_identity[[identity]], 1L, share)) {
      full <- TRUE
      break
    }
  }
  if (!full) {
    for (identity in identities) {
      rows <- ranked_by_identity[[identity]]
      if (take(rows, share + 1L, nrow(rows))) break
    }
  }
  .format_summary_markers(rbindlist(chosen), tag_identity)
}
#' Names listed in alternative_candidates -- the only identities the cooc columns may show.
#'
#' marker_evidence also keeps `rejected_subtype` rows and similar non-selected claims.
#' Those are audit trail, not co-occurring populations; putting them in
#' `cooccurring_markers` would invent a mix the delivery never declared.
.cooccurring_names <- function(alternative_candidates) {
  text <- trimws(as.character(alternative_candidates %||% "")[1])
  if (!nzchar(text) || identical(text, NOT_AVAILABLE)) {
    return(character(0))
  }
  parts <- trimws(strsplit(text, ";", fixed = TRUE)[[1]])
  unique(parts[nzchar(parts) & parts != NOT_AVAILABLE])
}
#' Table A -- concise reader-facing summary, one row per dataset and cluster.
scma_build_table_a <- function(tag, clustermap, marker_evidence) {
  out <- lapply(seq_len(nrow(clustermap)), function(i) {
    row <- clustermap[i]
    cluster <- as.character(row$cluster)
    group <- marker_evidence[as.character(cluster_id) == cluster]
    key <- list(markers = "", pmcids = "")
    cooc <- list(markers = "", pmcids = "")
    if (nrow(group)) {
      selected <- group[as.character(is_selected_annotation) == "True"]
      key <- .summary_marker_rows(
        if (nrow(selected)) selected else group, tag_identity = FALSE
      )
      names_wanted <- .cooccurring_names(row$alternative_candidates)
      if (length(names_wanted)) {
        cooc <- .summary_marker_rows(
          group[as.character(candidate_annotation) %in% names_wanted],
          tag_identity = TRUE
        )
      }
    }
    data.table(
      dataset = tag,
      cluster_id = cluster,
      n_cells = as.integer(row$n_cells),
      primary_annotation = row$primary_annotation,
      cell_type_annotation = row$cell_type_annotation,
      cell_subtype_annotation = row$cell_subtype_annotation,
      cell_state = row$cell_state,
      cell_ontology = row$cell_ontology,
      annotation_confidence = row$annotation_confidence,
      annotation_rationale = row$annotation_rationale,
      resolution_status = row$resolution_status,
      annotation_source = row$annotation_source,
      llm_status = row$llm_status,
      annotation_qc = row$annotation_qc,
      key_markers = na_display(key$markers),
      pmcid = na_display(key$pmcids),
      cooccurring_markers = na_display(cooc$markers),
      cooccurring_pmcid = na_display(cooc$pmcids),
      alternative_candidates = row$alternative_candidates,
      claim_warnings = row$claim_warnings,
      evidence_review_rounds = row$evidence_review_rounds,
      delivered_tier = row$delivered_tier,
      arbitrated = row$arbitrated
    )
  })
  table_a <- rbindlist(out)
  setcolorder(table_a, TABLE_A_COLS)
  table_a[]
}
#' Correctly-rounded decimal rounding of the double.
#'
#' R's own `round` is close but not the same function; the C library's `%.nf`
#' conversion IS correctly-rounded decimal rounding, so the round trip through
#' `sprintf` reproduces it digit for digit.
.round_decimal <- function(x, digits) {
  out <- suppressWarnings(as.numeric(sprintf(paste0("%.", digits, "f"), as.numeric(x))))
  out[!is.finite(as.numeric(x))] <- as.numeric(x)[!is.finite(as.numeric(x))]
  out
}
#' The POPULATION standard deviation, dividing by n. R's `sd` divides by n-1,
#' which would scale every dotplot z-value on small cluster counts.
.pop_sd <- function(x) {
  x <- as.numeric(x)
  sqrt(mean((x - mean(x))^2))
}
.JSON_EMPTY_OBJECT <- structure(list(), scma_json = "object")
.json_escape <- function(text) {
  text <- enc2utf8(as.character(text))
  text <- gsub("\\", "\\\\", text, fixed = TRUE)
  text <- gsub("\"", "\\\"", text, fixed = TRUE)
  text <- gsub("\n", "\\n", text, fixed = TRUE)
  text <- gsub("\r", "\\r", text, fixed = TRUE)
  text <- gsub("\t", "\\t", text, fixed = TRUE)
  if (grepl("[\x01-\x1f]", text, useBytes = TRUE)) {
    for (code in c(1:7, 11L, 14:31)) {
      text <- gsub(intToUtf8(code), sprintf("\\u%04x", code), text, fixed = TRUE)
    }
    text <- gsub("\b", "\\b", text, fixed = TRUE)
    text <- gsub("\f", "\\f", text, fixed = TRUE)
  }
  paste0("\"", text, "\"")
}
.json_scalar <- function(value) {
  if (is.null(value)) {
    return("null")
  }
  if (length(value) != 1L) {
    stop("json scalar of length != 1")
  }
  if (is.na(value) && !is.nan(value)) {
    return("null")
  }
  if (is.logical(value)) {
    return(if (isTRUE(value)) "true" else "false")
  }
  if (is.integer(value)) {
    return(format(value, scientific = FALSE, trim = TRUE))
  }
  if (is.numeric(value)) {
    if (!is.finite(value)) {
      return("null")
    }
    return(scma_shortest_double(value))
  }
  .json_escape(value)
}
#' Serialize an R structure to the canonical JSON every arm writes.
scma_json_canonical <- function(value) {
  marker <- attr(value, "scma_json", exact = TRUE)
  if (is.null(value)) {
    return("null")
  }
  if (is.list(value)) {
    is_object <- identical(marker, "object") || !is.null(names(value))
    if (is_object) {
      keys <- names(value) %||% character(0)
      keys <- sort(.cp(keys), method = "radix")
      parts <- vapply(
        keys,
        function(key) {
          paste0(.json_escape(key), ":", scma_json_canonical(value[[key]]))
        },
        ""
      )
      return(paste0("{", paste(parts, collapse = ","), "}"))
    }
    parts <- vapply(unname(value), scma_json_canonical, "")
    return(paste0("[", paste(parts, collapse = ","), "]"))
  }
  if (is.factor(value)) {
    value <- as.character(value)
  }
  if (is.atomic(value)) {
    if (length(value) == 0L) {
      return(if (identical(marker, "object")) "{}" else "[]")
    }
    if (length(value) == 1L && !identical(marker, "array")) {
      return(.json_scalar(value))
    }
    parts <- vapply(seq_along(value), function(i) .json_scalar(value[[i]]), "")
    return(paste0("[", paste(parts, collapse = ","), "]"))
  }
  stop(sprintf("cannot serialize object of class %s", paste(class(value), collapse = "/")))
}
#' Write a delivered table, with no row-index column.
#'
#' The parity comparator PARSES both files, so quoting style is free; what has to hold is
#' the parsed content. Logical columns are written as `True` / `False`, and any numeric
#' column whose full precision matters must already be a character column --
#' `scma_shortest_double` is how it gets there.
.write_csv <- function(table, path) {
  out <- copy(as.data.table(table))
  for (column in names(out)) {
    if (is.logical(out[[column]])) {
      set(out, j = column, value = ifelse(out[[column]], "True", "False"))
    }
  }
  fwrite(out, path, quote = "auto", na = "", eol = "\n")
  normalizePath(path)
}
#' Reference-free per-cell figure data: stable ID, cluster, and UMAP.
scma_build_cells <- function(fr) {
  meta <- fr$dm$meta
  cells <- data.table(
    cell = as.character(meta$cell),
    cluster = as.character(meta$cluster)
  )
  umap <- fr$dm$umap
  if (!is.null(umap)) {
    umap <- as.matrix(umap)
    cells[, umap_x := .round_decimal(umap[, 1], 6L)]
    cells[, umap_y := .round_decimal(umap[, 2], 6L)]
  } else {
    cells[, umap_x := NA_real_]
    cells[, umap_y := NA_real_]
  }
  cells[]
}
#' The markers that support each owner's assigned identity, one slot per gene.
#'
#' Genes are deduplicated across owners: the first owner in display order keeps the slot
#' and later owners read that same column. `owner_of` is a NAMED character vector from
#' cluster id to owning block; publication support orders the genes inside a block.
.supporting_panel <- function(marker_evidence, owner_of) {
  empty <- list(slots = list(), publications = numeric(0))
  if (!nrow(marker_evidence)) {
    return(empty)
  }
  panel <- marker_evidence[tolower(as.character(is_selected_annotation)) == "true"]
  if ("marker_polarity" %in% names(panel)) {
    panel <- panel[as.character(marker_polarity) != "negative"]
  }
  if (!nrow(panel)) {
    return(empty)
  }
  panel <- copy(panel)
  panel[, owner := unname(owner_of[as.character(cluster_id)])]
  panel <- panel[!is.na(owner)]
  panel[, gene_upper := toupper(as.character(gene))]
  panel[, support := {
    value <- suppressWarnings(as.numeric(publication_support_count))
    value[is.na(value)] <- 0
    value
  }]
  slots <- list()
  publications <- numeric(0)
  seen <- character(0)
  for (owner_name in unique(unname(owner_of))) {
    group <- panel[owner == owner_name]
    if (!nrow(group)) next
    ranked <- group[, .(support = max(support)), by = gene_upper]
    ranked <- ranked[order(-support, .cp(gene_upper), method = "radix")]
    for (i in seq_len(nrow(ranked))) {
      gene_name <- as.character(ranked$gene_upper[i])
      publications[[paste(owner_name, gene_name, sep = "\x1f")]] <-
        as.integer(ranked$support[i])
      if (gene_name %in% seen) next
      seen <- c(seen, gene_name)
      slots[[length(slots) + 1L]] <- list(owner = owner_name, gene = gene_name)
    }
  }
  list(slots = slots, publications = publications)
}
#' Shared core of the two dotplot builders: per-slot detection and scaled expression.
#'
#' `groups` is an ordered list of (name -> integer cell indices); expression statistics
#' are computed over each group's pooled cells.
.dotplot_stats <- function(norm_cells_by_gene, gene_column, groups) {
  expression <- as.numeric(norm_cells_by_gene[, gene_column])
  averages <- vapply(
    groups,
    function(idx) if (length(idx)) mean(expm1(expression[idx])) else 0,
    0
  )
  percentages <- vapply(
    groups,
    function(idx) if (length(idx)) 100 * mean(expression[idx] > 0) else 0,
    0
  )
  deviation <- .pop_sd(averages)
  scaled <- if (deviation > 0) (averages - mean(averages)) / deviation else averages * 0
  scaled <- pmin(pmax(scaled, -2.5), 2.5)
  list(pct = percentages, scaled = scaled)
}
#' Build a cluster-aligned dotplot from source-backed identity markers.
scma_build_identity_marker_dotplot <- function(fr, clustermap, marker_evidence) {
  menu <- as.character(rownames(fr$norm))
  order_map <- as.data.table(clustermap)[order(cluster_order)]
  cluster_order <- as.character(order_map$cluster)
  cluster_labels <- setNames(as.character(order_map$cluster_celltype), cluster_order)
  empty <- data.table(
    gene = character(0), marker_slot = character(0), gene_order = integer(0),
    gene_group = character(0), gene_group_order = integer(0), cluster = character(0),
    cluster_order = integer(0), cluster_celltype = character(0),
    pct_exp = numeric(0), avg_exp_scaled = numeric(0)
  )
  if (!nrow(marker_evidence)) {
    return(empty)
  }
  candidates <- copy(as.data.table(marker_evidence))
  candidates[, cluster_id := as.character(cluster_id)]
  candidates <- candidates[toupper(as.character(gene)) %in% toupper(menu)]
  panel <- .supporting_panel(
    candidates, setNames(cluster_order, cluster_order)
  )
  if (!length(panel$slots)) {
    return(empty)
  }
  meta_cluster <- as.character(fr$dm$meta$cluster)
  groups <- lapply(cluster_order, function(cluster) which(meta_cluster == cluster))
  names(groups) <- cluster_order
  norm_t <- Matrix::t(fr$norm[, as.character(fr$dm$meta$cell), drop = FALSE])
  gene_position <- setNames(seq_along(menu), toupper(menu))
  rows <- vector("list", length(panel$slots) * length(cluster_order))
  n <- 0L
  for (slot_index in seq_along(panel$slots)) {
    slot <- panel$slots[[slot_index]]
    stats <- .dotplot_stats(norm_t, gene_position[[slot$gene]], groups)
    owner_label <- cluster_labels[[slot$owner]]
    owner_rank <- match(slot$owner, cluster_order) - 1L
    for (row_order in seq_along(cluster_order)) {
      cluster <- cluster_order[row_order]
      n <- n + 1L
      rows[[n]] <- data.table(
        gene = slot$gene,
        marker_slot = paste0(slot$owner, "\x1f", slot$gene),
        gene_order = slot_index - 1L,
        gene_group = owner_label,
        gene_group_order = owner_rank,
        cluster = cluster,
        cluster_order = row_order - 1L,
        cluster_celltype = cluster_labels[[cluster]],
        pct_exp = .round_decimal(stats$pct[[row_order]], 3L),
        avg_exp_scaled = .round_decimal(stats$scaled[[row_order]], 4L)
      )
    }
  }
  rbindlist(rows)[]
}
CELLTYPE_DOTPLOT_COLS <- c(
  "gene", "marker_slot", "gene_order", "gene_group", "gene_group_order",
  "cell_type", "cell_type_order", "n_clusters", "n_cells", "pct_exp",
  "avg_exp_scaled", "n_pub"
)
#' The identity dotplot with cell types, not clusters, as rows.
#'
#' Clusters that resolved to the SAME cell type are one row, with detection and mean
#' expression recomputed over their pooled cells. `n_pub` is the publication support
#' `marker_evidence.csv` already reports for that (cell type, gene) pair.
scma_build_celltype_marker_dotplot <- function(fr, clustermap, marker_evidence) {
  menu <- as.character(rownames(fr$norm))
  order_map <- as.data.table(clustermap)[order(cluster_order)]
  label_of <- setNames(
    ifelse(
      nzchar(as.character(order_map$primary_annotation)),
      as.character(order_map$primary_annotation),
      as.character(order_map$cluster)
    ),
    as.character(order_map$cluster)
  )
  cell_types <- unique(unname(label_of))
  clusters_of <- lapply(
    cell_types,
    function(cell_type) names(label_of)[unname(label_of) == cell_type]
  )
  names(clusters_of) <- cell_types
  empty <- data.table(
    gene = character(0), marker_slot = character(0), gene_order = integer(0),
    gene_group = character(0), gene_group_order = integer(0), cell_type = character(0),
    cell_type_order = integer(0), n_clusters = integer(0), n_cells = integer(0),
    pct_exp = numeric(0), avg_exp_scaled = numeric(0), n_pub = integer(0)
  )
  if (!nrow(marker_evidence)) {
    return(empty)
  }
  candidates <- copy(as.data.table(marker_evidence))
  candidates[, cluster_id := as.character(cluster_id)]
  candidates <- candidates[toupper(as.character(gene)) %in% toupper(menu)]
  panel <- .supporting_panel(candidates, label_of)
  if (!length(panel$slots)) {
    return(empty)
  }
  meta_cluster <- as.character(fr$dm$meta$cluster)
  groups <- lapply(
    cell_types,
    function(cell_type) which(meta_cluster %in% clusters_of[[cell_type]])
  )
  names(groups) <- cell_types
  norm_t <- Matrix::t(fr$norm[, as.character(fr$dm$meta$cell), drop = FALSE])
  gene_position <- setNames(seq_along(menu), toupper(menu))
  rows <- vector("list", length(panel$slots) * length(cell_types))
  n <- 0L
  for (slot_index in seq_along(panel$slots)) {
    slot <- panel$slots[[slot_index]]
    stats <- .dotplot_stats(norm_t, gene_position[[slot$gene]], groups)
    owner_rank <- match(slot$owner, cell_types) - 1L
    for (row_order in seq_along(cell_types)) {
      cell_type <- cell_types[row_order]
      support <- unname(panel$publications[paste(cell_type, slot$gene, sep = "\x1f")])
      if (!length(support) || is.na(support)) support <- 0L
      n <- n + 1L
      rows[[n]] <- data.table(
        gene = slot$gene,
        marker_slot = paste0(slot$owner, "\x1f", slot$gene),
        gene_order = slot_index - 1L,
        gene_group = slot$owner,
        gene_group_order = owner_rank,
        cell_type = cell_type,
        cell_type_order = row_order - 1L,
        n_clusters = length(clusters_of[[cell_type]]),
        n_cells = length(groups[[row_order]]),
        pct_exp = .round_decimal(stats$pct[[row_order]], 3L),
        avg_exp_scaled = .round_decimal(stats$scaled[[row_order]], 4L),
        n_pub = as.integer(support)
      )
    }
  }
  rbindlist(rows)[]
}
MARKER_LIST_COLS <- c(
  "cluster", "cluster_celltype", "gene", "avg_log2FC", "pct_in", "pct_out", "auc",
  "mean_log_expression", "log_expression_difference", "pval", "padj",
  "in_marker_menu", "significant"
)
#' Per-cluster marker tables: every tested gene, and the significant subset.
#'
#' The source is the genome-wide FindAllMarkers-equivalent table written by
#' preprocessing, NOT the menu-restricted DE that drives annotation. `in_marker_menu`
#' marks the rows that also belong to the menu universe; `significant` applies the shared
#' marker gate to THIS table's genome-wide padj.
scma_build_marker_lists <- function(tag, clustermap, fr, cache = CACHE) {
  path <- file.path(cache, sprintf("%s_markers_all.rds", tag))
  if (!file.exists(path)) {
    stop(sprintf("%s: the genome-wide marker table is missing; rebuild prep", path),
      call. = FALSE
    )
  }
  markers <- as.data.table(readRDS(path))
  empty <- as.data.table(setNames(
    rep(list(character(0)), length(MARKER_LIST_COLS)), MARKER_LIST_COLS
  ))
  if (!nrow(markers)) {
    return(list(all = empty, significant = copy(empty)))
  }
  menu_genes <- toupper(as.character(rownames(fr$norm)))
  labels <- setNames(
    as.character(clustermap$cluster_celltype), as.character(clustermap$cluster)
  )
  out <- data.table(
    cluster = as.character(markers$group),
    gene = as.character(markers$feature),
    avg_log2FC = .round_decimal(as.numeric(markers$avg_log2FC), ND_LOGFC),
    pct_in = .round_decimal(as.numeric(markers$pct_in) / 100, ND_FRAC),
    pct_out = .round_decimal(as.numeric(markers$pct_out) / 100, ND_FRAC),
    auc = .round_decimal(as.numeric(markers$auc), ND_FRAC),
    mean_log_expression = .round_decimal(as.numeric(markers$avgExpr), ND_LOGFC),
    log_expression_difference = .round_decimal(as.numeric(markers$logFC), ND_LOGFC),
    pval = vapply(as.numeric(markers$pval), scma_shortest_double, ""),
    padj = vapply(as.numeric(markers$padj), scma_shortest_double, "")
  )
  out[, cluster_celltype := {
    mapped <- unname(labels[cluster])
    ifelse(is.na(mapped), cluster, mapped)
  }]
  out[, in_marker_menu := toupper(gene) %in% menu_genes]
  out[, significant := sig_pass(
    as.numeric(markers$avg_log2FC),
    as.numeric(markers$pct_in) / 100,
    as.numeric(markers$padj),
    as.numeric(markers$auc)
  )]
  cluster_levels <- unique(out$cluster)
  cluster_levels <- cluster_levels[.cluster_sort_key(cluster_levels)]
  ordering <- order(
    match(out$cluster, cluster_levels),
    -xtfrm(out$significant),
    -out$avg_log2FC,
    method = "radix"
  )
  out <- out[ordering]
  setcolorder(out, MARKER_LIST_COLS)
  list(all = out[], significant = out[significant == TRUE][])
}
#' Build the structured per-cluster audit sidecar kept out of the summary CSV.
#'
#' It carries what the summary cannot: every retrieved candidate with its three relative
#' scores and its complete measured panel, and the grouping the agent used to tell
#' "several names for one cell" from "several different cells".
scma_build_cluster_evidence <- function(tag, fr, marker_evidence) {
  evidence <- as.data.table(marker_evidence)
  marker_records <- list()
  if (nrow(evidence)) {
    for (cluster_id_value in unique(as.character(evidence$cluster_id))) {
      group <- evidence[as.character(cluster_id) == cluster_id_value]
      marker_records[[cluster_id_value]] <- lapply(
        seq_len(nrow(group)),
        function(i) {
          record <- as.list(group[i])
          raw <- record$publication_support_count
          if (grepl("^-?[0-9]+$", raw)) {
            record$publication_support_count <- as.integer(raw)
          }
          record
        }
      )
    }
  }
  scored <- as.data.table(fr$scoring$scored)
  cl_id <- fr$scoring$candidate_cl_id %||% character(0)
  top_candidates <- as.integer(fr$scoring$top_candidates)
  clusters <- names(fr$res)
  clusters <- clusters[.cluster_sort_key(clusters)]
  lapply(clusters, function(cluster) {
    cluster_id_value <- as.character(cluster)
    result <- fr$res[[cluster]]
    label <- as.character(result$annotation %||% "")[1]
    if (!nzchar(label)) label <- NOT_AVAILABLE
    subtype <- as.character(result$subtype %||% "")[1]
    pool <- scored[as.character(scored$cluster) == cluster_id_value]
    state <- fr$scoring$clusters[[cluster_id_value]] %||% list()
    shown <- head(pool, top_candidates)
    retrieval_order <- lapply(seq_len(nrow(shown)), function(i) {
      row <- shown[i]
      list(
        candidate = as.character(row$candidate),
        retrieval_rank = as.integer(row$retrieval_rank),
        retrieval_score = .round_decimal(as.numeric(row$retrieval_score), ND_CONF),
        marker_level = .round_decimal(as.numeric(row$marker_level), ND_CONF),
        cluster_level = .round_decimal(as.numeric(row$cluster_level), ND_CONF),
        single_cell_level = .round_decimal(as.numeric(row$single_cell_level), ND_CONF),
        significant_markers = as.integer(row$hits),
        measured_panel_size = as.integer(row$panel_size)
      )
    })
    review <- result$review %||% list()
    if (!length(review)) review <- .JSON_EMPTY_OBJECT
    arbitration <- result$arbitration %||% list()
    if (!length(arbitration)) arbitration <- .JSON_EMPTY_OBJECT
    program_reading <- result$program_reading %||% list()
    if (!length(program_reading)) program_reading <- .JSON_EMPTY_OBJECT
    definer_audit <- result$definer_audit %||% list()
    if (!length(definer_audit)) definer_audit <- .JSON_EMPTY_OBJECT
    confidence_cap <- result$confidence_cap %||% list()
    if (!length(confidence_cap)) confidence_cap <- .JSON_EMPTY_OBJECT
    consistency_note <- result$consistency_note %||% list()
    if (!length(consistency_note)) consistency_note <- .JSON_EMPTY_OBJECT
    reopened_by_consistency <- result$reopened_by_consistency %||% list()
    if (!length(reopened_by_consistency)) reopened_by_consistency <- .JSON_EMPTY_OBJECT
    list(
      schema_version = "scmarkeragent-cluster-evidence-v2",
      dataset = tag,
      cluster_id = cluster_id_value,
      annotation = list(
        primary_annotation = if (nzchar(subtype)) subtype else label,
        cell_type_annotation = label,
        cell_subtype_annotation = na_display(result$subtype),
        cell_state = na_display(result$state),
        co_occurring_identities = as.list(as.character(
          unlist(result$co_occurring_identities) %||% character(0)
        )),
        cell_ontology = na_display(.lookup(cl_id, label)),
        annotation_confidence = na_display(result$confidence),
        resolution_status = as.character(result$resolution_status %||% "")[1],
        annotation_source = as.character(result$annotation_source %||% "")[1],
        llm_status = as.character(result$llm_status %||% "")[1],
        support_markers = as.list(as.character(
          unlist(result$support_markers) %||% character(0)
        )),
        claim_warnings = as.list(as.character(
          unlist(result$claim_warnings) %||% character(0)
        ))
      ),
      evidence = list(
        annotation_rationale = na_display(result$rationale),
        key_markers = marker_records[[cluster_id_value]] %||% list(),
        candidates = unname(result$candidate_entries %||% list()),
        candidate_facts = result$candidate_facts %||% list()
      ),
      audit = list(
        annotator_turns = as.integer(result$turns %||% 0L),
        tool_calls = unname(result$tool_calls %||% list()),
        delivered_tier = as.integer(result$delivered_tier %||% 0L),
        tiers_available = as.integer(result$tiers_available %||% 0L),
        evidence_review = review,
        review_flags = unname(result$review_flags %||% list()),
        program_reading = program_reading,
        definer_audit = definer_audit,
        contested_rivals = unname(as.list(as.character(unlist(result$contested_rivals %||% list())))),
        contested = unname(as.list(as.character(unlist(result$contested %||% list())))),
        confidence_cap = confidence_cap,
        structural_notes = unname(as.list(as.character(unlist(result$structural_notes %||% list())))),
        consistency_note = consistency_note,
        reopened_by_consistency = reopened_by_consistency,
        arbitration = arbitration,
        resolution_detail = as.character(result$resolution_detail %||% "")[1],
        candidate_pool_size = nrow(pool),
        candidate_pool_size_before_hits_gate = as.integer(
          state$pool_size_before_hits_gate %||% nrow(pool)
        ),
        hits_threshold_applied = as.integer(state$hits_threshold_applied %||% 1L),
        candidates_shown = length(result$candidates %||% character(0)),
        retrieval_order = retrieval_order
      )
    )
  })
}
#' Write the audit sidecar as JSON Lines, one canonical record per line.
scma_write_cluster_evidence <- function(records, outdir) {
  path <- file.path(outdir, CLUSTER_EVIDENCE_FILE)
  lines <- vapply(records, scma_json_canonical, "")
  connection <- file(path, open = "wb")
  on.exit(close(connection))
  writeLines(enc2utf8(lines), connection, sep = "\n", useBytes = TRUE)
  normalizePath(path)
}
#' Write a concise, user-facing run and provenance summary.
scma_write_stats_note <- function(tag, outdir, clustermap, cells, table_a, table_b) {
  counts <- table(as.character(clustermap$resolution_status))
  status_names <- sort(.cp(names(counts)), method = "radix")
  lines <- c(
    sprintf("scMarkerAgent result summary: %s", tag),
    "",
    "Pipeline stages",
    "1. Preprocessing: QC, normalization, Leiden partitioning, genome-wide DE",
    "2. Candidate retrieval: marker-level, cluster-level and single-cell-level",
    "   relative evidence, combined into one retrieval order",
    "3. One annotating agent per cluster, opening with every candidate's complete",
    "   measured marker panel and free to query sources and cross-cluster expression",
    "4. Post-hoc audit: quoted measurements checked against the DE table, and claims",
    "   whose curated markers are unraised flagged without being altered",
    "",
    sprintf("Clusters: %d", nrow(clustermap)),
    sprintf("Cells scored: %d", nrow(cells)),
    sprintf("Cluster summary rows: %d", nrow(table_a)),
    sprintf("Marker evidence rows: %d", nrow(table_b)),
    paste0(
      "Resolution: ",
      paste(
        sprintf("%s %d", status_names, as.integer(counts[status_names])),
        collapse = ", "
      )
    ),
    "",
    sprintf("Structured cluster audit: %s", CLUSTER_EVIDENCE_FILE),
    "Every retrieved candidate, its three relative scores, its complete measured",
    "marker panel and the agent's identity grouping are kept in that JSON Lines",
    "sidecar rather than embedded in cluster_summary.csv.",
    "",
    "Marker evidence is measured in this dataset and carries its curated source. The",
    "retrieval order is a search order over evidence, not a probability and not a",
    "confidence. No author reference labels are used anywhere in this package.",
    sprintf("A field that does not apply to a row is written as %s.", NOT_AVAILABLE)
  )
  path <- file.path(outdir, STATS_FILE)
  connection <- file(path, open = "wb")
  on.exit(close(connection))
  writeLines(enc2utf8(lines), connection, sep = "\n", useBytes = TRUE)
  normalizePath(path)
}
FIGURE_CLUSTERMAP_COLUMNS <- c(
  "cluster", "cluster_order", "n_cells", "cell_type_annotation", "cell_ontology",
  "annotation_confidence", "resolution_status", "annotation_source", "cluster_celltype"
)
#' Build public annotation tables and an auditable evidence-data bundle.
#'
#' Same files and same fields every arm delivers, produced here from the saved `.rds`
#' run artifacts.
scma_generate_report <- function(tag, outdir, cache = CACHE, bundle_only = FALSE) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  bundle <- file.path(outdir, FIGDATA_DIR)
  dir.create(bundle, recursive = TRUE, showWarnings = FALSE)
  fr <- scma_load_run_artifacts(tag, cache)
  clustermap <- scma_build_clustermap(fr)
  table_b <- scma_build_marker_evidence(tag, fr)
  identity_markers <- scma_build_identity_marker_dotplot(fr, clustermap, table_b)
  celltype_markers <- scma_build_celltype_marker_dotplot(fr, clustermap, table_b)
  cells <- scma_build_cells(fr)
  paths <- list()
  write_out <- function(table, path) {
    paths[[basename(path)]] <<- .write_csv(table, path)
  }
  write_out(cells, file.path(bundle, "cells.csv"))
  write_out(
    as.data.table(clustermap)[, FIGURE_CLUSTERMAP_COLUMNS, with = FALSE],
    file.path(bundle, "clustermap.csv")
  )
  write_out(identity_markers, file.path(bundle, "dotplot_celltype_markers.csv"))
  write_out(
    celltype_markers, file.path(bundle, "dotplot_celltype_markers_by_celltype.csv")
  )
  if (!bundle_only) {
    table_a <- scma_build_table_a(tag, clustermap, table_b)
    cluster_evidence <- scma_build_cluster_evidence(tag, fr, table_b)
    marker_lists <- scma_build_marker_lists(tag, clustermap, fr, cache)
    write_out(table_a, file.path(outdir, "cluster_summary.csv"))
    write_out(table_b, file.path(outdir, "marker_evidence.csv"))
    write_out(marker_lists$all, file.path(outdir, "markers_all_by_cluster.csv"))
    write_out(
      marker_lists$significant,
      file.path(outdir, "markers_significant_by_cluster.csv")
    )
    evidence_path <- scma_write_cluster_evidence(cluster_evidence, outdir)
    paths[[basename(evidence_path)]] <- evidence_path
    note_path <- scma_write_stats_note(tag, outdir, clustermap, cells, table_a, table_b)
    paths[[basename(note_path)]] <- note_path
    if (!exists("scma_plot_dataset", mode = "function")) {
      source(file.path(.REP_SELF, "visualization.R"))
    }
    if (!exists("scma_build_dataset_viewer", mode = "function")) {
      source(file.path(.REP_SELF, "viewer.R"))
    }
    scma_plot_dataset(tag, outdir)
    viewer_paths <- scma_build_dataset_viewer(tag, outdir)
    for (name in names(viewer_paths)) paths[[name]] <- viewer_paths[[name]]
  }
  cat(sprintf("[report] %s: wrote reference-free result package -> %s\n", tag, outdir))
  invisible(paths)
}
.rep_invoked_directly <- local({
  match <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  length(match) == 1L &&
    identical(basename(sub("^--file=", "", match[1])), "reporting.R")
})
if (.rep_invoked_directly) {
  arguments <- commandArgs(trailingOnly = TRUE)
  if (!length(arguments)) {
    stop("usage: Rscript reporting.R <tag> [outdir]", call. = FALSE)
  }
  scma_generate_report(
    arguments[1],
    if (length(arguments) >= 2L) arguments[2] else file.path(RESULTS, arguments[1])
  )
}
