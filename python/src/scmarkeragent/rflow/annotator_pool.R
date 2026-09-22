#!/usr/bin/env Rscript
A_POOL <- CFG$cluster_annotation
SOURCES_PER_GENE <- as.integer(A_POOL$sources_per_marker)
SOURCE_GENES_PER_QUERY <- as.integer(A_POOL$source_genes_per_query)
SOURCE_MAX_CHARS <- as.integer(A_POOL$source_sentence_max_chars)
UNMEASURED_SHOWN <- as.integer(A_POOL$unmeasured_genes_shown)
WARNING_GENES_SHOWN <- as.integer(A_POOL$warning_genes_shown)
CANDIDATE_TIER_SIZE <- as.integer(A_POOL$candidate_tier_size)
CANDIDATE_TIERS <- as.integer(A_POOL$candidate_tiers)
DECISIVE_MIN_PUB <- as.integer(A_POOL$decisive_min_publications)
SINGLE_STUDY_PUB <- as.integer(A_POOL$single_study_publications)
SPARSE_LABEL_MAX_PAPERS <- as.integer(A_POOL$sparse_label_max_papers)
SPARSE_LABEL_PAPERS_MAX <- SPARSE_LABEL_MAX_PAPERS
EXCLUSION_MIN_PUB <- as.integer(A_POOL$exclusion_min_publications)
DM_GATE_TOP <- as.integer(A_POOL$defining_marker_gate_top)
DM_GATE_FLOOR_PCT <- as.numeric(A_POOL$defining_marker_gate_floor_pct)
DM_GATE_DOMINANT_PCT <- as.numeric(A_POOL$defining_marker_gate_dominant_pct)
GENOME_ENRICHED_ROWS <- as.integer(A_POOL$evidence_judge_enriched_rows)
GENOME_DEPLETED_ROWS <- as.integer(A_POOL$evidence_judge_depleted_rows)
REVIEW_OTHER_CANDIDATE_ROWS <- as.integer(A_POOL$review_other_candidate_rows)
REVIEW_SOURCES_PER_MARKER <- as.integer(A_POOL$review_sources_per_marker)
NEGATIVE_SOURCE_MIN_PCT_IN <- as.numeric(A_POOL$negative_source_min_pct_in)
CLUSTER_MARKER_ROWS <- as.integer(A_POOL$cluster_marker_rows)
DEFINERS_PER_CANDIDATE <- as.integer(A_POOL$definers_per_candidate %||% 8L)
DEFINER_TOP_CLUSTERS <- as.integer(A_POOL$definer_top_clusters %||% 3L)
GENOME_DETECTION_GAP_ROWS <- as.integer(A_POOL$genome_detection_gap_rows %||% 20L)
GENOME_CLAIMANTS_SHOWN <- as.integer(A_POOL$genome_claimants_shown %||% 3L)
PARTITION_OWN_SHOWN <- as.integer(A_POOL$partition_own_shown %||% 15L)
PARTITION_SHARED_SHOWN <- as.integer(A_POOL$partition_shared_shown %||% 8L)
TOOLS <- c("sources", "gene_across_clusters", "candidates_with_gene")
REVIEW_TOOLS <- c("sources", "gene_across_clusters", "candidates_with_gene")
PCT_DECIMALS <- 1L
LOGFC_DECIMALS <- 3L
SPECIFICITY_DECIMALS <- 3L
AUC_DECIMALS <- 3L
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
.round <- function(value, digits = 4L) {
  if (is.null(value) || !length(value)) {
    return(NULL)
  }
  number <- suppressWarnings(as.numeric(value[1]))
  if (!is.finite(number)) NULL else round(number, digits)
}
.clip <- function(text, limit) {
  value <- gsub("\\s+", " ", trimws(as.character(text %||% "")))
  if (nchar(value) <= limit) value else paste0(substr(value, 1, limit - 1), "\u2026")
}
.identifier <- function(value) {
  text <- trimws(as.character(value %||% "")[1])
  if (!nzchar(text) || tolower(text) %in% c("nan", "none", "na")) NOT_AVAILABLE else text
}
.num0 <- function(value) {
  number <- suppressWarnings(as.numeric(value %||% NA_real_)[1])
  if (!length(number) || is.na(number)) 0 else number
}
.by_npub <- function(rows) {
  if (!length(rows)) {
    return(rows)
  }
  rows[order(
    -vapply(rows, function(m) as.integer(m$n_pub %||% 0L), 0L),
    .cp(vapply(rows, function(m) as.character(m$gene), "")),
    method = "radix"
  )]
}
.positive_rows <- function(entry) {
  Filter(function(m) identical(m$polarity, "positive"), entry$markers %||% list())
}
.fmt_pct <- function(value) {
  if (is.null(value) || !length(value) || is.na(value[1])) {
    "None"
  } else {
    sprintf("%.1f", as.numeric(value[1]))
  }
}
.fmt_num <- function(value) sprintf("%.1f", .num0(value))
.fmt_int <- function(value) {
  if (is.null(value) || !length(value) || is.na(value[1])) {
    "None"
  } else {
    sprintf("%d", as.integer(value[1]))
  }
}
.py_list <- function(values) {
  items <- as.character(unlist(values %||% list()))
  if (!length(items)) {
    return("[]")
  }
  paste0("[", paste(sprintf("'%s'", items), collapse = ", "), "]")
}
.papers_table <- function(keys, values) {
  keys <- as.character(keys)
  values <- as.integer(values)
  keep <- !duplicated(keys, fromLast = TRUE)
  setNames(values[keep], keys[keep])
}
.LABEL_PAPERS <- new.env(parent = emptyenv())
scma_label_papers <- function(species, label) {
  if (!exists("table", envir = .LABEL_PAPERS, inherits = FALSE)) {
    path <- .env(
      "SCMA_LABEL_PAPERS",
      file.path(PACKAGE_ROOT, "resources", "label_papers_v4.csv")
    )
    table <- if (file.exists(path)) {
      rows <- fread(path, colClasses = "character", showProgress = FALSE)
      .papers_table(
        paste(as.character(rows$species), as.character(rows$cell_type), sep = "\x1f"),
        rows$n_papers
      )
    } else {
      integer(0)
    }
    assign("table", table, envir = .LABEL_PAPERS)
  }
  table <- get("table", envir = .LABEL_PAPERS, inherits = FALSE)
  value <- unname(table[paste(as.character(species), as.character(label), sep = "\x1f")])
  if (!length(value) || is.na(value)) NULL else as.integer(value)
}
.sentence_record <- function(record) {
  out <- list(
    pmid = .identifier(record$pmid),
    pmcid = .identifier(record$pmcid),
    sentence = .clip(record$sentence, SOURCE_MAX_CHARS)
  )
  n_pub <- record$n_pub
  if (!is.null(n_pub) && length(n_pub) && !is.na(n_pub[1])) {
    out$n_pub <- as.integer(n_pub[1])
    out$single_study <- out$n_pub <= SINGLE_STUDY_PUB
  }
  out
}
.opening_records <- function(sources, label, gene, k) {
  if (is.null(sources)) {
    return(list())
  }
  if (.is_source_server(sources)) {
    sources$opening(label, gene, k)
  } else {
    sources(label, gene, k)
  }
}
.is_raised <- function(marker) {
  value <- marker$avg_log2FC
  if (is.null(value) || !length(value) || !is.finite(suppressWarnings(as.numeric(value[1]))) ||
    as.numeric(value[1]) <= 0) {
    return(FALSE)
  }
  auc <- marker$auc
  !is.null(auc) && length(auc) && is.finite(suppressWarnings(as.numeric(auc[1]))) &&
    as.numeric(auc[1]) > SIG_AUC
}
.sorted_cluster_ids <- function(scoring) {
  ids <- names(scoring$clusters)
  ids[order(nchar(ids), ids)]
}
.panel_rows <- function(candidate_value, positive, negative, cluster, measured, curated,
                        native_of, specificity = NULL) {
  genes <- toupper(c(as.character(positive), as.character(negative)))
  polarity <- c(rep("positive", length(positive)), rep("negative", length(negative)))
  rows <- list()
  seen <- character(0)
  for (position in seq_along(genes)) {
    gene_value <- genes[position]
    if (gene_value %in% seen) next
    key <- paste(cluster, gene_value, sep = "\x1f")
    if (is.na(measured$pct_in[key]) || is.null(measured$pct_in[[key]])) next
    meta <- curated[.(candidate_value, gene_value), nomatch = NULL]
    curated_polarity <- if (nrow(meta)) {
      as.character(meta$marker_polarity[nrow(meta)])
    } else {
      polarity[position]
    }
    curated_npub <- as.integer(if (nrow(meta)) meta$n_pub[nrow(meta)] else 0L)
    if (identical(curated_polarity, "negative") && curated_npub < EXCLUSION_MIN_PUB) next
    seen <- c(seen, gene_value)
    native <- unname(native_of[gene_value])
    pct_in <- round(as.numeric(measured$pct_in[[key]]), PCT_DECIMALS)
    lfc <- round(as.numeric(measured$avg_log2FC[[key]]), LOGFC_DECIMALS)
    padj <- round(suppressWarnings(as.numeric(measured$padj[[key]] %||% NA_real_)), 6)
    auc_value <- round(suppressWarnings(as.numeric(measured$auc[[key]] %||% NA_real_)), AUC_DECIMALS)
    m_g <- suppressWarnings(as.numeric(
      if (gene_value %in% names(specificity)) specificity[[gene_value]] else NA_real_
    ))
    rows[[length(rows) + 1L]] <- list(
      gene = if (!length(native) || is.na(native)) gene_value else native,
      polarity = curated_polarity,
      n_pub = curated_npub,
      tier = if (nrow(meta)) .identifier(meta$tier[nrow(meta)]) else NOT_AVAILABLE,
      pct_in = pct_in,
      pct_out = round(as.numeric(measured$pct_out[[key]]), PCT_DECIMALS),
      avg_log2FC = lfc,
      auc = if (is.na(auc_value)) NULL else auc_value,
      significant = isTRUE(!is.na(padj) && sig_pass(lfc, pct_in / 100, padj, auc_value)),
      specificity = if (length(m_g) && !is.na(m_g)) round(m_g, SPECIFICITY_DECIMALS) else NULL,
      recommended = FALSE
    )
  }
  .by_npub(rows)
}
.exclusion_sources <- function(candidate_value, markers, sources) {
  found <- list()
  for (marker in markers) {
    if (!identical(marker$polarity, "negative")) next
    pct_in <- marker$pct_in
    if (is.null(pct_in) || is.na(pct_in) || pct_in < NEGATIVE_SOURCE_MIN_PCT_IN) next
    gene <- toupper(as.character(marker$gene))
    found[[gene]] <- lapply(
      .opening_records(sources, candidate_value, gene, SOURCES_PER_GENE),
      .sentence_record
    )
  }
  found
}
.build_pool <- function(scoring, de, native_of, sources = NULL, markers_all = NULL,
                        prep_genes = NULL) {
  keys <- paste(de$group, de$gene_key, sep = "\x1f")
  measured <- list(
    pct_in = setNames(de$pct_in, keys),
    pct_out = setNames(de$pct_out, keys),
    avg_log2FC = setNames(de$avg_log2FC, keys),
    auc = setNames(de$auc, keys),
    padj = setNames(de$padj, keys)
  )
  gene_clusters <- list()
  for (gene in unique(de$gene_key)) {
    rows <- de[gene_key == gene]
    gene_clusters[[gene]] <- lapply(seq_len(nrow(rows)), function(i) {
      list(
        cluster_id = as.character(rows$group[i]),
        pct_in = round(as.numeric(rows$pct_in[i]), PCT_DECIMALS),
        pct_out = round(as.numeric(rows$pct_out[i]), PCT_DECIMALS)
      )
    })
  }
  genome_stats <- list(pct_in = c(), pct_out = c(), avg_log2FC = c(), auc = c(), padj = c())
  if (!is.null(markers_all) && nrow(markers_all)) {
    ga <- as.data.table(markers_all)
    ga[, group := as.character(group)]
    ga[, gene_key := toupper(as.character(feature))]
    ga_keys <- paste(ga$group, ga$gene_key, sep = "\x1f")
    is_new <- !(ga_keys %in% keys)
    if (any(is_new)) {
      sub <- ga[is_new]
      sub_keys <- ga_keys[is_new]
      genome_stats <- list(
        pct_in = setNames(sub$pct_in, sub_keys), pct_out = setNames(sub$pct_out, sub_keys),
        avg_log2FC = setNames(sub$avg_log2FC, sub_keys), auc = setNames(sub$auc, sub_keys),
        padj = setNames(sub$padj, sub_keys)
      )
      for (gene in unique(sub$gene_key)) {
        rows <- sub[gene_key == gene]
        gene_clusters[[gene]] <- c(gene_clusters[[gene]], lapply(seq_len(nrow(rows)), function(i) {
          list(
            cluster_id = as.character(rows$group[i]),
            pct_in = round(as.numeric(rows$pct_in[i]), PCT_DECIMALS),
            pct_out = round(as.numeric(rows$pct_out[i]), PCT_DECIMALS)
          )
        }))
      }
      for (gene in unique(sub$gene_key)) {
        if (is.na(unname(native_of[gene]))) {
          native_of[gene] <- as.character(sub[gene_key == gene]$feature[1])
        }
      }
    }
  }
  measured_genes <- if (!is.null(prep_genes)) unique(toupper(as.character(prep_genes))) else character(0)
  curated <- as.data.table(scoring$panel_records)
  curated[, gene_key := toupper(as.character(gene_key))]
  setkey(curated, candidate, gene_key)
  gene_carriers <- list()
  curated_positive <- list()
  for (i in seq_len(nrow(curated))) {
    gene <- curated$gene_key[i]
    row <- list(
      cell_type = as.character(curated$candidate[i]),
      n_pub = as.integer(curated$n_pub[i] %||% 0L),
      tier = .identifier(curated$tier[i]),
      polarity = as.character(curated$marker_polarity[i])
    )
    gene_carriers[[gene]] <- c(gene_carriers[[gene]] %||% list(), list(row))
    if (identical(row$polarity, "positive")) {
      curated_positive[[row$cell_type]] <- unique(c(
        curated_positive[[row$cell_type]] %||% character(0), gene
      ))
    }
  }
  tissue_contexts <- as.character(scoring$tissue_contexts %||% character(0))
  candidate_contexts <- scoring$candidate_tissue_contexts %||% list()
  retrieval_context_of <- scoring$candidate_retrieval_context %||% character(0)
  .retrieval_context <- function(name) {
    value <- retrieval_context_of[name]
    if (!length(value) || is.na(value)) "" else unname(as.character(value))
  }
  scored <- as.data.table(scoring$scored)
  clusters <- list()
  for (cluster in .sorted_cluster_ids(scoring)) {
    state <- scoring$clusters[[cluster]]
    cluster_value <- cluster
    frame <- scored[cluster == cluster_value]
    entries <- list()
    for (name in as.character(state$candidates %||% character(0))) {
      candidate_value <- name
      row <- frame[candidate == candidate_value]
      if (!nrow(row)) next
      markers <- .panel_rows(
        name,
        scoring$measured_panels[[name]] %||% character(0),
        scoring$negative_panels[[name]] %||% character(0),
        cluster, measured, curated, native_of, scoring$marker_specificity
      )
      shown_genes <- toupper(vapply(markers, function(m) m$gene, ""))
      unmeasured <- sort(setdiff(curated_positive[[name]] %||% character(0), shown_genes))
      shown <- unmeasured[seq_len(min(UNMEASURED_SHOWN, length(unmeasured)))]
      entries[[length(entries) + 1L]] <- list(
        cell_type = name,
        retrieval_rank = as.integer(row$retrieval_rank[1]),
        tissue_context = as.list(as.character(candidate_contexts[[name]] %||% character(0))),
        retrieval_context = .retrieval_context(name),
        markers = markers,
        exclusion_sources = .exclusion_sources(name, markers, sources),
        unmeasured_curated_genes = list(
          count = length(unmeasured),
          genes = as.list(unname(vapply(shown, function(g) {
            native <- unname(native_of[g])
            if (!length(native) || is.na(native)) g else native
          }, "")))
        ),
        program = list(
          median_in = .round(row$program_in_median[1]),
          median_out = .round(row$program_out_median[1])
        )
      )
    }
    clusters[[cluster]] <- list(
      cluster_id = cluster,
      status = state$status,
      n_cells = as.integer(state$n_cells),
      candidates = entries
    )
  }
  context <- scoring$context
  list(
    context = list(
      species = context$species,
      tissue = context$tissue,
      disease = as.character(context$disease),
      development_stage = if (nzchar(context$development_stage %||% "")) {
        context$development_stage
      } else {
        NOT_AVAILABLE
      },
      clusters_in_dataset = as.integer(context$n_clusters)
    ),
    clusters = clusters,
    tissue_contexts = as.list(tissue_contexts),
    candidate_tissue_contexts = lapply(candidate_contexts, function(v) as.list(as.character(v))),
    candidate_retrieval_context = as.list(retrieval_context_of),
    retrieval_semantics = as.character(scoring$retrieval_semantics %||% ""),
    gene_clusters = gene_clusters,
    gene_carriers = gene_carriers,
    native_gene = as.list(native_of),
    native_menu = as.list(sort(.cp(names(curated_positive)), method = "radix")),
    native_menu_genes = as.list(sort(.cp(unique(unlist(curated_positive, use.names = FALSE))),
      method = "radix"
    )),
    de_stats = measured,
    genome_stats = genome_stats,
    measured_genes = measured_genes,
    marker_specificity = scoring$marker_specificity %||% list()
  )
}
scma_assign_tiers <- function(pool, scoring) {
  for (cluster in names(pool$clusters)) {
    scored_state <- (scoring$clusters %||% list())[[cluster]] %||% list()
    extras <- character(0)
    for (block in scored_state$extra_contexts %||% list()) {
      extras <- c(extras, as.character(unlist(block$candidates %||% list())))
    }
    extras <- unique(extras)
    entries <- pool$clusters[[cluster]]$candidates
    seen_primary <- 0L
    for (index in seq_along(entries)) {
      entry <- entries[[index]]
      if (as.character(entry$cell_type) %in% extras || !is.null(entry$borrowed_context)) {
        entries[[index]]$tier <- 1L
        next
      }
      entries[[index]]$tier <- as.integer(seen_primary %/% CANDIDATE_TIER_SIZE + 1L)
      seen_primary <- seen_primary + 1L
    }
    pool$clusters[[cluster]]$candidates <- entries
  }
  pool
}
.tier_of <- function(entry) as.integer(entry$tier %||% 1L)
.tiers_available <- function(pool, cluster) {
  entries <- pool$clusters[[as.character(cluster)]]$candidates
  highest <- if (length(entries)) max(vapply(entries, .tier_of, 1L)) else 1L
  if (highest > CANDIDATE_TIERS) {
    stop(sprintf(
      "cluster %s: candidates stamped tier %d above the %d-tier ceiling; scma_assign_tiers and retrieval.top_candidates disagree",
      as.character(cluster), as.integer(highest), as.integer(CANDIDATE_TIERS)
    ), call. = FALSE)
  }
  min(CANDIDATE_TIERS, highest)
}
.candidate_names <- function(pool, cluster, tier = NULL) {
  entries <- pool$clusters[[as.character(cluster)]]$candidates
  if (!is.null(tier)) {
    entries <- Filter(function(e) .tier_of(e) <= as.integer(tier), entries)
  }
  vapply(entries, function(e) as.character(e$cell_type), "")
}
.find_candidate <- function(pool, cluster, name) {
  for (entry in pool$clusters[[as.character(cluster)]]$candidates) {
    if (identical(entry$cell_type, as.character(name))) {
      return(entry)
    }
  }
  NULL
}
PACKET_MARKER_FIELDS <- c(
  "gene", "polarity", "n_pub", "recommended", "specificity",
  "pct_in", "pct_out", "avg_log2FC", "auc", "significant",
  "shared_on_page"
)
.packet_marker <- function(row) {
  row[intersect(PACKET_MARKER_FIELDS, names(row))]
}
MARKER_ROW_FORMAT <- paste0(
  "GENE n<publications>[R] <pct_in>/<pct_out>[!] fc<avg_log2FC> auc<auc> ",
  "sp<specificity>[ sig][ neg][ ~Name1; Name2]: R the resource's recommended flag; ",
  "! not raised here (avg_log2FC at or below 0, or auc at or below 0.5); sig the row ",
  "clears the pipeline's significance gate (avg_log2FC above 0.25, pct_in above 10, ",
  "auc above 0.5, adjusted p below 0.05); neg the resource curates the gene as a ",
  "NEGATIVE marker of this identity, one its cells do not express; ~ the other ",
  "candidates of this page the species' literature also ties the gene to, full names ",
  "separated by '; '. A '-' where a number belongs means not measured."
)
.SHARED_SEP <- "; "
.compact_fmt_num <- function(value) {
  if (is.null(value) || !length(value) || is.na(value[1])) {
    return("-")
  }
  if (is.logical(value)) {
    return(if (isTRUE(value[1])) "1" else "0")
  }
  number <- suppressWarnings(as.numeric(value[1]))
  if (!is.finite(number)) {
    return("-")
  }
  text <- as.character(number)
  if (!grepl("[.eE]", text)) text <- paste0(text, ".0")
  text
}
compact_marker_row <- function(row) {
  n_pub <- row$n_pub
  tag <- sprintf("%s n%s", row$gene, if (is.null(n_pub) || is.na(n_pub)) "-" else as.character(as.integer(n_pub)))
  if (isTRUE(row$recommended)) tag <- paste0(tag, "R")
  tag <- paste0(tag, sprintf(" %s/%s", .compact_fmt_num(row$pct_in), .compact_fmt_num(row$pct_out)))
  if (!.is_raised(row)) tag <- paste0(tag, "!")
  tag <- paste0(tag, sprintf(" fc%s auc%s", .compact_fmt_num(row$avg_log2FC), .compact_fmt_num(row$auc)))
  tag <- paste0(tag, sprintf(" sp%s", .compact_fmt_num(row$specificity)))
  if (isTRUE(row$significant)) tag <- paste0(tag, " sig")
  if (identical(as.character(row$polarity %||% "positive"), "negative")) tag <- paste0(tag, " neg")
  shared <- Filter(nzchar, as.character(unlist(row$shared_on_page %||% list())))
  if (length(shared)) tag <- paste0(tag, " ~", paste(shared, collapse = .SHARED_SEP))
  tag
}
.COMPACT_ROW_RE <- paste0(
  "^(\\S+) n(\\d+|-)(R?) (\\S+)/(\\S+)(!?) fc(\\S+) auc(\\S+) sp(\\S+)",
  "((?: sig| neg)*)(?: ~(.*))?$"
)
compact_panel <- function(rows) {
  vapply(rows %||% list(), compact_marker_row, "")
}
compact_row_genes <- function(rows) {
  out <- character(0)
  for (row in rows %||% list()) {
    gene <- if (is.list(row)) {
      row$gene
    } else {
      m <- regmatches(as.character(row), regexpr("^\\S+", as.character(row)))
      if (length(m) && nzchar(m)) m else NULL
    }
    if (!is.null(gene) && nzchar(gene)) out <- c(out, toupper(as.character(gene)))
  }
  unique(out)
}
compact_definer_rows <- function(rows) {
  vapply(rows %||% list(), function(r) {
    tag <- sprintf("%s n%s%s", r$gene, as.integer(r$n_pub), if (isTRUE(r$recommended)) "R" else "")
    tag <- paste0(tag, if (isTRUE(r$in_this_tissue_panel)) "" else "*")
    if (!is.null(r$pct_in)) {
      tag <- paste0(tag, sprintf(" %s/%s", r$pct_in, r$pct_out))
      if (!isTRUE(r$raised)) tag <- paste0(tag, "!")
      peaks <- as.character(unlist(r$peak_clusters %||% list()))
      if (length(peaks)) tag <- paste0(tag, sprintf(" peak=%s", peaks[1]))
    } else {
      tag <- paste0(tag, sprintf(" %s", if (nzchar(as.character(r$measurement %||% ""))) r$measurement else "unmeasured"))
    }
    shared <- Filter(nzchar, as.character(unlist(r$shared_with_page %||% list())))
    if (length(shared)) {
      shown <- shared[seq_len(min(3L, length(shared)))]
      tag <- paste0(
        tag, " ~", paste(shown, collapse = ","),
        if (length(shared) > 3L) sprintf(",+%d", length(shared) - 3L) else ""
      )
    }
    tag
  }, "")
}
.cluster_marker_table <- function(pool, cluster, tier = NULL) {
  rows <- list()
  for (entry in pool$clusters[[as.character(cluster)]]$candidates) {
    if (!is.null(tier) && .tier_of(entry) > as.integer(tier)) next
    for (marker in entry$markers) {
      if (!identical(marker$polarity, "positive") || !isTRUE(marker$significant)) next
      gene <- as.character(marker$gene)
      if (is.null(rows[[gene]])) {
        rows[[gene]] <- list(
          gene = gene,
          pct_in = marker$pct_in,
          pct_out = marker$pct_out,
          specificity = marker$specificity,
          claimed_by = list()
        )
      }
      rows[[gene]]$claimed_by <- c(
        rows[[gene]]$claimed_by,
        list(list(cell_type = entry$cell_type, n_pub = as.integer(marker$n_pub)))
      )
    }
  }
  if (!length(rows)) {
    return(list())
  }
  rows <- lapply(rows, function(row) {
    claims <- row$claimed_by
    order_by <- order(
      -vapply(claims, function(c) c$n_pub, 0L),
      vapply(claims, function(c) c$cell_type, "")
    )
    row$claimed_by <- claims[order_by]
    row
  })
  contrast <- vapply(rows, function(r) {
    as.numeric(r$pct_in %||% 0) - as.numeric(r$pct_out %||% 0)
  }, 0)
  gene <- vapply(rows, function(r) r$gene, "")
  ordered <- unname(rows[order(-contrast, gene)])
  ordered[seq_len(min(CLUSTER_MARKER_ROWS, length(ordered)))]
}
.candidate_block <- function(pool, entry) {
  species <- as.character(pool$context$species %||% "")
  block <- list(
    cell_type = entry$cell_type,
    markers = as.list(compact_panel(entry$markers)),
    exclusion_sources = entry$exclusion_sources,
    unmeasured_curated_genes = entry$unmeasured_curated_genes,
    program = entry$program
  )
  if (!is.null(entry$borrowed_context)) {
    block$borrowed_context <- entry$borrowed_context
    block$flat_or_undetected_curated_genes <- entry$flat_or_undetected_curated_genes
  }
  papers <- scma_label_papers(species, as.character(entry$cell_type))
  if (!is.null(papers)) {
    block$n_papers <- as.integer(papers)
  }
  if (!is.null(entry$definers)) {
    block$species_definers <- as.list(compact_definer_rows(entry$definers))
  }
  block
}
.packet_facts <- function(facts, entry) {
  positives <- .positive_rows(entry)
  by_pub <- .by_npub(positives)
  facts$recommended_markers <- as.list(vapply(
    Filter(function(m) isTRUE(m$recommended), positives), function(m) as.character(m$gene), ""
  ))
  facts$best_published_markers <- as.list(vapply(
    by_pub[seq_len(min(2L, length(by_pub)))], function(m) as.character(m$gene), ""
  ))
  facts$detected_exclusions <- as.list(compact_panel(.detected_exclusion_markers(entry)))
  facts
}
.cluster_packet <- function(pool, cluster, tier = 1L) {
  state <- pool$clusters[[as.character(cluster)]]
  tier <- as.integer(tier)
  entries <- Filter(function(e) .tier_of(e) <= tier, state$candidates)
  ordered <- entries[order(
    .cp(vapply(entries, function(e) as.character(e$cell_type), "")), method = "radix"
  )]
  c(
    list(
      query = c(
        pool$context,
        list(cluster_id = as.character(cluster), cells_in_cluster = as.integer(state$n_cells))
      )
    ),
    genome_pages(pool, cluster, tier),
    list(
      cluster_markers = .cluster_marker_table(pool, cluster, tier = tier),
      row_format = MARKER_ROW_FORMAT,
      candidate_facts = lapply(entries, function(entry) {
        .packet_facts(.candidate_facts(pool, cluster, entry), entry)
      }),
      candidates = lapply(ordered, function(entry) .candidate_block(pool, entry))
    )
  )
}
.tier_packet <- function(pool, cluster, tier) {
  tier <- as.integer(tier)
  entries <- Filter(
    function(e) .tier_of(e) == tier, pool$clusters[[as.character(cluster)]]$candidates
  )
  ordered <- entries[order(
    .cp(vapply(entries, function(e) as.character(e$cell_type), "")), method = "radix"
  )]
  list(
    additional_candidates_tier = tier,
    cluster_markers = .cluster_marker_table(pool, cluster, tier = tier),
    row_format = MARKER_ROW_FORMAT,
    candidate_facts = lapply(entries, function(entry) {
      .packet_facts(.candidate_facts(pool, cluster, entry), entry)
    }),
    candidates = lapply(ordered, function(entry) .candidate_block(pool, entry))
  )
}
.review_panel <- function(entry) {
  markers <- entry$markers %||% list()
  positives <- Filter(function(r) identical(r$polarity, "positive"), markers)
  negatives <- Filter(function(r) identical(r$polarity, "negative"), markers)
  lapply(c(positives, negatives), .packet_marker)
}
.review_sources <- function(label, panel, sources) {
  found <- list()
  if (is.null(sources)) {
    return(found)
  }
  for (row in panel) {
    gene <- toupper(as.character(row$gene))
    if (!is.null(found[[gene]])) next
    records <- .opening_records(sources, label, gene, REVIEW_SOURCES_PER_MARKER)
    if (length(records)) found[[gene]] <- lapply(records, .sentence_record)
  }
  found
}
.other_candidate_blocks <- function(pool, cluster, claimed, tier) {
  blocks <- list()
  species <- as.character(pool$context$species %||% "")
  for (entry in pool$clusters[[as.character(cluster)]]$candidates) {
    name <- as.character(entry$cell_type)
    if (name %in% claimed || .tier_of(entry) > as.integer(tier)) next
    rows <- .review_panel(entry)
    positives <- Filter(function(r) identical(r$polarity, "positive"), rows)
    positives <- positives[seq_len(min(REVIEW_OTHER_CANDIDATE_ROWS, length(positives)))]
    negatives <- Filter(function(r) identical(r$polarity, "negative"), rows)
    block <- list(cell_type = name, panel = as.list(compact_panel(c(positives, negatives))))
    if (!is.null(entry$borrowed_context)) block$borrowed_context <- entry$borrowed_context
    papers <- scma_label_papers(species, name)
    if (!is.null(papers)) block$n_papers <- as.integer(papers)
    blocks[[length(blocks) + 1L]] <- block
  }
  blocks
}
.review_packet <- function(pool, cluster, final, sources = NULL, tier = 1L,
                           unknown_token = "Unknown") {
  state <- pool$clusters[[as.character(cluster)]]
  selected <- as.character(final$selected %||% "")
  subtype <- as.character(final$subtype %||% "")
  others <- as.character(unlist(final$co_occurring_identities %||% list()))
  claims <- list()
  claim_names <- character(0)
  if (nzchar(selected) && !identical(selected, unknown_token)) {
    claims[[length(claims) + 1L]] <- list(role = "selected", name = selected)
    claim_names <- c(claim_names, selected)
  }
  if (nzchar(subtype) && !identical(subtype, selected)) {
    claims[[length(claims) + 1L]] <- list(role = "finer_name", name = subtype)
    claim_names <- c(claim_names, subtype)
  }
  for (name in others) {
    if (nzchar(name) && !(name %in% claim_names)) {
      claims[[length(claims) + 1L]] <- list(role = "co_occurring", name = name)
      claim_names <- c(claim_names, name)
    }
  }
  .claimed_block <- function(name, role) {
    entry <- .find_candidate(pool, cluster, name)
    if (is.null(entry)) return(NULL)
    panel <- .review_panel(entry)
    block <- list(
      cell_type = name, role = role, panel = as.list(compact_panel(panel)),
      sources = .review_sources(name, panel, sources),
      detected_exclusions = as.list(compact_panel(.detected_exclusion_markers(entry))),
      unmeasured_curated_genes = entry$unmeasured_curated_genes,
      facts = .packet_facts(.candidate_facts(pool, cluster, entry), entry)
    )
    if (!is.null(entry$definers)) block$species_definers <- entry$definers
    if (!is.null(entry$borrowed_context)) block$borrowed_context <- entry$borrowed_context
    audit <- (final$definer_audit %||% list())[[name]]
    if (!is.null(audit)) block$definer_audit <- audit
    block
  }
  claimed <- Filter(Negate(is.null), lapply(claims, function(c) .claimed_block(c$name, c$role)))
  contested_names <- Filter(function(n) nzchar(n) && !(n %in% claim_names),
                            as.character(unlist(final$contested_rivals %||% list())))
  contested_blocks <- Filter(Negate(is.null), lapply(contested_names, function(n) {
    .claimed_block(n, "contested_rival")
  }))
  covered <- unique(c(claim_names, vapply(contested_blocks, function(b) b$cell_type, "")))
  packet <- c(
    list(
      query = c(
        pool$context,
        list(cluster_id = as.character(cluster), cells_in_cluster = as.integer(state$n_cells))
      ),
      delivered = list(
        selected = selected, subtype = subtype, co_occurring_identities = as.list(others),
        state = as.character(final$state %||% ""), reason = as.character(final$reason %||% "")
      ),
      row_format = MARKER_ROW_FORMAT,
      claimed_identities = claimed
    ),
    genome_pages(pool, cluster, tier),
    list(candidates_not_claimed = .other_candidate_blocks(pool, cluster, covered, tier))
  )
  if (length(contested_blocks)) {
    packet$contested_identities <- contested_blocks
    packet$shared_with_contested <- .shared_with(
      pool, cluster, selected, vapply(contested_blocks, function(b) b$cell_type, "")
    )
    sel_entry <- .find_candidate(pool, cluster, selected)
    if (!is.null(sel_entry)) {
      packet$pair_partitions <- lapply(contested_blocks, function(b) {
        rentry <- .find_candidate(pool, cluster, b$cell_type)
        if (is.null(rentry)) NULL else pair_partition(sel_entry, rentry)
      })
      packet$pair_partitions <- Filter(Negate(is.null), packet$pair_partitions)
    }
  }
  if (length(final$definer_audit %||% list())) {
    packet$definer_audit <- final$definer_audit
  }
  screened <- .screened_out_claimants(pool, cluster)
  if (length(screened)) {
    packet$screened_out_claimants <- screened
  }
  packet
}
SCREENED_OUT_SHOWN <- 5L
SCREENED_OUT_ANCHORS_SHOWN <- 12L
.screened_out_claimants <- function(pool, cluster) {
  log <- (pool$borrowed %||% list())[[as.character(cluster)]] %||% list()
  native <- pool$native_gene %||% list()
  prefix <- "tissue screen: "
  .native <- function(genes) {
    genes <- as.character(unlist(genes))
    as.list(unname(vapply(genes, function(g) as.character(native[[g]] %||% g), "")))
  }
  out <- list()
  for (item in log$considered %||% list()) {
    rejected <- as.character(item$rejected %||% "")
    if (!startsWith(rejected, prefix)) next
    anchors <- as.character(unlist(item$anchor_genes %||% list()))
    out[[length(out) + 1L]] <- list(
      cell_type = as.character(item$cell_type %||% ""),
      borrow_rank = as.integer(item$borrow_rank %||% 0L),
      anchors_covered = as.integer(item$anchors_covered %||% length(anchors)),
      anchor_genes = .native(anchors[seq_len(min(SCREENED_OUT_ANCHORS_SHOWN, length(anchors)))]),
      recommended_anchors = .native(item$recommended_anchors %||% list()),
      screen_reason = substr(rejected, nchar(prefix) + 1L, nchar(rejected))
    )
  }
  if (!length(out)) return(list())
  ranks <- vapply(out, function(r) if (r$borrow_rank > 0L) r$borrow_rank else 1000000L, 0L)
  covered <- vapply(out, function(r) r$anchors_covered, 0L)
  names_ <- vapply(out, function(r) r$cell_type, "")
  ord <- order(ranks, -covered, .cp(names_), method = "radix")
  out <- out[ord]
  out[seq_len(min(SCREENED_OUT_SHOWN, length(out)))]
}
.shared_with <- function(pool, cluster, identity, others) {
  entry <- .find_candidate(pool, cluster, identity)
  if (is.null(entry) || !length(others)) return(list())
  genes <- unique(toupper(vapply(
    Filter(function(m) identical(m$polarity, "positive"), entry$markers), function(m) as.character(m$gene), ""
  )))
  genes <- unique(c(genes, toupper(vapply(entry$definers %||% list(), function(d) as.character(d$gene), ""))))
  genes <- sort(.cp(genes), method = "radix")
  claimants <- pool$species_claimants %||% list()
  out <- list()
  for (gene in genes) {
    also <- Filter(function(n) n %in% others, vapply(claimants[[gene]] %||% list(), function(c) c[[1]], ""))
    if (length(also)) out[[gene]] <- as.list(also)
  }
  out
}
.carried <- function(marker) {
  pct_in <- .num0(marker$pct_in)
  raised <- .is_raised(marker) && isTRUE(marker$significant) && pct_in >= DM_GATE_FLOOR_PCT
  raised || pct_in >= DM_GATE_DOMINANT_PCT
}
scma_definers_absent <- function(rows, k) {
  top <- .by_npub(rows)
  top <- top[seq_len(min(as.integer(k), length(top)))]
  length(top) == as.integer(k) &&
    all(vapply(top, function(m) .num0(m$pct_in) < DM_GATE_FLOOR_PCT, TRUE))
}
scma_genome_de_summary <- function(markers_all) {
  out <- list()
  if (is.null(markers_all) || !nrow(markers_all)) {
    return(out)
  }
  frame <- as.data.table(markers_all)
  frame[, group := as.character(group)]
  .rows <- function(sub, fields) {
    lapply(seq_len(nrow(sub)), function(i) {
      row <- list(
        gene = as.character(sub$feature[i]),
        log2FC = .round(sub$avg_log2FC[i], 2L),
        pct_in = .round(sub$pct_in[i], 1L),
        pct_out = .round(sub$pct_out[i], 1L)
      )
      if ("auc" %in% fields) row$auc <- .round(sub$auc[i], 3L)
      row
    })
  }
  for (cluster in unique(frame$group)) {
    cluster_value <- cluster
    sub <- frame[group == cluster_value]
    up <- sub[padj < 0.05 & avg_log2FC > 0 & pct_in >= 10]
    up <- head(up[order(-avg_log2FC)], GENOME_ENRICHED_ROWS)
    down <- sub[avg_log2FC < 0 & pct_out >= 20]
    down <- head(down[order(avg_log2FC)], GENOME_DEPLETED_ROWS)
    anchors <- sub[padj < 0.05 & avg_log2FC > 0 & pct_in >= DM_GATE_DOMINANT_PCT]
    anchors <- head(anchors[order(-avg_log2FC)], GENOME_ENRICHED_ROWS)
    out[[as.character(cluster)]] <- list(
      top_enriched = .rows(up, c("auc")),
      top_depleted = .rows(down, character(0)),
      anchor_genes = .rows(anchors, character(0))
    )
  }
  out
}
.measurement_for <- function(pool, cluster, gene) {
  key <- paste(as.character(cluster), toupper(as.character(gene)), sep = "\x1f")
  pct_in <- unname((pool$de_stats$pct_in %||% list())[key])
  if (length(pct_in) && !is.na(pct_in)) {
    return(list(row = list(
      pct_in = pct_in, pct_out = unname(pool$de_stats$pct_out[key]),
      avg_log2FC = unname(pool$de_stats$avg_log2FC[key]), auc = unname(pool$de_stats$auc[key])
    ), how = "de_table"))
  }
  gpct_in <- unname((pool$genome_stats$pct_in %||% list())[key])
  if (length(gpct_in) && !is.na(gpct_in)) {
    return(list(row = list(
      pct_in = gpct_in, pct_out = unname(pool$genome_stats$pct_out[key]),
      avg_log2FC = unname(pool$genome_stats$avg_log2FC[key]), auc = unname(pool$genome_stats$auc[key])
    ), how = "genome_table"))
  }
  if (toupper(as.character(gene)) %in% (pool$measured_genes %||% character(0))) {
    return(list(row = NULL, how = "below_reporting_screen"))
  }
  list(row = NULL, how = "not_measured")
}
.top_clusters_for <- function(pool, gene, k) {
  rows <- pool$gene_clusters[[toupper(as.character(gene))]] %||% list()
  ranked <- rows[order(-vapply(rows, function(r) r$pct_in %||% 0, 0))]
  ranked[seq_len(min(as.integer(k), length(ranked)))]
}
attach_definers <- function(pool, definers_by_name) {
  n <- 0L
  for (cluster in names(pool$clusters)) {
    state <- pool$clusters[[cluster]]
    tissue_claimants <- list()
    for (e in state$candidates) {
      e_name <- as.character(e$cell_type)
      for (m in e$markers) {
        if (identical(m$polarity, "positive")) {
          gene_u <- toupper(as.character(m$gene))
          tissue_claimants[[gene_u]] <- unique(c(tissue_claimants[[gene_u]], e_name))
        }
      }
    }
    for (index in seq_along(state$candidates)) {
      entry <- state$candidates[[index]]
      name <- as.character(entry$cell_type)
      panel_genes <- unique(toupper(vapply(
        Filter(function(m) identical(m$polarity, "positive"), entry$markers), function(m) as.character(m$gene), ""
      )))
      rows <- list()
      for (d in definers_by_name[[name]] %||% list()) {
        gene_u <- toupper(as.character(d$gene))
        measurement <- .measurement_for(pool, cluster, gene_u)
        in_panel <- gene_u %in% panel_genes
        row <- list(
          gene = unname(pool$native_gene[[gene_u]]) %||% d$gene,
          n_pub = as.integer(d$n_pub), recommended = isTRUE(d$recommended),
          in_this_tissue_panel = in_panel
        )
        if (!in_panel) {
          others <- setdiff(tissue_claimants[[gene_u]] %||% character(0), name)
          if (length(others)) row$shared_with_page <- sort(.cp(others), method = "radix")
        }
        if (!is.null(measurement$row)) {
          row$pct_in <- measurement$row$pct_in
          row$pct_out <- measurement$row$pct_out
          row$auc <- measurement$row$auc
          row$raised <- .is_raised(measurement$row)
          peaks <- .top_clusters_for(pool, gene_u, DEFINER_TOP_CLUSTERS)
          row$peak_clusters <- as.list(vapply(peaks, function(c) sprintf("%s:%s", c$cluster_id, c$pct_in), ""))
        } else {
          row$measurement <- measurement$how
        }
        rows[[length(rows) + 1L]] <- row
      }
      entry$definers <- rows
      state$candidates[[index]] <- entry
      n <- n + as.integer(length(rows) > 0L)
    }
    pool$clusters[[cluster]] <- state
  }
  list(pool = pool, n = n)
}
attach_page_sharing <- function(pool, species_claimants) {
  n <- 0L
  for (cluster in names(pool$clusters)) {
    state <- pool$clusters[[cluster]]
    here <- unique(as.character(vapply(state$candidates, function(e) as.character(e$cell_type), "")))
    cache <- new.env(parent = emptyenv())
    others_for <- function(gene_u) {
      if (!exists(gene_u, envir = cache, inherits = FALSE)) {
        claimants <- vapply(species_claimants[[gene_u]] %||% list(), function(c) c[[1]], "")
        assign(gene_u, claimants[claimants %in% here], envir = cache)
      }
      get(gene_u, envir = cache, inherits = FALSE)
    }
    for (index in seq_along(state$candidates)) {
      entry <- state$candidates[[index]]
      name <- as.character(entry$cell_type)
      for (mi in seq_along(entry$markers)) {
        marker <- entry$markers[[mi]]
        if (!identical(marker$polarity, "positive")) {
          entry$markers[[mi]]$shared_on_page <- NULL
          next
        }
        shared <- setdiff(others_for(toupper(as.character(marker$gene))), name)
        if (length(shared)) {
          entry$markers[[mi]]$shared_on_page <- shared
          n <- n + 1L
        } else {
          entry$markers[[mi]]$shared_on_page <- NULL
        }
      }
      for (di in seq_along(entry$definers)) {
        row <- entry$definers[[di]]
        shared <- setdiff(others_for(toupper(as.character(row$gene))), name)
        prior <- setdiff(as.character(unlist(row$shared_with_page %||% list())), name)
        merged <- sort(.cp(union(shared, prior)), method = "radix")
        if (length(merged)) {
          entry$definers[[di]]$shared_with_page <- merged
          n <- n + 1L
        } else {
          entry$definers[[di]]$shared_with_page <- NULL
        }
      }
      state$candidates[[index]] <- entry
    }
    pool$clusters[[cluster]] <- state
  }
  list(pool = pool, n = n)
}
.partition_string <- function(row, with_sharers = FALSE) {
  tag <- sprintf("%s n%d", row$gene, as.integer(row$n_pub %||% 0L))
  if (!is.null(row$pct_in)) {
    tag <- paste0(tag, sprintf(" %s/%s", row$pct_in, row$pct_out), if (isTRUE(row$raised)) "" else "!")
  } else {
    tag <- paste0(tag, sprintf(" %s", if (nzchar(as.character(row$measurement %||% ""))) row$measurement else "unmeasured"))
  }
  shared <- as.character(unlist(row$shared %||% list()))
  if (isTRUE(with_sharers) && length(shared)) {
    shown <- shared[seq_len(min(3L, length(shared)))]
    tag <- paste0(tag, " ~", paste(shown, collapse = ","),
                  if (length(shared) > 3L) sprintf(",+%d", length(shared) - 3L) else "")
  }
  tag
}
.gene_rows <- function(entry) {
  by_gene <- list()
  for (marker in entry$markers %||% list()) {
    if (!identical(marker$polarity, "positive")) next
    gene_u <- toupper(as.character(marker$gene))
    by_gene[[gene_u]] <- list(
      gene = marker$gene, n_pub = marker$n_pub, pct_in = marker$pct_in, pct_out = marker$pct_out,
      auc = marker$auc, raised = .is_raised(marker), shared = marker$shared_on_page %||% list()
    )
  }
  for (row in entry$definers %||% list()) {
    gene_u <- toupper(as.character(row$gene))
    if (!is.null(by_gene[[gene_u]])) {
      if (as.integer(row$n_pub %||% 0L) > as.integer(by_gene[[gene_u]]$n_pub %||% 0L)) {
        by_gene[[gene_u]]$n_pub <- row$n_pub
      }
      next
    }
    by_gene[[gene_u]] <- list(
      gene = row$gene, n_pub = row$n_pub, pct_in = row$pct_in, pct_out = row$pct_out,
      auc = row$auc, raised = isTRUE(row$raised), measurement = row$measurement,
      shared = row$shared_with_page %||% list()
    )
  }
  if (!length(by_gene)) return(list())
  has_auc <- vapply(by_gene, function(r) !is.null(r$auc), TRUE)
  auc_val <- vapply(by_gene, function(r) if (is.null(r$auc)) 0 else as.numeric(r$auc), 0)
  ord <- order(as.integer(!has_auc), -auc_val, .cp(vapply(by_gene, function(r) as.character(r$gene), "")),
               method = "radix")
  unname(by_gene[ord])
}
lineage_partition <- function(entry) {
  rows <- .gene_rows(entry)
  own <- Filter(function(r) !length(r$shared), rows)
  shared_rows <- Filter(function(r) length(r$shared) > 0, rows)
  list(
    own = as.list(vapply(own[seq_len(min(PARTITION_OWN_SHOWN, length(own)))], .partition_string, "")),
    own_total = length(own),
    shared = as.list(vapply(
      shared_rows[seq_len(min(PARTITION_OWN_SHOWN, length(shared_rows)))],
      function(r) .partition_string(r, with_sharers = TRUE), ""
    )),
    shared_total = length(shared_rows)
  )
}
pair_partition <- function(entry_a, entry_b) {
  name_a <- as.character(entry_a$cell_type)
  name_b <- as.character(entry_b$cell_type)
  rows_a <- .gene_rows(entry_a)
  names(rows_a) <- vapply(rows_a, function(r) toupper(as.character(r$gene)), "")
  rows_b <- .gene_rows(entry_b)
  names(rows_b) <- vapply(rows_b, function(r) toupper(as.character(r$gene)), "")
  both_keys <- character(0)
  for (g in names(rows_a)) {
    if (g %in% names(rows_b) || name_b %in% as.character(unlist(rows_a[[g]]$shared %||% list()))) {
      both_keys <- c(both_keys, g)
    }
  }
  for (g in names(rows_b)) {
    if (g %in% names(rows_a) || name_a %in% as.character(unlist(rows_b[[g]]$shared %||% list()))) {
      both_keys <- c(both_keys, g)
    }
  }
  both_keys <- unique(both_keys)
  only_a <- rows_a[setdiff(names(rows_a), both_keys)]
  only_b <- rows_b[setdiff(names(rows_b), both_keys)]
  both <- lapply(both_keys, function(g) rows_a[[g]] %||% rows_b[[g]])
  .order_rows <- function(rows) {
    if (!length(rows)) return(rows)
    has_auc <- vapply(rows, function(r) !is.null(r$auc), TRUE)
    auc_val <- vapply(rows, function(r) if (is.null(r$auc)) 0 else as.numeric(r$auc), 0)
    ord <- order(as.integer(!has_auc), -auc_val, .cp(vapply(rows, function(r) as.character(r$gene), "")),
                 method = "radix")
    unname(rows[ord])
  }
  out <- list(pair = list(name_a, name_b))
  out[[paste0("only_", name_a)]] <- as.list(vapply(
    .order_rows(only_a)[seq_len(min(PARTITION_OWN_SHOWN, length(only_a)))], .partition_string, ""
  ))
  out[[paste0("only_", name_b)]] <- as.list(vapply(
    .order_rows(only_b)[seq_len(min(PARTITION_OWN_SHOWN, length(only_b)))], .partition_string, ""
  ))
  out$both <- as.list(vapply(
    .order_rows(both)[seq_len(min(PARTITION_OWN_SHOWN, length(both)))], .partition_string, ""
  ))
  totals <- list(length(only_a), length(only_b), length(both))
  names(totals) <- c(name_a, name_b, "both")
  out$totals <- totals
  out
}
attach_genome_claims <- function(pool, claimants, known_genes, markers_all) {
  genome <- pool$genome_de %||% list()
  gaps <- list()
  if (!is.null(markers_all) && nrow(markers_all)) {
    ma <- as.data.table(markers_all)
    ma[, group := as.character(group)]
    for (cluster in unique(ma$group)) {
      cluster_value <- cluster
      dom <- ma[group == cluster_value & pct_in >= DM_GATE_DOMINANT_PCT]
      dom[, gap := pct_in - pct_out]
      dom <- head(dom[order(-gap)], GENOME_DETECTION_GAP_ROWS * 2L)
      gaps[[as.character(cluster)]] <- lapply(seq_len(nrow(dom)), function(i) {
        list(
          gene = as.character(dom$feature[i]), pct_in = .round(dom$pct_in[i], 1L),
          pct_out = .round(dom$pct_out[i], 1L), log2FC = .round(dom$avg_log2FC[i], 2L),
          auc = .round(dom$auc[i], 3L)
        )
      })
    }
  }
  .claims <- function(gene) {
    rows <- claimants[[toupper(as.character(gene))]] %||% list()
    rows <- rows[seq_len(min(GENOME_CLAIMANTS_SHOWN, length(rows)))]
    lapply(rows, function(r) list(cell_type = r[[1]], n_pub = as.integer(r[[2]])))
  }
  .filter_rows <- function(rows, cap) {
    kept <- list()
    dropped <- 0L
    for (row in rows) {
      if (!(toupper(as.character(row$gene)) %in% known_genes)) {
        dropped <- dropped + 1L
        next
      }
      row$claimed_by_species <- .claims(row$gene)
      kept[[length(kept) + 1L]] <- row
      if (length(kept) >= cap) break
    }
    list(kept = kept, dropped = dropped)
  }
  for (cluster in names(pool$clusters)) {
    block <- genome[[cluster]] %||% list()
    fe <- .filter_rows(block$top_enriched %||% list(), GENOME_ENRICHED_ROWS)
    fd <- .filter_rows(block$top_depleted %||% list(), GENOME_DEPLETED_ROWS)
    fg <- .filter_rows(gaps[[cluster]] %||% list(), GENOME_DETECTION_GAP_ROWS)
    block$top_enriched <- fe$kept
    block$top_depleted <- fd$kept
    block$top_detection_gap <- fg$kept
    block$rows_not_curated_anywhere <- list(
      enriched = fe$dropped, depleted = fd$dropped, detection_gap = fg$dropped
    )
    genome[[cluster]] <- block
  }
  pool$genome_de <- genome
  pool
}
.mark_on_page <- function(rows, on_page) {
  lapply(rows %||% list(), function(row) {
    row$claimed_by_species <- lapply(row$claimed_by_species %||% list(), function(c) {
      c$on_page <- as.character(c$cell_type) %in% on_page
      c
    })
    row
  })
}
genome_pages <- function(pool, cluster, tier) {
  block <- (pool$genome_de %||% list())[[as.character(cluster)]] %||% list()
  on_page <- .candidate_names(pool, cluster, tier)
  list(
    cluster_top_enriched = .mark_on_page(block$top_enriched, on_page),
    cluster_top_detection_gap = .mark_on_page(block$top_detection_gap, on_page),
    cluster_top_depleted = .mark_on_page(block$top_depleted, on_page),
    genome_rows_not_curated_anywhere = block$rows_not_curated_anywhere %||% list()
  )
}
definer_block <- function(entry) {
  entry$definers %||% list()
}
scma_attach_recommended <- function(pool, recommended_pairs) {
  pairs <- unique(as.character(recommended_pairs))
  marked <- 0L
  for (cluster in names(pool$clusters)) {
    entries <- pool$clusters[[cluster]]$candidates
    for (index in seq_along(entries)) {
      name <- as.character(entries[[index]]$cell_type)
      markers <- entries[[index]]$markers
      if (!length(markers)) next
      keys <- paste(name, toupper(vapply(markers, function(m) {
        as.character(m$gene)
      }, "")), sep = "\x1f")
      flags <- keys %chin% pairs
      for (position in seq_along(markers)) {
        markers[[position]]$recommended <- flags[position]
      }
      marked <- marked + sum(flags)
      entries[[index]]$markers <- markers
    }
    pool$clusters[[cluster]]$candidates <- entries
  }
  list(pool = pool, marked = as.integer(marked))
}
.measured_row <- function(marker) {
  list(
    gene = marker$gene,
    n_pub = as.integer(marker$n_pub %||% 0L),
    recommended = isTRUE(marker$recommended),
    pct_in = marker$pct_in,
    pct_out = marker$pct_out,
    avg_log2FC = marker$avg_log2FC,
    auc = marker$auc,
    significant = isTRUE(marker$significant)
  )
}
.detected_exclusions <- function(entry) {
  rows <- list()
  for (marker in entry$markers %||% list()) {
    if (!identical(marker$polarity, "negative")) next
    pct_in <- marker$pct_in
    if (is.null(pct_in) || !length(pct_in) || is.na(pct_in[1]) ||
      as.numeric(pct_in[1]) < NEGATIVE_SOURCE_MIN_PCT_IN) next
    rows[[length(rows) + 1L]] <- .measured_row(marker)
  }
  if (!length(rows)) {
    return(rows)
  }
  rows[order(-vapply(rows, function(r) .num0(r$pct_in), 0))]
}
.detected_exclusion_markers <- function(entry) {
  rows <- Filter(function(marker) {
    identical(marker$polarity, "negative") && !is.null(marker$pct_in) &&
      length(marker$pct_in) && !is.na(marker$pct_in[1]) &&
      as.numeric(marker$pct_in[1]) >= NEGATIVE_SOURCE_MIN_PCT_IN
  }, entry$markers %||% list())
  if (!length(rows)) {
    return(rows)
  }
  rows[order(-vapply(rows, function(r) .num0(r$pct_in), 0))]
}
.candidate_facts <- function(pool, cluster, entry) {
  species <- as.character(pool$context$species %||% "")
  name <- as.character(entry$cell_type)
  positives <- .positive_rows(entry)
  by_pub <- .by_npub(positives)
  raised <- Filter(function(m) .is_raised(m) && isTRUE(m$significant), positives)
  papers <- scma_label_papers(species, name)
  list(
    cell_type = name,
    retrieval_rank = as.integer(entry$retrieval_rank %||% 0L),
    tier = .tier_of(entry),
    n_papers = if (is.null(papers)) NULL else as.integer(papers),
    recommended_markers = lapply(
      Filter(function(m) isTRUE(m$recommended), positives), .measured_row
    ),
    best_published_markers = lapply(by_pub[seq_len(min(2L, length(by_pub)))], .measured_row),
    positive_markers_measured = length(positives),
    positive_markers_raised_and_significant = length(raised),
    detected_exclusions = .detected_exclusions(entry),
    own_vs_shared = lineage_partition(entry)
  )
}
.candidate_facts_table <- function(pool, cluster, tier = NULL) {
  entries <- pool$clusters[[as.character(cluster)]]$candidates
  if (!is.null(tier)) {
    entries <- Filter(function(e) .tier_of(e) <= as.integer(tier), entries)
  }
  lapply(entries, function(entry) .candidate_facts(pool, cluster, entry))
}
.parent_own_marker_sets <- function(parent_entry) {
  positives <- .positive_rows(parent_entry)
  rec <- unique(vapply(
    Filter(function(m) isTRUE(m$recommended), positives),
    function(m) as.character(m$gene), ""
  ))
  top <- .by_npub(positives)
  top2 <- vapply(top[seq_len(min(2L, length(top)))], function(m) as.character(m$gene), "")
  list(recommended = as.character(rec), top2 = as.character(top2))
}
.finer_vs_parent_facts <- function(entry, parent_entry) {
  parent_sets <- .parent_own_marker_sets(parent_entry)
  parent_own <- toupper(c(parent_sets$recommended, parent_sets$top2))
  positives <- .positive_rows(entry)
  by_pub <- .by_npub(positives)
  own_rows <- Filter(function(m) !(toupper(as.character(m$gene)) %in% parent_own), positives)
  own_by_pub <- .by_npub(own_rows)
  own_raised <- Filter(function(m) .is_raised(m) && isTRUE(m$significant), own_by_pub)
  parent_raised_n <- length(Filter(function(m) {
    identical(m$polarity, "positive") && .is_raised(m) && isTRUE(m$significant)
  }, parent_entry$markers %||% list()))
  list(
    finer_name = as.character(entry$cell_type),
    parent = as.character(parent_entry$cell_type),
    parent_own_markers = list(
      recommended = as.list(
        parent_sets$recommended[order(.cp(toupper(parent_sets$recommended)), method = "radix")]
      ),
      best_published = as.list(parent_sets$top2)
    ),
    own_markers = lapply(own_by_pub, .measured_row),
    own_markers_raised_and_significant = lapply(own_raised, .measured_row),
    finer_best_published = lapply(by_pub[seq_len(min(2L, length(by_pub)))], .measured_row),
    parent_positive_markers_raised_and_significant = as.integer(parent_raised_n)
  )
}
.borrow_admission_row <- function(pool, cluster, entry) {
  name <- as.character(entry$cell_type)
  positives <- .positive_rows(entry)
  rec <- Filter(function(m) isTRUE(m$recommended), positives)
  rec_carried <- Filter(.carried, rec)
  ordered <- .by_npub(positives)
  top2 <- ordered[seq_len(min(2L, length(ordered)))]
  top2_carried <- Filter(.carried, top2)
  if (length(rec)) {
    gate <- length(rec_carried) >= 1L
    gate_basis <- "recommended"
  } else if (length(top2) >= 2L) {
    gate <- length(top2_carried) >= 1L
    gate_basis <- "top2_n_pub"
  } else {
    gate <- TRUE
    gate_basis <- "too_thin_to_gate"
  }
  if (gate && scma_definers_absent(positives, DM_GATE_TOP)) {
    gate <- FALSE
    gate_basis <- "top2_absent"
  }
  .genes <- function(rows) as.list(unname(vapply(rows, function(m) as.character(m$gene), "")))
  list(
    cell_type = name,
    recommended_measured = .genes(rec),
    recommended_carried = .genes(rec_carried),
    top2_n_pub = .genes(top2),
    top2_carried = .genes(top2_carried),
    identity_gate = if (isTRUE(gate)) "pass" else "fail",
    identity_gate_basis = gate_basis
  )
}
.run_tool <- function(pool, cluster, tool, args, sources = NULL) {
  if (!is.list(args)) args <- list()
  if (identical(tool, "sources")) {
    return(.tool_sources(pool, cluster, args, sources))
  }
  if (identical(tool, "gene_across_clusters")) {
    return(.tool_gene_across_clusters(pool, args))
  }
  if (identical(tool, "candidates_with_gene")) {
    return(.tool_candidates_with_gene(pool, cluster, args))
  }
  list(
    tool = as.character(tool)[1],
    error = sprintf("unknown tool; available: %s", paste(TOOLS, collapse = ", "))
  )
}
.run_review_tool <- function(pool, cluster, labels, tool, args, sources = NULL) {
  if (!is.list(args)) args <- list()
  if (!(as.character(tool)[1] %in% REVIEW_TOOLS)) {
    return(list(
      tool = as.character(tool)[1],
      error = sprintf("unknown tool; available: %s", paste(REVIEW_TOOLS, collapse = ", "))
    ))
  }
  if (identical(tool, "gene_across_clusters")) {
    return(.tool_gene_across_clusters(pool, args))
  }
  if (identical(tool, "candidates_with_gene")) {
    return(.tool_candidates_with_gene(pool, cluster, args))
  }
  asked <- as.character(args$label %||% args$candidate %||% "")[1]
  if (nzchar(asked) && !(asked %in% labels)) {
    return(list(
      tool = "sources",
      error = sprintf(
        "'%s' is not under test here; available: %s", asked, paste(labels, collapse = ", ")
      )
    ))
  }
  label <- if (nzchar(asked)) asked else (labels[1] %||% "")
  args$candidate <- label
  .tool_sources(pool, cluster, args, sources)
}
.is_source_server <- function(sources) {
  is.environment(sources) && is.function(sources$take)
}
.tool_sources <- function(pool, cluster, args, sources) {
  candidate <- as.character(args$candidate %||% "")[1]
  if (is.null(.find_candidate(pool, cluster, candidate))) {
    return(list(
      tool = "sources",
      error = sprintf("'%s' is not a candidate of this cluster", candidate)
    ))
  }
  genes <- unlist(args$genes %||% list())
  genes <- trimws(as.character(genes))
  genes <- genes[nzchar(genes)]
  truncated <- length(genes) > SOURCE_GENES_PER_QUERY
  genes <- genes[seq_len(min(SOURCE_GENES_PER_QUERY, length(genes)))]
  found <- list()
  left <- list()
  for (gene in genes) {
    key <- toupper(gene)
    if (is.null(sources)) {
      found[[key]] <- list()
      next
    }
    records <- if (.is_source_server(sources)) {
      answer <- sources$take(candidate, key)
      if (answer$remaining > 0L || isTRUE(answer$limit_reached)) {
        left[[key]] <- if (isTRUE(answer$limit_reached)) {
          "limit reached"
        } else {
          sprintf("%d more", answer$remaining)
        }
      }
      answer$sources
    } else {
      sources(candidate, key, SOURCES_PER_GENE)
    }
    found[[key]] <- lapply(records, .sentence_record)
  }
  out <- list(
    tool = "sources", candidate = candidate, sources = found, truncated = truncated
  )
  if (length(left)) out$not_yet_shown <- left
  out
}
scma_register_packet_sources <- function(server, pool, cluster) {
  if (!.is_source_server(server) || !is.function(server$note_delivered)) {
    return(invisible(NULL))
  }
  for (entry in pool$clusters[[as.character(cluster)]]$candidates) {
    candidate <- as.character(entry$cell_type)
    records_by_gene <- entry$exclusion_sources %||% list()
    for (gene in names(records_by_gene)) {
      server$note_delivered(candidate, gene, records_by_gene[[gene]])
    }
  }
  invisible(NULL)
}
.tool_gene_across_clusters <- function(pool, args) {
  gene <- toupper(trimws(as.character(args$gene %||% "")[1]))
  rows <- pool$gene_clusters[[gene]]
  if (is.null(rows) || !length(rows)) {
    return(list(
      tool = "gene_across_clusters", gene = gene,
      error = "not measured in this dataset"
    ))
  }
  rows <- rows[order(-vapply(rows, function(r) r$pct_in %||% 0, 0))]
  list(
    tool = "gene_across_clusters",
    gene = pool$native_gene[[gene]] %||% gene,
    clusters = rows
  )
}
.tool_candidates_with_gene <- function(pool, cluster, args) {
  gene <- toupper(trimws(as.character(args$gene %||% "")[1]))
  here <- .candidate_names(pool, cluster)
  rows <- Filter(
    function(row) row$cell_type %in% here,
    pool$gene_carriers[[gene]] %||% list()
  )
  if (length(rows)) {
    rows <- rows[order(
      -vapply(rows, function(r) as.numeric(r$n_pub), 0),
      vapply(rows, function(r) r$cell_type, "")
    )]
  }
  list(
    tool = "candidates_with_gene",
    gene = pool$native_gene[[gene]] %||% gene,
    candidates = rows
  )
}
.claim_warnings <- function(pool, cluster, claims) {
  lines <- list()
  for (claim in claims) {
    role <- as.character(claim[[1]])
    name <- as.character(claim[[2]])
    entry <- .find_candidate(pool, cluster, name)
    if (is.null(entry)) next
    positive <- Filter(function(m) identical(m$polarity, "positive"), entry$markers)
    unraised <- Filter(function(m) !.is_raised(m), positive)
    if (!length(positive) || !length(unraised)) next
    shown <- unraised[seq_len(min(WARNING_GENES_SHOWN, length(unraised)))]
    listed <- paste(vapply(shown, function(m) {
      sprintf(
        "%s(%s/%s)", m$gene,
        formatC(m$pct_in, format = "f", digits = PCT_DECIMALS),
        formatC(m$pct_out, format = "f", digits = PCT_DECIMALS)
      )
    }, ""), collapse = ",")
    lines[[length(lines) + 1L]] <- list(name, sprintf(
      "%s %s: not_raised %d/%d | top_n_pub: %s | +%d more",
      role, name, length(unraised), length(positive), listed,
      length(unraised) - length(shown)
    ))
  }
  lines
}
