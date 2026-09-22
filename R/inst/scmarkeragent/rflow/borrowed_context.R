suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})
B_BORROW <- CFG$borrowed_context
BORROW_MIN_ANCHORS <- as.integer(B_BORROW$min_anchor_genes)
BORROW_MIN_PUB <- DECISIVE_MIN_PUB
BORROW_MAX_PER_CLUSTER <- as.integer(B_BORROW$max_per_cluster)
BORROW_ANCHOR_ROWS <- GENOME_ENRICHED_ROWS
BORROW_ANCHOR_MIN_PCT <- DM_GATE_DOMINANT_PCT
BORROW_OWN_SIGNIFICANT_MIN <- as.integer(B_BORROW$own_significant_markers_min)
BORROW_QUALIFIER_PROMISCUITY_MIN <- as.integer(B_BORROW$qualifier_promiscuity_min)
GENERIC_HEAD_NOUNS <- as.character(B_BORROW$generic_head_nouns)
eps <- as.numeric(CFG$retrieval$epsilon)
SCREEN_FOREIGN <- "foreign"
SCREEN_PLAUSIBLE <- "plausible"
SCREEN_VERDICTS <- c(SCREEN_PLAUSIBLE, SCREEN_FOREIGN)
SCREEN_REJECTED_PREFIX <- "tissue screen: "
scma_borrow_norm_name <- function(value) {
  text <- tolower(trimws(as.character(value)))
  text <- gsub(",", " ", text, fixed = TRUE)
  trimws(gsub("[[:space:]]+", " ", text))
}
.LABEL_PAPERS_NORM <- new.env(parent = emptyenv())
.label_papers_normalized <- function() {
  if (!exists("table", envir = .LABEL_PAPERS_NORM, inherits = FALSE)) {
    path <- .env(
      "SCMA_LABEL_PAPERS",
      file.path(PACKAGE_ROOT, "resources", "label_papers_v4.csv")
    )
    table <- if (file.exists(path)) {
      rows <- fread(path, colClasses = "character", showProgress = FALSE)
      .papers_table(
        paste(
          as.character(rows$species),
          scma_borrow_norm_name(rows$cell_type),
          sep = "\x1f"
        ),
        rows$n_papers
      )
    } else {
      integer(0)
    }
    assign("table", table, envir = .LABEL_PAPERS_NORM)
  }
  get("table", envir = .LABEL_PAPERS_NORM, inherits = FALSE)
}
.lineage_classes <- function(parents) {
  items <- sort(.cp(unique(as.character(parents))), method = "radix")
  n <- length(items)
  if (!n) {
    return(0L)
  }
  state <- new.env(parent = emptyenv())
  state$root <- seq_len(n)
  find <- function(a) {
    while (state$root[a] != a) {
      state$root[a] <- state$root[state$root[a]]
      a <- state$root[a]
    }
    a
  }
  tokens <- lapply(items, function(p) {
    setdiff(strsplit(p, " ", fixed = TRUE)[[1]], GENERIC_HEAD_NOUNS)
  })
  if (n > 1L) {
    for (i in seq_len(n - 1L)) {
      for (j in seq(i + 1L, n)) {
        if (length(intersect(tokens[[i]], tokens[[j]]))) {
          ri <- find(i)
          rj <- find(j)
          if (ri != rj) state$root[ri] <- rj
        }
      }
    }
  }
  length(unique(vapply(seq_len(n), find, 0L)))
}
scma_qualifier_lexicon <- function(db, cutoff = BORROW_QUALIFIER_PROMISCUITY_MIN) {
  papers <- .label_papers_normalized()
  frame <- unique(as.data.table(db)[
    !is.na(species) & !is.na(cell_type),
    .(species = as.character(species), cell_type = as.character(cell_type))
  ])
  frame <- frame[nzchar(species) & nzchar(cell_type)]
  frame[, norm := scma_borrow_norm_name(cell_type)]
  raw <- list()
  for (one in unique(frame$species)) {
    names_here <- unique(frame[species == one]$norm)
    present <- rep(TRUE, length(names_here))
    names(present) <- names_here
    for (name in names_here) {
      tokens <- strsplit(name, " ", fixed = TRUE)[[1]]
      if (length(tokens) < 2L) next
      parent <- paste(tokens[-1], collapse = " ")
      if (is.na(present[parent]) || parent %in% GENERIC_HEAD_NOUNS) next
      n_papers <- unname(papers[paste(one, parent, sep = "\x1f")])
      if (!length(n_papers) || is.na(n_papers)) n_papers <- 0L
      if (n_papers <= SPARSE_LABEL_PAPERS_MAX) next
      raw[[tokens[1]]] <- unique(c(raw[[tokens[1]]], parent))
    }
  }
  keep <- vapply(raw, function(parents) {
    .lineage_classes(sort(.cp(parents), method = "radix")) >= as.integer(cutoff)
  }, TRUE)
  tokens <- sort(.cp(names(raw)[keep]), method = "radix")
  digest_value <- digest(
    enc2utf8(paste(tokens, collapse = "\n")),
    algo = "sha256", serialize = FALSE
  )
  list(
    tokens = as.list(tokens), size = length(tokens), sha256 = digest_value,
    cutoff = as.integer(cutoff), vocabulary_pairs = nrow(frame),
    modifier_tokens = length(raw)
  )
}
scma_qualified_form_of <- function(claimant, eligible, lexicon) {
  tokens <- strsplit(scma_borrow_norm_name(claimant), " ", fixed = TRUE)[[1]]
  eligible <- as.character(eligible)
  local_norm <- scma_borrow_norm_name(eligible)
  keep <- !(local_norm %in% GENERIC_HEAD_NOUNS)
  local_norm <- local_norm[keep]
  local_name <- eligible[keep]
  if (length(tokens) < 2L) {
    return(NULL)
  }
  for (k in seq_len(length(tokens) - 1L)) {
    extra <- tokens[seq_len(k)]
    rest <- paste(tokens[seq(k + 1L, length(tokens))], collapse = " ")
    hit <- which(local_norm == rest)
    if (length(hit) && all(extra %in% lexicon)) {
      return(list(
        eligible_candidate = local_name[hit[length(hit)]],
        modifier = paste(extra, collapse = " ")
      ))
    }
  }
  NULL
}
scma_borrow_resource_slice <- function(db, species, disease) {
  terms <- .normalize_disease_query(disease)
  frame <- as.data.table(db)
  wanted_species <- as.character(species)
  keep <- frame$species == wanted_species &
    frame$marker_polarity %in% c("positive", "negative")
  frame <- frame[keep]
  qc_ok <- frame$gene_qc_pass == "TRUE" | is.na(frame$gene_qc_pass)
  mask <- qc_ok & !frame$is_in_vitro &
    !is.na(frame$gene_symbol) & nzchar(frame$gene_symbol) &
    !is.na(frame$cell_type) & nzchar(frame$cell_type)
  if (!is.null(terms)) {
    mask <- mask & (.normalize_disease_name(frame$disease_normalized) %in% terms)
  }
  out <- frame[mask, .(
    tissue_type = as.character(tissue_type), gene_symbol = as.character(gene_symbol),
    cell_type = as.character(cell_type), marker_polarity = as.character(marker_polarity),
    n_pub_support = as.integer(n_pub_support), confidence_tier = as.character(confidence_tier),
    is_recommended_marker = as.character(is_recommended_marker)
  )]
  out[, gu := toupper(gene_symbol)]
  out[, recommended := toupper(is_recommended_marker) == "TRUE"]
  out[]
}
scma_translate_resource <- function(resource, gene_map) {
  out <- copy(resource)
  mapped <- unname(gene_map[out$gu])
  out <- out[!is.na(mapped)]
  mapped <- mapped[!is.na(mapped)]
  out[, donor_gene := gene_symbol]
  out[, gene_symbol := mapped]
  out[, gu := toupper(as.character(gene_symbol))]
  out[]
}
.borrow_n_tissues <- function(tissues_key, name) {
  value <- unname(tissues_key[as.character(name)])
  if (!length(value) || is.na(value)) 0L else as.integer(value)
}
.borrow_panel_meta <- function(rows) {
  tier_rank <- c(high = 0L, medium = 1L, low = 2L)
  rows <- copy(as.data.table(rows))
  rank <- tier_rank[tolower(as.character(rows$confidence_tier))]
  rank[is.na(rank)] <- 3L
  rows[, tr := as.integer(rank)]
  rows <- rows[order(-n_pub_support, tr, seq_len(nrow(rows)))]
  name <- as.character(rows$cell_type[1])
  meta <- list()
  positive <- character(0)
  negative <- character(0)
  tissues <- list()
  for (i in seq_len(nrow(rows))) {
    gene <- as.character(rows$gu[i])
    pol <- as.character(rows$marker_polarity[i])
    tissues[[gene]] <- unique(c(tissues[[gene]], as.character(rows$tissue_type[i])))
    key <- paste(name, gene, sep = "\x1f")
    if (!is.null(meta[[key]])) {
      next
    }
    meta[[key]] <- list(
      n_pub = as.integer(rows$n_pub_support[i]),
      tier = as.character(rows$confidence_tier[i]), polarity = pol
    )
    if (identical(pol, "positive")) {
      positive <- c(positive, gene)
    } else {
      negative <- c(negative, gene)
    }
  }
  list(meta = meta, positive = positive, negative = negative, tissues = tissues)
}
.borrow_curated_table <- function(name, meta) {
  if (!length(meta)) {
    table <- data.table(
      candidate = character(0), gene_key = character(0),
      marker_polarity = character(0), n_pub = integer(0), tier = character(0)
    )
  } else {
    table <- data.table(
      candidate = name,
      gene_key = vapply(strsplit(names(meta), "\x1f", fixed = TRUE), `[`, "", 2L),
      marker_polarity = vapply(meta, function(m) as.character(m$polarity), ""),
      n_pub = vapply(meta, function(m) as.integer(m$n_pub), 0L),
      tier = vapply(meta, function(m) as.character(m$tier), "")
    )
  }
  setkey(table, candidate, gene_key)
  table
}
.borrow_percentile_rank <- function(values) {
  if (length(values) <= 1L) {
    return(rep(1, length(values)))
  }
  frank(round(values, 10L), ties.method = "average") / length(values)
}
.borrow_specificity <- function(universe) {
  if (!nrow(universe)) {
    return(setNames(numeric(0), character(0)))
  }
  total <- uniqueN(universe$cell)
  holders <- universe[, .(n = uniqueN(cell)), by = gene_key]
  setNames(log((total + 1) / (holders$n + 1)) / log(total + 1), holders$gene_key)
}
scma_borrow_candidates <- function(pool, resource, measured_genes, markers_all,
                                  sources = NULL, screen = NULL,
                                  species = NULL, cross_species_resource = NULL,
                                  ortho_dir = NULL) {
  native_genes <- as.character(unlist(pool$native_menu_genes %||% list()))
  all_positive <- resource[marker_polarity == "positive"]
  positive <- all_positive[n_pub_support >= BORROW_MIN_PUB]
  by_type_genes <- lapply(split(positive$gu, positive$cell_type), unique)
  rec_rows <- positive[recommended == TRUE]
  rec_by_type <- lapply(split(rec_rows$gu, rec_rows$cell_type), unique)
  npub_by_pair <- positive[, .(n = max(n_pub_support)), by = .(cell_type, gu)]
  npub_key <- setNames(npub_by_pair$n, paste(npub_by_pair$cell_type, npub_by_pair$gu, sep = "\x1f"))
  panel_genes_by_type <- lapply(split(all_positive$gu, all_positive$cell_type), unique)
  tissues_by_type <- resource[, .(n = uniqueN(tissue_type)), by = cell_type]
  tissues_key <- setNames(tissues_by_type$n, tissues_by_type$cell_type)
  cross_tables <- list()
  if (!is.null(species) && length(cross_species_resource) && !is.null(ortho_dir)) {
    for (other in names(cross_species_resource)) {
      other_resource <- cross_species_resource[[other]]
      gene_map <- ortho_one_to_one(other, species, ortho_dir)
      donor_positive <- other_resource[
        marker_polarity == "positive" & n_pub_support >= BORROW_MIN_PUB
      ]
      translated <- scma_translate_resource(donor_positive, gene_map)
      translated <- translated[!(cell_type %in% names(by_type_genes))]
      by_type_genes_x <- lapply(split(translated$gu, translated$cell_type), unique)
      rec_rows_x <- translated[recommended == TRUE]
      rec_by_type_x <- lapply(split(rec_rows_x$gu, rec_rows_x$cell_type), unique)
      npub_by_pair_x <- translated[, .(n = max(n_pub_support)), by = .(cell_type, gu)]
      npub_key_x <- setNames(
        npub_by_pair_x$n, paste(npub_by_pair_x$cell_type, npub_by_pair_x$gu, sep = "\x1f")
      )
      tissues_by_type_x <- other_resource[, .(n = uniqueN(tissue_type)), by = cell_type]
      tissues_key_x <- setNames(tissues_by_type_x$n, tissues_by_type_x$cell_type)
      cross_tables[[other]] <- list(
        translated = translated,
        by_type_genes = by_type_genes_x,
        rec_by_type = rec_by_type_x,
        npub_key = npub_key_x,
        tissues_key = tissues_key_x,
        ortholog_pairs = length(gene_map)
      )
    }
  }
  ma <- as.data.table(markers_all)
  ma[, group := as.character(group)]
  ma[, gu := toupper(as.character(feature))]
  sig_ok <- sig_pass(ma$avg_log2FC, ma$pct_in / 100, ma$padj, ma$auc)
  significant <- split(ma$gu[sig_ok], ma$group[sig_ok])
  significant <- lapply(significant, function(g) unique(g))
  wide <- dcast(ma, gu ~ group, value.var = "avg_log2FC", fun.aggregate = function(x) x[1])
  rel_genes <- wide$gu
  cluster_cols <- setdiff(names(wide), "gu")
  rel_mat <- as.matrix(wide[, ..cluster_cols])
  relative <- apply(rel_mat, 2, function(col) {
    frank(col, ties.method = "average", na.last = "keep") / sum(!is.na(col))
  })
  dimnames(relative) <- list(rel_genes, cluster_cols)
  universe_parts <- c(
    list(all_positive[, .(cell = cell_type, gene_key = gu)]),
    lapply(cross_tables, function(tab) tab$translated[, .(cell = cell_type, gene_key = gu)])
  )
  universe <- rbindlist(universe_parts, use.names = TRUE)
  claimant_specificity <- .borrow_specificity(universe)
  display_specificity <- pool$marker_specificity %||% list()
  events <- list()
  pending <- list()
  for (cluster in names(pool$clusters)) {
    state <- pool$clusters[[cluster]]
    genome <- (pool$genome_de %||% list())[[cluster]] %||% list()
    anchors <- toupper(vapply(genome$anchor_genes %||% list(), function(a) {
      as.character(a$gene)
    }, ""))
    anchor_set <- unique(anchors)
    log <- list(
      anchor_genes = as.list(anchors),
      unclaimed_anchor_genes = list(),
      native_best = list(cell_type = "", anchors_covered = 0L),
      considered = list(), borrowed = list(), name_gated = list(),
      qualifier_gate = "removed"
    )
    events[[cluster]] <- log
    if (!length(anchor_set)) next
    unclaimed <- anchor_set[!(anchor_set %in% native_genes)]
    native_best <- 0L
    native_best_name <- ""
    for (entry in state$candidates) {
      covered <- sum(vapply(entry$markers, function(m) {
        identical(m$polarity, "positive") && toupper(as.character(m$gene)) %in% anchor_set
      }, TRUE))
      if (covered > native_best) {
        native_best <- as.integer(covered)
        native_best_name <- as.character(entry$cell_type)
      }
    }
    admitted <- unique(as.character(vapply(state$candidates, function(e) {
      as.character(e$cell_type)
    }, "")))
    log$unclaimed_anchor_genes <- as.list(unclaimed)
    log$native_best <- list(cell_type = native_best_name, anchors_covered = native_best)
    events[[cluster]] <- log
    if (length(unclaimed) < BORROW_MIN_ANCHORS) next
    unclaimed_set <- unclaimed
    hits_here <- significant[[cluster]] %||% character(0)
    relative_col <- if (cluster %in% colnames(relative)) relative[, cluster] else NULL
    .scored_claimant <- function(name, genes, rec_set, npub_key_set, panel_set, donor) {
      if (name %in% admitted) return(NULL)
      hit <- sort(.cp(intersect(genes, unclaimed_set)), method = "radix")
      if (length(hit) < BORROW_MIN_ANCHORS) return(NULL)
      panel <- sort(.cp(intersect(panel_set, measured_genes)), method = "radix")
      sig_hits <- intersect(panel, hits_here)
      if (!length(panel) || !length(sig_hits)) return(NULL)
      w_hit <- unname(claimant_specificity[sig_hits])
      w_hit[is.na(w_hit)] <- 0
      w_panel <- unname(claimant_specificity[panel])
      w_panel[is.na(w_panel)] <- 0
      perc <- if (is.null(relative_col)) {
        rep(0, length(panel))
      } else {
        value <- unname(relative_col[panel])
        value[is.na(value)] <- 0
        value
      }
      marker_level <- sum(w_hit) / sqrt(length(panel))
      cluster_level <- if (sum(w_panel) > 0) {
        sum(w_panel * perc) / max(sum(w_panel), eps)
      } else {
        mean(perc)
      }
      list(
        cell_type = name, anchors_covered = length(hit), anchor_genes = hit,
        recommended_anchors = sort(.cp(intersect(hit, rec_set[[name]] %||% character(0))), method = "radix"),
        anchor_n_pub_sum = as.integer(sum(vapply(hit, function(g) {
          value <- unname(npub_key_set[paste(name, g, sep = "\x1f")])
          if (!length(value) || is.na(value)) 0L else as.integer(value)
        }, 0L))),
        donor_species = donor,
        marker_level = marker_level, cluster_level = cluster_level
      )
    }
    scored <- list()
    for (name in names(by_type_genes)) {
      genes <- by_type_genes[[name]]
      panel_set <- panel_genes_by_type[[name]] %||% character(0)
      row <- .scored_claimant(name, genes, rec_by_type, npub_key, panel_set, NULL)
      if (!is.null(row)) scored[[length(scored) + 1L]] <- row
    }
    for (other in names(cross_tables)) {
      tab <- cross_tables[[other]]
      for (name in names(tab$by_type_genes)) {
        genes <- tab$by_type_genes[[name]]
        row <- .scored_claimant(name, genes, tab$rec_by_type, tab$npub_key, genes, other)
        if (!is.null(row)) scored[[length(scored) + 1L]] <- row
      }
    }
    if (!length(scored)) next
    rank_m <- .borrow_percentile_rank(vapply(scored, function(s) s$marker_level, 0))
    rank_c <- .borrow_percentile_rank(vapply(scored, function(s) s$cluster_level, 0))
    for (i in seq_along(scored)) {
      scored[[i]]$borrow_score <- sqrt(rank_m[i] * rank_c[i])
    }
    ord <- order(
      -vapply(scored, function(s) s$borrow_score, 0),
      -vapply(scored, function(s) s$cluster_level, 0),
      -vapply(scored, function(s) s$marker_level, 0),
      .cp(vapply(scored, function(s) s$cell_type, "")),
      method = "radix"
    )
    scored <- scored[ord]
    pending[[cluster]] <- list(
      state = state, log = log, scored = scored,
      native_best = native_best, native_best_name = native_best_name
    )
  }
  screened <- list()
  if (!is.null(screen) && length(pending)) {
    all_names <- sort(.cp(unique(unlist(lapply(pending, function(p) {
      vapply(p$scored, function(s) s$cell_type, "")
    })))), method = "radix")
    screened <- screen(all_names) %||% list()
  }
  for (cluster in names(pending)) {
    p <- pending[[cluster]]
    state <- p$state
    log <- p$log
    native_best <- p$native_best
    native_best_name <- p$native_best_name
    for (item in p$scored) {
      name <- item$cell_type
      hit <- item$anchor_genes
      donor_species <- item$donor_species
      source_frame <- if (is.null(donor_species)) resource else cross_tables[[donor_species]]$translated
      source_tissues_key <- if (is.null(donor_species)) tissues_key else cross_tables[[donor_species]]$tissues_key
      source_rec_by_type <- if (is.null(donor_species)) rec_by_type else cross_tables[[donor_species]]$rec_by_type
      item$tissue_contexts <- .borrow_n_tissues(source_tissues_key, name)
      verdict <- screened[[name]] %||% list()
      if (identical(as.character(verdict$verdict %||% ""), SCREEN_FOREIGN)) {
        item$rejected <- paste0(SCREEN_REJECTED_PREFIX, trimws(as.character(verdict$reason %||% "")))
        log$considered[[length(log$considered) + 1L]] <- item
        next
      }
      if (length(log$borrowed) >= BORROW_MAX_PER_CLUSTER) {
        item$rejected <- sprintf("cap of %d borrowed candidates reached", BORROW_MAX_PER_CLUSTER)
        log$considered[[length(log$considered) + 1L]] <- item
        next
      }
      name_value <- name
      rows <- source_frame[cell_type == name_value]
      panel_meta <- .borrow_panel_meta(rows)
      meta <- panel_meta$meta
      pos_genes <- panel_meta$positive
      neg_genes <- panel_meta$negative
      cluster_value <- cluster
      here <- ma[group == cluster_value]
      keys <- paste(cluster, here$gu, sep = "\x1f")
      stats <- list(
        pct_in = setNames(as.numeric(here$pct_in), keys),
        pct_out = setNames(as.numeric(here$pct_out), keys),
        avg_log2FC = setNames(as.numeric(here$avg_log2FC), keys),
        auc = setNames(as.numeric(here$auc), keys),
        padj = setNames(as.numeric(here$padj), keys)
      )
      native_sym <- unlist(pool$native_gene)
      extra_sym <- setNames(as.character(here$feature), here$gu)
      extra_sym <- extra_sym[!(names(extra_sym) %in% names(native_sym))]
      extra_sym <- extra_sym[!duplicated(names(extra_sym))]
      native_sym <- c(native_sym, extra_sym)
      for (field in c("pct_in", "pct_out", "avg_log2FC", "auc", "padj")) {
        native_field <- pool$de_stats[[field]]
        mine <- native_field[startsWith(names(native_field), paste0(cluster, "\x1f"))]
        stats[[field]][names(mine)] <- mine
      }
      curated <- .borrow_curated_table(name, meta)
      markers <- .panel_rows(
        name, pos_genes[pos_genes %in% measured_genes],
        neg_genes[neg_genes %in% measured_genes],
        cluster, stats, curated, native_sym, display_specificity
      )
      in_table <- toupper(vapply(markers, function(m) as.character(m$gene), ""))
      unmeasured <- sort(.cp(pos_genes[!(pos_genes %in% measured_genes)]), method = "radix")
      flat <- pos_genes[pos_genes %in% measured_genes & !(pos_genes %in% in_table)]
      if (length(flat)) {
        flat_npub <- vapply(flat, function(g) {
          as.integer(meta[[paste(name, g, sep = "\x1f")]]$n_pub)
        }, 0L)
        flat <- flat[order(-flat_npub, .cp(flat), method = "radix")]
      }
      .shown <- function(genes) {
        genes <- genes[seq_len(min(UNMEASURED_SHOWN, length(genes)))]
        as.list(unname(vapply(genes, function(g) {
          value <- unname(native_sym[g])
          if (!length(value) || is.na(value)) g else value
        }, "")))
      }
      entry <- list(
        cell_type = name,
        retrieval_rank = length(state$candidates) + 1L,
        markers = markers,
        exclusion_sources = .exclusion_sources(name, markers, sources),
        unmeasured_curated_genes = list(count = length(unmeasured), genes = .shown(unmeasured)),
        flat_or_undetected_curated_genes = list(count = length(flat), genes = .shown(flat)),
        program = list(median_in = NULL, median_out = NULL),
        borrowed_context = list(
          source_tissues = as.list(sort(.cp(unique(unlist(
            lapply(hit, function(g) panel_meta$tissues[[g]] %||% character(0))
          ))), method = "radix")),
          tissue_contexts_documenting_type = item$tissue_contexts,
          donor_species = donor_species,
          anchor_genes = as.list(hit),
          recommended_anchors = as.list(item$recommended_anchors),
          native_best = list(cell_type = native_best_name, anchors_covered = native_best)
        )
      )
      rec_here <- source_rec_by_type[[name]] %||% character(0)
      for (position in seq_along(markers)) {
        gene <- toupper(as.character(markers[[position]]$gene))
        has_meta <- !is.null(meta[[paste(name, gene, sep = "\x1f")]]$n_pub) &&
          as.integer(meta[[paste(name, gene, sep = "\x1f")]]$n_pub) != 0L
        markers[[position]]$recommended <- isTRUE(has_meta && gene %in% rec_here)
      }
      entry$markers <- markers
      own_sig <- Filter(function(m) identical(m$polarity, "positive") && .carried(m), markers)
      if (length(own_sig) < BORROW_OWN_SIGNIFICANT_MIN) {
        n_pos <- length(Filter(function(m) identical(m$polarity, "positive"), markers))
        item$rejected <- sprintf(
          "fewer than %d of its own positive markers carried here (%d of %d measured)",
          BORROW_OWN_SIGNIFICANT_MIN, length(own_sig), n_pos
        )
        log$considered[[length(log$considered) + 1L]] <- item
        next
      }
      item$own_markers_significant <- as.list(vapply(own_sig, function(m) as.character(m$gene), ""))
      item$panel_measured <- length(markers)
      log$borrowed[[length(log$borrowed) + 1L]] <- item
      state$candidates <- c(state$candidates, list(entry))
      for (m in markers) {
        gene <- toupper(as.character(m$gene))
        carriers <- pool$gene_carriers[[gene]] %||% list()
        if (!any(vapply(carriers, function(c) identical(c[[1]], name), TRUE))) {
          carriers[[length(carriers) + 1L]] <- list(
            name, as.integer(m$n_pub %||% 0L), as.character(m$tier %||% NOT_AVAILABLE), m$polarity
          )
          pool$gene_carriers[[gene]] <- carriers
        }
        if (is.null(pool$gene_clusters[[gene]])) {
          gene_value <- gene
          rows_g <- ma[gu == gene_value]
          pool$gene_clusters[[gene]] <- lapply(seq_len(nrow(rows_g)), function(i) {
            list(
              cluster_id = as.character(rows_g$group[i]),
              pct_in = .round(rows_g$pct_in[i], 1L),
              pct_out = .round(rows_g$pct_out[i], 1L)
            )
          })
        }
        if (is.null(pool$native_gene[[gene]])) {
          value <- unname(native_sym[gene])
          pool$native_gene[[gene]] <- if (!length(value) || is.na(value)) gene else value
        }
      }
    }
    pool$clusters[[cluster]] <- state
    events[[cluster]] <- log
  }
  list(pool = pool, events = events)
}
scma_borrowed_names <- function(events) {
  names_found <- character(0)
  for (log in events %||% list()) {
    for (item in log$borrowed %||% list()) {
      names_found <- c(names_found, as.character(item$cell_type))
    }
  }
  sort(.cp(unique(names_found)), method = "radix")
}
#' For each borrowed name, the donor species it actually crossed from, if any.
#'
#' A same-species (cross-tissue) borrow contributes nothing here; a name is a key only
#' if at least one cluster admitted it via `features.cross_species_borrow`, so callers
#' can tell the two borrowed shapes apart and fetch sentences accordingly. Named list:
#' name -> unique character vector of donor species.
scma_borrowed_donor_species <- function(events) {
  out <- list()
  for (log in events %||% list()) {
    for (item in log$borrowed %||% list()) {
      donor <- item$donor_species
      if (!is.null(donor) && nzchar(as.character(donor))) {
        nm <- as.character(item$cell_type)
        out[[nm]] <- union(out[[nm]], as.character(donor))
      }
    }
  }
  out
}
