#!/usr/bin/env Rscript
# =============================================================================
# annotator_pool.R -- the measured evidence one annotating agent reasons over,
# and the tools it may query, field for field and ordering for ordering.
#
# Everything here is a lookup over artifacts the pipeline already produced: the curated
# marker panels, the one-versus-rest DE table and the curated source sentences. Nothing is
# re-decided, nothing is scored and no threshold removes a candidate. The agent compares
# identities; this file only supplies facts and, afterwards, checks the numbers it quoted.
#
# Two properties of the opening packet carry the design.
#
# Each candidate's panel arrives WHOLE, ordered by curated publication support. Showing
# only the up-regulated markers hides the two facts that decide a call. A claimed identity
# whose best-corroborated markers are measured near zero is refuted by exactly those rows,
# and a cell type that dominates the dataset has its own defining markers driven to a flat
# fold change by the one-versus-rest comparison, so its cluster reads as having no
# evidence for it at all unless the whole panel is on the page.
#
# Publication support is the ordering, not a gate. It is the curated database's own record
# of which gene defines the cell type, and sorting on it puts that gene first; whether the
# cluster raises it is a separate question the measurements answer on the same row.
#
# Sourced by cluster_annotation.R; defines no side effects of its own.
# =============================================================================

A_POOL <- CFG$cluster_annotation
SOURCES_PER_GENE <- as.integer(A_POOL$sources_per_marker)
SOURCE_GENES_PER_QUERY <- as.integer(A_POOL$source_genes_per_query)
SOURCE_MAX_CHARS <- as.integer(A_POOL$source_sentence_max_chars)
UNMEASURED_SHOWN <- as.integer(A_POOL$unmeasured_genes_shown)
WARNING_GENES_SHOWN <- as.integer(A_POOL$warning_genes_shown)
## An exclusion the cluster actually detects is the one curated claim whose weight turns
## entirely on what the literature says, and it is not the kind of row an agent thinks to
## ask about. Those sentences ship with the packet; every other source is still fetched by
## name. The threshold selects which rows carry provenance and nothing else -- every
## negative marker appears in the panel either way, and the reading stays with the agent.
NEGATIVE_SOURCE_MIN_PCT_IN <- as.numeric(A_POOL$negative_source_min_pct_in)
## How many of a judged label's genes arrive with their sentences already attached. It caps
## the SENTENCES, not the panel: the panel is whole, because a truncated panel let 10% of
## the genes an answer was actually built on -- C1QB for Kupffer cell, VSIG4 for alveolar
## macrophage -- reach the judge nowhere at all, and a check cannot weigh a row it was never
## shown. Sentences are the expensive half, so the best-published ones ship and the judge
## asks for the rest.
JUDGE_SOURCES_PUSHED <- as.integer(A_POOL$judge_sources_pushed)
JUDGE_SOURCES_PER_GENE <- as.integer(A_POOL$judge_sources_per_gene)
## How many cell types one `pool_search` answer names. Enough to show whether the retrieved
## menu is missing an identity, short enough that it cannot become a second candidate pool.
POOL_SEARCH_ROWS <- as.integer(A_POOL$pool_search_rows)
## How many rows of the cluster's own ranked evidence travel with the packet. A reading
## order, not a shortlist: every row is already on some candidate's panel, and the panels
## remain the complete record.
CLUSTER_MARKER_ROWS <- as.integer(A_POOL$cluster_marker_rows)

## The names the annotator may use. `sources` is the one piece of evidence deliberately
## kept out of the opening packet: the curated sentences run to thousands of rows for one
## dataset and only a handful ever decide a cluster. The others answer questions the packet
## cannot: how a gene behaves in the OTHER clusters, which candidates claim it, and -- since
## retrieval hands over 15 of a median 67 admitted cell types and nothing downstream could
## see past that -- which cell types OUTSIDE the menu claim a gene this cluster raises.
TOOLS <- c("sources", "gene_across_clusters", "candidates_with_gene", "pool_search")

## The tools a judge may call, and deliberately fewer. `candidates_with_gene` and
## `pool_search` both name identities other than the one under test, which is exactly what a
## judge must not be given: handed a menu it could prefer a different label, and that
## disagreement says nothing about whether THIS label is carried. What it may ask for is
## more of its own label's sentences, and how a gene behaves in the other clusters.
JUDGE_TOOLS <- c("sources", "gene_across_clusters")

## Detection is carried on the DE table's own scale, percent, everywhere the agent sees it
## and everywhere the audit compares it. One scale, one place.
PCT_DECIMALS <- 1L
LOGFC_DECIMALS <- 3L
SPECIFICITY_DECIMALS <- 3L

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

## Sentences shipped with a packet, registered as delivered when a server is serving.
##
## Registering matters: the first thing an agent asks about is usually a gene it can
## already see sentences for, and a request that replied with those same sentences would
## read as a refusal to look further. Mirrors `annotator_pool._opening_records`.
.opening_records <- function(sources, label, gene, k) {
  if (is.null(sources)) {
    return(list())
  }
  records <- if (.is_source_server(sources)) {
    sources$opening(label, gene, k)
  } else {
    sources(label, gene, k)
  }
  lapply(records, function(record) {
    list(
      pmid = .identifier(record$pmid),
      pmcid = .identifier(record$pmcid),
      sentence = .clip(record$sentence, SOURCE_MAX_CHARS)
    )
  })
}

## Whether the cluster raises this marker: up-regulated one-versus-rest. One predicate,
## one meaning, used by the packet's own wording and by the post-hoc audit. It is
## deliberately not a defensibility test -- a gene can be raised and uninformative, or
## unraised and decisive -- so nothing here removes a candidate.
.is_raised <- function(marker) {
  value <- marker$avg_log2FC
  !is.null(value) && length(value) && is.finite(as.numeric(value[1])) &&
    as.numeric(value[1]) > 0
}

## ---- building ------------------------------------------------------------------
## Cluster identifiers in one deterministic order, so the packet, the record and the
## printed summary agree on it whichever of them is read first.
.sorted_cluster_ids <- function(scoring) {
  ids <- names(scoring$clusters)
  ids[order(nchar(ids), ids)]
}

## One candidate's curated panel in this cluster, best-corroborated marker first.
## Positive and negative markers are one list with a polarity column rather than two
## blocks, so a reader meets the identity's strongest curated evidence first whichever
## direction it points in. Only genes this dataset measured appear; the rest are counted
## separately, because a gene with no measurement supports nothing and refutes nothing.
## `candidate_value` / `gene_value`, not `candidate` / `gene`: inside `.(...)` data.table
## evaluates a bare name against the table's own columns first, so a local sharing a
## column's name silently becomes that whole column and the lookup turns into a
## cartesian join whose first row is then read as every marker's polarity and n_pub.
## `significant` and `specificity` carry the two quantities retrieval already decided the
## row on -- `sig_pass` admitted it, `m_g` weighed it -- neither of which `n_pub` stands
## in for. They state no conclusion; they are the run's own numbers arriving where the
## decision is made.
.panel_rows <- function(candidate_value, positive, negative, cluster, measured, curated,
                        native_of, specificity = NULL) {
  genes <- c(positive, negative)
  polarity <- c(rep("positive", length(positive)), rep("negative", length(negative)))
  keep <- !duplicated(genes)
  genes <- genes[keep]
  polarity <- polarity[keep]
  rows <- list()
  for (position in seq_along(genes)) {
    gene_value <- toupper(genes[position])
    key <- paste(cluster, gene_value, sep = "\x1f")
    if (is.na(measured$pct_in[key]) || is.null(measured$pct_in[[key]])) next
    meta <- curated[.(candidate_value, gene_value), nomatch = NULL]
    native <- unname(native_of[gene_value])
    pct_in <- round(as.numeric(measured$pct_in[[key]]), PCT_DECIMALS)
    lfc <- round(as.numeric(measured$avg_log2FC[[key]]), LOGFC_DECIMALS)
    padj <- suppressWarnings(as.numeric(measured$padj[[key]] %||% NA_real_))
    ## Specificity is derived from the positive panels alone, but this loop also walks the
    ## negative ones, so a gene that was never scored does reach here. `[[` raises on a name
    ## the table does not carry, before %||% can supply the default; every arm tolerates
    ## that miss, which is what keeps them in step.
    m_g <- suppressWarnings(as.numeric(
      if (gene_value %in% names(specificity)) specificity[[gene_value]] else NA_real_
    ))
    rows[[length(rows) + 1L]] <- list(
      gene = if (!length(native) || is.na(native)) gene_value else native,
      polarity = if (nrow(meta)) as.character(meta$marker_polarity[1]) else polarity[position],
      n_pub = as.integer(if (nrow(meta)) meta$n_pub[1] else 0L),
      tier = if (nrow(meta)) .identifier(meta$tier[1]) else NOT_AVAILABLE,
      pct_in = pct_in,
      pct_out = round(as.numeric(measured$pct_out[[key]]), PCT_DECIMALS),
      avg_log2FC = lfc,
      significant = isTRUE(!is.na(padj) && sig_pass(lfc, pct_in / 100, padj)),
      specificity = if (length(m_g) && !is.na(m_g)) round(m_g, SPECIFICITY_DECIMALS) else NULL
    )
  }
  if (length(rows)) {
    rows <- rows[order(
      -vapply(rows, function(r) r$n_pub, 0L),
      vapply(rows, function(r) toupper(r$gene), "")
    )]
  }
  rows
}

## Every cluster's candidates with their complete measured panels. `de` is the
## one-versus-rest table keyed by (group, gene_key); `native_of` maps an upper-cased gene
## key back to the symbol the dataset spells. The result is a plain nested structure so
## every arm builds the identical thing and they can be compared field for field.
## `sources` is the same callable `.run_tool` takes; without it a detected exclusion
## simply carries no provenance and everything else is unchanged.
.exclusion_sources <- function(candidate_value, markers, sources) {
  found <- list()
  for (marker in markers) {
    if (!identical(marker$polarity, "negative")) next
    pct_in <- marker$pct_in
    if (is.null(pct_in) || is.na(pct_in) || pct_in < NEGATIVE_SOURCE_MIN_PCT_IN) next
    gene <- toupper(as.character(marker$gene))
    found[[gene]] <- .opening_records(
      sources, candidate_value, gene, SOURCES_PER_GENE
    )
  }
  found
}

.build_pool <- function(scoring, de, native_of, sources = NULL) {
  keys <- paste(de$group, de$gene_key, sep = "\x1f")
  measured <- list(
    pct_in = setNames(de$pct_in, keys),
    pct_out = setNames(de$pct_out, keys),
    avg_log2FC = setNames(de$avg_log2FC, keys),
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
    gene_clusters = gene_clusters,
    gene_carriers = gene_carriers,
    native_gene = as.list(native_of)
  )
}

.candidate_names <- function(pool, cluster) {
  vapply(pool$clusters[[as.character(cluster)]]$candidates, function(e) e$cell_type, "")
}

.find_candidate <- function(pool, cluster, name) {
  for (entry in pool$clusters[[as.character(cluster)]]$candidates) {
    if (identical(entry$cell_type, as.character(name))) {
      return(entry)
    }
  }
  NULL
}

## ---- the packet the agent opens with -----------------------------------------
## Source sentences are the largest thing available and the smallest fraction of it that
## ever matters, so they are fetched by name -- except behind a detected exclusion, where
## the sentences are the evidence and arrive with the panel.
## What the agent is shown of a marker row: the run's record minus what it cannot use.
## `tier` is omitted because it carries nothing `n_pub` does not -- the resource derives
## it from the publication count by a fixed rule, with no exception anywhere in a panel,
## and the packet's own field list never described it. The row keeps it for
## `marker_evidence.csv`, where a reader asked for a tier rather than a count.
PACKET_MARKER_FIELDS <- c(
  "gene", "polarity", "n_pub", "specificity",
  "pct_in", "pct_out", "avg_log2FC", "significant"
)

.packet_marker <- function(row) {
  row[intersect(PACKET_MARKER_FIELDS, names(row))]
}

## The cluster's measured evidence keyed by gene instead of by identity.
##
## The candidate blocks answer "is this identity's curated definition met here". They
## cannot answer "which identity does this cluster's evidence belong to", because that
## question is about a gene's claimants and each block shows one claimant. Ordering by
## contrast rather than by `n_pub` is the other half: publication support is the right
## order for reading one definition and measurably the wrong one for comparing two, the
## two correlating at -0.08 within a panel across the shipped library.
##
## Nothing is computed that is not already on a panel row, and no row is a conclusion.
## Significant rows only, on the pipeline's own gate.
.cluster_marker_table <- function(pool, cluster) {
  rows <- list()
  for (entry in pool$clusters[[as.character(cluster)]]$candidates) {
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

.cluster_packet <- function(pool, cluster) {
  state <- pool$clusters[[as.character(cluster)]]
  entries <- state$candidates
  entries <- entries[order(vapply(entries, function(e) e$cell_type, ""))]
  list(
    query = c(
      pool$context,
      list(
        cluster_id = as.character(cluster),
        cells_in_cluster = as.integer(state$n_cells)
      )
    ),
    cluster_markers = .cluster_marker_table(pool, cluster),
    candidates = lapply(entries, function(entry) {
      list(
        cell_type = entry$cell_type,
        markers = lapply(entry$markers, .packet_marker),
        exclusion_sources = entry$exclusion_sources,
        unmeasured_curated_genes = entry$unmeasured_curated_genes,
        program = entry$program
      )
    })
  )
}

## ---- quality control ---------------------------------------------------------
## The identity's WHOLE curated definition as this cluster measures it.
##
## Ordered by publication support and NOT filtered to what the cluster raises, which is
## the whole difference between a check and a rubber stamp. A gene the literature ties
## most strongly to this identity, sitting near zero here, is the evidence AGAINST it;
## keep only the raised rows and every label arrives looking supported, because the only
## rows left are the ones that supported it. Each row carries its own measurement, so
## which of them the cluster actually raises is visible without pre-selecting them.
##
## Not cut to a length either, for the same reason one step further on. Cut at the fourteen
## best-published rows, 61 of the genes the annotator's own answers rested on -- 10.3% of
## them, over 29 of 104 clusters -- were absent from the judge's page entirely: C1QB and
## CSF1R under Kupffer cell, VSIG4 under alveolar macrophage, PKHD1 under cholangiocyte. A
## judge cannot weigh a row it was never shown, and it had no way to ask for one. Rows are
## cheap; it is the sentences that are not, and those are what it asks for.
##
## The detected exclusions come too: a curated negative marker the cluster expresses is
## the other half of the case against.
.judge_panel <- function(entry) {
  positives <- Filter(function(row) {
    identical(row$polarity, "positive")
  }, entry$markers)
  detected <- Filter(function(row) {
    pct_in <- suppressWarnings(as.numeric(row$pct_in %||% NA_real_))
    identical(row$polarity, "negative") &&
      is.finite(pct_in) && pct_in >= NEGATIVE_SOURCE_MIN_PCT_IN
  }, entry$markers)
  lapply(c(positives, detected), .packet_marker)
}

## The sentences that ship with the panel: the best-published, and every exclusion.
##
## A whole panel's sentences run to hundreds of rows and the judge reads a handful, so the
## packet opens with the ones a verdict usually turns on and the `sources` tool answers for
## any other gene on the panel. Every DETECTED EXCLUSION is pushed whatever its rank -- an
## exclusion's whole weight is what its sentence does, so a judge holding the row without
## the sentence holds the half that cannot be read.
.judge_sources <- function(label, panel, sources) {
  found <- list()
  if (is.null(sources)) {
    return(found)
  }
  positives <- 0L
  for (row in panel) {
    gene <- toupper(as.character(row$gene))
    if (identical(row$polarity, "positive")) {
      if (positives >= JUDGE_SOURCES_PUSHED) next
      positives <- positives + 1L
    }
    if (!is.null(found[[gene]])) next
    records <- .opening_records(sources, label, gene, JUDGE_SOURCES_PER_GENE)
    if (length(records)) found[[gene]] <- records
  }
  found
}

## The evidence one judgment turns on, and deliberately nothing else.
##
## A judge handed the whole pool would be a second annotator: it could prefer a different
## candidate, and that disagreement says nothing about whether THIS label is carried by
## the measurements. So the packet narrows to the label under test -- the rows the cluster
## raises for it, the exclusions it detects against it, and the curated sentences behind
## those genes.
##
## Three things are absent by construction. The annotator's rationale, confidence and
## candidate ranking are the opinion under test, not evidence for it. The retrieval order
## would say which label the pipeline preferred. And no reference or author label exists
## anywhere in this package to leak in the first place.
##
## With `parent`, the question becomes whether the finer name is separable rather than
## whether it is plausible, so the parent's panel travels with it: a refinement whose
## raised markers all belong to the parent is the parent read at a narrower name.
.judge_packet <- function(pool, cluster, label, sources = NULL, parent = "") {
  entry <- .find_candidate(pool, cluster, label)
  if (is.null(entry)) {
    return(list())
  }
  state <- pool$clusters[[as.character(cluster)]]
  panel <- .judge_panel(entry)
  packet <- list(
    query = c(
      pool$context,
      list(
        cluster_id = as.character(cluster),
        cells_in_cluster = as.integer(state$n_cells)
      )
    ),
    label_under_test = as.character(label),
    panel = panel,
    unmeasured_curated_genes = entry$unmeasured_curated_genes,
    sources = .judge_sources(as.character(label), panel, sources)
  )
  parent_entry <- if (nzchar(parent)) .find_candidate(pool, cluster, parent) else NULL
  if (!is.null(parent_entry)) {
    packet$parent_label <- as.character(parent)
    packet$parent_panel <- .judge_panel(parent_entry)
    packet$raised_and_absent_from_parent_panel <- .subtype_exclusive_markers(
      pool, cluster, as.character(parent), as.character(label)
    )
  }
  packet
}

## ---- tools -------------------------------------------------------------------
## One tool call, answered from the artifacts. Read-only, and short by construction. An
## unknown tool or an unusable argument is answered with an error observation rather than
## an exception: the agent is mid-conversation, and telling it what went wrong lets it
## correct the call, where stopping would lose the whole cluster.
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
  if (identical(tool, "pool_search")) {
    return(.tool_pool_search(pool, cluster, args))
  }
  list(
    tool = as.character(tool)[1],
    error = sprintf("unknown tool; available: %s", paste(TOOLS, collapse = ", "))
  )
}

## One tool call from a judge, answerable only about the labels it is judging.
##
## The same read-only answers the annotator gets, through a gate that keeps the judge from
## becoming a second annotator: `sources` only for a label under test, and no tool at all
## that names another identity. Without this a judge could ask for a rival's sentences and
## start choosing between them, and its verdict would no longer be about whether THIS label
## is carried.
.run_judge_tool <- function(pool, cluster, labels, tool, args, sources = NULL) {
  if (!is.list(args)) args <- list()
  if (!(as.character(tool)[1] %in% JUDGE_TOOLS)) {
    return(list(
      tool = as.character(tool)[1],
      error = sprintf("unknown tool; available: %s", paste(JUDGE_TOOLS, collapse = ", "))
    ))
  }
  if (identical(tool, "gene_across_clusters")) {
    return(.tool_gene_across_clusters(pool, args))
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

## Whether `sources` is a streaming server rather than a plain fetch callable. Both are
## accepted because the pool is built once for every cluster, before any per-cluster server
## exists, while a conversation is served by its own. Mirrors `hasattr(sources, "take")`.
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
      ## Only the genes with something still behind them are reported, so the note stays
      ## short and says something when it appears.
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
    found[[key]] <- lapply(records, function(record) {
      list(
        pmid = .identifier(record$pmid),
        pmcid = .identifier(record$pmcid),
        sentence = .clip(record$sentence, SOURCE_MAX_CHARS)
      )
    })
  }
  out <- list(
    tool = "sources", candidate = candidate, sources = found, truncated = truncated
  )
  if (length(left)) out$not_yet_shown <- left
  out
}

## Tell a cluster's server which sentences its opening packet already shows.
##
## `exclusion_sources` are resolved once, when the pool is built for every cluster at once,
## so the per-cluster server cannot have served them and would otherwise offer them again as
## though they were new. The packet is the first thing the agent reads; a request answered
## with what it just read looks like a refusal to go further. Mirrors
## `annotator_pool.register_packet_sources`.
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

## The curated cell types claiming this gene that are NOT among the candidates.
##
## The one way anything downstream can see past retrieval. Retrieval admits a median of 67
## cell types per cluster and hands over 15, and the agent picked outside the top three in
## 45% of clusters and from the last four slots in 11% -- an answer rate that has not decayed
## at the boundary, so the rows beyond it are not empty, only invisible. Reading them is what
## lets a missing menu be reported instead of silently absorbed.
##
## What comes back is a name and a publication count, deliberately not a panel: it is
## evidence that the menu may be short, not a sixteenth candidate. Nothing here can be
## answered with -- `selected` is still checked against the supplied candidates -- and the
## agent reports it in `unlisted_identity`.
.tool_pool_search <- function(pool, cluster, args) {
  gene <- toupper(trimws(as.character(args$gene %||% "")[1]))
  here <- .candidate_names(pool, cluster)
  rows <- Filter(
    function(row) !(row$cell_type %in% here),
    pool$gene_carriers[[gene]] %||% list()
  )
  total <- length(rows)
  if (total) {
    rows <- rows[order(
      -vapply(rows, function(r) as.numeric(r$n_pub), 0),
      vapply(rows, function(r) r$cell_type, "")
    )]
    rows <- head(rows, POOL_SEARCH_ROWS)
  }
  list(
    tool = "pool_search",
    gene = pool$native_gene[[gene]] %||% gene,
    not_in_candidates = rows,
    total = total
  )
}

## ---- warnings ----------------------------------------------------------------
## The raised markers that belong to the subtype's definition and not the parent's.
##
## A finer name is only a finer CALL when the cluster raises a marker the parent's own
## curated definition does not contain. A gene curated for both identities is the parent
## program read at a narrower name, however wide its contrast, so the comparison is
## against the parent's WHOLE curated panel rather than the part of it this cluster
## raises.
.subtype_exclusive_markers <- function(pool, cluster, selected, subtype) {
  parent <- .find_candidate(pool, cluster, selected)
  finer <- .find_candidate(pool, cluster, subtype)
  if (is.null(parent) || is.null(finer)) {
    return(character(0))
  }
  positive <- function(entry) {
    Filter(function(marker) identical(marker$polarity, "positive"), entry$markers)
  }
  gene_of <- function(markers) {
    toupper(vapply(markers, function(marker) as.character(marker$gene), ""))
  }
  raised <- Filter(.is_raised, positive(finer))
  sort(unique(setdiff(gene_of(raised), gene_of(positive(parent)))))
}


## The candidates the answer put in one group with `name`, `name` included.
.group_of <- function(identity_groups, name) {
  group <- as.character(name)
  for (entry in identity_groups %||% list()) {
    members <- as.character(unlist(entry$candidates %||% list()))
    if (as.character(name) %in% members) {
      group <- union(group, members)
    }
  }
  group
}


.detected_exclusions <- function(pool, cluster, group) {
  found <- list()
  for (name in sort(group)) {
    entry <- .find_candidate(pool, cluster, name)
    if (is.null(entry)) next
    for (marker in entry$markers) {
      if (!identical(marker$polarity, "negative")) next
      pct_in <- suppressWarnings(as.numeric(marker$pct_in %||% NA_real_))
      if (!is.finite(pct_in) || pct_in < NEGATIVE_SOURCE_MIN_PCT_IN) next
      gene <- as.character(marker$gene)
      key <- toupper(gene)
      previous <- found[[key]]
      if (is.null(previous) || pct_in > previous$pct_in) {
        found[[key]] <- list(
          gene = gene,
          curated_negative_for = name,
          pct_in = pct_in,
          pct_out = marker$pct_out,
          n_pub = as.integer(marker$n_pub)
        )
      }
    }
  }
  if (!length(found)) {
    return(list())
  }
  unname(found[sort(names(found))])
}


## Detected exclusions bound to the identity the answer selected AND to each it reported
## as merely co-occurring, one block per identity.
##
## `.claim_warnings` performs the mirror of this check on the positive half of the panel,
## and both read the same packet the same way. The difference is when: that one runs after
## the answer is fixed and is filed for a human reader, while an exclusion on the identity
## the agent is about to assert is still answerable.
##
## Both sides are returned because one side alone cannot be read. A block whose
## `detected_exclusions` is empty says so explicitly rather than by absence: the fact that
## decides a cluster like a phagocytic macrophage called hepatocyte is not that the pick
## carries exclusions, it is that the identity it demoted carries none. Ambient RNA,
## engulfed cargo and doublets all raise a positive marker, and none of them can make a
## cell express the gene whose absence defines it.
##
## Facts only: the gene, the candidate the resource curates it under, and the measurement.
## Nothing here says a detected exclusion refutes anything, or ranks one against another;
## that reading stays with the agent, answered against the `exclusion_sources` it holds.
.binding_exclusions <- function(pool, cluster, selected, identity_groups,
                                co_occurring = NULL) {
  selected <- as.character(selected)
  blocks <- list(list(
    identity = selected,
    role = "selected",
    detected_exclusions = .detected_exclusions(
      pool, cluster, .group_of(identity_groups, selected)
    )
  ))
  seen <- selected
  for (name in as.character(unlist(co_occurring %||% list()))) {
    if (!nzchar(name) || name %in% seen) next
    seen <- c(seen, name)
    blocks[[length(blocks) + 1L]] <- list(
      identity = name,
      role = "co_occurring",
      detected_exclusions = .detected_exclusions(
        pool, cluster, .group_of(identity_groups, name)
      )
    )
  }
  carries <- vapply(blocks, function(block) length(block$detected_exclusions) > 0L, TRUE)
  if (any(carries)) blocks else list()
}


## One short line per claimed identity whose curated markers are not all raised.
##
## The warning is computed over the identity's WHOLE measured positive panel, not over
## the gene the agent chose to quote. That is the point of it: a claim can be written on
## a single incidental marker with a wide contrast while the genes the database says
## define that cell sit flat, and a check that only looked at the quoted gene would see
## nothing wrong. It removes nothing and renames nothing; publication support only orders
## which unraised genes are worth naming, and being unraised is the whole trigger.
##
## Each line travels WITH its identity (a list(name, line) pair) rather than the
## identity being parsed back out of the line: 212 curated cell-type names carry a
## colon of their own (`BCR::ABL1 positive primitive cell`), so cutting the line at
## its first colon filed the warning against an identity that does not exist.
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
    ## One fixed decimal, the percentage rounded to PCT_DECIMALS, so every arm produces
    ## the same warning string byte for byte.
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
