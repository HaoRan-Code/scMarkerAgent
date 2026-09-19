#!/usr/bin/env Rscript
STAGE_MODEL <- as.character(A_POOL$stage_model %||% A_POOL$review_model)
STAGE_EFFORT <- as.character(A_POOL$stage_reasoning_effort %||% A_POOL$review_reasoning_effort)
AUDIT_PANEL_ROWS <- as.integer(A_POOL$audit_panel_rows %||% 15L)
AUDIT_SENTENCES_PER_GENE <- as.integer(A_POOL$audit_sentences_per_gene %||% 2L)
AUDIT_RIVALS_MAX <- as.integer(A_POOL$audit_rivals_max %||% 4L)
AUDIT_RIVAL_PANEL_ROWS <- as.integer(A_POOL$audit_rival_panel_rows %||% 10L)
AUDIT_GENOME_ROWS <- as.integer(A_POOL$audit_genome_rows %||% 15L)
CONSISTENCY_GROUP_MAX <- as.integer(A_POOL$consistency_group_max %||% 8L)
OWN_PROGRAM_VALUES <- c("dominant", "minority", "absent", "single_study_only")
CONSISTENCY_VALUES <- c("consistent", "weaker_member", "conflict")
LINEAGE_SPECIFICITY_VALUES <- c("one_lineage", "several_lineages", "uncertain")
RELATION_VALUES <- c("coarser_parent", "finer_child", "sibling_or_pole", "different_lineage",
                     "state_label")
OCCUPANCY_VALUES <- c("identity_only", "rival_only", "both_cluster_level", "both_minority",
                      "neither")
POLE_RELATIONS <- c("sibling_or_pole", "different_lineage")
CONTESTED_OCCUPANCY <- c("rival_only", "both_cluster_level")
.stage_call_many <- function(clusters, prompts, trace_ids, turn, validate, retry_note,
                             api_url, api_key, schema_retries) {
  out <- setNames(vector("list", length(clusters)), clusters)
  active <- clusters
  attempt <- setNames(rep(0L, length(clusters)), clusters)
  while (length(active)) {
    responses <- scma_cached_call_llm_many(
      lapply(active, function(cl) prompts[[cl]]), api_url, api_key,
      reasoning_effort = STAGE_EFFORT, model = STAGE_MODEL,
      trace_ids = lapply(active, function(cl) trace_ids[[cl]]),
      turn_indexes = lapply(active, function(cl) as.integer(turn))
    )
    still <- character(0)
    for (i in seq_along(active)) {
      cl <- active[i]
      content <- responses[[i]][[1]]
      parsed <- if (is.null(content)) NULL else scma_parse_json(content)
      problems <- if (is.null(parsed)) "return one JSON object" else validate(cl, parsed)
      if (!length(problems)) {
        out[[cl]] <- parsed
        next
      }
      if (attempt[[cl]] >= schema_retries) {
        out[[cl]] <- NULL
        next
      }
      attempt[[cl]] <- attempt[[cl]] + 1L
      prompts[[cl]] <- paste0(
        prompts[[cl]], sprintf(
          "\n\n# Retry %d: %s -- %s.\n", attempt[[cl]], retry_note,
          paste(head(problems, 4L), collapse = "; ")
        )
      )
      still <- c(still, cl)
    }
    active <- still
  }
  out
}
compact_definers <- compact_definer_rows
program_packet <- function(pool, cluster, tier) {
  state <- pool$clusters[[as.character(cluster)]]
  c(
    list(query = c(
      pool$context,
      list(cluster_id = as.character(cluster), cells_in_cluster = as.integer(state$n_cells))
    )),
    genome_pages(pool, cluster, tier)
  )
}
.page_genes <- function(packet) {
  genes <- character(0)
  for (key in c("cluster_top_enriched", "cluster_top_detection_gap", "cluster_top_depleted")) {
    genes <- c(genes, toupper(vapply(packet[[key]] %||% list(), function(r) as.character(r$gene), "")))
  }
  unique(genes)
}
.mentions_candidate <- function(text, names) {
  low <- tolower(as.character(text %||% "")[1])
  for (name in names) {
    n <- tolower(trimws(name))
    if (nchar(n) >= 6L && grepl(n, low, fixed = TRUE)) return(name)
  }
  NULL
}
validate_program <- function(value, page_genes, names) {
  if (!is.list(value)) return("return one JSON object")
  problems <- character(0)
  programs <- value$programs
  if (!is.list(programs)) {
    problems <- c(problems, "programs must be a list")
    programs <- list()
  }
  for (block in programs) {
    if (!is.list(block)) {
      problems <- c(problems, "each program must be an object")
      next
    }
    for (g in as.character(unlist(block$genes %||% list()))) {
      if (!(toupper(g) %in% page_genes)) {
        problems <- c(problems, sprintf("gene %s is not on the genome-wide pages", g))
        break
      }
    }
    hit <- .mentions_candidate(block$`function` %||% "", names)
    if (!is.null(hit)) {
      problems <- c(problems, sprintf("'function' names a candidate cell type (%s); describe the function", hit))
    }
    if (!(as.character(block$lineage_specificity %||% "") %in% LINEAGE_SPECIFICITY_VALUES)) {
      problems <- c(problems, sprintf(
        "each program needs lineage_specificity, one of %s", .py_list(LINEAGE_SPECIFICITY_VALUES)
      ))
    }
    hit <- .mentions_candidate(block$lineage_note %||% "", names)
    if (!is.null(hit)) {
      problems <- c(problems, sprintf("'lineage_note' names a candidate cell type (%s)", hit))
    }
  }
  hit <- .mentions_candidate(value$mixture_note %||% "", names)
  if (!is.null(hit)) {
    problems <- c(problems, sprintf("mixture_note names a candidate cell type (%s)", hit))
  }
  if (!is.numeric(suppressWarnings(as.integer(value$lineage_programs_count)))) {
    problems <- c(problems, "lineage_programs_count must be an integer")
  }
  problems
}
program_reading_many <- function(pool, clusters, tiers, template, headers, api_key,
                                 api_url, trace_ids, turn, schema_retries) {
  prompts <- setNames(lapply(clusters, function(cl) {
    packet <- program_packet(pool, cl, tiers[[cl]])
    paste0(headers[[cl]], .render(template, packet))
  }), clusters)
  page_genes_of <- setNames(lapply(clusters, function(cl) .page_genes(program_packet(pool, cl, tiers[[cl]]))), clusters)
  names_of <- setNames(lapply(clusters, function(cl) .candidate_names(pool, cl)), clusters)
  .stage_call_many(
    clusters, prompts, trace_ids, turn,
    function(cl, v) validate_program(v, page_genes_of[[cl]], names_of[[cl]]),
    "the reading was rejected", api_url, api_key, schema_retries
  )
}
.sentences_for <- function(sources, identity, genes) {
  out <- list()
  if (is.null(sources)) return(out)
  for (gene in genes) {
    key <- toupper(as.character(gene))
    if (!is.null(out[[key]])) next
    records <- .opening_records(sources, identity, key, AUDIT_SENTENCES_PER_GENE)
    if (length(records)) out[[key]] <- lapply(records, .sentence_record)
  }
  out
}
.rival_names <- function(pool, cluster, identity, claimed, tier) {
  rivals <- setdiff(claimed, identity)
  genome <- (pool$genome_de %||% list())[[as.character(cluster)]] %||% list()
  on_page <- .candidate_names(pool, cluster, tier)
  votes <- list()
  rows <- (genome$top_enriched %||% list())[seq_len(min(AUDIT_GENOME_ROWS, length(genome$top_enriched %||% list())))]
  for (row in rows) {
    for (c in row$claimed_by_species %||% list()) {
      name <- as.character(c$cell_type)
      if (name %in% on_page && name != identity && !(name %in% rivals)) {
        votes[[name]] <- (votes[[name]] %||% 0L) + 1L
      }
    }
  }
  if (length(votes)) {
    vote_names <- names(votes)
    vote_counts <- vapply(votes, function(v) as.integer(v), 0L)
    ord <- order(-vote_counts, .cp(vote_names), method = "radix")
    for (name in vote_names[ord]) {
      if (length(rivals) >= AUDIT_RIVALS_MAX) break
      if (!(name %in% rivals)) rivals <- c(rivals, name)
    }
  }
  rivals[seq_len(min(AUDIT_RIVALS_MAX, length(rivals)))]
}
audit_packet <- function(pool, cluster, identity, role, claimed, tier, sources, program) {
  entry <- .find_candidate(pool, cluster, identity)
  if (is.null(entry)) return(NULL)
  state <- pool$clusters[[as.character(cluster)]]
  definers <- entry$definers %||% list()
  panel <- .review_panel(entry)
  positives <- Filter(function(r) identical(r$polarity, "positive"), panel)
  positives <- positives[seq_len(min(AUDIT_PANEL_ROWS, length(positives)))]
  negatives <- Filter(function(r) identical(r$polarity, "negative"), panel)
  species_sources <- pool$species_sources %||% sources
  sentences <- .sentences_for(species_sources, identity, vapply(definers, function(d) as.character(d$gene), ""))
  rival_names <- .rival_names(pool, cluster, identity, claimed, tier)
  rivals <- list()
  for (name in rival_names) {
    rentry <- .find_candidate(pool, cluster, name)
    if (is.null(rentry)) next
    rdefiners <- rentry$definers %||% list()
    rpanel <- Filter(function(r) identical(r$polarity, "positive"), .review_panel(rentry))
    rival <- list(
      cell_type = name, claimed_here = name %in% claimed,
      species_definers = as.list(compact_definers(rdefiners)),
      panel = as.list(compact_panel(rpanel[seq_len(min(AUDIT_RIVAL_PANEL_ROWS, length(rpanel)))])),
      sentences = .sentences_for(species_sources, name, vapply(rdefiners, function(d) as.character(d$gene), "")),
      versus_identity = pair_partition(entry, rentry)
    )
    if (!is.null(rentry$borrowed_context)) rival$borrowed_context <- rentry$borrowed_context
    rivals[[length(rivals) + 1L]] <- rival
  }
  genome <- genome_pages(pool, cluster, tier)
  packet <- list(
    query = c(
      pool$context,
      list(cluster_id = as.character(cluster), cells_in_cluster = as.integer(state$n_cells))
    ),
    identity = identity, role = role, row_format = MARKER_ROW_FORMAT,
    species_definers = definers,
    own_vs_shared = lineage_partition(entry),
    panel = as.list(compact_panel(c(positives, negatives))),
    sentences = sentences,
    detected_exclusions = as.list(compact_panel(.detected_exclusion_markers(entry))),
    rivals = rivals,
    program_reading = program %||% list(),
    program_genes_that_are_not_identity = program_non_identity_genes(program),
    cluster_top_enriched = (genome$cluster_top_enriched %||% list())[
      seq_len(min(AUDIT_GENOME_ROWS, length(genome$cluster_top_enriched %||% list())))
    ]
  )
  if (!is.null(entry$borrowed_context)) packet$borrowed_context <- entry$borrowed_context
  packet
}
program_non_identity_genes <- function(program) {
  out <- list(state = list(), several_lineages = list())
  if (!is.list(program)) return(out)
  seen <- character(0)
  for (block in program$state_programs %||% list()) {
    for (g in as.character(unlist(block$genes %||% list()))) {
      if (!(toupper(g) %in% seen)) {
        seen <- c(seen, toupper(g))
        out$state[[length(out$state) + 1L]] <- g
      }
    }
  }
  for (block in program$programs %||% list()) {
    if (!identical(as.character(block$lineage_specificity %||% ""), "several_lineages")) next
    for (g in as.character(unlist(block$genes %||% list()))) {
      if (!(toupper(g) %in% seen)) {
        seen <- c(seen, toupper(g))
        out$several_lineages[[length(out$several_lineages) + 1L]] <- g
      }
    }
  }
  out
}
.GENE_LIST_KEYS <- c("genes", "state", "several_lineages", "background_or_shared",
                    "unplaced_top_genes")
.COMPACT_FIRST_RE <- "^(\\S+) n(?:[0-9]+|-)"
.genes_in <- function(obj, key = NULL) {
  out <- character(0)
  if (is.list(obj)) {
    is_named <- !is.null(names(obj)) && any(nzchar(names(obj)))
    gene <- obj$gene
    if (is.character(gene) && length(gene) == 1L) out <- c(out, toupper(gene))
    if (identical(key, "sentences") || identical(key, "sources") || identical(key, "exclusion_sources")) {
      out <- c(out, toupper(names(obj)))
    }
    if (is_named) {
      for (k in names(obj)) out <- c(out, .genes_in(obj[[k]], k))
    } else {
      for (item in obj) {
        if (is.character(item) && length(item) == 1L && !is.null(key) && key %in% .GENE_LIST_KEYS &&
            !grepl(" ", item, fixed = TRUE)) {
          out <- c(out, toupper(item))
        } else {
          out <- c(out, .genes_in(item, key))
        }
      }
    }
    return(out)
  }
  if (is.character(obj) && length(obj) == 1L) {
    m <- regmatches(obj, regexpr(.COMPACT_FIRST_RE, obj, perl = TRUE))
    if (length(m) && nzchar(m)) {
      g <- sub(" n.*$", "", m, perl = TRUE)
      out <- c(out, toupper(g))
    }
  }
  out
}
packet_genes <- function(packet) {
  unique(.genes_in(packet))
}
identity_side_genes <- function(packet) {
  identity <- as.character(packet$identity %||% "")
  outer <- packet[setdiff(names(packet), "rivals")]
  out <- unique(.genes_in(outer))
  for (rival in packet$rivals %||% list()) {
    versus <- rival$versus_identity %||% list()
    out <- c(out, .genes_in(versus[[paste0("only_", identity)]] %||% list(), "genes"))
    out <- c(out, .genes_in(versus$both %||% list(), "genes"))
  }
  unique(out)
}
validate_audit <- function(value, packet) {
  if (!is.list(value)) return("return one JSON object")
  problems <- character(0)
  if (!(as.character(value$own_program %||% "") %in% OWN_PROGRAM_VALUES)) {
    problems <- c(problems, sprintf("own_program must be one of %s", .py_list(OWN_PROGRAM_VALUES)))
  }
  known <- identity_side_genes(packet)
  rival_known <- packet_genes(packet)
  for (key in c("definers_present", "definers_absent")) {
    genes <- value[[key]]
    if (!is.list(genes) && !is.character(genes)) {
      problems <- c(problems, sprintf("%s must be a list", key))
      next
    }
    genes <- as.character(unlist(genes %||% list()))
    bad <- genes[!(toupper(genes) %in% known)]
    if (length(bad)) {
      problems <- c(problems, sprintf("%s names genes not in this packet: %s", key, .py_list(head(bad, 4L))))
    }
  }
  rival_names <- vapply(packet$rivals %||% list(), function(r) as.character(r$cell_type), "")
  read <- value$rivals_read
  if (!is.list(read)) {
    problems <- c(problems, "rivals_read must be a list with one entry per rival")
    read <- list()
  }
  seen <- character(0)
  for (entry in read) {
    if (!is.list(entry)) {
      problems <- c(problems, "each rivals_read entry must be an object")
      next
    }
    name <- as.character(entry$rival %||% "")
    if (!(name %in% rival_names)) {
      problems <- c(problems, sprintf("rivals_read names '%s', which is not a rival in this packet", name))
      next
    }
    seen <- c(seen, name)
    if (!(as.character(entry$relation %||% "") %in% RELATION_VALUES)) {
      problems <- c(problems, sprintf("%s: relation must be one of %s", name, .py_list(RELATION_VALUES)))
    }
    if (!(as.character(entry$whose_program_occupies %||% "") %in% OCCUPANCY_VALUES)) {
      problems <- c(problems, sprintf("%s: whose_program_occupies must be one of %s", name, .py_list(OCCUPANCY_VALUES)))
    }
    for (key in c("separating_genes", "rival_separating_genes")) {
      genes <- entry[[key]]
      if (!is.list(genes) && !is.character(genes)) {
        problems <- c(problems, sprintf("%s: %s must be a list", name, key))
        next
      }
      genes <- as.character(unlist(genes %||% list()))
      bad <- genes[!(toupper(genes) %in% rival_known)]
      if (length(bad)) {
        problems <- c(problems, sprintf("%s: %s names genes not in this packet: %s", name, key, .py_list(head(bad, 4L))))
      }
    }
    if (!nzchar(trimws(as.character(entry$note %||% "")))) {
      problems <- c(problems, sprintf("%s: note must quote the measured percentages of both sides", name))
    }
  }
  missing <- setdiff(rival_names, seen)
  if (length(missing)) {
    problems <- c(problems, sprintf("every rival in this packet needs one rivals_read entry; missing %s", .py_list(head(missing, 4L))))
  }
  if (!nzchar(trimws(as.character(value$reason %||% "")))) {
    problems <- c(problems, "reason must quote the measured percentages")
  }
  problems
}
audit_many <- function(pool, jobs, sources, programs, template, headers, api_key, api_url,
                       trace_ids, base_turn, schema_retries) {
  if (!length(jobs)) return(list())
  keys <- vapply(jobs, function(j) paste(j$cluster, j$identity, sep = "\x1f"), "")
  packets <- setNames(lapply(jobs, function(j) {
    audit_packet(pool, j$cluster, j$identity, j$role, j$claimed, j$tier, sources, programs[[j$cluster]])
  }), keys)
  prompts <- setNames(lapply(jobs, function(j) {
    key <- paste(j$cluster, j$identity, sep = "\x1f")
    paste0(headers[[j$cluster]], .render(template, packets[[key]]))
  }), keys)
  turn_of <- setNames(lapply(seq_along(jobs), function(i) as.integer(base_turn[[jobs[[i]]$cluster]] * 100 + i)), keys)
  trace_of <- setNames(lapply(jobs, function(j) trace_ids[[j$cluster]]), keys)
  out <- setNames(vector("list", length(keys)), keys)
  active <- keys
  attempt <- setNames(rep(0L, length(keys)), keys)
  while (length(active)) {
    responses <- scma_cached_call_llm_many(
      lapply(active, function(k) prompts[[k]]), api_url, api_key,
      reasoning_effort = STAGE_EFFORT, model = STAGE_MODEL,
      trace_ids = lapply(active, function(k) trace_of[[k]]),
      turn_indexes = lapply(active, function(k) turn_of[[k]])
    )
    still <- character(0)
    for (i in seq_along(active)) {
      k <- active[i]
      content <- responses[[i]][[1]]
      parsed <- if (is.null(content)) NULL else scma_parse_json(content)
      problems <- if (is.null(parsed)) "return one JSON object" else validate_audit(parsed, packets[[k]])
      if (!length(problems)) {
        parsed$identity <- packets[[k]]$identity
        parsed$role <- packets[[k]]$role
        parsed$species_definers <- as.list(compact_definers(packets[[k]]$species_definers))
        out[[k]] <- parsed
        next
      }
      if (attempt[[k]] >= schema_retries) {
        out[[k]] <- NULL
        next
      }
      attempt[[k]] <- attempt[[k]] + 1L
      prompts[[k]] <- paste0(
        prompts[[k]], sprintf(
          "\n\n# Retry %d: the audit was rejected -- %s.\n", attempt[[k]],
          paste(head(problems, 4L), collapse = "; ")
        )
      )
      still <- c(still, k)
    }
    active <- still
  }
  out
}
cached_audit <- function(cache, tier, identity, need_rivals) {
  if (is.null(cache)) return(NULL)
  key <- paste(as.integer(tier), identity, sep = "\x1f")
  if (!exists(key, envir = cache, inherits = FALSE)) return(NULL)
  audit <- get(key, envir = cache, inherits = FALSE)
  read <- vapply(rivals_read(audit), function(e) as.character(e$rival %||% ""), "")
  need <- setdiff(need_rivals, identity)
  if (any(!(need %in% read))) return(NULL)
  audit
}
rivals_read <- function(audit) {
  Filter(is.list, audit$rivals_read %||% list())
}
contested_rival_names <- function(final, audits, claimed_names) {
  selected <- as.character(final$selected %||% "")
  claimed <- unique(claimed_names)
  out <- character(0)
  for (entry in rivals_read(audits[[selected]])) {
    name <- as.character(entry$rival %||% "")
    if (!nzchar(name) || name %in% claimed || identical(name, selected)) next
    if (as.character(entry$relation %||% "") %in% POLE_RELATIONS && !(name %in% out)) {
      out <- c(out, name)
    }
  }
  out
}
.pair_reading <- function(audit, other) {
  for (entry in rivals_read(audit)) {
    if (identical(as.character(entry$rival %||% ""), other)) return(entry)
  }
  NULL
}
delivery_contested <- function(final, audits, contested_names) {
  selected <- as.character(final$selected %||% "")
  audit <- audits[[selected]]
  if (is.null(audit)) return(character(0))
  reasons <- character(0)
  own <- as.character(audit$own_program %||% "")
  if (nzchar(own) && !identical(own, "dominant")) {
    absent <- as.character(unlist(audit$definers_absent %||% list()))
    reasons <- c(reasons, sprintf(
      "the definer audit of '%s' reads its own defining program as %s here%s",
      selected, own,
      if (length(absent)) sprintf(" (%s)", substr(paste(absent, collapse = ", "), 1, 160)) else ""
    ))
  }
  for (name in contested_names) {
    entry <- .pair_reading(audit, name)
    if (is.null(entry)) next
    relation <- as.character(entry$relation %||% "")
    occupancy <- as.character(entry$whose_program_occupies %||% "")
    if (identical(occupancy, "rival_only")) {
      reasons <- c(reasons, sprintf(
        "its own audit reads the pair '%s' / '%s' (%s) as %s here: %s",
        selected, name, relation, occupancy, substr(as.character(entry$note %||% ""), 1, 200)
      ))
      next
    }
    if (identical(occupancy, "neither")) {
      reasons <- c(reasons, sprintf(
        "its own audit reads the pair '%s' / '%s' (%s) as neither: no gene in the packet tells them apart, so the delivered pole is not established against '%s' and the level the rows reach is their common parent where this page holds it: %s",
        selected, name, relation, name, substr(as.character(entry$note %||% ""), 1, 200)
      ))
      next
    }
  }
  reasons
}
POPULATION_VERDICTS <- c("established_second_population", "shared_program_of_same_population", "neither_separates", "not_present")
ESTABLISHED <- "established_second_population"
NEITHER_SEPARATES <- "neither_separates"
population_names <- function(packet) {
  delivered <- packet$delivered %||% list()
  out <- character(0)
  for (name in as.character(unlist(delivered$co_occurring_identities %||% list()))) {
    if (nzchar(name) && !(name %in% out)) out <- c(out, name)
  }
  for (block in packet$contested_identities %||% list()) {
    name <- as.character(block$cell_type %||% "")
    if (nzchar(name) && !(name %in% out)) out <- c(out, name)
  }
  out
}
population_conflicts <- function(final, review) {
  if (is.null(review)) return(character(0))
  claimed <- Filter(nzchar, as.character(unlist(final$co_occurring_identities %||% list())))
  out <- character(0)
  for (entry in review$population_verdicts %||% list()) {
    if (!is.list(entry)) next
    name <- as.character(entry$identity %||% "")
    verdict <- as.character(entry$verdict %||% "")
    why <- substr(as.character(entry$reason %||% ""), 1, 220)
    if (name %in% claimed && !identical(verdict, ESTABLISHED)) {
      out <- c(out, sprintf(
        "'%s' is reported as a co-occurring identity, and the review reads it as %s: %s",
        name, verdict, why
      ))
    } else if (!(name %in% claimed) && identical(verdict, ESTABLISHED)) {
      out <- c(out, sprintf(
        "'%s' is not reported, and the review reads it as an established second population here: %s",
        name, why
      ))
    }
  }
  out
}
CONF_ORDER <- c(high = 3L, medium = 2L, low = 1L)
confidence_cap <- function(final, audits) {
  selected <- as.character(final$selected %||% "")
  audit <- audits[[selected]]
  if (is.null(audit)) return(list(cap = NULL, why = ""))
  if (identical(as.character(audit$own_program %||% ""), "single_study_only")) {
    return(list(cap = "low", why = sprintf("'%s' is carried only by single-study genes in its definer audit", selected)))
  }
  if (!identical(as.character(audit$own_program %||% ""), "dominant")) {
    return(list(cap = "medium", why = sprintf(
      "the definer audit reads '%s' as %s, not dominant", selected, audit$own_program
    )))
  }
  for (entry in rivals_read(audit)) {
    if (as.character(entry$relation %||% "") %in% c("coarser_parent", "finer_child") &&
        identical(as.character(entry$whose_program_occupies %||% ""), "rival_only")) {
      return(list(cap = "medium", why = sprintf(
        "the audit reads the distinction '%s' draws against '%s' as carried by '%s' here",
        selected, entry$rival, entry$rival
      )))
    }
  }
  list(cap = NULL, why = "")
}
apply_cap <- function(final, cap) {
  if (is.null(cap)) return(list(final = final, applied = FALSE))
  current <- as.character(final$confidence %||% "low")
  if ((CONF_ORDER[current] %||% 1L) > CONF_ORDER[[cap]]) {
    final$confidence <- cap
    return(list(final = final, applied = TRUE))
  }
  list(final = final, applied = FALSE)
}
.batched <- function(items, size) {
  if (length(items) <= size) return(list(items))
  starts <- seq(1L, length(items), by = size)
  lapply(starts, function(s) items[s:min(s + size - 1L, length(items))])
}
consistency_groups <- function(results, unknown_token) {
  by_name <- list()
  for (cluster in names(results)) {
    name <- as.character(results[[cluster]]$annotation %||% "")
    if (!nzchar(name) || identical(name, unknown_token)) next
    by_name[[name]] <- c(by_name[[name]] %||% character(0), cluster)
  }
  groups <- list()
  for (name in sort(.cp(names(by_name)), method = "radix")) {
    clusters <- by_name[[name]]
    if (length(clusters) < 2L) next
    ordered <- clusters[order(nchar(clusters), clusters)]
    batches <- .batched(ordered, CONSISTENCY_GROUP_MAX)
    for (i in seq_along(batches)) {
      suffix <- if (length(batches) > 1L) sprintf(" (batch %d of %d)", i, length(batches)) else ""
      groups[[length(groups) + 1L]] <- list(
        why_grouped = sprintf("same delivered name: %s%s", name, suffix),
        clusters = sort(batches[[i]], method = "radix")
      )
    }
  }
  seen_pairs <- character(0)
  cluster_order <- names(results)[order(nchar(names(results)), names(results))]
  for (cluster in cluster_order) {
    record <- results[[cluster]]
    name <- as.character(record$annotation %||% "")
    audit <- (record$definer_audit %||% list())[[name]]
    for (entry in rivals_read(audit)) {
      rival <- as.character(entry$rival %||% "")
      if (!(as.character(entry$relation %||% "") %in% POLE_RELATIONS) ||
          !(as.character(entry$whose_program_occupies %||% "") %in% CONTESTED_OCCUPANCY)) next
      others <- Filter(function(c) identical(as.character(results[[c]]$annotation %||% ""), rival), names(results))
      others <- others[order(nchar(others), others)]
      if (!length(others)) next
      members <- sort(unique(c(cluster, others)))
      members <- members[order(nchar(members), members)]
      key <- paste(members, collapse = "\x1f")
      if (length(members) < 2L || key %in% seen_pairs) next
      seen_pairs <- c(seen_pairs, key)
      for (batch in .batched(members, CONSISTENCY_GROUP_MAX)) {
        if (length(batch) >= 2L) {
          groups[[length(groups) + 1L]] <- list(
            why_grouped = sprintf(
              "contested pair: c%s was delivered '%s' while its own audit read '%s' as %s here",
              cluster, name, rival, entry$whose_program_occupies
            ),
            clusters = batch
          )
        }
      }
    }
  }
  groups
}
consistency_packet <- function(pool, group, results) {
  names_list <- character(0)
  for (cluster in group$clusters) {
    name <- as.character(results[[cluster]]$annotation %||% "")
    if (nzchar(name) && !(name %in% names_list)) names_list <- c(names_list, name)
  }
  members <- lapply(group$clusters, function(cluster) {
    record <- results[[cluster]]
    name <- as.character(record$annotation %||% "")
    genome <- (pool$genome_de %||% list())[[as.character(cluster)]] %||% list()
    audits <- record$definer_audit %||% list()
    definers_here <- list()
    for (other in names_list) {
      entry <- .find_candidate(pool, cluster, other)
      if (!is.null(entry)) definers_here[[other]] <- as.list(compact_definers(entry$definers %||% list()))
    }
    audit <- audits[[name]] %||% list()
    top <- (genome$top_enriched %||% list())[seq_len(min(12L, length(genome$top_enriched %||% list())))]
    list(
      cluster_id = as.character(cluster),
      n_cells = as.integer(pool$clusters[[as.character(cluster)]]$n_cells),
      delivered = list(
        selected = name,
        co_occurring_identities = as.list(record$co_occurring_identities %||% list()),
        confidence = as.character(record$confidence %||% "")
      ),
      species_definers_of_each_delivered_name = definers_here,
      definer_audit = audit[intersect(names(audit), c("own_program", "definers_present", "definers_absent", "rivals_read"))],
      top_enriched = as.list(vapply(top, function(r) sprintf("%s %s/%s", r$gene, r$pct_in, r$pct_out), ""))
    )
  })
  list(
    query = pool$context,
    group = list(why_grouped = group$why_grouped, delivered_names = as.list(names_list)),
    members = members
  )
}
validate_consistency <- function(value, clusters) {
  if (!is.list(value)) return("return one JSON object")
  problems <- character(0)
  members <- value$members
  if (!is.list(members)) return("members must be a list")
  seen <- character(0)
  for (m in members) {
    if (!is.list(m)) {
      problems <- c(problems, "each member must be an object")
      next
    }
    cid <- as.character(m$cluster_id %||% "")
    if (!(cid %in% clusters)) {
      problems <- c(problems, sprintf("cluster_id '%s' is not in this group", cid))
    }
    seen <- c(seen, cid)
    if (!(as.character(m$status %||% "") %in% CONSISTENCY_VALUES)) {
      problems <- c(problems, sprintf("status must be one of %s", .py_list(CONSISTENCY_VALUES)))
    } else if (!identical(as.character(m$status), "consistent") &&
               !nzchar(trimws(as.character(m$note_for_review %||% "")))) {
      problems <- c(problems, sprintf("cluster %s: a marked member needs note_for_review", cid))
    }
  }
  missing <- setdiff(clusters, seen)
  if (length(missing)) {
    problems <- c(problems, sprintf("every member needs a status; missing %s", .py_list(head(missing, 4L))))
  }
  problems
}
dataset_consistency_many <- function(pool, results, template, header, api_key, api_url,
                                     trace_id, schema_retries) {
  groups <- consistency_groups(results, as.character(A_POOL$unknown_token %||% "Unknown"))
  if (!length(groups)) return(list())
  packets <- lapply(groups, function(g) consistency_packet(pool, g, results))
  prompts <- lapply(packets, function(p) paste0(header, .render(template, p)))
  active <- seq_along(groups)
  attempt <- rep(0L, length(groups))
  out <- vector("list", length(groups))
  while (length(active)) {
    responses <- scma_cached_call_llm_many(
      lapply(active, function(i) prompts[[i]]), api_url, api_key,
      reasoning_effort = STAGE_EFFORT, model = STAGE_MODEL,
      trace_ids = lapply(active, function(i) trace_id),
      turn_indexes = lapply(active, function(i) as.integer(900L + i))
    )
    still <- integer(0)
    for (pos in seq_along(active)) {
      i <- active[pos]
      content <- responses[[pos]][[1]]
      parsed <- if (is.null(content)) NULL else scma_parse_json(content)
      problems <- if (is.null(parsed)) "return one JSON object" else validate_consistency(parsed, groups[[i]]$clusters)
      if (!length(problems)) {
        out[[i]] <- list(group = groups[[i]], reading = parsed)
        next
      }
      if (attempt[i] >= schema_retries) {
        out[[i]] <- list(group = groups[[i]], reading = NULL)
        next
      }
      attempt[i] <- attempt[i] + 1L
      prompts[[i]] <- paste0(
        prompts[[i]], sprintf(
          "\n\n# Retry %d: the reading was rejected -- %s.\n", attempt[i],
          paste(head(problems, 4L), collapse = "; ")
        )
      )
      still <- c(still, i)
    }
    active <- still
  }
  out
}
