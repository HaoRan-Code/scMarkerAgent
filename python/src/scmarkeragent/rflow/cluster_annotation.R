#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})
.script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(.script_argument) != 1L) {
  stop("cluster_annotation.R must be executed with Rscript")
}
.sd <- dirname(normalizePath(sub("^--file=", "", .script_argument)))
source(file.path(.sd, "config.R"))
source(file.path(.sd, "llm_client.R"))
sys.source(file.path(.sd, "evidence_gate.R"), envir = environment())
source(file.path(.sd, "uberon_ontology.R"))
source(file.path(.sd, "marker_database.R"))
source(file.path(.sd, "marker_sources.R"))
source(file.path(.sd, "ortho_map.R"))
sys.source(file.path(.sd, "annotator_pool.R"), envir = environment())
sys.source(file.path(.sd, "borrowed_context.R"), envir = environment())
sys.source(file.path(.sd, "stages.R"), envir = environment())
args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args) >= 1 && nzchar(args[1])) args[1] else INPUT_TAG
A <- CFG$cluster_annotation
ANNOTATOR_MODEL <- as.character(A$annotator_model)
ANNOTATOR_EFFORT <- as.character(A$annotator_reasoning_effort)
ANNOTATOR_SCHEMA <- as.character(A$annotator_schema_version)
MAX_TURNS <- as.integer(A$max_turns)
SCHEMA_RETRIES <- as.integer(A$schema_retries)
SOURCES_PER_MARKER <- as.integer(A$sources_per_marker)
SOURCE_BATCHES_PER_MARKER <- as.integer(A$source_batches_per_marker)
UNKNOWN <- as.character(A$unknown_token)
REVIEW_MODEL <- as.character(A$review_model)
REVIEW_EFFORT <- as.character(A$review_reasoning_effort)
REVIEW_MAX_ROUNDS_PER_TIER <- as.integer(A$review_max_rounds_per_tier)
REVIEW_QUERY_TURNS <- as.integer(A$review_query_turns)
ARBITER_MODEL <- as.character(A$arbiter_model)
BORROW_SCREEN_MODEL <- as.character(A$borrow_screen_model %||% A$review_model)
BORROW_SCREEN_EFFORT <- as.character(A$borrow_screen_reasoning_effort %||% A$review_reasoning_effort)
ARBITER_EFFORT <- as.character(A$arbiter_reasoning_effort)
REVIEW_CITATIONS_REQUIRED <- 3L
CONFIDENCE_VALUES <- c("high", "medium", "low")
QUOTE_TOLERANCE <- 0.15
MAX_REPEATED_QUERIES <- 3L
RESOLVED <- "resolved"
MIXED <- "mixed"
UNRESOLVED <- "unresolved"
UNSUPPORTED <- "unsupported_empty_candidate_pool"
IDENTITY_PASS <- "accept"
IDENTITY_VALUES <- c("accept", "reject", "insufficient_evidence")
SUBTYPE_PASS <- "separable"
SUBTYPE_VALUES <- c("separable", "not_separable", "")
QC_VALUES <- c(QC_PASSED, QC_REVISED, QC_ARBITRATED, QC_FAILED, QC_UNCHECKED)
.prompt <- function(name) {
  text <- paste(
    readLines(file.path(PROMPT_DIR, paste0(name, ".txt")), warn = FALSE),
    collapse = "\n"
  )
  if (grepl(PANEL_READING_PLACEHOLDER, text, fixed = TRUE)) {
    shared <- paste(
      readLines(file.path(PROMPT_DIR, paste0(PANEL_READING, ".txt")), warn = FALSE),
      collapse = "\n"
    )
    text <- sub(PANEL_READING_PLACEHOLDER, trimws(shared, "both"), text, fixed = TRUE)
  }
  if (isTRUE(CROSS_SPECIES)) {
    note <- trimws(paste(
      readLines(file.path(PROMPT_DIR, "cross_species_note.txt"), warn = FALSE),
      collapse = "\n"
    ), "both")
    if (grepl("\n# Input\n", text, fixed = TRUE)) {
      text <- sub("\n# Input\n", paste0("\n", note, "\n\n# Input\n"), text, fixed = TRUE)
    } else if (grepl("EVIDENCE_PACKET_JSON=", text, fixed = TRUE)) {
      text <- sub(
        "EVIDENCE_PACKET_JSON=", paste0(note, "\n\nEVIDENCE_PACKET_JSON="),
        text,
        fixed = TRUE
      )
    }
  }
  text
}
.render <- function(template, packet) {
  sub(
    "{{EVIDENCE_PACKET_JSON}}",
    toJSON(packet, auto_unbox = TRUE, null = "null", digits = NA),
    template,
    fixed = TRUE
  )
}
.compact_json <- function(value) {
  toJSON(value, auto_unbox = TRUE, null = "null", digits = NA)
}
.context_header <- function(pool, cluster) {
  ctx <- pool$context %||% list()
  state <- pool$clusters[[as.character(cluster)]]
  disease <- paste(as.character(unlist(ctx$disease %||% "")), collapse = ", ")
  tissue <- paste(as.character(unlist(ctx$tissue %||% "")), collapse = ", ")
  sprintf(
    paste0(
      "# Dataset context\n",
      "species: %s | tissue: %s | disease: %s | development stage: %s | ",
      "cluster %s of %s clusters | %s cells in this cluster\n\n"
    ),
    as.character(ctx$species %||% ""), tissue, disease,
    as.character(ctx$development_stage %||% ""), as.character(cluster),
    as.character(ctx$clusters_in_dataset %||% ""),
    formatC(as.integer(state$n_cells), format = "d", big.mark = ",")
  )
}
.observation_block <- function(index, observation) {
  sprintf("\n\n# Observation %d\n%s\n", index, .compact_json(observation))
}
.valid_query <- function(value) {
  is.list(value) && identical(value$action %||% "", "query") &&
    is.character(value$tool %||% NULL) && length(value$tool) == 1L
}
.claimed_identities <- function(value) {
  claims <- list()
  selected <- as.character(value$selected %||% "")
  if (nzchar(selected) && !identical(selected, UNKNOWN)) {
    claims[[length(claims) + 1L]] <- list("selected", selected)
  }
  subtype <- as.character(value$subtype %||% "")
  if (nzchar(subtype) && !identical(subtype, selected)) {
    claims[[length(claims) + 1L]] <- list("subtype", subtype)
  }
  for (name in value$co_occurring_identities %||% list()) {
    claims[[length(claims) + 1L]] <- list("cooc", as.character(name))
  }
  claims
}
.listed_genes <- function(pool, cluster, selected, subtype) {
  listed <- character(0)
  for (name in unique(c(selected, subtype))) {
    if (!nzchar(name)) next
    entry <- .find_candidate(pool, cluster, name)
    if (is.null(entry)) next
    listed <- c(listed, toupper(vapply(entry$markers, function(m) as.character(m$gene), "")))
  }
  unique(listed)
}
.sanitize_support_markers <- function(value, pool, cluster) {
  if (!is.list(value) || !identical(as.character(value$action %||% "")[1], "final")) {
    return(list(value = value, dropped = character(0)))
  }
  support <- as.character(unlist(value$support_markers %||% list()))
  if (!length(support)) {
    return(list(value = value, dropped = character(0)))
  }
  selected <- as.character(value$selected %||% "")[1]
  subtype <- as.character(value$subtype %||% "")[1]
  listed <- .listed_genes(pool, cluster, selected, subtype)
  keep <- support[toupper(support) %in% listed]
  dropped <- support[!(toupper(support) %in% listed)]
  if (!length(dropped)) {
    return(list(value = value, dropped = character(0)))
  }
  value$support_markers <- as.list(keep)
  list(value = value, dropped = dropped)
}
.quotes_match <- function(value, pool, cluster) {
  entries <- value$claim_evidence %||% list()
  if (!is.list(entries)) {
    return(FALSE)
  }
  quoted <- character(0)
  for (entry in entries) {
    if (!is.list(entry)) {
      return(FALSE)
    }
    identity <- as.character(entry$identity %||% "")[1]
    gene <- toupper(as.character(entry$decisive_gene %||% "")[1])
    candidate <- .find_candidate(pool, cluster, identity)
    if (is.null(candidate) || !nzchar(gene)) {
      return(FALSE)
    }
    marker <- NULL
    for (row in candidate$markers) {
      if (identical(toupper(as.character(row$gene)), gene)) marker <- row
    }
    if (is.null(marker)) {
      return(FALSE)
    }
    for (field in c("pct_in", "pct_out")) {
      stated <- suppressWarnings(as.numeric(entry[[field]] %||% NA))
      measured <- suppressWarnings(as.numeric(marker[[field]] %||% NA))
      if (!length(stated) || !is.finite(stated) || !is.finite(measured) ||
        abs(stated - measured) > QUOTE_TOLERANCE) {
        return(FALSE)
      }
    }
    quoted <- c(quoted, identity)
  }
  for (claim in .claimed_identities(value)) {
    if (!(as.character(claim[[2]]) %in% quoted)) {
      return(FALSE)
    }
  }
  TRUE
}
.program_fit_problems <- function(value, program) {
  programs <- vapply(program$programs %||% list(), function(b) as.character(b$`function` %||% ""), "")
  entries <- value$program_fit
  if (!is.null(entries) && !is.list(entries) && !is.character(entries)) {
    return("program_fit must be a list, one entry per identity-bearing program of program_reading")
  }
  selected <- as.character(value$selected %||% "")[1]
  if (!length(programs) || identical(selected, UNKNOWN)) {
    return(character(0))
  }
  allowed <- c(unique(vapply(.claimed_identities(value), function(c) as.character(c[[2]]), "")), "unexplained")
  problems <- character(0)
  seen <- character(0)
  entries <- if (is.list(entries)) entries else list()
  for (entry in entries) {
    if (!is.list(entry)) {
      problems <- c(problems, "each program_fit entry must be an object")
      next
    }
    name <- as.character(entry$program %||% "")[1]
    if (!(name %in% programs)) {
      problems <- c(problems, sprintf("program_fit names a program that program_reading did not: '%s'", name))
      next
    }
    if (!(as.character(entry$accounted_for_by %||% "")[1] %in% allowed)) {
      problems <- c(problems, sprintf(
        "'%s': accounted_for_by must be one of %s exactly as written", name, .py_list(sort(.cp(allowed), method = "radix"))
      ))
    }
    seen <- c(seen, name)
  }
  missing <- setdiff(programs, seen)
  if (length(missing)) {
    problems <- c(problems, sprintf(
      "every identity-bearing program needs one program_fit entry; missing %s", .py_list(head(missing, 3L))
    ))
  }
  problems
}
.valid_final <- function(value, pool, cluster, tier = 1L, program = NULL) {
  if (!is.list(value) || !identical(value$action %||% "", "final")) {
    return(FALSE)
  }
  if (!identical(value$schema_version %||% "", ANNOTATOR_SCHEMA)) {
    return(FALSE)
  }
  if (!((value$confidence %||% "") %in% CONFIDENCE_VALUES)) {
    return(FALSE)
  }
  if (!is.character(value$reason %||% "") || !is.character(value$state %||% "")) {
    return(FALSE)
  }
  names_available <- as.character(.candidate_names(pool, cluster, tier = tier))
  selected <- value$selected
  if (!is.character(selected) || !nzchar(trimws(as.character(selected %||% "")[1]))) {
    return(FALSE)
  }
  selected <- as.character(selected)[1]
  if (!identical(selected, UNKNOWN) && !(selected %in% names_available)) {
    return(FALSE)
  }
  subtype <- value$subtype %||% ""
  if (!is.character(subtype)) {
    return(FALSE)
  }
  subtype <- as.character(subtype)[1]
  if (nzchar(subtype) && (identical(selected, UNKNOWN) || identical(subtype, selected) ||
    !(subtype %in% names_available))) {
      return(FALSE)
    }
  lineage <- value$lineage %||% ""
  if (!is.character(lineage)) {
    return(FALSE)
  }
  lineage <- as.character(lineage)[1]
  if (nzchar(lineage) && !identical(lineage, selected)) {
    if (identical(selected, UNKNOWN) || !(lineage %in% names_available)) {
      return(FALSE)
    }
  }
  others <- value$co_occurring_identities
  if (!is.null(others) && !is.list(others) && !is.character(others)) {
      return(FALSE)
    }
  others <- as.character(unlist(others %||% list()))
  if (any(!(others %in% names_available)) || any(others == selected)) {
      return(FALSE)
    }
  if (identical(selected, UNKNOWN) && length(others)) {
    return(FALSE)
  }
  refinements <- value$possible_refinements
  if (!is.null(refinements) && !is.list(refinements) && !is.character(refinements)) {
    return(FALSE)
  }
  refinements <- as.character(unlist(refinements %||% list()))
  if (any(!(refinements %in% names_available))) {
    return(FALSE)
  }
  listed <- .listed_genes(pool, cluster, selected, subtype)
  support <- as.character(unlist(value$support_markers %||% list()))
  if (length(support) && any(!(toupper(support) %in% listed))) {
      return(FALSE)
    }
  if (length(.program_fit_problems(value, program))) {
    return(FALSE)
  }
  .quotes_match(value, pool, cluster)
}
.final_problems <- function(value, pool, cluster, tier = 1L, program = NULL) {
  problems <- character(0)
  if (!is.list(value)) {
    return(problems)
  }
  names_available <- as.character(.candidate_names(pool, cluster, tier = tier))
  selected <- as.character(value$selected %||% "")[1]
  entries <- value$claim_evidence
  if (is.null(entries) || !is.list(entries)) {
    problems <- c(problems, "claim_evidence must be a list with one entry per claimed identity")
    entries <- list()
  }
  for (entry in entries) {
    if (!is.list(entry)) next
    identity <- as.character(entry$identity %||% "")[1]
    gene <- as.character(entry$decisive_gene %||% "")[1]
    candidate <- .find_candidate(pool, cluster, identity)
    if (is.null(candidate)) {
      problems <- c(problems, sprintf("'%s' is not a supplied candidate", identity))
      next
    }
    marker <- NULL
    for (row in candidate$markers) {
      if (identical(toupper(as.character(row$gene)), toupper(gene))) marker <- row
    }
    if (is.null(marker)) {
      for (row in candidate$definers %||% list()) {
        if (is.list(row) && !is.null(row$pct_in) &&
            identical(toupper(as.character(row$gene %||% "")), toupper(gene))) marker <- row
      }
    }
    if (is.null(marker)) {
      problems <- c(problems, sprintf("%s is not on the curated panel of '%s'", gene, identity))
      next
    }
    for (field in c("pct_in", "pct_out")) {
      stated <- suppressWarnings(as.numeric(entry[[field]] %||% NA))
      if (!length(stated) || !is.finite(stated)) {
        problems <- c(problems, sprintf("%s for %s must be the number shown", field, gene))
        next
      }
      measured <- suppressWarnings(as.numeric(marker[[field]] %||% NA))
      if (!is.finite(measured) || abs(stated - measured) > QUOTE_TOLERANCE) {
        problems <- c(problems, sprintf(
          "%s for %s under '%s' must be copied exactly (%s shown)",
          field, gene, identity, .fmt_pct(marker[[field]])
        ))
      }
    }
  }
  claimed <- unique(vapply(.claimed_identities(value), function(c) as.character(c[[2]]), ""))
  quoted <- unique(vapply(
    Filter(is.list, entries), function(e) as.character(e$identity %||% "")[1], ""
  ))
  for (name in sort(.cp(setdiff(claimed, quoted)), method = "radix")) {
    problems <- c(problems, sprintf("'%s' is claimed but has no claim_evidence entry", name))
  }
  problems <- c(problems, .program_fit_problems(value, program))
  if (!nzchar(selected)) {
    problems <- c(problems, "selected must be a supplied candidate or Unknown")
  } else if (!identical(selected, UNKNOWN) && !(selected %in% names_available)) {
    problems <- c(problems, sprintf("selected '%s' is not a supplied candidate", selected))
  }
  if (!identical(value$schema_version %||% "", ANNOTATOR_SCHEMA)) {
    problems <- c(problems, sprintf("schema_version must be '%s'", ANNOTATOR_SCHEMA))
  }
  if (!((value$confidence %||% "") %in% CONFIDENCE_VALUES)) {
    problems <- c(problems, sprintf("confidence must be one of %s", .py_list(CONFIDENCE_VALUES)))
  }
  for (field in c("reason", "state", "subtype", "lineage")) {
    if (!is.character(value[[field]] %||% "")) {
      problems <- c(problems, sprintf("%s must be a string", field))
    }
  }
  subtype_name <- as.character(value$subtype %||% "")[1]
  if (nzchar(subtype_name) && identical(subtype_name, selected)) {
    problems <- c(problems, paste0(
      "subtype must differ from selected (leave it empty when there is no finer name)"
    ))
  }
  for (field in c("subtype", "lineage")) {
    name <- as.character(value[[field]] %||% "")[1]
    if (nzchar(name) && !identical(name, selected)) {
      if (identical(selected, UNKNOWN)) {
        problems <- c(problems, sprintf("%s cannot be given with selected Unknown", field))
      } else if (!(name %in% names_available)) {
        problems <- c(problems, sprintf("%s '%s' is not a supplied candidate", field, name))
      }
    }
  }
  others <- value$co_occurring_identities
  if (!is.null(others) && !is.list(others) && !is.character(others)) {
    problems <- c(problems, "co_occurring_identities must be a list")
  }
  for (name in as.character(unlist(others %||% list()))) {
    if (!(name %in% names_available)) {
      problems <- c(problems, sprintf("co-occurring '%s' is not a supplied candidate", name))
    } else if (identical(name, selected)) {
      problems <- c(problems, sprintf("co-occurring '%s' is the selected identity itself", name))
    } else if (identical(selected, UNKNOWN)) {
      problems <- c(problems, "co_occurring_identities must be empty with selected Unknown")
    }
  }
  refinements <- value$possible_refinements
  if (!is.null(refinements) && !is.list(refinements) && !is.character(refinements)) {
    problems <- c(problems, "possible_refinements must be a list")
  }
  for (name in as.character(unlist(refinements %||% list()))) {
    if (!(name %in% names_available)) {
      problems <- c(problems, sprintf(
        "possible_refinements '%s' is not a supplied candidate; only candidate names can go there",
        name
      ))
    }
  }
  listed_genes <- .listed_genes(pool, cluster, selected, subtype_name)
  bad_support <- as.character(unlist(value$support_markers %||% list()))
  bad_support <- bad_support[!(toupper(bad_support) %in% listed_genes)]
  if (length(bad_support)) {
    where <- paste0(
      sprintf("'%s'", selected),
      if (nzchar(subtype_name) && !identical(subtype_name, selected)) {
        sprintf(" or '%s'", subtype_name)
      } else {
        ""
      }
    )
    problems <- c(problems, sprintf(
      "support_markers not on %s panel: %s", where, .py_list(head(bad_support, 6L))
    ))
  }
  problems
}
.packet_citations <- function(packet) {
  out <- character(0)
  for (block_key in c("claimed_identities", "contested_identities")) {
    for (block in packet[[block_key]] %||% list()) {
      identity <- as.character(block$cell_type %||% "")[1]
      sources_by_gene <- block$sources %||% list()
      for (gene in names(sources_by_gene)) {
        for (record in sources_by_gene[[gene]] %||% list()) {
          for (key in c("pmcid", "pmid")) {
            value <- trimws(as.character(record[[key]] %||% "")[1])
            if (nzchar(value) && !identical(value, NOT_AVAILABLE)) {
              out <- c(out, paste(identity, toupper(gene), value, sep = "\x1f"))
            }
          }
        }
      }
    }
  }
  unique(out)
}
.packet_sentence_count <- function(packet) {
  count <- 0L
  for (block_key in c("claimed_identities", "contested_identities")) {
    for (block in packet[[block_key]] %||% list()) {
      sources_by_gene <- block$sources %||% list()
      for (gene in names(sources_by_gene)) {
        for (record in sources_by_gene[[gene]] %||% list()) {
          pmcid <- trimws(as.character(record$pmcid %||% "")[1])
          pmid <- trimws(as.character(record$pmid %||% "")[1])
          if ((nzchar(pmcid) && !identical(pmcid, NOT_AVAILABLE)) ||
              (nzchar(pmid) && !identical(pmid, NOT_AVAILABLE))) {
            count <- count + 1L
          }
        }
      }
    }
  }
  count
}
.citations_required <- function(packet) {
  min(REVIEW_CITATIONS_REQUIRED, .packet_sentence_count(packet))
}
.cited_ok <- function(verdict, packet) {
  available <- .packet_citations(packet)
  required <- .citations_required(packet)
  entries <- verdict$evidence_cited
  if (!is.list(entries)) {
    return(list(ok = FALSE, problems = "evidence_cited must be a list of {identity, gene, citation}"))
  }
  seen <- character(0)
  problems <- character(0)
  for (entry in entries) {
    if (!is.list(entry)) {
      problems <- c(problems, "each evidence_cited entry must be an object")
      next
    }
    identity <- as.character(entry$identity %||% "")[1]
    gene <- toupper(as.character(entry$gene %||% "")[1])
    citation <- trimws(as.character(entry$citation %||% "")[1])
    key <- paste(identity, gene, citation, sep = "\x1f")
    if (key %in% available) {
      seen <- union(seen, key)
    } else {
      problems <- c(problems, sprintf(
        "evidence_cited %s for %s under '%s' is not a sentence this packet supplied",
        if (nzchar(citation)) citation else "(empty)",
        if (nzchar(gene)) gene else "(no gene)", identity
      ))
    }
  }
  if (length(seen) < required) {
    problems <- c(problems, sprintf(
      paste0(
        "cite at least %d distinct sentences from `sources` in this packet (identity, gene ",
        "and the citation exactly as shown); %d of them were real"
      ), required, length(seen)
    ))
  }
  list(ok = !length(problems), problems = problems)
}
.absent_definers <- function(packet) {
  selected <- as.character((packet$delivered %||% list())$selected %||% "")
  audit <- (packet$definer_audit %||% list())[[selected]] %||% list()
  Filter(nzchar, as.character(unlist(audit$definers_absent %||% list())))
}
ABSENT_DISPOSITIONS <- c("dropout_single_gene", "tissue_context_absent", "refutes")
.disposition_ok <- function(verdict, packet) {
  absent <- toupper(.absent_definers(packet))
  if (!length(absent)) return(character(0))
  entries <- verdict$absent_definers_disposition
  if (!is.list(entries)) {
    return(sprintf(
      "absent_definers_disposition must be a list, one entry per gene the definer audit reads as absent here: %s",
      .py_list(sort(.cp(absent), method = "radix")[seq_len(min(6L, length(absent)))])
    ))
  }
  problems <- character(0)
  seen <- character(0)
  for (entry in entries) {
    if (!is.list(entry)) {
      problems <- c(problems, "each absent_definers_disposition entry must be an object")
      next
    }
    gene <- toupper(as.character(entry$gene %||% "")[1])
    if (!(gene %in% absent)) {
      problems <- c(problems, sprintf("absent_definers_disposition names '%s', which the audit does not read as absent here", gene))
      next
    }
    seen <- c(seen, gene)
    if (!(as.character(entry$disposition %||% "")[1] %in% ABSENT_DISPOSITIONS)) {
      problems <- c(problems, sprintf("%s: disposition must be one of %s", gene, .py_list(ABSENT_DISPOSITIONS)))
    }
    if (!nzchar(trimws(as.character(entry$reason %||% "")[1]))) {
      problems <- c(problems, sprintf("%s: reason must say what the rows and the sentences show", gene))
    }
  }
  missing <- setdiff(absent, seen)
  if (length(missing)) {
    problems <- c(problems, sprintf(
      "every gene the audit reads as absent needs a disposition; missing %s",
      .py_list(sort(.cp(missing), method = "radix")[seq_len(min(6L, length(missing)))])
    ))
  }
  if (identical(as.character(verdict$identity_verdict %||% "")[1], IDENTITY_PASS) &&
      any(vapply(entries, function(e) is.list(e) && identical(as.character(e$disposition %||% ""), "refutes"), TRUE))) {
    problems <- c(problems, paste0(
      "a definer marked `refutes` and an `accept` are the verdict disagreeing with itself: ",
      "return 'reject' (with the candidate this packet holds in better_candidate), or say ",
      "what the rows show that makes that gene a dropout or absent from this tissue's ",
      "context rather than a refutation"
    ))
  }
  problems
}
.detected_exclusions <- function(packet) {
  selected <- as.character((packet$delivered %||% list())$selected %||% "")[1]
  for (block in packet$claimed_identities %||% list()) {
    if (identical(as.character(block$cell_type %||% ""), selected) &&
        identical(as.character(block$role %||% ""), "selected")) {
      genes <- character(0)
      for (row in block$detected_exclusions %||% list()) {
        gene <- if (is.list(row)) as.character(row$gene %||% "") else sub("^(\\S+).*$", "\\1", as.character(row))
        if (nzchar(gene)) genes <- c(genes, toupper(gene))
      }
      return(sort(unique(genes), method = "radix"))
    }
  }
  character(0)
}
EXCLUSION_DISPOSITIONS <- c("background", "minority", "refutes")
.exclusion_disposition_ok <- function(verdict, packet) {
  detected <- .detected_exclusions(packet)
  if (!length(detected)) return(character(0))
  entries <- verdict$detected_exclusions_disposition
  if (!is.list(entries)) {
    return(sprintf(
      "detected_exclusions_disposition must be a list, one entry per curated exclusion of the delivered identity this cluster detects: %s",
      .py_list(detected[seq_len(min(6L, length(detected)))])
    ))
  }
  problems <- character(0)
  seen <- character(0)
  for (entry in entries) {
    if (!is.list(entry)) {
      problems <- c(problems, "each detected_exclusions_disposition entry must be an object")
      next
    }
    gene <- toupper(as.character(entry$gene %||% "")[1])
    if (!(gene %in% detected)) {
      problems <- c(problems, sprintf("detected_exclusions_disposition names '%s', which is not a detected exclusion of the delivered identity here", gene))
      next
    }
    seen <- c(seen, gene)
    if (!(as.character(entry$disposition %||% "")[1] %in% EXCLUSION_DISPOSITIONS)) {
      problems <- c(problems, sprintf("%s: disposition must be one of %s", gene, .py_list(EXCLUSION_DISPOSITIONS)))
    }
    if (!nzchar(trimws(as.character(entry$reason %||% "")[1]))) {
      problems <- c(problems, sprintf("%s: reason must say what the rows and the sentences show", gene))
    }
  }
  missing <- setdiff(detected, seen)
  if (length(missing)) {
    problems <- c(problems, sprintf(
      "every detected exclusion of the delivered identity needs a disposition; missing %s",
      .py_list(missing[seq_len(min(6L, length(missing)))])
    ))
  }
  if (identical(as.character(verdict$identity_verdict %||% "")[1], IDENTITY_PASS) &&
      any(vapply(entries, function(e) is.list(e) && identical(as.character(e$disposition %||% ""), "refutes"), TRUE))) {
    problems <- c(problems, paste0(
      "a detected exclusion marked `refutes` and an `accept` are the verdict disagreeing with itself: ",
      "return 'reject' (with the candidate the evidence carries in better_candidate, where this ",
      "packet holds it), or say what the rows show that makes that gene the tissue's background ",
      "or a minority here rather than a refutation"
    ))
  }
  problems
}
.separating_sentences <- function(packet) {
  selected <- as.character((packet$delivered %||% list())$selected %||% "")[1]
  shared <- toupper(names(packet$shared_with_contested %||% list()))
  out <- character(0)
  for (key in .packet_citations(packet)) {
    parts <- strsplit(key, "\x1f", fixed = TRUE)[[1]]
    if (length(parts) >= 2L && identical(parts[1], selected) && !(parts[2] %in% shared)) {
      out <- c(out, parts[2])
    }
  }
  sort(unique(out))
}
.distinguishing_citation <- function(verdict, packet) {
  contested <- vapply(packet$contested_identities %||% list(), function(b) as.character(b$cell_type), "")
  if (!length(contested)) return(character(0))
  shared <- toupper(names(packet$shared_with_contested %||% list()))
  selected <- as.character((packet$delivered %||% list())$selected %||% "")
  own <- character(0)
  for (entry in verdict$evidence_cited %||% list()) {
    if (is.list(entry) && identical(as.character(entry$identity %||% ""), selected)) {
      own <- c(own, toupper(as.character(entry$gene %||% "")))
    }
  }
  if (any(nzchar(own) & !(own %in% shared))) return(character(0))
  ## Two names whose curated definitions the resource shares gene for gene leave nothing
  ## to cite, and the demand would then refuse every verdict. What the evidence reaches
  ## for such a pair is their common parent or a `mixed` answer, which is the reviewer's
  ## reading to give.
  if (!length(.separating_sentences(packet))) return(character(0))
  sprintf(
    paste0(
      "every sentence cited for '%s' is for a gene that %s also curates, so the citations do not ",
      "reach the distinction this verdict settles: cite at least one sentence for a gene of '%s' ",
      "that the contested identity does not claim"
    ), selected, paste(sprintf("'%s'", contested), collapse = ", "), selected
  )
}
.population_ok <- function(verdict, packet) {
  required <- population_names(packet)
  entries <- verdict$population_verdicts
  if (is.null(entries)) {
    if (length(required)) {
      return(sprintf("population_verdicts is required: one entry for each of %s", .py_list(head(required, 6L))))
    }
    entries <- list()
  }
  if (!is.list(entries)) {
    return("population_verdicts must be a list")
  }
  problems <- character(0)
  known <- packet_genes(packet)
  seen <- character(0)
  for (entry in entries) {
    if (!is.list(entry)) {
      problems <- c(problems, "each population_verdicts entry must be an object")
      next
    }
    name <- as.character(entry$identity %||% "")[1]
    if (!(name %in% required)) {
      problems <- c(problems, sprintf(
        "population_verdicts names '%s', which is neither a co-occurring identity of the answer nor a contested identity of this packet", name
      ))
      next
    }
    seen <- c(seen, name)
    if (!(as.character(entry$verdict %||% "")[1] %in% POPULATION_VERDICTS)) {
      problems <- c(problems, sprintf("%s: verdict must be one of %s", name, .py_list(POPULATION_VERDICTS)))
    }
    genes <- entry$genes
    if ((!is.list(genes) && !is.character(genes)) || !length(genes)) {
      problems <- c(problems, sprintf("%s: genes must list the measured genes the verdict rests on", name))
    } else {
      genes <- as.character(unlist(genes))
      bad <- genes[!(toupper(genes) %in% known)]
      if (length(bad)) {
        problems <- c(problems, sprintf("%s: genes names genes not in this packet: %s", name, .py_list(head(bad, 4L))))
      }
    }
    if (!nzchar(trimws(as.character(entry$reason %||% "")[1]))) {
      problems <- c(problems, sprintf("%s: reason must say what the rows and the sentences show", name))
    }
  }
  missing <- setdiff(required, seen)
  if (length(missing)) {
    problems <- c(problems, sprintf("every name the packet asks about needs a population verdict; missing %s", .py_list(head(missing, 4L))))
  }
  if (identical(as.character(verdict$identity_verdict %||% "")[1], IDENTITY_PASS)) {
    undivided <- character(0)
    for (entry in entries) {
      if (is.list(entry) && identical(as.character(entry$verdict %||% ""), NEITHER_SEPARATES)) {
        undivided <- c(undivided, as.character(entry$identity %||% "")[1])
      }
    }
    if (length(undivided)) {
      problems <- c(problems, sprintf(paste0(
        "'%s' for %s and an `accept` are the verdict disagreeing with itself: where no gene in this ",
        "packet tells the delivered identity from that name, neither is established alone -- return ",
        "'reject' with their common parent in better_candidate where this packet holds it (empty ",
        "otherwise), or read the pair as one of the other verdicts"
      ), NEITHER_SEPARATES, .py_list(head(undivided, 3L))))
    }
  }
  problems
}
.valid_review <- function(verdict, packet, subtype) {
  if (!is.list(verdict)) {
    return(list(ok = FALSE, problems = "return one JSON object with the verdict"))
  }
  problems <- character(0)
  identity_verdict <- as.character(verdict$identity_verdict %||% "")[1]
  if (!(identity_verdict %in% IDENTITY_VALUES)) {
    problems <- c(problems, sprintf("identity_verdict must be one of %s", .py_list(IDENTITY_VALUES)))
  } else if (identical(identity_verdict, IDENTITY_PASS) &&
    nzchar(trimws(as.character(verdict$better_candidate %||% "")[1]))) {
    problems <- c(problems, sprintf(paste0(
      "identity_verdict is '%s' but better_candidate is not empty -- naming a better ",
      "candidate is a rejection of the delivered identity: return 'reject' with that ",
      "candidate in better_candidate, or clear better_candidate if the delivered identity ",
      "is in fact the one to accept"
    ), IDENTITY_PASS))
  }
  subtype_verdict <- as.character(verdict$subtype_verdict %||% "")[1]
  if (!(subtype_verdict %in% SUBTYPE_VALUES)) {
    problems <- c(problems, sprintf("subtype_verdict must be one of %s", .py_list(SUBTYPE_VALUES)))
  } else if (nzchar(subtype) && !nzchar(subtype_verdict)) {
    problems <- c(problems, sprintf(
      "a finer name ('%s') was delivered, so subtype_verdict must be 'separable' or 'not_separable'",
      subtype
    ))
  }
  if (!is.character(verdict$reason %||% "") || !nzchar(trimws(as.character(verdict$reason %||% "")[1]))) {
    problems <- c(problems, "reason must name genes with their measured percentages")
  }
  cited <- .cited_ok(verdict, packet)
  problems <- c(problems, cited$problems)
  if (isTRUE(cited$ok)) {
    problems <- c(problems, .distinguishing_citation(verdict, packet))
  }
  problems <- c(problems, .disposition_ok(verdict, packet))
  problems <- c(problems, .exclusion_disposition_ok(verdict, packet))
  problems <- c(problems, .population_ok(verdict, packet))
  list(ok = !length(problems), problems = problems)
}
.review_passed <- function(verdict, subtype) {
  if (is.null(verdict) || !length(verdict)) {
    return(FALSE)
  }
  if (!identical(as.character(verdict$identity_verdict %||% ""), IDENTITY_PASS)) {
    return(FALSE)
  }
  !nzchar(subtype) || identical(as.character(verdict$subtype_verdict %||% ""), SUBTYPE_PASS)
}
.review_many <- function(pool, clusters, finals, source_servers, conversation, api_url, api_key) {
  out <- setNames(vector("list", length(clusters)), clusters)
  packets <- list()
  labels_of <- list()
  prompts <- list()
  review_servers <- list()
  seen <- list()
  unusable <- list()
  asked <- list()
  ## The last well-formed verdict that did not clear the citation requirement. Kept so a
  ## reviewer that answered every time is recorded as having answered and decided by
  ## arbitration, rather than reported as a reviewer nobody could reach.
  last_verdict <- list()
  active <- character(0)
  for (cluster in clusters) {
    final <- finals[[cluster]]
    annotator_server <- source_servers[[cluster]]
    server <- if (.is_source_server(annotator_server)) {
      scma_source_server(
        annotator_server$context, batch = SOURCES_PER_MARKER, max_batches = SOURCE_BATCHES_PER_MARKER
      )
    } else {
      annotator_server
    }
    review_servers[[cluster]] <- server
    packet <- .review_packet(
      pool, cluster, final, server, tier = conversation[[cluster]]$tier, unknown_token = UNKNOWN
    )
    if (length(packet$claimed_identities)) {
      required <- population_names(packet)
      if (length(required)) packet$population_verdicts_required <- as.list(required)
      if (length(packet$contested_identities %||% list())) {
        separating <- .separating_sentences(packet)
        if (length(separating)) packet$separating_genes <- as.list(separating)
      }
    }
    packets[[cluster]] <- packet
    if (!length(packet$claimed_identities)) {
      out[[cluster]] <- list(verdict = NULL, packet = packet)
      next
    }
    here <- unique(as.character(vapply(
      packet$claimed_identities, function(b) as.character(b$cell_type %||% ""), ""
    )))
    labels_of[[cluster]] <- here[nzchar(here)]
    prompts[[cluster]] <- paste0(
      .context_header(pool, cluster), .render(.prompt("evidence_review"), packet)
    )
    seen[[cluster]] <- list()
    unusable[[cluster]] <- 0L
    asked[[cluster]] <- 0L
    active <- c(active, cluster)
  }
  if (!length(active)) {
    return(out)
  }
  for (round_index in seq_len(REVIEW_QUERY_TURNS + 1L + SCHEMA_RETRIES)) {
    responses <- scma_cached_call_llm_many(
      lapply(active, function(cl) prompts[[cl]]), api_url, api_key,
      reasoning_effort = REVIEW_EFFORT, model = REVIEW_MODEL,
      trace_ids = lapply(active, function(cl) conversation[[cl]]$trace_id),
      turn_indexes = lapply(active, function(cl) as.integer(conversation[[cl]]$turn))
    )
    still <- character(0)
    for (position in seq_along(active)) {
      cluster <- active[position]
      subtype <- as.character(finals[[cluster]]$subtype %||% "")[1]
      content <- responses[[position]][[1]]
      parsed <- if (is.null(content)) NULL else scma_parse_json(content)
      if (is.list(parsed) && nzchar(as.character(parsed$identity_verdict %||% ""))) {
        outcome <- .valid_review(parsed, packets[[cluster]], subtype)
        if (isTRUE(outcome$ok)) {
          parsed$citations_required <- .citations_required(packets[[cluster]])
          if (length(packets[[cluster]]$contested_identities %||% list())) {
            parsed$separating_sentences <- as.list(.separating_sentences(packets[[cluster]]))
          }
          out[[cluster]] <- list(verdict = parsed, packet = packets[[cluster]], status = "ok")
          next
        }
        if (as.character(parsed$identity_verdict %||% "")[1] %in% IDENTITY_VALUES) {
          last_verdict[[cluster]] <- parsed
        }
        if (unusable[[cluster]] >= SCHEMA_RETRIES) {
          if (!is.null(last_verdict[[cluster]])) {
            verdict <- last_verdict[[cluster]]
            verdict$citations_required <- .citations_required(packets[[cluster]])
            if (length(packets[[cluster]]$contested_identities %||% list())) {
              verdict$separating_sentences <- as.list(.separating_sentences(packets[[cluster]]))
            }
            verdict$uncited <- TRUE
            verdict$uncited_problems <- as.list(head(outcome$problems, 4L))
            out[[cluster]] <- list(verdict = verdict, packet = packets[[cluster]], status = "uncited")
            next
          }
          out[[cluster]] <- list(verdict = NULL, packet = packets[[cluster]])
          next
        }
        unusable[[cluster]] <- unusable[[cluster]] + 1L
        prompts[[cluster]] <- paste0(
          prompts[[cluster]], sprintf(
            "\n\n# Retry %d: the verdict was rejected -- %s. Return one valid verdict object.\n",
            unusable[[cluster]], paste(head(outcome$problems, 4L), collapse = "; ")
          )
        )
        still <- c(still, cluster)
        next
      }
      answerable <- is.list(parsed) && .valid_query(parsed) && asked[[cluster]] < REVIEW_QUERY_TURNS
      if (!answerable) {
        if (unusable[[cluster]] >= SCHEMA_RETRIES) {
          out[[cluster]] <- list(verdict = NULL, packet = packets[[cluster]])
          next
        }
        unusable[[cluster]] <- unusable[[cluster]] + 1L
        prompts[[cluster]] <- paste0(
          prompts[[cluster]],
          sprintf("\n\n# Retry %d: return the verdict as one JSON object.\n", unusable[[cluster]])
        )
        still <- c(still, cluster)
        next
      }
      call <- .compact_json(list(tool = as.character(parsed$tool)[1], args = parsed$args %||% list()))
      streaming <- identical(as.character(parsed$tool)[1], "sources")
      observation <- if (streaming) NULL else seen[[cluster]][[call]]
      if (is.null(observation)) {
        observation <- .run_review_tool(
          pool, cluster, labels_of[[cluster]], parsed$tool, parsed$args %||% list(),
          review_servers[[cluster]]
        )
        seen[[cluster]][[call]] <- observation
      } else {
        observation <- c(observation, list(duplicate_query = TRUE))
      }
      asked[[cluster]] <- asked[[cluster]] + 1L
      prompts[[cluster]] <- paste0(
        prompts[[cluster]], .observation_block(asked[[cluster]], observation),
        if (asked[[cluster]] >= REVIEW_QUERY_TURNS) {
          "\n# No further queries. Return the verdict now.\n"
        } else {
          ""
        }
      )
      still <- c(still, cluster)
    }
    active <- still
    if (!length(active)) break
  }
  for (cluster in active) {
    out[[cluster]] <- list(verdict = NULL, packet = packets[[cluster]])
  }
  out
}
.pair_tables <- function(pool, cluster, selected, rivals) {
  sel <- .find_candidate(pool, cluster, selected)
  if (is.null(sel)) return(list())
  out <- list()
  for (name in rivals) {
    rentry <- .find_candidate(pool, cluster, as.character(name))
    if (!is.null(rentry)) out[[length(out) + 1L]] <- pair_partition(sel, rentry)
  }
  out
}
.contested_feedback <- function(verdict, reasons, tier, rounds_left, tiers_left,
                                pair_tables_list = NULL, population = NULL) {
  tier <- as.integer(tier)
  rounds_left <- as.integer(rounds_left)
  tiers_left <- as.integer(tiers_left)
  grow_text <- if (tiers_left > 0L) {
    sprintf(" Only where that decision is Unknown are the next %d candidates added.", tiers_left * 15L)
  } else {
    ""
  }
  block <- list(
    evidence_review = list(
      identity_verdict = verdict$identity_verdict, subtype_verdict = verdict$subtype_verdict %||% "",
      supporting_markers = verdict$supporting_markers %||% list(), conflicting_markers = verdict$conflicting_markers %||% list(),
      evidence_cited = verdict$evidence_cited %||% list(), population_verdicts = verdict$population_verdicts %||% list(),
      reason = verdict$reason, flags = verdict$flags %||% list()
    ),
    contested = as.list(reasons), population_conflicts = as.list(population %||% list()),
    note = paste0(
      "the evidence review accepted your identity. What did not close is listed above: under ",
      "`contested`, the definer audit of the name you delivered, made on its OWN rows, ",
      "contradicting that delivery; under `population_conflicts`, the reviewer's reading of each ",
      "co-occurring identity you reported and of each contested identity you did not ",
      "(`population_verdicts`) where it disagrees with your list. Answer again from the ",
      "candidates you already hold, by [R8] and [A5]; `pair_partitions` below lists, for each ",
      "disputed pair, what each side curates that the other does not and what both do, as ",
      "measured here. ",
      sprintf("Candidates %d of the retrieval order are on this page. ", tier * 15L),
      if (rounds_left > 0L) {
        sprintf("%d further attempt(s) at this tier, then the answer is decided from the rounds so far.%s", rounds_left, grow_text)
      } else {
        paste0("This is the last attempt at this tier; after it the answer is decided from the rounds so far.", grow_text)
      }
    )
  )
  if (length(pair_tables_list)) block$pair_partitions <- pair_tables_list
  block
}
.review_feedback <- function(verdict, tier, rounds_left, tiers_left, audits = NULL,
                             better_audit = NULL, pair_tables_list = NULL) {
  tier <- as.integer(tier)
  rounds_left <- as.integer(rounds_left)
  tiers_left <- as.integer(tiers_left)
  tail_text <- if (rounds_left > 0L) {
    paste0(
      sprintf("%d further attempt(s) at this tier", rounds_left),
      if (tiers_left > 0L) {
        sprintf(", then the next %d candidates are added.", tiers_left * 15L)
      } else {
        ", then the answer is decided from the rounds so far."
      }
    )
  } else if (tiers_left > 0L) {
    sprintf("The next %d candidates are added after this answer.", tiers_left * 15L)
  } else {
    "This is the last attempt; after it the answer is decided from the rounds so far."
  }
  note <- paste0(
    "your answer was checked against the curated source sentences behind the whole panel ",
    "of every identity you claimed -- positive rows and curated exclusions alike -- and it ",
    "did not survive. The reviewer read the sentences it cites above; it did not see your ",
    "candidate ranking or your confidence. Answer again against the candidates you already ",
    "hold: keep the call only where the panel and its sentences still carry it, otherwise ",
    "name the candidate they do carry, or drop a finer name that is not separable. ",
    sprintf("Candidates %d of the retrieval order are on this page. ", tier * 15L),
    tail_text,
    " If no candidate on this page is carried by its own evidence, answer Unknown and say ",
    "which evidence nothing accounts for."
  )
  has_better <- nzchar(trimws(as.character(verdict$better_candidate %||% "")[1]))
  legal <- paste0(
    "A refusal is not an instruction to move to another leaf of the page. The answers the ",
    "evidence can carry are: (1) the name you gave, where its own lineage genes stand in ",
    "these cells and you can answer the reviewer's objection from the rows and the ",
    "sentences -- say which rows; ",
    if (has_better) {
      "(2) the reviewer's better_candidate, only under the conditions in better_candidate_note; "
    } else {
      "(2) -- the reviewer named no better candidate, so none is on offer; "
    },
    "(3) another candidate whose OWN lineage genes occupy these cells at cluster level -- it ",
    "is audited on the same terms before any review, and a name whose own program reads ",
    "minority or absent there comes back refused; (4) their common parent, where this page ",
    "holds it and neither separates; (5) a finer name dropped for the name it was delivered ",
    "under, where the finer one is not separable; or (6) Unknown, where nothing on this page ",
    "is carried by its own evidence -- say which measurements nothing accounts for. Moving to ",
    "a candidate whose own program is weaker here than the one refused is none of these, and ",
    "a co-occurring identity is reported on its own lineage genes in its own share of the ",
    "cells, not as a way of answering a refusal. "
  )
  block <- list(
    evidence_review = list(
      identity_verdict = verdict$identity_verdict,
      subtype_verdict = verdict$subtype_verdict %||% "",
      better_candidate = as.character(verdict$better_candidate %||% "")[1],
      supporting_markers = verdict$supporting_markers %||% list(),
      conflicting_markers = verdict$conflicting_markers %||% list(),
      evidence_cited = verdict$evidence_cited %||% list(),
      population_verdicts = verdict$population_verdicts %||% list(),
      reason = verdict$reason,
      flags = verdict$flags %||% list()
    ),
    note = paste0(
      "your answer was checked against the curated source sentences behind the whole panel ",
      "of every identity you claimed -- positive rows and curated exclusions alike -- and it ",
      "did not survive. The reviewer read the sentences it cites above; it did not see your ",
      "candidate ranking or your confidence. Answer again against the candidates you already ",
      "hold. ", legal, tail_text
    )
  )
  if (length(audits)) {
    block$definer_audit_of_your_claims <- lapply(audits, function(a) a[setdiff(names(a), "species_definers")])
  }
  if (length(pair_tables_list)) block$pair_partitions <- pair_tables_list
  if (length(better_audit)) {
    block$definer_audit_of_better_candidate <- better_audit[setdiff(names(better_audit), "species_definers")]
    block$better_candidate_note <- paste0(
      "the reviewer's better_candidate was put through the same definer audit. A name the ",
      "reviewer prefers is a name to check, not to adopt: take it as `selected` only where ",
      "that audit reads its own program as `dominant` here on its OWN lineage genes. ",
      "`minority` means at most a co-occurring identity beside a `selected` whose own ",
      "program does occupy these cells; `absent` means the page does not carry that name at ",
      "all, and the reviewer's objection to your answer still has to be met from the rows -- ",
      "keep your name where its own lineage genes carry it, move to the common parent where ",
      "neither name separates, or answer Unknown"
    )
  }
  block
}
.arbiter_packet <- function(pool, cluster, attempts, packet, tier, contested = NULL,
                            audits = NULL, pair_tables_list = NULL) {
  claimed <- list()
  for (block in packet$claimed_identities %||% list()) {
    claimed[[as.character(block$cell_type)]] <- block
  }
  for (item in attempts) {
    for (claim in .claimed_identities(item$final)) {
      name <- as.character(claim[[2]])
      if (!is.null(claimed[[name]])) next
      entry <- .find_candidate(pool, cluster, name)
      if (is.null(entry)) next
      claimed[[name]] <- list(
        cell_type = name,
        panel = as.list(compact_panel(.review_panel(entry))),
        facts = .packet_facts(.candidate_facts(pool, cluster, entry), entry),
        detected_exclusions = as.list(compact_panel(.detected_exclusion_markers(entry)))
      )
    }
  }
  answer_fields <- c(
    "selected", "subtype", "lineage", "state", "co_occurring_identities",
    "possible_refinements", "support_markers", "claim_evidence", "confidence", "reason"
  )
  review_fields <- c(
    "identity_verdict", "subtype_verdict", "better_candidate",
    "supporting_markers", "conflicting_markers", "evidence_cited", "population_verdicts",
    "reason", "flags", "uncited"
  )
  audits_with_rows <- audits %||% list()
  rounds <- lapply(seq_along(attempts), function(index) {
    item <- attempts[[index]]
    review <- if (length(item$review)) {
      item$review[intersect(review_fields, names(item$review))]
    } else {
      list(identity_verdict = "not_checked")
    }
    for (name in names(item$final$definer_audit %||% list())) {
      if (is.null(audits_with_rows[[name]])) audits_with_rows[[name]] <- item$final$definer_audit[[name]]
    }
    list(
      round = index, tier = as.integer(item$tier),
      answer = item$final[intersect(answer_fields, names(item$final))],
      review = review,
      definer_audit = lapply(item$final$definer_audit %||% list(), function(a) a[setdiff(names(a), "species_definers")]),
      structural_note = as.character(item$structural_note %||% "")
    )
  })
  packet_out <- list(
    query = packet$query %||% list(),
    rounds = rounds,
    row_format = MARKER_ROW_FORMAT,
    identities_claimed = unname(claimed),
    definer_audits = lapply(audits_with_rows, function(a) a[setdiff(names(a), "species_definers")]),
    program_reading = packet$program_reading %||% list(),
    candidates_not_claimed = packet$candidates_not_claimed %||% list(),
    cluster_top_enriched = packet$cluster_top_enriched %||% list(),
    cluster_top_detection_gap = packet$cluster_top_detection_gap %||% list(),
    cluster_top_depleted = packet$cluster_top_depleted %||% list(),
    candidates_available = as.list(.candidate_names(pool, cluster, tier = tier))
  )
  screened <- .screened_out_claimants(pool, cluster)
  if (length(screened)) packet_out$screened_out_claimants <- screened
  if (length(contested)) packet_out$contested <- as.list(contested)
  if (length(pair_tables_list)) packet_out$pair_partitions <- pair_tables_list
  packet_out
}
.population_verdicts_by_name <- function(attempts) {
  out <- list()
  for (item in attempts %||% list()) {
    for (entry in (item$review %||% list())$population_verdicts %||% list()) {
      if (is.list(entry) && nzchar(as.character(entry$identity %||% "")[1])) {
        name <- as.character(entry$identity)[1]
        out[[name]] <- c(out[[name]] %||% character(0), as.character(entry$verdict %||% "")[1])
      }
    }
  }
  out
}
.arbiter_problems <- function(value, names_available, verdicts_of = NULL) {
  if (!is.list(value)) {
    return("return one JSON object with the decision")
  }
  problems <- character(0)
  selected <- value$selected
  selected_value <- ""
  if (!is.character(selected) || !nzchar(trimws(as.character(selected %||% "")[1]))) {
    problems <- c(problems, "selected must be a candidate name or Unknown")
  } else {
    selected_value <- as.character(selected)[1]
    if (!identical(selected_value, UNKNOWN) && !(selected_value %in% names_available)) {
      problems <- c(problems, sprintf(
        "selected '%s' is not one of the candidates offered", selected_value
      ))
    }
  }
  for (field in c("subtype", "lineage")) {
    name <- value[[field]]
    if (!is.null(name) && !is.character(name)) {
      problems <- c(problems, sprintf("%s must be a string", field))
    } else {
      name <- as.character(name %||% "")[1]
      if (nzchar(name) && !identical(name, selected_value) && !(name %in% names_available)) {
        problems <- c(problems, sprintf(
          "%s '%s' is not one of the candidates offered", field, name
        ))
      }
    }
  }
  others <- value$co_occurring_identities
  if (!is.null(others) && !is.list(others) && !is.character(others)) {
    problems <- c(problems, "co_occurring_identities must be a list")
  } else {
    for (name in as.character(unlist(others %||% list()))) {
      if (!(name %in% names_available)) {
        problems <- c(problems, sprintf(
          "co-occurring '%s' is not one of the candidates offered", name
        ))
        next
      }
      read <- (verdicts_of %||% list())[[name]] %||% character(0)
      if (length(read) && !(ESTABLISHED %in% read)) {
        problems <- c(problems, sprintf(paste0(
          "co-occurring '%s' was read by the rounds' reviews as %s and by none as an established ",
          "second population; a co-occurring identity is delivered only on a review that read it ",
          "so ([F2c]) -- leave it out"
        ), name, .py_list(sort(unique(read), method = "radix"))))
      }
    }
  }
  if (!((value$confidence %||% "") %in% CONFIDENCE_VALUES)) {
    problems <- c(problems, sprintf("confidence must be one of %s", .py_list(CONFIDENCE_VALUES)))
  }
  if (!is.character(value$reason %||% "") || !nzchar(trimws(as.character(value$reason %||% "")[1]))) {
    problems <- c(problems, "reason must say which evidence decided it")
  }
  problems
}
.arbiter_final <- function(base, decision) {
  selected <- as.character(decision$selected %||% "")[1]
  subtype <- as.character(decision$subtype %||% "")[1]
  keep <- unique(c(
    Filter(function(n) nzchar(n) && !identical(n, UNKNOWN), c(selected, subtype)),
    as.character(unlist(decision$co_occurring_identities %||% list()))
  ))
  out <- base
  out$selected <- selected
  out$subtype <- subtype
  lineage_value <- as.character(decision$lineage %||% "")[1]
  out$lineage <- if (nzchar(lineage_value)) lineage_value else selected
  out$co_occurring_identities <- as.list(as.character(unlist(
    decision$co_occurring_identities %||% list()
  )))
  state_value <- as.character(decision$state %||% "")[1]
  out$state <- if (nzchar(state_value)) state_value else as.character(base$state %||% "")
  confidence_value <- as.character(decision$confidence %||% "")[1]
  out$confidence <- if (nzchar(confidence_value)) {
    confidence_value
  } else {
    base_confidence <- as.character(base$confidence %||% "")[1]
    if (nzchar(base_confidence)) base_confidence else "low"
  }
  out$reason <- as.character(decision$reason %||% "")[1]
  out$claim_evidence <- Filter(function(item) {
    is.list(item) && as.character(item$identity %||% "")[1] %in% keep
  }, base$claim_evidence %||% list())
  out$support_markers <- base$support_markers %||% list()
  out
}
.arbitrate_many <- function(pool, clusters, conversation, api_url, api_key,
                           contested_of = NULL, audits_of = NULL, pair_tables_of = NULL,
                           turn_offset = 0L) {
  out <- setNames(vector("list", length(clusters)), clusters)
  prompts <- list()
  names_avail <- list()
  active <- character(0)
  for (cluster in clusters) {
    conv <- conversation[[cluster]]
    if (!length(conv$attempts)) {
      out[[cluster]] <- NULL
      next
    }
    arbiter_packet <- .arbiter_packet(
      pool, cluster, conv$attempts, conv$last_packet, conv$tier,
      contested = (contested_of %||% list())[[cluster]],
      audits = (audits_of %||% list())[[cluster]],
      pair_tables_list = (pair_tables_of %||% list())[[cluster]]
    )
    prompts[[cluster]] <- paste0(
      .context_header(pool, cluster), .render(.prompt("final_arbiter"), arbiter_packet)
    )
    names_avail[[cluster]] <- unique(as.character(unlist(arbiter_packet$candidates_available)))
    active <- c(active, cluster)
  }
  if (!length(active)) {
    return(out)
  }
  for (attempt in seq_len(SCHEMA_RETRIES + 1L)) {
    responses <- scma_cached_call_llm_many(
      lapply(active, function(cl) prompts[[cl]]), api_url, api_key,
      reasoning_effort = ARBITER_EFFORT, model = ARBITER_MODEL,
      trace_ids = lapply(active, function(cl) conversation[[cl]]$trace_id),
      turn_indexes = lapply(active, function(cl) as.integer(conversation[[cl]]$turn + turn_offset))
    )
    still <- character(0)
    for (i in seq_along(active)) {
      cluster <- active[i]
      content <- responses[[i]][[1]]
      parsed <- if (is.null(content)) NULL else scma_parse_json(content)
      problems <- .arbiter_problems(parsed, names_avail[[cluster]],
                                    .population_verdicts_by_name(conversation[[cluster]]$attempts))
      if (!length(problems)) {
        out[[cluster]] <- parsed
        next
      }
      if (attempt >= SCHEMA_RETRIES + 1L) {
        out[[cluster]] <- NULL
        next
      }
      prompts[[cluster]] <- paste0(
        prompts[[cluster]],
        sprintf(
          "\n\n# Retry %d: the decision was rejected -- %s. Return one valid decision object.\n",
          attempt, paste(head(problems, 4L), collapse = "; ")
        )
      )
      still <- c(still, cluster)
    }
    active <- still
    if (!length(active)) break
  }
  out
}
.answer_key <- function(final) {
  paste(c(
    as.character(final$selected %||% "")[1], as.character(final$subtype %||% "")[1],
    sort(as.character(unlist(final$co_occurring_identities %||% list())))
  ), collapse = "\x1f")
}
.audit_delivery_many <- function(pool, clusters, finals, tiers, programs, audit_template,
                                 headers, api_key, api_url, trace_ids, turn_of, caches) {
  claimed_of <- setNames(lapply(clusters, function(cl) {
    unique(vapply(.claimed_identities(finals[[cl]]), function(c) as.character(c[[2]]), ""))
  }), clusters)
  role_of <- setNames(lapply(clusters, function(cl) {
    roles <- list()
    for (c in .claimed_identities(finals[[cl]])) roles[[as.character(c[[2]])]] <- as.character(c[[1]])
    roles
  }), clusters)
  .run_phase <- function(names_of) {
    jobs <- list()
    for (cl in clusters) {
      for (name in names_of[[cl]] %||% character(0)) {
        need <- .rival_names(pool, cl, name, claimed_of[[cl]], tiers[[cl]])
        hit <- cached_audit(caches[[cl]], tiers[[cl]], name, need)
        if (!is.null(hit)) {
          finals[[cl]]$definer_audit[[name]] <<- c(hit, list(role = role_of[[cl]][[name]] %||% "contested_rival"))
          next
        }
        jobs[[length(jobs) + 1L]] <- list(
          cluster = cl, identity = name, role = role_of[[cl]][[name]] %||% "contested_rival",
          claimed = claimed_of[[cl]], tier = tiers[[cl]]
        )
      }
    }
    if (!length(jobs)) return(invisible(NULL))
    results <- audit_many(pool, jobs, NULL, programs, audit_template, headers, api_key,
                          api_url, trace_ids, turn_of, SCHEMA_RETRIES)
    for (job in jobs) {
      key <- paste(job$cluster, job$identity, sep = "\x1f")
      audit <- results[[key]]
      if (is.null(audit)) next
      finals[[job$cluster]]$definer_audit[[job$identity]] <<- audit
      key2 <- paste(as.integer(job$tier), job$identity, sep = "\x1f")
      assign(key2, audit, envir = caches[[job$cluster]])
    }
    invisible(NULL)
  }
  for (cl in clusters) if (is.null(finals[[cl]]$definer_audit)) finals[[cl]]$definer_audit <- list()
  .run_phase(claimed_of)
  contested_of <- setNames(lapply(clusters, function(cl) {
    contested_rival_names(finals[[cl]], finals[[cl]]$definer_audit, claimed_of[[cl]])
  }), clusters)
  for (cl in clusters) finals[[cl]]$contested_rivals <- contested_of[[cl]]
  out <- setNames(lapply(clusters, function(cl) {
    delivery_contested(finals[[cl]], finals[[cl]]$definer_audit, contested_of[[cl]])
  }), clusters)
  list(finals = finals, contested = out)
}
.escalate <- function(conv, pool, cluster, reason) {
  if (conv$tier >= conv$tiers) {
    return(list(conv = conv, escalated = FALSE))
  }
  conv$tier <- as.integer(conv$tier + 1L)
  conv$rounds_in_tier <- 0L
  block <- .tier_packet(pool, cluster, conv$tier)
  block$note <- sprintf(
    paste0(
      "%s Candidates ranked %d and below in the retrieval order are added here, with ",
      "their whole measured panels and facts. Everything already on this page is still ",
      "available: the earlier candidates were not withdrawn, and one of them may still be ",
      "the answer. Answer again."
    ),
    reason, as.integer((conv$tier - 1L) * CANDIDATE_TIER_SIZE + 1L)
  )
  conv$prompt <- paste0(
    conv$prompt, .observation_block(length(conv$transcript) + length(conv$attempts) + 1L, block)
  )
  list(conv = conv, escalated = TRUE)
}
.deliver <- function(conv, final, arbitration, contested = character(0)) {
  capped <- NULL
  cap_why <- ""
  if (!identical(as.character(final$selected %||% "")[1], UNKNOWN)) {
    cap_res <- confidence_cap(final, final$definer_audit %||% list())
    step <- apply_cap(final, cap_res$cap)
    final <- step$final
    if (isTRUE(step$applied)) {
      capped <- cap_res$cap
      cap_why <- cap_res$why
    }
  }
  last_item <- if (length(conv$attempts)) conv$attempts[[length(conv$attempts)]] else NULL
  if (!is.null(last_item) && !is.null(last_item$review) &&
      identical(as.character(final$selected %||% "")[1], as.character(last_item$final$selected %||% "")[1])) {
    minority <- character(0)
    for (e in last_item$review$detected_exclusions_disposition %||% list()) {
      if (is.list(e) && identical(as.character(e$disposition %||% ""), "minority")) {
        minority <- c(minority, as.character(e$gene %||% "")[1])
      }
    }
    if (length(minority)) {
      step <- apply_cap(final, "medium")
      final <- step$final
      if (isTRUE(step$applied)) {
        capped <- "medium"
        cap_why <- sprintf(
          "the review reads %s as detected exclusions of '%s' that are real in a minority of these cells",
          .py_list(head(minority, 4L)), as.character(final$selected %||% "")[1]
        )
      }
    }
  }
  if (isTRUE(conv$review_unreachable %||% FALSE)) {
    final <- apply_cap(final, "low")$final
  }
  if (length(contested)) {
    final <- apply_cap(final, "low")$final
  }
  agent_final <- if (length(conv$attempts)) {
    conv$attempts[[length(conv$attempts)]]$final
  } else {
    final
  }
  rounds <- lapply(seq_along(conv$attempts), function(index) {
    item <- conv$attempts[[index]]
    list(
      round = index, tier = as.integer(item$tier),
      selected = as.character(item$final$selected %||% "")[1],
      subtype = as.character(item$final$subtype %||% "")[1],
      verdict = item$review
    )
  })
  checked <- length(conv$attempts) > 0L &&
    any(vapply(conv$attempts, function(item) length(item$review) > 0L, TRUE))
  last <- if (length(conv$attempts)) conv$attempts[[length(conv$attempts)]] else NULL
  passed <- !is.null(last) && !is.null(last$review) && !isTRUE(last$review$uncited) &&
    isTRUE(.review_passed(last$review, as.character(last$final$subtype %||% "")[1])) &&
    is.null(arbitration) && !nzchar(as.character(last$structural_note %||% ""))
  list(
    final = final,
    agent_selected = as.character(agent_final$selected %||% "")[1],
    arbitration = arbitration,
    turns = conv$turn,
    transcript = conv$transcript,
    trace_id = conv$trace_id,
    turn_budget_exhausted = conv$forced,
    delivered_tier = conv$tier,
    tiers_available = conv$tiers,
    program_reading = conv$program %||% list(),
    definer_audit = final$definer_audit %||% list(),
    contested_rivals = as.character(unlist(final$contested_rivals %||% list())),
    contested = as.character(contested),
    confidence_cap = if (!is.null(capped)) list(cap = capped, why = cap_why) else NULL,
    review_unreachable = isTRUE(conv$review_unreachable %||% FALSE),
    structural_notes = conv$structural_notes %||% list(),
    consistency_note = conv$consistency_note %||% NULL,
    review = list(rounds = rounds, checked = isTRUE(checked), passed = isTRUE(passed))
  )
}
.measurements <- function(entry, cluster, limit = 30L) {
  positives <- list()
  negatives <- list()
  for (marker in entry$markers) {
    pct_in <- marker$pct_in %||% 0
    if (identical(marker$polarity, "positive")) {
      if (!.is_raised(marker)) next
      if (length(positives) >= limit) next
    } else if (identical(marker$polarity, "negative")) {
      if (!is.finite(pct_in) || pct_in < NEGATIVE_SOURCE_MIN_PCT_IN) next
    } else {
      next
    }
    gene <- toupper(marker$gene)
    row <- list(
      gene = marker$gene,
      polarity = marker$polarity,
      detection_fraction_in = .round(pct_in / 100),
      detection_fraction_out = .round(marker$pct_out / 100),
      avg_log2FC = marker$avg_log2FC,
      auc = marker$auc,
      cross_cluster_percentile = .round(
        if (gene %in% rownames(rel_pct)) rel_pct[gene, cluster] else NA_real_
      ),
      publication_support = marker$n_pub,
      evidence_tier = marker$tier
    )
    if (identical(marker$polarity, "positive")) {
      positives[[length(positives) + 1L]] <- row
    } else {
      negatives[[length(negatives) + 1L]] <- row
    }
  }
  c(positives, negatives)
}
.candidate_entry <- function(entry, cluster, role, evidence, warning) {
  list(
    cell_type = entry$cell_type,
    retrieval_rank = entry$retrieval_rank,
    tier = .tier_of(entry),
    borrowed_context = entry$borrowed_context,
    tissue_context = entry$tissue_context %||% list(),
    retrieval_context = entry$retrieval_context %||% "",
    claim_role = role,
    claim_evidence = evidence %||% structure(list(), names = character(0)),
    claim_warning = warning,
    panel = entry$markers,
    unmeasured_curated_genes = entry$unmeasured_curated_genes,
    single_cell_program = list(
      in_cluster_median = entry$program$median_in,
      out_of_cluster_median = entry$program$median_out
    ),
    decisive_marker_measurements = .measurements(entry, cluster)
  )
}
.plain_entries <- function(cluster) {
  lapply(pool$clusters[[cluster]]$candidates, function(entry) {
    .candidate_entry(entry, cluster, "", NULL, "")
  })
}
.fallback <- function(cluster, frame, reason, detail) {
  if (!nrow(frame)) {
    return(list(
      annotation_qc = QC_UNCHECKED, review = list(),
      cluster_id = cluster, annotation = UNKNOWN, subtype = "", lineage = "", state = "",
      co_occurring_identities = list(), confidence = NOT_AVAILABLE,
      rationale = paste(
        "no candidate's curated positive markers are significantly",
        "up-regulated in this cluster"
      ),
      resolution_status = UNRESOLVED, resolution_detail = "empty_candidate_pool",
      annotation_source = "no_candidate", llm_status = reason,
      support_markers = list(), claim_warnings = list()
    ))
  }
  top <- frame[order(retrieval_rank)][1]
  list(
    annotation_qc = QC_UNCHECKED, review = list(),
    cluster_id = cluster, annotation = as.character(top$candidate),
    subtype = "", lineage = as.character(top$candidate), state = "",
    co_occurring_identities = list(), confidence = NOT_AVAILABLE,
    rationale = sprintf(
      paste(
        "no model judgement was available; this candidate leads the joint retrieval",
        "order over marker-level, cluster-level and single-cell-level evidence with",
        "%d of %d measured positive markers significantly up-regulated here.",
        "The retrieval order is not a confidence."
      ),
      as.integer(top$hits), as.integer(top$panel_size)
    ),
    resolution_status = RESOLVED, resolution_detail = detail,
    annotation_source = "relative_score_fallback", llm_status = reason,
    support_markers = list(), claim_warnings = list()
  )
}
.head_chars <- function(value, limit) {
  text <- as.character(value %||% "")[1]
  if (is.na(text)) text <- ""
  substr(text, 1L, as.integer(limit))
}
.qc_status <- function(outcome) {
  review <- outcome$review %||% list()
  rounds <- review$rounds %||% list()
  if (!isTRUE(review$checked)) {
    return(QC_UNCHECKED)
  }
  if (!is.null(outcome$arbitration) && length(outcome$arbitration)) {
    return(QC_ARBITRATED)
  }
  if (!isTRUE(review$passed)) {
    return(QC_FAILED)
  }
  if (length(rounds) > 1L) QC_REVISED else QC_PASSED
}
.delivered_facts <- function(pool, cluster, names_wanted) {
  precomputed <- (pool$facts %||% list())[[cluster]] %||% list()
  out <- list()
  for (name in unique(names_wanted)) {
    if (!nzchar(name) || identical(name, UNKNOWN)) next
    row <- NULL
    for (candidate_row in precomputed) {
      if (identical(as.character(candidate_row$cell_type), name)) {
        row <- candidate_row
        break
      }
    }
    if (is.null(row)) {
      entry <- .find_candidate(pool, cluster, name)
      if (!is.null(entry)) row <- .candidate_facts(pool, cluster, entry)
    }
    if (!is.null(row)) out[[name]] <- row
  }
  out
}
.result <- function(cluster, pool, outcome) {
  names_available <- as.character(.candidate_names(pool, cluster))
  final <- outcome$final
  selected <- as.character(final$selected)[1]
  subtype <- as.character(final$subtype %||% "")[1]
  others <- as.character(unlist(final$co_occurring_identities %||% list()))
  lineage_value <- as.character(final$lineage %||% "")[1]
  lineage <- if (nzchar(lineage_value)) lineage_value else selected
  claims <- .claimed_identities(final)
  warned <- .claim_warnings(pool, cluster, claims)
  arbitration <- outcome$arbitration
  if (!is.null(arbitration) && length(arbitration)) {
    agent <- as.character(outcome$agent_selected %||% "")[1]
    if (nzchar(agent) && !identical(agent, selected)) {
      warned[[length(warned) + 1L]] <- list(agent, sprintf(
        paste0(
          "arbitration %s -> %s: the evidence review refused the annotator's answer and ",
          "the arbitration turn delivered this one | %s"
        ),
        agent, selected, .head_chars(arbitration$reason, 600L)
      ))
    }
  }
  warning_of <- list()
  for (item in warned) {
    name <- as.character(item[[1]])
    line <- as.character(item[[2]])
    warning_of[[name]] <- if (is.null(warning_of[[name]])) {
      line
    } else {
      paste(warning_of[[name]], line, sep = " || ")
    }
  }
  warnings <- if (length(warned)) vapply(warned, function(item) as.character(item[[2]]), "") else character(0)
  role_of <- list()
  for (claim in claims) role_of[[as.character(claim[[2]])]] <- as.character(claim[[1]])
  evidence_of <- list()
  for (item in final$claim_evidence %||% list()) {
    if (is.list(item)) evidence_of[[as.character(item$identity %||% "")]] <- item
  }
  entries <- lapply(pool$clusters[[cluster]]$candidates, function(entry) {
    .candidate_entry(
      entry, cluster,
      role_of[[entry$cell_type]] %||% "",
      evidence_of[[entry$cell_type]],
      warning_of[[entry$cell_type]] %||% ""
    )
  })
  contested <- as.character(outcome$contested %||% character(0))
  if (identical(selected, UNKNOWN)) {
    status <- UNRESOLVED
    detail <- "no_candidate_carried_by_its_evidence"
  } else if (isTRUE(outcome$review_unreachable %||% FALSE)) {
    status <- UNRESOLVED
    detail <- "evidence_review_unreachable"
  } else if (length(contested)) {
    status <- UNRESOLVED
    detail <- "contested_after_rounds"
  } else if (length(others)) {
    status <- MIXED
    detail <- "agent_majority_of_several_identities"
  } else {
    status <- RESOLVED
    detail <- "agent_selected"
  }
  cap <- outcome$confidence_cap %||% NULL
  if (!is.null(cap)) {
    warnings <- c(warnings, sprintf("confidence capped at %s: %s", cap$cap, cap$why))
  }
  for (note in as.character(unlist(outcome$structural_notes %||% list()))) {
    warnings <- c(warnings, sprintf("structural: %s", note))
  }
  review <- outcome$review %||% list()
  review_flags <- list()
  for (item in review$rounds %||% list()) {
    for (text in (item$verdict %||% list())$flags %||% list()) {
      text_value <- trimws(as.character(text %||% "")[1])
      if (nzchar(text_value)) {
        review_flags[[length(review_flags) + 1L]] <- list(
          round = item$round, label = item$selected, text = .head_chars(text_value, 300L)
        )
      }
    }
  }
  list(
    program_reading = outcome$program_reading %||% list(),
    definer_audit = outcome$definer_audit %||% list(),
    contested = as.list(contested),
    contested_rivals = as.list(as.character(unlist(outcome$contested_rivals %||% list()))),
    confidence_cap = cap,
    structural_notes = as.list(outcome$structural_notes %||% list()),
    consistency_note = outcome$consistency_note %||% NULL,
    cluster_id = cluster,
    annotation = selected,
    subtype = subtype,
    lineage = lineage,
    agent_selected = as.character(outcome$agent_selected %||% final$selected)[1],
    arbitration = arbitration,
    possible_refinements = as.list(sort(.cp(unique(Filter(
      nzchar, as.character(unlist(final$possible_refinements %||% list()))
    ))), method = "radix")),
    candidate_facts = .delivered_facts(pool, cluster, c(selected, subtype)),
    borrowed_context = (pool$borrowed %||% list())[[cluster]] %||% list(),
    delivered_borrowed = if (identical(selected, UNKNOWN)) {
      FALSE
    } else {
      !is.null((.find_candidate(
        pool, cluster, if (nzchar(subtype)) subtype else selected
      ) %||% list())$borrowed_context)
    },
    state = as.character(final$state %||% ""),
    co_occurring_identities = as.list(others),
    confidence = as.character(final$confidence),
    rationale = as.character(final$reason %||% ""),
    resolution_status = status,
    resolution_detail = detail,
    annotation_source = "cluster_annotation",
    llm_status = "annotated",
    support_markers = as.list(Filter(
      nzchar, as.character(unlist(final$support_markers %||% list()))
    )),
    claim_warnings = as.list(warnings),
    candidates = names_available,
    candidate_entries = entries,
    delivered_tier = as.integer(outcome$delivered_tier %||% 1L),
    tiers_available = as.integer(outcome$tiers_available %||% 1L),
    turns = as.integer(outcome$turns %||% 0L),
    tool_calls = lapply(outcome$transcript %||% list(), function(item) item$call),
    turn_budget_exhausted = isTRUE(outcome$turn_budget_exhausted),
    annotation_qc = .qc_status(outcome),
    review = review,
    review_flags = review_flags
  )
}
scoring <- readRDS(file.path(CACHE, sprintf("%s_candidate_scoring.rds", tag)))
cc <- readRDS(file.path(CACHE, sprintf("%s_de_meta.rds", tag)))
de <- as.data.table(cc$de)
de[, group := as.character(group)]
de[, gene_key := toupper(as.character(feature))]
setkey(de, group, gene_key)
native_of <- setNames(as.character(cc$menu_genes), toupper(as.character(cc$menu_genes)))
clusters_sorted <- .sorted_cluster_ids(scoring)
context <- scoring$context
scored <- as.data.table(scoring$scored)
wide <- dcast(de, gene_key ~ group, value.var = "avg_log2FC")
rel_genes <- wide$gene_key
cluster_columns <- setdiff(names(wide), "gene_key")
rel <- as.matrix(wide[, ..cluster_columns])
rel_pct <- t(matrix(
  apply(rel, 1, function(row) frank(row, ties.method = "average") / length(row)),
  nrow = length(cluster_columns), ncol = length(rel_genes)
))
dimnames(rel_pct) <- list(rel_genes, cluster_columns)
UB_SOURCES <- uberon_load(OBO_UBERON)
source_context <- scma_source_context(
  UB_SOURCES, DB_SOURCES, context$species, context$tissue, context$disease
)
.fetch_from <- function(context) {
  function(candidate, gene, k) scma_source_head(context, candidate, gene, k)
}
.sources_for <- .fetch_from(source_context)
.screen_problems <- function(value, names) {
  if (!is.list(value) || (!is.list(value$screen) && !is.character(value$screen))) {
    return("return one JSON object with a `screen` list")
  }
  entries <- value$screen
  problems <- character(0)
  seen <- character(0)
  for (entry in entries) {
    if (!is.list(entry)) {
      problems <- c(problems, "each screen entry must be an object")
      next
    }
    name <- as.character(entry$name %||% "")[1]
    if (!(name %in% names)) {
      problems <- c(problems, sprintf("'%s' is not one of the names given", name))
      next
    }
    seen <- c(seen, name)
    if (!(as.character(entry$verdict %||% "")[1] %in% SCREEN_VERDICTS)) {
      problems <- c(problems, sprintf("%s: verdict must be one of %s", name, .py_list(SCREEN_VERDICTS)))
    }
  }
  missing <- setdiff(names, seen)
  if (length(missing)) {
    problems <- c(problems, sprintf("every name needs a verdict; missing %s", .py_list(head(missing, 6L))))
  }
  problems
}
.make_borrow_screen <- function(context, api_key, api_url) {
  record <- new.env(parent = emptyenv())
  record$model <- scma_resolve_model(BORROW_SCREEN_MODEL)
  record$reasoning_effort <- BORROW_SCREEN_EFFORT
  record$status <- "not_run"
  record$names <- list()
  record$foreign <- list()
  record$verdicts <- list()
  template <- paste(readLines(file.path(PROMPT_DIR, "borrow_tissue_screen.txt"), warn = FALSE), collapse = "\n")
  disease <- context$disease
  disease <- if (is.list(disease) || length(disease) > 1L) paste(as.character(unlist(disease)), collapse = ", ") else as.character(disease)
  tissue <- context$tissue
  tissue <- if (is.list(tissue) || length(tissue) > 1L) paste(as.character(unlist(tissue)), collapse = ", ") else as.character(tissue)
  dataset_context <- list(
    species = as.character(context$species), tissue = tissue, disease = disease,
    development_stage = as.character(context$development_stage %||% "")
  )
  screen_fn <- function(names) {
    names <- as.character(names)
    record$names <- as.list(names)
    if (!length(names)) {
      record$status <- "nothing_to_screen"
      return(list())
    }
    prompt <- .render(template, list(dataset_context = dataset_context, names = as.list(names)))
    trace_id <- substr(paste0(format(Sys.time(), "%H%M%S"), basename(tempfile(""))), 1, 16)
    parsed <- NULL
    for (attempt in seq_len(SCHEMA_RETRIES + 1L)) {
      pair <- scma_cached_call_llm(
        prompt, api_url, api_key, reasoning_effort = BORROW_SCREEN_EFFORT, model = BORROW_SCREEN_MODEL,
        trace_id = trace_id, turn_index = as.integer(attempt - 1L)
      )
      content <- pair[[1]]
      candidate <- if (is.null(content)) NULL else scma_parse_json(content)
      problems <- .screen_problems(candidate, names)
      if (!length(problems)) {
        parsed <- candidate
        break
      }
      if (attempt >= SCHEMA_RETRIES + 1L) {
        record$status <- sprintf("unavailable: %s", paste(head(problems, 3L), collapse = "; "))
        break
      }
      prompt <- paste0(prompt, sprintf(
        "\n\n# Retry %d: the answer was rejected -- %s. Return one valid screen object.\n",
        attempt, paste(head(problems, 4L), collapse = "; ")
      ))
    }
    if (is.null(parsed)) return(list())
    verdicts <- list()
    for (entry in parsed$screen) {
      verdicts[[as.character(entry$name)]] <- list(
        verdict = as.character(entry$verdict), reason = as.character(entry$reason %||% "")
      )
    }
    record$status <- "ok"
    record$verdicts <- verdicts
    record$foreign <- as.list(sort(.cp(Filter(
      function(n) identical(as.character(verdicts[[n]]$verdict), SCREEN_FOREIGN), names(verdicts)
    )), method = "radix"))
    verdicts
  }
  list(screen = screen_fn, record = record)
}
scma_prepare_pool <- function(scoring, de, native_of, markers_all, source_context,
                             prep_genes, quiet = FALSE, borrow_screen = NULL) {
  context <- scoring$context
  .say <- function(...) if (!isTRUE(quiet)) cat(...)
  native_sources <- source_context
  pool <- .build_pool(scoring, de, native_of, .fetch_from(native_sources),
                      markers_all = markers_all, prep_genes = prep_genes)
  if (length(pool$tissue_contexts %||% list()) > 1L) {
    by_retrieval <- table(as.character(unlist(pool$candidate_retrieval_context)))
    .say(sprintf(
      "  tissue contexts (operator list): %s; %s; eligible candidates by the context that admits them: %s\n",
      paste(as.character(unlist(pool$tissue_contexts)), collapse = ", "),
      as.character(pool$retrieval_semantics %||% ""),
      paste(sprintf("%s=%d", names(by_retrieval), as.integer(by_retrieval)), collapse = ", ")
    ))
  }
  pool$genome_de <- if (!is.null(markers_all)) {
    scma_genome_de_summary(markers_all)
  } else {
    .say("  [warn] markers_all not found; the evidence review runs without genome-wide DE and nothing is borrowed\n")
    list()
  }
  dbf <- scma_load_marker_db(DB_CSV)
  sources <- native_sources
  pool$borrowed <- list()
  donor_species_by_name <- list()
  if (!is.null(markers_all)) {
    measured <- unique(toupper(as.character(prep_genes)))
    resource <- scma_borrow_resource_slice(dbf, context$species, context$disease)
    cross_species_args <- list()
    if (isTRUE(CROSS_SPECIES_BORROW)) {
      other_species <- setdiff(ORTHO_SPECIES, as.character(context$species))
      cross_species_resource <- setNames(
        lapply(other_species, function(sp) scma_borrow_resource_slice(dbf, sp, context$disease)),
        other_species
      )
      cross_species_args <- list(
        species = as.character(context$species),
        cross_species_resource = cross_species_resource,
        ortho_dir = ORTHO_DIR
      )
    }
    borrowed <- do.call(scma_borrow_candidates, c(
      list(pool, resource, measured, markers_all, sources = NULL, screen = borrow_screen),
      cross_species_args
    ))
    pool <- borrowed$pool
    pool$borrowed <- borrowed$events
    screened_out <- sort(.cp(unique(unlist(lapply(pool$borrowed, function(log) {
      vapply(Filter(function(item) startsWith(as.character(item$rejected %||% ""), SCREEN_REJECTED_PREFIX),
                    log$considered %||% list()), function(item) as.character(item$cell_type), "")
    })))), method = "radix")
    if (length(screened_out)) {
      .say(sprintf(
        "  borrowed-name tissue screen: %d name(s) read foreign to %s and not admitted (%s)\n",
        length(screened_out), paste(as.character(unlist(context$tissue)), collapse = ", "),
        paste(screened_out, collapse = ", ")
      ))
    }
    names_borrowed <- scma_borrowed_names(pool$borrowed)
    donor_species_by_name <- scma_borrowed_donor_species(pool$borrowed)
    if (length(names_borrowed)) {
      extra <- scma_source_types_across_tissues(
        DB_SOURCES, context$species, names_borrowed, context$disease,
        donor_species_by_name = donor_species_by_name, ortho_dir = ORTHO_DIR
      )
      sources <- scma_composite_sources(native_sources, extra, names_borrowed)
      fetch <- .fetch_from(sources)
      for (cluster in names(pool$clusters)) {
        entries <- pool$clusters[[cluster]]$candidates
        for (index in seq_along(entries)) {
          if (!is.null(entries[[index]]$borrowed_context)) {
            entries[[index]]$exclusion_sources <- .exclusion_sources(
              entries[[index]]$cell_type, entries[[index]]$markers, fetch
            )
          }
        }
        pool$clusters[[cluster]]$candidates <- entries
      }
    }
    n_events <- sum(vapply(pool$borrowed, function(v) length(v$borrowed %||% list()), 0L))
    .say(sprintf(
      "  borrowed-context retrieval: %d candidates borrowed in %d clusters (%s)\n",
      n_events,
      sum(vapply(pool$borrowed, function(v) length(v$borrowed %||% list()) > 0L, TRUE)),
      if (length(names_borrowed)) paste(names_borrowed, collapse = ", ") else "none"
    ))
  }
  wanted_species <- as.character(context$species)
  recf <- dbf[dbf$species == wanted_species & dbf$marker_polarity == "positive" &
    toupper(as.character(dbf$is_recommended_marker)) == "TRUE"]
  rec_pairs <- unique(paste(
    as.character(recf$cell_type), toupper(as.character(recf$gene_symbol)), sep = "\x1f"
  ))
  attached <- scma_attach_recommended(pool, rec_pairs)
  pool <- attached$pool
  species_rows <- scma_species_positive_rows(dbf, wanted_species)
  all_names <- sort(.cp(unique(unlist(lapply(pool$clusters, function(s) {
    vapply(s$candidates, function(e) as.character(e$cell_type), "")
  })))), method = "radix")
  definers <- scma_species_definers(species_rows, all_names, k = DEFINERS_PER_CANDIDATE)
  res_def <- attach_definers(pool, definers)
  pool <- res_def$pool
  species_source_ctx <- scma_composite_sources(
    native_sources, scma_source_types_across_tissues(
      DB_SOURCES, wanted_species, all_names, context$disease,
      donor_species_by_name = donor_species_by_name, ortho_dir = ORTHO_DIR
    ),
    all_names
  )
  pool$species_sources <- .fetch_from(species_source_ctx)
  pool$species_claimants <- scma_species_gene_claimants(species_rows)
  res_share <- attach_page_sharing(pool, pool$species_claimants)
  pool <- res_share$pool
  pool <- attach_genome_claims(pool, pool$species_claimants, scma_known_genes(dbf), markers_all)
  .say(sprintf(
    "  evidence: species-wide definers on %d candidate entries (%d names); %d rows marked as shared with another page candidate\n",
    res_def$n, length(definers), res_share$n
  ))
  pool <- scma_assign_tiers(pool, scoring)
  pool$facts <- setNames(
    lapply(names(pool$clusters), function(cluster) .candidate_facts_table(pool, cluster)),
    names(pool$clusters)
  )
  per_tier <- list()
  for (cluster in names(pool$clusters)) {
    for (entry in pool$clusters[[cluster]]$candidates) {
      key <- as.character(.tier_of(entry))
      per_tier[[key]] <- (per_tier[[key]] %||% 0L) + 1L
    }
  }
  tier_keys <- names(per_tier)
  tier_order <- order(as.integer(tier_keys))
  .say(sprintf(
    "  facts layer: %d recommended (cell type, gene) pairs for %s; %d panel rows marked; candidates by delivery tier %s\n",
    length(rec_pairs), wanted_species, attached$marked,
    paste(sprintf(
      "%s=%d", tier_keys[tier_order], as.integer(unlist(per_tier))[tier_order]
    ), collapse = ", ")
  ))
  list(pool = pool, sources = sources)
}
markers_all_path <- file.path(CACHE, sprintf("%s_markers_all.rds", tag))
markers_all <- if (file.exists(markers_all_path)) readRDS(markers_all_path) else NULL
credentials <- .resolve_credentials()
api_key <- credentials$key
api_url <- credentials$url
llm_status <- if (!isTRUE(CLUSTER_ANNOTATION_ENABLED)) {
  "disabled"
} else if (!nzchar(api_key) || !nzchar(api_url)) {
  "skipped_no_credentials"
} else {
  "enabled"
}
borrow_screen_fn <- NULL
borrow_screen_record <- list(status = "not_run")
if (identical(llm_status, "enabled")) {
  screen_bundle <- .make_borrow_screen(scoring$context, api_key, api_url)
  borrow_screen_fn <- screen_bundle$screen
  borrow_screen_record_env <- screen_bundle$record
}
.prepared <- scma_prepare_pool(
  scoring, de, native_of, markers_all, source_context, cc$genes, borrow_screen = borrow_screen_fn
)
pool <- .prepared$pool
sources_context <- .prepared$sources
if (exists("borrow_screen_record_env")) {
  borrow_screen_record <- as.list(borrow_screen_record_env)
}
pool$borrow_screen <- borrow_screen_record
cat(sprintf(
  "==== annotate(r) %s: %d clusters, up to %d candidates each in %d tiers of %d, llm=%s ====\n",
  tag, length(clusters_sorted), as.integer(scoring$top_candidates), CANDIDATE_TIERS,
  CANDIDATE_TIER_SIZE, llm_status
))
cat(sprintf(
  "  annotator: %s (%s); review: %s (%s), %d rounds per tier; arbiter: %s (%s); borrow screen: %s (%s), %s\n",
  scma_resolve_model(ANNOTATOR_MODEL), ANNOTATOR_EFFORT,
  scma_resolve_model(REVIEW_MODEL), REVIEW_EFFORT, REVIEW_MAX_ROUNDS_PER_TIER,
  scma_resolve_model(ARBITER_MODEL), ARBITER_EFFORT,
  scma_resolve_model(BORROW_SCREEN_MODEL), BORROW_SCREEN_EFFORT, borrow_screen_record$status %||% "not_run"
))
template <- .prompt("cluster_annotator")
.subset_raw <- trimws(Sys.getenv("SCMA_CLUSTER_SUBSET", ""))
cluster_subset <- if (nzchar(.subset_raw)) {
  Filter(nzchar, trimws(strsplit(.subset_raw, ",", fixed = TRUE)[[1]]))
} else {
  NULL
}
runnable <- character(0)
results <- list()
for (cluster in clusters_sorted) {
  cluster_value <- cluster
  frame <- scored[cluster == cluster_value]
  if (identical(pool$clusters[[cluster]]$status, UNSUPPORTED)) {
    results[[cluster]] <- .fallback(cluster, frame[0], llm_status, "empty_candidate_pool")
    next
  }
  if (!is.null(cluster_subset) && !(as.character(cluster) %in% cluster_subset)) {
    results[[cluster]] <- .fallback(cluster, frame, "not_run_in_subset", "subset_skipped")
    next
  }
  if (!identical(llm_status, "enabled")) {
    record <- .fallback(cluster, frame, llm_status, "relative_score_fallback_top1")
    record$candidates <- as.character(.candidate_names(pool, cluster))
    record$candidate_entries <- .plain_entries(cluster)
    results[[cluster]] <- record
    next
  }
  runnable <- c(runnable, cluster)
}
.open_conversations <- function(pool, cluster_ids, sources_context, api_key, api_url,
                                headers_of = NULL, consistency_notes = list()) {
  if (!length(cluster_ids)) {
    return(list(conversation = list(), source_servers = list(), headers = list()))
  }
  source_servers <- list()
  for (cluster in cluster_ids) {
    server <- scma_source_server(
      sources_context,
      batch = SOURCES_PER_MARKER, max_batches = SOURCE_BATCHES_PER_MARKER
    )
    scma_register_packet_sources(server, pool, cluster)
    source_servers[[cluster]] <- server
  }
  headers <- headers_of %||% setNames(
    lapply(cluster_ids, function(cl) .context_header(pool, cl)), cluster_ids
  )
  trace_ids <- setNames(lapply(cluster_ids, function(cl) {
    substr(paste0(format(Sys.time(), "%H%M%S"), basename(tempfile(""))), 1, 16)
  }), cluster_ids)
  programs <- program_reading_many(
    pool, cluster_ids, setNames(as.list(rep(1L, length(cluster_ids))), cluster_ids),
    .prompt("program_reading"), headers, api_key, api_url, trace_ids, 0L, SCHEMA_RETRIES
  )
  conversation <- list()
  for (cluster in cluster_ids) {
    opening <- .cluster_packet(pool, cluster, tier = 1L)
    opening$program_reading <- programs[[cluster]] %||% list()
    note <- consistency_notes[[cluster]]
    if (!is.null(note)) {
      opening$dataset_consistency_note <- note
    }
    conversation[[cluster]] <- list(
      prompt = paste0(headers[[cluster]], .render(template, opening)),
      trace_id = trace_ids[[cluster]],
      turn = 0L, schema_failures = 0L, repeated = 0L, forced = FALSE,
      seen = list(), transcript = list(), outcome = NULL,
      tier = 1L, tiers = .tiers_available(pool, cluster), rounds_in_tier = 0L,
      attempts = list(), last_packet = list(),
      program = programs[[cluster]] %||% list(),
      audit_cache = new.env(parent = emptyenv()),
      review_unreachable = FALSE, structural_notes = list(), consistency_note = note
    )
  }
  list(conversation = conversation, source_servers = source_servers, headers = headers)
}
audit_template <- .prompt("definer_audit")
.opened <- .open_conversations(pool, runnable, sources_context, api_key, api_url)
conversation <- .opened$conversation
source_servers <- .opened$source_servers
headers <- .opened$headers
.run_conversation_rounds <- function(pool, conversation, source_servers, headers,
                                     api_key, api_url) {
active <- names(conversation)
while (length(active)) {
  prompts <- vapply(active, function(cluster) conversation[[cluster]]$prompt, "")
  turns <- vapply(active, function(cluster) conversation[[cluster]]$turn + 1L, 0L)
  traces <- vapply(active, function(cluster) conversation[[cluster]]$trace_id, "")
  responses <- scma_cached_call_llm_many(
    as.list(unname(prompts)), api_url, api_key,
    reasoning_effort = ANNOTATOR_EFFORT, model = ANNOTATOR_MODEL,
    trace_ids = as.list(unname(traces)), turn_indexes = as.list(as.integer(unname(turns)))
  )
  still <- character(0)
  to_review <- character(0)
  redelivered <- character(0)
  pending_final <- list()
  for (position in seq_along(active)) {
    cluster <- active[position]
    conv <- conversation[[cluster]]
    conv$turn <- conv$turn + 1L
    pair <- responses[[position]]
    content <- pair[[1]]
    error <- pair[[2]]
    if (is.null(content)) {
      if (identical(error %||% "", "cache_miss_no_credentials")) {
        conv$outcome <- list(error = error, turns = conv$turn, transcript = conv$transcript)
        conversation[[cluster]] <- conv
        next
      }
      conv$schema_failures <- conv$schema_failures + 1L
      if (conv$schema_failures > SCHEMA_RETRIES) {
        conv$outcome <- list(
          error = error %||% "no response", turns = conv$turn, transcript = conv$transcript
        )
        conversation[[cluster]] <- conv
        next
      }
      conv$prompt <- paste0(conv$prompt, sprintf("\n\n# Retry %d\n", conv$schema_failures))
      conversation[[cluster]] <- conv
      still <- c(still, cluster)
      next
    }
    parsed <- scma_parse_json(content)
    if (!is.null(parsed)) {
      sanitized <- .sanitize_support_markers(parsed, pool, cluster)
      parsed <- sanitized$value
      if (length(sanitized$dropped)) {
        message(sprintf(
          "  [warn] cluster %s: dropped support_markers not on %s panel: %s",
          cluster, as.character(parsed$selected %||% "")[1],
          paste(sanitized$dropped, collapse = ", ")
        ))
      }
    }
    if (!is.null(parsed) && .valid_final(parsed, pool, cluster, conv$tier, conv$program)) {
      conv$schema_failures <- 0L
      if (identical(as.character(parsed$selected %||% "")[1], UNKNOWN)) {
        step <- .escalate(
          conv, pool, cluster, "You reported that no candidate on the page is carried."
        )
        conv <- step$conv
        if (isTRUE(step$escalated)) {
        conversation[[cluster]] <- conv
        still <- c(still, cluster)
        next
      }
        conv$outcome <- .deliver(conv, parsed, NULL)
        conversation[[cluster]] <- conv
        next
      }
      if (length(conv$attempts) && identical(.answer_key(parsed), conv$last_refused_answer %||% "")) {
        conversation[[cluster]] <- conv
        redelivered <- c(redelivered, cluster)
        pending_final[[cluster]] <- parsed
        next
      }
      conversation[[cluster]] <- conv
      to_review <- c(to_review, cluster)
      pending_final[[cluster]] <- parsed
      next
    }
    if (!is.null(parsed) && .valid_query(parsed) && !isTRUE(conv$forced)) {
      conv$schema_failures <- 0L
      call <- .compact_json(list(
        tool = as.character(parsed$tool)[1], args = parsed$args %||% list()
      ))
      streaming <- identical(as.character(parsed$tool)[1], "sources")
      if (!is.null(conv$seen[[call]]) && !streaming) {
        conv$repeated <- conv$repeated + 1L
        observation <- c(conv$seen[[call]], list(duplicate_query = TRUE))
      } else {
        conv$repeated <- 0L
        observation <- .run_tool(
          pool, cluster, as.character(parsed$tool)[1], parsed$args %||% list(),
          source_servers[[cluster]]
        )
        conv$seen[[call]] <- observation
      }
      conv$transcript[[length(conv$transcript) + 1L]] <- list(
        turn = conv$turn,
        call = list(tool = as.character(parsed$tool)[1], args = parsed$args %||% list()),
        duplicate = conv$repeated > 0L
      )
      conv$prompt <- paste0(
        conv$prompt,
        .observation_block(length(conv$transcript) + length(conv$attempts), observation)
      )
      if (conv$repeated >= MAX_REPEATED_QUERIES ||
        (MAX_TURNS > 0L && conv$turn >= MAX_TURNS)) {
        conv$forced <- TRUE
        conv$prompt <- paste0(
          conv$prompt, "\n# No further queries. Return the final answer now.\n"
        )
      }
      conversation[[cluster]] <- conv
      still <- c(still, cluster)
      next
    }
    conv$schema_failures <- conv$schema_failures + 1L
    if (conv$schema_failures > SCHEMA_RETRIES) {
      conv$outcome <- list(
        error = "schema violation", turns = conv$turn, transcript = conv$transcript,
        trace_id = conv$trace_id
      )
      conversation[[cluster]] <- conv
      next
    }
    problems <- if (is.list(parsed) && identical(parsed$action %||% "", "final")) {
      .final_problems(parsed, pool, cluster, conv$tier, conv$program)
    } else {
      character(0)
    }
    conv$prompt <- paste0(
      conv$prompt,
      if (length(problems)) {
        sprintf(
          "\n\n# Retry %d: the final answer was rejected -- %s. Return one valid final JSON object.\n",
          conv$schema_failures, paste(head(problems, 6L), collapse = "; ")
        )
      } else {
        sprintf("\n\n# Retry %d: return one valid JSON object.\n", conv$schema_failures)
      }
    )
    conversation[[cluster]] <- conv
    still <- c(still, cluster)
  }
  audited <- c(to_review, redelivered)
  if (length(audited)) {
    audit_step <- .audit_delivery_many(
      pool, audited, pending_final,
      setNames(lapply(audited, function(cl) conversation[[cl]]$tier), audited),
      setNames(lapply(audited, function(cl) conversation[[cl]]$program), audited),
      audit_template, setNames(lapply(audited, function(cl) headers[[cl]] %||% .context_header(pool, cl)), audited),
      api_key, api_url,
      setNames(lapply(audited, function(cl) conversation[[cl]]$trace_id), audited),
      setNames(lapply(audited, function(cl) conversation[[cl]]$turn), audited),
      setNames(lapply(audited, function(cl) conversation[[cl]]$audit_cache), audited)
    )
    pending_final <- audit_step$finals
    reviews <- if (length(to_review)) {
      .review_many(pool, to_review, pending_final, source_servers, conversation, api_url, api_key)
    } else {
      list()
    }
    to_arbitrate <- character(0)
    contested_of_arb <- list()
    for (cluster in redelivered) {
      conv <- conversation[[cluster]]
      final <- pending_final[[cluster]]
      note <- sprintf(
        "'%s' was re-delivered unchanged after a refusal; rounds ended, decided by arbitration",
        as.character(final$selected %||% "")[1]
      )
      conv$structural_notes <- c(conv$structural_notes, list(note))
      conv$attempts <- c(conv$attempts, list(list(tier = conv$tier, final = final, review = NULL, structural_note = note)))
      conversation[[cluster]] <- conv
      to_arbitrate <- c(to_arbitrate, cluster)
    }
    for (cluster in to_review) {
      conv <- conversation[[cluster]]
      final <- pending_final[[cluster]]
      selected_now <- as.character(final$selected %||% "")[1]
      review_pair <- reviews[[cluster]]
      review <- review_pair$verdict
      status <- review_pair$status %||% (if (is.null(review)) "unreachable" else "ok")
      structural_note <- ""
      conv$attempts <- c(
        conv$attempts, list(list(tier = conv$tier, final = final, review = review))
      )
      conv$last_packet <- review_pair$packet %||% list()
      subtype <- as.character(final$subtype %||% "")[1]
      if (identical(status, "unreachable") || is.null(review)) {
        conv$review_unreachable <- TRUE
        conv$structural_notes <- c(conv$structural_notes, list("evidence review unreachable; delivered unresolved at low confidence"))
        conv$outcome <- .deliver(conv, final, NULL)
        conversation[[cluster]] <- conv
        next
      }
      if (identical(status, "uncited")) {
        note <- sprintf(
          "review verdict '%s' could not be validated (citations short after retries); decided by arbitration",
          as.character(review$identity_verdict %||% "")
        )
        conv$structural_notes <- c(conv$structural_notes, list(note))
        conv$attempts[[length(conv$attempts)]]$structural_note <- note
        conversation[[cluster]] <- conv
        to_arbitrate <- c(to_arbitrate, cluster)
        next
      }
      if (.review_passed(review, subtype)) {
        contested <- delivery_contested(final, final$definer_audit %||% list(), final$contested_rivals %||% character(0))
        population <- population_conflicts(final, review)
        if (!length(contested) && !length(population)) {
          conv$outcome <- .deliver(conv, final, NULL)
          conversation[[cluster]] <- conv
          next
        }
        note <- paste0("review accepted, but ", paste(c(contested, population), collapse = "; "), "; sent back to the annotator")
        conv$structural_notes <- c(conv$structural_notes, list(note))
        conv$attempts[[length(conv$attempts)]]$structural_note <- note
        conv$rounds_in_tier <- conv$rounds_in_tier + 1L
        rounds_left <- max(0L, REVIEW_MAX_ROUNDS_PER_TIER - conv$rounds_in_tier)
        if (rounds_left > 0L && !isTRUE(conv$forced)) {
          conv$last_refused_answer <- .answer_key(final)
          conv$prompt <- paste0(
            conv$prompt,
            .observation_block(
              length(conv$transcript) + length(conv$attempts) + 1L,
              .contested_feedback(
                review, contested, conv$tier, rounds_left, conv$tiers - conv$tier,
                .pair_tables(pool, cluster, selected_now, final$contested_rivals %||% character(0)),
                population
              )
            )
          )
          conversation[[cluster]] <- conv
          still <- c(still, cluster)
          next
        }
        conversation[[cluster]] <- conv
        to_arbitrate <- c(to_arbitrate, cluster)
        contested_of_arb[[cluster]] <- contested
        next
      }
      conv$last_refused_answer <- .answer_key(final)
      conv$rounds_in_tier <- conv$rounds_in_tier + 1L
      rounds_left <- max(0L, REVIEW_MAX_ROUNDS_PER_TIER - conv$rounds_in_tier)
      if (rounds_left > 0L && !isTRUE(conv$forced)) {
        better <- as.character(review$better_candidate %||% "")[1]
        better_audit <- NULL
        if (nzchar(better) && !identical(better, UNKNOWN) && is.null((final$definer_audit %||% list())[[better]])) {
          hit <- audit_many(
            pool, list(list(cluster = cluster, identity = better, role = "better_candidate",
                            claimed = unique(c(vapply(.claimed_identities(final), function(c) c[[2]], ""), better)),
                            tier = conv$tier)),
            NULL, list(as.list(conv$program)), audit_template, list(headers[[cluster]] %||% .context_header(pool, cluster)),
            api_key, api_url, list(conv$trace_id), list(as.integer(conv$turn + 70L)), SCHEMA_RETRIES
          )
          better_audit <- hit[[paste(cluster, better, sep = "\x1f")]]
        }
        versus <- setdiff(c(if (nzchar(better) && !identical(better, UNKNOWN)) better,
                            as.character(unlist(final$contested_rivals %||% list()))), selected_now)
        conv$prompt <- paste0(
          conv$prompt,
          .observation_block(
            length(conv$transcript) + length(conv$attempts) + 1L,
            .review_feedback(review, conv$tier, rounds_left, conv$tiers - conv$tier,
                             final$definer_audit, better_audit, .pair_tables(pool, cluster, selected_now, versus))
          )
        )
        conversation[[cluster]] <- conv
        still <- c(still, cluster)
        next
      }
      conversation[[cluster]] <- conv
      to_arbitrate <- c(to_arbitrate, cluster)
    }
    if (length(to_arbitrate)) {
      decisions <- .arbitrate_many(
        pool, to_arbitrate, conversation, api_url, api_key,
        contested_of = contested_of_arb,
        audits_of = setNames(lapply(to_arbitrate, function(cl) pending_final[[cl]]$definer_audit), to_arbitrate),
        pair_tables_of = setNames(lapply(to_arbitrate, function(cl) {
          .pair_tables(pool, cl, as.character(pending_final[[cl]]$selected %||% "")[1],
                       as.character(unlist(pending_final[[cl]]$contested_rivals %||% list())))
        }), to_arbitrate)
      )
      for (cluster in to_arbitrate) {
        conv <- conversation[[cluster]]
        final <- pending_final[[cluster]]
        decision <- decisions[[cluster]]
        if (is.null(decision)) {
          conv$outcome <- .deliver(conv, final, NULL, contested_of_arb[[cluster]] %||% character(0))
          conversation[[cluster]] <- conv
          next
        }
        if (identical(as.character(decision$selected %||% "")[1], UNKNOWN) && !isTRUE(conv$forced) && conv$tier < conv$tiers) {
          why <- trimws(as.character(decision$reason %||% ""))
          note <- sprintf(
            "arbitration at tier %d answered Unknown: no candidate on this page is carried by its own evidence; the next tier of candidates is added",
            conv$tier
          )
          conv$structural_notes <- c(conv$structural_notes, list(note))
          conv$attempts <- c(conv$attempts, list(list(
            tier = conv$tier, review = NULL, structural_note = note,
            final = list(selected = UNKNOWN, subtype = "", lineage = "", co_occurring_identities = list(),
                        claim_evidence = list(), support_markers = list(), confidence = "low", reason = why)
          )))
          step <- .escalate(conv, pool, cluster, paste0(
            "The rounds at this tier were decided by arbitration, and the arbitration found no candidate on this page carried by its own evidence",
            if (nzchar(why)) sprintf(" (%s)", substr(why, 1, 400)) else "", "."
          ))
          conv <- step$conv
          conv$last_refused_answer <- ""
          conversation[[cluster]] <- conv
          still <- c(still, cluster)
          next
        }
        final_answer <- if (identical(as.character(decision$selected %||% "")[1], UNKNOWN)) {
          out <- final
          out$selected <- UNKNOWN
          out$subtype <- ""
          out$lineage <- ""
          out$co_occurring_identities <- list()
          out$claim_evidence <- list()
          out$support_markers <- list()
          out$reason <- as.character(decision$reason %||% "")[1]
          confidence_value <- as.character(decision$confidence %||% "")[1]
          out$confidence <- if (nzchar(confidence_value)) confidence_value else "low"
          out
    } else {
          base <- final
          decision_selected <- as.character(decision$selected %||% "")[1]
          for (item in rev(conv$attempts)) {
            if (identical(as.character(item$final$selected %||% "")[1], decision_selected)) {
              base <- item$final
              break
            }
          }
          .arbiter_final(base, decision)
        }
        recheck <- if (!identical(as.character(final_answer$selected %||% "")[1], UNKNOWN)) {
          re_audit <- .audit_delivery_many(
            pool, cluster, setNames(list(final_answer), cluster),
            setNames(list(conv$tier), cluster), setNames(list(conv$program), cluster),
            audit_template, setNames(list(headers[[cluster]] %||% .context_header(pool, cluster)), cluster),
            api_key, api_url, setNames(list(conv$trace_id), cluster), setNames(list(as.integer(conv$turn + 50L)), cluster),
            setNames(list(conv$audit_cache), cluster)
          )
          list(final = re_audit$finals[[cluster]], contested = re_audit$contested[[cluster]])
        } else {
          list(final = final_answer, contested = character(0))
        }
        final_answer <- recheck$final
        contested2 <- recheck$contested
        if (length(contested2)) {
          note <- paste0("the arbitration was read again because ", paste(contested2, collapse = "; "))
          conv$structural_notes <- c(conv$structural_notes, list(note))
          again <- .arbitrate_many(
            pool, cluster, conversation, api_url, api_key,
            contested_of = setNames(list(contested2), cluster),
            audits_of = setNames(list(final_answer$definer_audit), cluster),
            pair_tables_of = setNames(list(.pair_tables(
              pool, cluster, as.character(final_answer$selected %||% "")[1],
              as.character(unlist(final_answer$contested_rivals %||% list()))
            )), cluster),
            turn_offset = 60L
          )[[cluster]]
          if (!is.null(again)) {
            decision <- again
            final_answer <- if (identical(as.character(decision$selected %||% "")[1], UNKNOWN)) {
              out <- final
              out$selected <- UNKNOWN; out$subtype <- ""; out$lineage <- ""
              out$co_occurring_identities <- list(); out$claim_evidence <- list(); out$support_markers <- list()
              out$reason <- as.character(decision$reason %||% "")[1]
              confval <- as.character(decision$confidence %||% "")[1]
              out$confidence <- if (nzchar(confval)) confval else "low"
              out
            } else {
              .arbiter_final(final, decision)
            }
            re_audit2 <- if (!identical(as.character(final_answer$selected %||% "")[1], UNKNOWN)) {
              .audit_delivery_many(
                pool, cluster, setNames(list(final_answer), cluster),
                setNames(list(conv$tier), cluster), setNames(list(conv$program), cluster),
                audit_template, setNames(list(headers[[cluster]] %||% .context_header(pool, cluster)), cluster),
                api_key, api_url, setNames(list(conv$trace_id), cluster), setNames(list(as.integer(conv$turn + 70L)), cluster),
                setNames(list(conv$audit_cache), cluster)
              )
            } else {
              list(finals = setNames(list(final_answer), cluster), contested = setNames(list(character(0)), cluster))
            }
            final_answer <- re_audit2$finals[[cluster]]
            contested2 <- re_audit2$contested[[cluster]]
          }
          if (length(contested2)) {
            conv$structural_notes <- c(conv$structural_notes, list(paste0("delivered unresolved: ", paste(contested2, collapse = "; "))))
          }
        }
        conv$outcome <- .deliver(conv, final_answer, decision, contested2)
        conversation[[cluster]] <- conv
      }
    }
  }
  active <- still
}
conversation
}
conversation <- .run_conversation_rounds(pool, conversation, source_servers, headers, api_key, api_url)
for (cluster in runnable) {
  conv <- conversation[[cluster]]
  outcome <- conv$outcome %||% list(
    error = "no answer", turns = conv$turn, transcript = conv$transcript
  )
  if (is.null(outcome$final)) {
    cluster_value <- cluster
    frame <- scored[cluster == cluster_value]
    record <- .fallback(
      cluster, frame, sprintf("failed:%s", outcome$error %||% ""), "annotation_failed"
    )
    record$candidates <- as.character(.candidate_names(pool, cluster))
    record$candidate_entries <- .plain_entries(cluster)
    record$turns <- as.integer(outcome$turns %||% 0L)
    results[[cluster]] <- record
    next
  }
  results[[cluster]] <- .result(cluster, pool, outcome)
}
dataset_consistency <- list()
reopened_by_consistency <- list()
if (length(runnable) && identical(llm_status, "enabled") && !is.null(cluster_subset)) {
  cat(sprintf("  dataset consistency skipped: SCMA_CLUSTER_SUBSET annotates %d of %d clusters\n",
              length(runnable), length(results)))
}
if (length(runnable) && identical(llm_status, "enabled") && is.null(cluster_subset)) {
  dataset_consistency <- dataset_consistency_many(
    pool, results, .prompt("dataset_consistency"),
    sprintf(
      "# Dataset context\nspecies: %s | tissue: %s | disease: %s\n\n",
      as.character(context$species), paste(as.character(unlist(context$tissue)), collapse = ", "),
      paste(as.character(unlist(context$disease)), collapse = ", ")
    ),
    api_key, api_url, substr(paste0(format(Sys.time(), "%H%M%S"), basename(tempfile())), 1, 16),
    SCHEMA_RETRIES
  )
  to_reopen <- list()
  for (item in dataset_consistency) {
    reading <- item$reading %||% list()
    for (member in reading$members %||% list()) {
      if (as.character(member$status %||% "") %in% c("weaker_member", "conflict")) {
        cid <- as.character(member$cluster_id)
        if (!is.null(results[[cid]]) && cid %in% runnable) {
          to_reopen[[cid]] <- list(
            status = member$status, why_grouped = item$group$why_grouped,
            note = as.character(member$note_for_review %||% ""),
            group_reason = as.character(reading$reason %||% "")
          )
        }
      }
    }
  }
  if (length(to_reopen)) {
    reopen_ids <- names(to_reopen)
    cat(sprintf(
      "  dataset consistency: %d group reading(s); re-opening %d cluster(s): %s\n",
      length(dataset_consistency), length(to_reopen), paste(sort(reopen_ids), collapse = ", ")
    ))
    reopened_state <- .open_conversations(
      pool, reopen_ids, sources_context, api_key, api_url,
      consistency_notes = to_reopen
    )
    reopened_state$conversation <- .run_conversation_rounds(
      pool, reopened_state$conversation, reopened_state$source_servers,
      reopened_state$headers, api_key, api_url
    )
    for (cid in reopen_ids) {
      conv <- reopened_state$conversation[[cid]]
      outcome <- conv$outcome
      if (is.null(outcome) || is.null(outcome$final)) {
        next
      }
      before <- results[[cid]]
      after <- .result(cid, pool, outcome)
      after$reopened_by_consistency <- c(
        to_reopen[[cid]],
        list(before = list(
          annotation = before$annotation, confidence = before$confidence,
          co_occurring_identities = before$co_occurring_identities
        ))
      )
      results[[cid]] <- after
      reopened_by_consistency[[cid]] <- after$reopened_by_consistency
      conversation[[cid]] <- conv
    }
  }
}
results <- results[clusters_sorted]
scma_flush_response_cache()
saveRDS(
  list(
    tag = tag, context = context, llm_status = llm_status,
    dataset_consistency = dataset_consistency,
    reopened_by_consistency = reopened_by_consistency,
    borrow_screen = pool$borrow_screen %||% list(status = "not_run"),
    models = list(
      annotator = scma_resolve_model(ANNOTATOR_MODEL),
      annotator_reasoning_effort = ANNOTATOR_EFFORT,
      review = scma_resolve_model(REVIEW_MODEL),
      review_reasoning_effort = REVIEW_EFFORT,
      arbiter = scma_resolve_model(ARBITER_MODEL),
      arbiter_reasoning_effort = ARBITER_EFFORT,
      borrow_screen = scma_resolve_model(BORROW_SCREEN_MODEL),
      borrow_screen_reasoning_effort = BORROW_SCREEN_EFFORT
    ),
    results = results,
    transcripts = lapply(conversation, function(conv) conv$transcript),
    borrow_qualifier_lexicon = pool$borrow_qualifier_lexicon %||% list(),
    tissue_contexts = pool$tissue_contexts %||% list(),
    candidate_tissue_contexts = pool$candidate_tissue_contexts %||% list(),
    candidate_retrieval_context = pool$candidate_retrieval_context %||% list(),
    retrieval_semantics = pool$retrieval_semantics %||% "",
    candidate_tier_size = CANDIDATE_TIER_SIZE,
    candidate_tiers = CANDIDATE_TIERS
  ),
  file.path(CACHE, sprintf("%s_annotations.rds", tag))
)
statuses <- vapply(results, function(r) r$resolution_status, "")
turn_counts <- vapply(results, function(r) as.integer(r$turns %||% 0L), 0L)
warned <- sum(vapply(results, function(r) length(r$claim_warnings %||% list()) > 0L, TRUE))
qc_counts <- setNames(rep(0L, length(QC_VALUES)), QC_VALUES)
for (r in results) {
  key <- as.character(r$annotation_qc %||% QC_UNCHECKED)
  qc_counts[key] <- as.integer(qc_counts[key]) + 1L
}
tier_counts <- list()
for (r in results) {
  key <- suppressWarnings(as.integer(r$delivered_tier %||% 0L))
  if (length(key) && !is.na(key) && key > 0L) {
    key_name <- as.character(key)
    tier_counts[[key_name]] <- (tier_counts[[key_name]] %||% 0L) + 1L
  }
}
tier_names <- names(tier_counts)
tier_order <- if (length(tier_names)) order(as.integer(tier_names)) else integer(0)
not_model_decided <- sum(vapply(
  results, function(r) !identical(r$annotation_source, "cluster_annotation"), TRUE
))
cat(sprintf(
  paste0(
    "[done] %s: %d resolved, %d mixed, %d unresolved of %d clusters; turns median %d, ",
    "max %d; evidence review %s; delivered from tier %s; %d clusters carry a claim ",
    "warning; %d not model-decided\n"
  ),
  tag, sum(statuses == RESOLVED), sum(statuses == MIXED), sum(statuses == UNRESOLVED),
  length(results),
  as.integer(stats::median(turn_counts[turn_counts > 0L] %||% 0L)),
  max(turn_counts, 0L),
  paste(sprintf("%s=%d", names(qc_counts), as.integer(qc_counts)), collapse = ", "),
  if (length(tier_names)) {
    paste(sprintf(
      "%s=%d", tier_names[tier_order], as.integer(unlist(tier_counts))[tier_order]
    ), collapse = ", ")
  } else {
    "none"
  },
  warned, not_model_decided
))
