#!/usr/bin/env Rscript
# =============================================================================
# cluster_annotation.R -- one agent, able to query the data, decides each cluster.
# Renders the shared prompt template over the shared packet fields.
#
# The agent opens with every retrieved candidate's COMPLETE curated panel as measured
# here, ordered by publication support, and may then call read-only tools for the curated
# source sentences, for a gene's behaviour in the other clusters, and for which candidates
# claim a gene. It answers when it is ready.
#
# Nothing is removed before it answers: the candidates are compared together, once, by the
# only participant that can see the alternatives. A single-candidate reading cannot
# foreclose the comparison that is the whole judgement.
#
# The code decides nothing about biology. Afterwards it checks the numbers the agent
# quoted against the DE table and flags claims whose curated markers are not raised,
# without deleting or renaming anything.
#
# Output: cache/<tag>_annotations.rds
#
#   Usage: Rscript cluster_annotation.R <tag>
# =============================================================================
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
source(file.path(.sd, "marker_sources.R"))
sys.source(file.path(.sd, "annotator_pool.R"), envir = environment())

args <- commandArgs(trailingOnly = TRUE)
tag <- if (length(args) >= 1 && nzchar(args[1])) args[1] else INPUT_TAG

A <- CFG$cluster_annotation
ANNOTATOR_MODEL <- as.character(A$annotator_model)
ANNOTATOR_EFFORT <- as.character(A$annotator_reasoning_effort)
ANNOTATOR_SCHEMA <- as.character(A$annotator_schema_version)
MAX_TURNS <- as.integer(A$max_turns)
SCHEMA_RETRIES <- as.integer(A$schema_retries)
UNKNOWN <- as.character(A$unknown_token)
MAX_RATIONALE <- as.integer(A$max_rationale_chars)
SUBTYPE_EXCLUSIVE_REQUIRED <- as.integer(A$subtype_exclusive_markers_required)
JUDGE_MODEL <- as.character(A$judge_model)
JUDGE_EFFORT <- as.character(A$judge_reasoning_effort)
JUDGE_MAX_ROUNDS <- as.integer(A$judge_max_rounds)
## How many times one judgment may ask for evidence the packet did not carry.
JUDGE_QUERY_TURNS <- as.integer(A$judge_query_turns)

CONFIDENCE_VALUES <- c("high", "medium", "low")
## A faithfully copied percentage differs from the shown one by less than a rounding step.
QUOTE_TOLERANCE <- 0.15
## Re-asking an identical question cannot make progress, so that is what the loop counts
## rather than a turn budget that would cut off a cluster still learning something.
MAX_REPEATED_QUERIES <- 3L

RESOLVED <- "resolved"
MIXED <- "mixed"
UNRESOLVED <- "unresolved"
UNSUPPORTED <- "unsupported_empty_candidate_pool"

## One template, with the shared panel-reading rules spliced in. How a marker row is read
## -- what counts as raised, what an unraised one means, how publication counts weigh --
## is one question, and every template that reads a panel has to answer it the same way.
## Held in one file and included, so the annotator and the judges cannot drift into two
## readings of the same number.
##
## Under the cross-species switch the pooled evidence deliberately crosses species, and
## every template that reads it is told so in one shared note (cross_species_note.txt,
## the same file the Python arm splices): a sentence may name another species, and that
## is not grounds against a label. OFF leaves every template byte-identical.
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
## `fixed = TRUE` is load-bearing: R reads backslashes in a regex REPLACEMENT, so a
## curated sentence carrying a quote reached the model as `the "inflammatory" subset`
## with the JSON escape eaten, and the packet was no longer parseable JSON.
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

## One tool result appended to the running prompt. Only this observation is appended:
## never the opening packet again, and never a result the agent did not ask for. The fixed
## head plus the packet therefore stay byte-identical from turn to turn, which is the part
## a provider can serve from a cached prefix.
.observation_block <- function(index, observation) {
  sprintf("\n\n# Observation %d\n%s\n", index, .compact_json(observation))
}

## ---- quality control ---------------------------------------------------------
SUBTYPE_PASS <- "separable"
ANNOTATION_PASS <- "supported"
## The verdicts a second answer could plausibly change, and only those. A label pitched
## finer than the evidence reaches can be answered from candidates the agent already
## holds, and a contradiction names the gene to answer against. `insufficient_evidence` is
## not there on purpose: the panel too thin to decide is the same panel the next answer
## would read, so asking again buys the same reply at the price of a turn. Nor is a
## refinement that was merely demoted -- the parent it fell back to is a real answer, and
## re-opening a call the evidence already carries is how a check turns into a haggle.
##
## There is no `too_broad` here because this judge cannot reach it: saying a label is too
## broad means naming a narrower identity, and this judge is handed one label and told to
## use nothing else. An inverted refinement is caught upstream by the subtype judge's
## direction rule, which is given both labels and can compare them.
RETRYABLE_VERDICTS <- c("too_specific", "contradicted")

## How good an answer was, for choosing between the versions the rounds produced. A label
## only pitched too fine still identified the lineage; a contradicted label did not, and a
## panel too thin to decide settled nothing either way.
.verdict_rank <- function(verdict) {
  switch(as.character(verdict %||% ""),
    supported = 3L,
    too_specific = 2L,
    insufficient_evidence = 1L,
    contradicted = 0L,
    # An unreadable verdict ranks with a refuted one: neither is an answer to deliver.
    0L
  )
}

## One batch of quality-control judgments, with the turns they ask for.
##
## Its own model and its own prompt, so the cache key differs from the annotator's and a
## judgment is never served the annotator's reply.
##
## A judge used to get one call and whatever evidence code had decided to send. That made
## code the last word on what evidence exists, and code was wrong about it: a fourteen-row
## cut hid 10% of the genes the answers rested on. It now reads the whole panel and may ask
## for what the packet did not carry -- the same growing prompt the annotator uses, so the
## head stays byte-identical and a provider can serve it from cache. Batched by round rather
## than one judgment at a time, because 38 clusters asked in series would be unusable.
.judge_many <- function(pool, clusters, packets, template, conversation,
                        servers, api_url, api_key) {
  out <- vector("list", length(packets))
  active <- which(vapply(packets, function(packet) length(packet) > 0L, TRUE))
  if (!length(active)) {
    return(out)
  }
  prompts <- list()
  labels <- list()
  seen <- list()
  unusable <- list()
  asked <- list()
  for (index in active) {
    prompts[[as.character(index)]] <- .render(template, packets[[index]])
    packet <- packets[[index]]
    labels[[as.character(index)]] <- as.character(unlist(
      packet[intersect(c("label_under_test", "parent_label"), names(packet))]
    ))
    seen[[as.character(index)]] <- list()
    unusable[[as.character(index)]] <- 0L
    asked[[as.character(index)]] <- 0L
  }
  ## One pass per question the budget allows, plus the pass that answers, plus room for a
  ## reply that was neither so it can be asked again instead of being filed as unchecked.
  for (round in seq_len(JUDGE_QUERY_TURNS + 1L + SCHEMA_RETRIES)) {
    responses <- scma_cached_call_llm_many(
      lapply(active, function(index) prompts[[as.character(index)]]),
      api_url, api_key,
      reasoning_effort = JUDGE_EFFORT, model = JUDGE_MODEL,
      trace_ids = lapply(active, function(index) conversation[[clusters[index]]]$trace_id),
      turn_indexes = lapply(
        active, function(index) as.integer(conversation[[clusters[index]]]$turn)
      )
    )
    still <- integer(0)
    for (position in seq_along(active)) {
      index <- active[position]
      key <- as.character(index)
      content <- responses[[position]][[1]]
      parsed <- if (is.null(content)) NULL else scma_parse_json(content)
      if (is.list(parsed) && !is.null(parsed$verdict)) {
        out[[index]] <- parsed
        next
      }
      answerable <- .valid_query(parsed) && asked[[key]] < JUDGE_QUERY_TURNS
      if (!answerable) {
        ## Neither a verdict nor a question that can still be answered, so it is asked
        ## again rather than filed as unchecked. `not_checked` has to keep meaning "no judge
        ## was reachable": reporting a label as never checked because one reply was badly
        ## formatted tells a reader the check did not run, which is a different fact.
        if (unusable[[key]] >= SCHEMA_RETRIES) next
        unusable[[key]] <- unusable[[key]] + 1L
        prompts[[key]] <- paste0(
          prompts[[key]],
          sprintf("\n\n# Retry %d: return the verdict as one JSON object.\n", unusable[[key]])
        )
        still <- c(still, index)
        next
      }
      call <- .compact_json(list(
        tool = as.character(parsed$tool)[1], args = parsed$args %||% list()
      ))
      ## A repeated call returns the same answer, so it is answered once and marked. The
      ## budget is what stops the loop; this only stops it buying the same row twice.
      ## `sources` is exempt: it serves the next sentences rather than the same ones, so
      ## asking again is a new answer. Its own per-pair limit is the brake instead.
      streaming <- identical(as.character(parsed$tool)[1], "sources")
      observation <- if (streaming) NULL else seen[[key]][[call]]
      if (is.null(observation)) {
        observation <- .run_judge_tool(
          pool, clusters[index], labels[[key]], parsed$tool,
          parsed$args %||% list(), servers[[clusters[index]]]
        )
        seen[[key]][[call]] <- observation
      } else {
        observation <- c(observation, list(duplicate_query = TRUE))
      }
      asked[[key]] <- asked[[key]] + 1L
      prompts[[key]] <- paste0(
        prompts[[key]], .observation_block(asked[[key]], observation),
        if (asked[[key]] >= JUDGE_QUERY_TURNS) {
          "\n# No further queries. Return the verdict now.\n"
        } else {
          ""
        }
      )
      still <- c(still, index)
    }
    active <- still
    if (!length(active)) break
  }
  out
}

## Both quality-control questions asked of a batch of answers.
##
## Two judges rather than one, because they are different questions asked of different
## evidence: whether the finer name is separable from ITS PARENT, and whether the
## delivered label is carried by its own panel. A single verdict would have to mean both
## at once, and a refinement dropped for being unseparable would be indistinguishable from
## a lineage call that was simply wrong.
##
## A refinement is demoted whenever it is not affirmatively separable, so an unchecked
## refinement costs granularity rather than passing by default. Where no judge answered at
## all, nothing is demoted and nothing is retried: silently coarsening every label of an
## offline replay would be a worse failure than leaving the answer unchecked.
.review_many <- function(pool, clusters, finals, servers, conversation, api_url, api_key) {
  reviews <- list()
  selected_of <- list()
  subtype_of <- list()
  for (cluster in clusters) {
    final <- finals[[cluster]]
    selected_of[[cluster]] <- as.character(final$selected %||% "")[1]
    subtype_of[[cluster]] <- as.character(final$subtype %||% "")[1]
    reviews[[cluster]] <- list(
      subtype_verdict = NULL, annotation_verdict = NULL, demoted_subtype = "",
      judged_label = "", checked = FALSE, passed = FALSE, retry = FALSE, rank = 0L
    )
  }
  ask <- Filter(function(cluster) {
    nzchar(subtype_of[[cluster]]) && nzchar(selected_of[[cluster]]) &&
      !identical(selected_of[[cluster]], UNKNOWN)
  }, clusters)
  if (length(ask)) {
    verdicts <- .judge_many(
      pool, ask,
      lapply(ask, function(cluster) {
        .judge_packet(
          pool, cluster, subtype_of[[cluster]], servers[[cluster]],
          parent = selected_of[[cluster]]
        )
      }),
      .prompt("subtype_judge"), conversation, servers, api_url, api_key
    )
    for (position in seq_along(ask)) {
      cluster <- ask[position]
      verdict <- verdicts[[position]]
      reviews[[cluster]]$subtype_verdict <- verdict
      if (!is.null(verdict)) {
        reviews[[cluster]]$checked <- TRUE
        if (!identical(as.character(verdict$verdict %||% "")[1], SUBTYPE_PASS)) {
          reviews[[cluster]]$demoted_subtype <- subtype_of[[cluster]]
        }
      }
    }
  }
  for (cluster in clusters) {
    subtype <- subtype_of[[cluster]]
    reviews[[cluster]]$judged_label <- if (
      nzchar(subtype) && !nzchar(reviews[[cluster]]$demoted_subtype)
    ) {
      subtype
    } else {
      selected_of[[cluster]]
    }
  }
  ask <- Filter(function(cluster) {
    label <- reviews[[cluster]]$judged_label
    nzchar(label) && !identical(label, UNKNOWN)
  }, clusters)
  ## `Unknown` claims no identity, so there is nothing for a judge to carry or to refute,
  ## and the answer stands as given. That counts as checked: `not_checked` is reserved for
  ## a judge that could not be reached, and an abstention reported the same way would be
  ## indistinguishable from a run whose model was down.
  for (cluster in setdiff(clusters, ask)) {
    reviews[[cluster]]$checked <- TRUE
    reviews[[cluster]]$passed <- TRUE
    reviews[[cluster]]$rank <- .verdict_rank(ANNOTATION_PASS)
  }
  if (length(ask)) {
    verdicts <- .judge_many(
      pool, ask,
      lapply(ask, function(cluster) {
        .judge_packet(pool, cluster, reviews[[cluster]]$judged_label, servers[[cluster]])
      }),
      .prompt("annotation_judge"), conversation, servers, api_url, api_key
    )
    for (position in seq_along(ask)) {
      cluster <- ask[position]
      verdict <- verdicts[[position]]
      reviews[[cluster]]$annotation_verdict <- verdict
      if (!is.null(verdict)) {
        answer <- as.character(verdict$verdict %||% "")[1]
        reviews[[cluster]]$checked <- TRUE
        reviews[[cluster]]$rank <- .verdict_rank(answer)
        reviews[[cluster]]$passed <- identical(answer, ANNOTATION_PASS)
        reviews[[cluster]]$retry <- answer %in% RETRYABLE_VERDICTS
      }
    }
  }
  reviews
}

## What the annotator is told, and only that.
##
## The packet is already at the top of this conversation, so none of it is sent again:
## what comes back is the verdict, the genes it turned on, and the direction the level was
## wrong in. That is the whole point of asking -- an error found and not carried to the
## next step is an error nobody acted on -- but it is also the whole of it. Which candidate
## to name is the agent's decision, made against the candidates it already holds, and
## re-deciding it here would put two annotators on one cluster.
.qc_feedback <- function(review) {
  block <- list()
  verdict <- review$subtype_verdict
  if (!is.null(verdict)) {
    block$subtype <- list(
      name = review$demoted_subtype,
      verdict = verdict$verdict,
      conflicting_markers = verdict$conflicting_markers %||% list(),
      reason = verdict$reason
    )
  }
  verdict <- review$annotation_verdict
  if (!is.null(verdict)) {
    block$annotation <- list(
      label = review$judged_label,
      verdict = verdict$verdict,
      conflicting_markers = verdict$conflicting_markers %||% list(),
      reason = verdict$reason
    )
  }
  block$note <- paste0(
    "quality control on the answer above, judged from that identity's own panel and ",
    "its curated sentences alone. A refinement that is not separable has already been ",
    "dropped to the parent. Answer again against the candidates you already hold: keep ",
    "the call where the evidence still carries it, or name the candidate the evidence ",
    "does carry. Query only where the verdict turns on something not already on this page."
  )
  block
}

## How the delivered label left quality control, in one word.
##
## `not_checked` is its own answer rather than a pass: a replay with no reachable judge
## would otherwise report every cluster as having survived a check nobody ran.
.qc_status <- function(review, demoted) {
  if (!length(review) || !isTRUE(review$checked)) {
    return(QC_UNCHECKED)
  }
  label <- as.character(review$judged_label %||% "")
  verdict <- as.character((review$annotation_verdict %||% list())$verdict %||% "")
  carried <- identical(verdict, ANNOTATION_PASS) || !nzchar(label) ||
    identical(label, UNKNOWN)
  if (!carried) {
    return(QC_FAILED)
  }
  if (isTRUE(demoted)) {
    return(QC_DEMOTED)
  }
  if (as.integer(review$selected_round %||% 1L) > 1L) QC_REVISED else QC_PASSED
}

## ---- answer shape ------------------------------------------------------------
.valid_query <- function(value) {
  is.list(value) && identical(value$action %||% "", "query") &&
    is.character(value$tool %||% NULL) && length(value$tool) == 1L
}

## (role, identity) for every identity the answer asserts is present. Selected, subtype
## and each co-occurring name are one kind of thing -- a claim that these cells carry that
## identity -- so the evidence requirement and the warning apply to all of them alike.
.claimed_identities <- function(value) {
  claims <- list()
  selected <- as.character(value$selected %||% "")
  if (nzchar(selected) && !identical(selected, UNKNOWN)) {
    claims[[length(claims) + 1L]] <- list("selected", selected)
  }
  subtype <- as.character(value$subtype %||% "")
  if (nzchar(subtype)) claims[[length(claims) + 1L]] <- list("subtype", subtype)
  for (name in value$co_occurring_identities %||% list()) {
    claims[[length(claims) + 1L]] <- list("cooc", as.character(name))
  }
  claims
}

## The subtype to drop and the line saying why, or two empty strings.
##
## Shared lineage markers cannot establish a finer identity, so the refinement is kept
## only when the cluster raises enough markers that the parent's curated definition does
## not contain. The line records the count against the requirement and names whichever
## exclusive markers there were, so a demotion can be checked rather than trusted.
.subtype_demotion <- function(pool, cluster, selected, subtype) {
  empty <- list(rejected = "", line = "")
  if (!nzchar(subtype) || !nzchar(selected) || identical(selected, UNKNOWN)) {
    return(empty)
  }
  exclusive <- .subtype_exclusive_markers(pool, cluster, selected, subtype)
  if (length(exclusive) >= SUBTYPE_EXCLUSIVE_REQUIRED) {
    return(empty)
  }
  found <- if (length(exclusive)) paste(exclusive, collapse = ",") else "none"
  list(
    rejected = subtype,
    line = sprintf(
      paste0(
        "rejected_subtype %s: exclusive_defining_markers %d/%d",
        " | raised panel is shared with %s | exclusive: %s"
      ),
      subtype, length(exclusive), SUBTYPE_EXCLUSIVE_REQUIRED, selected, found
    )
  )
}

## Every claim names a gene from its own panel and copies its measurement. The gene must
## belong to the claimed identity's curated panel, not merely to the dataset: a claim
## supported by a gene the resource never tied to that cell type is not a claim the
## curated evidence can carry.
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
      if (identical(toupper(row$gene), gene)) marker <- row
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

## Whether one answer is answerable from the evidence this cluster was given. Three kinds
## of check, none of them a threshold: the answer has to be structurally complete, every
## name in it has to be a supplied candidate, and every claim has to quote a gene from its
## own identity's curated panel with the measurement the packet showed. Whether that
## measurement is good enough for the claim is the agent's judgement and stays so.
.valid_final <- function(value, pool, cluster) {
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
  names_available <- .candidate_names(pool, cluster)
  groups <- value$identity_groups %||% list()
  if (!length(groups)) {
    return(FALSE)
  }
  grouped <- character(0)
  members_of <- list()
  for (group in groups) {
    if (!is.list(group) || !nzchar(trimws(as.character(group$identity %||% "")))) {
      return(FALSE)
    }
    members <- group$candidates %||% list()
    if (!length(members)) {
      return(FALSE)
    }
    entries <- vapply(members, as.character, "")
    members_of[[length(members_of) + 1L]] <- entries
    grouped <- c(grouped, entries)
  }
  if (!identical(sort(grouped), sort(as.character(names_available)))) {
    return(FALSE)
  }
  selected <- as.character(value$selected %||% "")
  if (!length(selected) || !nzchar(selected)) {
    return(FALSE)
  }
  if (!identical(selected, UNKNOWN) && !(selected %in% names_available)) {
    return(FALSE)
  }
  home <- character(0)
  for (entries in members_of) {
    if (selected %in% entries) home <- entries
  }
  subtype <- as.character(value$subtype %||% "")
  if (nzchar(subtype)) {
    if (identical(selected, UNKNOWN) || identical(subtype, selected) ||
      !(subtype %in% names_available) || !(subtype %in% home)) {
      return(FALSE)
    }
  }
  others <- vapply(value$co_occurring_identities %||% list(), as.character, "")
  if (length(others)) {
    if (identical(selected, UNKNOWN)) {
      return(FALSE)
    }
    ## A co-occurring identity has to be a DIFFERENT identity. Another name from the
    ## selected candidate's own group is the same cell read twice, not a second population.
    if (!all(others %in% names_available) || any(others %in% home) ||
      any(others == selected)) {
      return(FALSE)
    }
  }
  ## An identity retrieval never offered is reportable but not answerable. Reportable
  ## because the cut is real -- 15 of a median 67 admitted cell types reach the agent, and
  ## until now a menu that missed the answer looked exactly like a menu that did not. Not
  ## answerable because a name with no supplied panel has no `claim_evidence` to quote, so
  ## asserting it would put an unverifiable label in the same column as verified ones.
  ## Naming a candidate here is a contradiction: that one could just be answered with.
  unlisted <- value$unlisted_identity %||% ""
  if (!is.character(unlisted) || length(unlisted) != 1L ||
    as.character(unlisted) %in% names_available) {
    return(FALSE)
  }
  listed <- character(0)
  entry <- .find_candidate(pool, cluster, selected)
  if (!is.null(entry)) {
    listed <- toupper(vapply(entry$markers, function(m) m$gene, ""))
  }
  for (gene in value$support_markers %||% list()) {
    if (!is.character(gene) || !(toupper(as.character(gene)) %in% listed)) {
      return(FALSE)
    }
  }
  .quotes_match(value, pool, cluster)
}

## ---- everything below runs the stage -----------------------------------------
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

## Where each cluster ranks among all clusters for one gene's fold change. Reported in the
## marker table only: it saturates at 1.0 for most of the markers shown, so as evidence it
## separates almost nothing and never reaches the agent.
wide <- dcast(de, gene_key ~ group, value.var = "avg_log2FC")
rel_genes <- wide$gene_key
cluster_columns <- setdiff(names(wide), "gene_key")
rel <- as.matrix(wide[, ..cluster_columns])
rel_pct <- t(matrix(
  apply(rel, 1, function(row) frank(row, ties.method = "average") / length(row)),
  nrow = length(cluster_columns), ncol = length(rel_genes)
))
dimnames(rel_pct) <- list(rel_genes, cluster_columns)

## Source sentences are resolved on this arm. The order is computed from the seed and
## the hash every arm uses, so all of them serve identical sentences without one asking
## another.
UB_SOURCES <- uberon_load(OBO_UBERON)
source_context <- scma_source_context(
  UB_SOURCES, DB_SOURCES, context$species, context$tissue, context$disease
)

.sources_for <- function(candidate, gene, k) {
  records <- scma_source_head(source_context, candidate, gene, k)
  lapply(records, function(record) {
    list(
      pmid = .identifier(record$pmid),
      pmcid = .identifier(record$pmcid),
      sentence = .clip(record$sentence, SOURCE_MAX_CHARS)
    )
  })
}

pool <- .build_pool(scoring, de, native_of, .sources_for)

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
cat(sprintf(
  "==== annotate(r) %s: %d clusters, top %d candidates each, llm=%s ====\n",
  tag, length(clusters_sorted), scoring$top_candidates, llm_status
))
cat(sprintf(
  "  annotator: %s (%s)\n", scma_resolve_model(ANNOTATOR_MODEL), ANNOTATOR_EFFORT
))

## ---- records ----------------------------------------------------------------
.fallback <- function(cluster, frame, reason, detail) {
  if (!nrow(frame)) {
    return(list(
      cluster_id = cluster, annotation = UNKNOWN, subtype = "", state = "",
      co_occurring_identities = list(), confidence = NOT_AVAILABLE,
      rationale = paste(
        "no candidate's curated positive markers are significantly",
        "up-regulated in this cluster"
      ),
      resolution_status = UNRESOLVED, resolution_detail = "empty_candidate_pool",
      annotation_source = "no_candidate", llm_status = reason,
      support_markers = list(), claim_warnings = list(), unlisted_identity = "",
      annotation_qc = QC_UNCHECKED, quality_control = list()
    ))
  }
  top <- frame[order(retrieval_rank)][1]
  list(
    cluster_id = cluster, annotation = as.character(top$candidate),
    subtype = "", state = "", co_occurring_identities = list(),
    confidence = NOT_AVAILABLE,
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
    support_markers = list(), claim_warnings = list(), unlisted_identity = "",
    annotation_qc = QC_UNCHECKED, quality_control = list()
  )
}

## The measured markers reported for one claimed identity: the raised ones,
## best-corroborated first, on the fraction scale the published tables use. The scale
## changes exactly here, at the boundary between what the agent read and what a reader
## reads, and nowhere else.
.measurements <- function(entry, cluster, limit = 30L) {
  ## Raised positives capped at `limit`, plus every detected
  ## exclusion on its own budget. The reason argues from those exclusions, so they belong
  ## in the evidence table beside the rows it names.
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

## ---- the conversation, run for every cluster in lock step --------------------
## Each cluster is one sequential conversation, and the clusters advance together: every
## cluster still talking contributes one prompt to a round, and the round goes out
## through the shared parallel transport. That gives R the same concurrency a thread
## pool would, without a per-cluster thread.
template <- .prompt("cluster_annotator")
runnable <- character(0)
results <- list()
for (cluster in clusters_sorted) {
  cluster_value <- cluster
  frame <- scored[cluster == cluster_value]
  if (identical(pool$clusters[[cluster]]$status, UNSUPPORTED)) {
    results[[cluster]] <- .fallback(cluster, frame[0], llm_status, "empty_candidate_pool")
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

## One server per cluster: "already shown" is a fact about one conversation, and sharing it
## across clusters would let one cluster's reading exhaust another's. The judge that checks
## a cluster is served by the same one, as it is on every arm, so a check is never
## answered with the sentences its own annotation just read.
source_servers <- list()
for (cluster in runnable) {
  server <- scma_source_server(
    source_context,
    batch = SOURCES_PER_MARKER, max_batches = SOURCE_BATCHES_PER_MARKER
  )
  scma_register_packet_sources(server, pool, cluster)
  source_servers[[cluster]] <- server
}

conversation <- list()
for (cluster in runnable) {
  conversation[[cluster]] <- list(
    prompt = .render(template, .cluster_packet(pool, cluster)),
    # Correlates one cluster's turns in the cold-call log and nothing else, so any
    # per-run unique string will do.
    trace_id = substr(
      paste0(format(Sys.time(), "%H%M%S"), basename(tempfile(""))), 1, 16
    ),
    turn = 0L, schema_failures = 0L, repeated = 0L, forced = FALSE,
    seen = list(), transcript = list(), outcome = NULL, audited = NULL,
    attempts = list()
  )
}

active <- runnable
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
        conv$outcome <- list(error = error, turns = conv$turn)
        conversation[[cluster]] <- conv
        next
      }
      conv$schema_failures <- conv$schema_failures + 1L
      if (conv$schema_failures > SCHEMA_RETRIES) {
        conv$outcome <- list(error = error %||% "no response", turns = conv$turn)
        conversation[[cluster]] <- conv
        next
      }
      conv$prompt <- paste0(conv$prompt, sprintf("\n\n# Retry %d\n", conv$schema_failures))
      conversation[[cluster]] <- conv
      still <- c(still, cluster)
      next
    }
    parsed <- scma_parse_json(content)
    if (!is.null(parsed) && .valid_final(parsed, pool, cluster)) {
      # The audit the run already performs, returned in time to be answered. It asserts
      # nothing about the identity: the agent reads the exclusion against the sources it
      # already holds and answers again, keeping or changing the call. Asked once, so a
      # cluster costs at most one extra turn.
      pending <- if (!is.null(conv$audited) || isTRUE(conv$forced)) {
        list()
      } else {
        .binding_exclusions(
          pool, cluster, as.character(parsed$selected %||% "")[1],
          parsed$identity_groups, parsed$co_occurring_identities
        )
      }
      if (length(pending)) {
        conv$audited <- pending
        conv$prompt <- paste0(conv$prompt, .observation_block(
          length(conv$transcript) + 1L,
          list(
            exclusion_audit = pending,
            note = paste0(
              "curated exclusions the cluster detects on the identity you selected ",
              "and on each you reported as co-occurring; an empty list is a measured ",
              "absence, not a gap"
            )
          )
        ))
        conversation[[cluster]] <- conv
        still <- c(still, cluster)
        next
      }
      # Quality control is asked of the answer, not of the conversation, so the clusters
      # that reached one are collected and judged in a batch rather than one at a time.
      conv$audited <- conv$audited %||% list()
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
      ## `sources` is the one tool whose answer changes when it is asked again: it serves
      ## the next sentences rather than the same ones. Short-circuiting it as a repeat is
      ## what used to make "these three did not settle it" an unanswerable request. Its own
      ## per-pair limit is the brake instead.
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
        conv$prompt, .observation_block(length(conv$transcript), observation)
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
        error = "schema violation", turns = conv$turn, transcript = conv$transcript
      )
      conversation[[cluster]] <- conv
      next
    }
    conv$prompt <- paste0(
      conv$prompt,
      sprintf("\n\n# Retry %d: return one valid JSON object.\n", conv$schema_failures)
    )
    conversation[[cluster]] <- conv
    still <- c(still, cluster)
  }
  if (length(to_review)) {
    reviews <- .review_many(
      pool, to_review, pending_final, source_servers, conversation, api_url, api_key
    )
    for (cluster in to_review) {
      conv <- conversation[[cluster]]
      review <- reviews[[cluster]]
      conv$attempts <- c(
        conv$attempts, list(list(final = pending_final[[cluster]], review = review))
      )
      if (isTRUE(review$retry) && !isTRUE(conv$forced) &&
        length(conv$attempts) < JUDGE_MAX_ROUNDS) {
        conv$prompt <- paste0(conv$prompt, .observation_block(
          length(conv$transcript) + 1L, list(quality_control = .qc_feedback(review))
        ))
        conversation[[cluster]] <- conv
        still <- c(still, cluster)
        next
      }
      # Every round produced a complete answer, so the one delivered is the one the
      # judges rated highest; a tie goes to the latest, which saw the most feedback.
      ranks <- vapply(conv$attempts, function(attempt) {
        as.integer(attempt$review$rank %||% 0L)
      }, 0L)
      best <- max(which(ranks == max(ranks)))
      chosen <- conv$attempts[[best]]
      conv$outcome <- list(
        final = chosen$final, turns = conv$turn, transcript = conv$transcript,
        trace_id = conv$trace_id, turn_budget_exhausted = conv$forced,
        exclusion_audit = conv$audited %||% list(),
        quality_control = c(
          list(rounds = length(conv$attempts), selected_round = best), chosen$review
        )
      )
      conversation[[cluster]] <- conv
    }
  }
  active <- still
}

## ---- one post-hoc audit per answer ------------------------------------------
for (cluster in runnable) {
  conv <- conversation[[cluster]]
  outcome <- conv$outcome %||% list(error = "no answer", turns = conv$turn)
  cluster_value <- cluster
  frame <- scored[cluster == cluster_value]
  if (is.null(outcome$final)) {
    record <- .fallback(
      cluster, frame, sprintf("failed:%s", outcome$error %||% ""), "annotation_failed"
    )
    record$candidates <- as.character(.candidate_names(pool, cluster))
    record$candidate_entries <- .plain_entries(cluster)
    record$turns <- as.integer(outcome$turns %||% 0L)
    results[[cluster]] <- record
    next
  }
  final <- outcome$final
  selected <- as.character(final$selected)
  others <- vapply(final$co_occurring_identities %||% list(), as.character, "")
  ## The finer level has to earn its own name. A subtype whose raised markers all belong
  ## to the parent's curated definition is the parent program under a narrower label, so
  ## the reported identity falls back to the parent rather than claiming a refinement the
  ## measurements do not separate. The rejected name keeps its evidence and its role: why
  ## the cluster is NOT the finer type is as auditable as why it is the parent.
  subtype <- as.character(final$subtype %||% "")
  demotion <- .subtype_demotion(pool, cluster, selected, subtype)
  ## The free gate and the judge answer the same question from different evidence -- one
  ## from whether the parent's panel already contains the genes, the other from what the
  ## curated sentences make of them -- and a refinement has to survive both. Either verdict
  ## alone is enough to drop it, and whichever spoke first keeps its own line so the record
  ## says which check the refinement failed.
  review <- outcome$quality_control %||% list()
  judged_out <- as.character(review$demoted_subtype %||% "")
  if (nzchar(judged_out) && !nzchar(demotion$rejected)) {
    verdict <- review$subtype_verdict %||% list()
    demotion <- list(
      rejected = judged_out,
      line = trimws(sprintf(
        "rejected_subtype %s: quality_control %s | %s", judged_out,
        as.character(verdict$verdict %||% "unavailable"),
        as.character(verdict$reason %||% "")
      ))
    )
  }
  if (nzchar(demotion$rejected)) subtype <- ""
  final_for_claims <- final
  final_for_claims$subtype <- subtype
  claims <- .claimed_identities(final_for_claims)
  if (nzchar(demotion$rejected)) {
    claims[[length(claims) + 1L]] <- list("rejected_subtype", demotion$rejected)
  }
  warned <- .claim_warnings(pool, cluster, claims)
  if (nzchar(demotion$line)) {
    warned[[length(warned) + 1L]] <- list(demotion$rejected, demotion$line)
  }
  warnings <- vapply(warned, function(item) as.character(item[[2]]), "")
  role_of <- list()
  for (claim in claims) role_of[[as.character(claim[[2]])]] <- as.character(claim[[1]])
  evidence_of <- list()
  for (item in final$claim_evidence %||% list()) {
    if (is.list(item)) evidence_of[[as.character(item$identity %||% "")]] <- item
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
  entries <- lapply(pool$clusters[[cluster]]$candidates, function(entry) {
    .candidate_entry(
      entry, cluster,
      role_of[[entry$cell_type]] %||% "",
      evidence_of[[entry$cell_type]],
      warning_of[[entry$cell_type]] %||% ""
    )
  })
  if (identical(selected, UNKNOWN)) {
    status <- UNRESOLVED
    detail <- "agent_established_no_identity"
  } else if (length(others)) {
    status <- MIXED
    detail <- "agent_majority_of_several_identities"
  } else {
    status <- RESOLVED
    detail <- "agent_selected"
  }
  results[[cluster]] <- list(
    cluster_id = cluster, annotation = selected,
    subtype = subtype,
    state = as.character(final$state %||% ""),
    co_occurring_identities = as.list(others),
    confidence = as.character(final$confidence),
    rationale = .clip(final$reason, MAX_RATIONALE),
    resolution_status = status, resolution_detail = detail,
    annotation_source = "cluster_annotation", llm_status = "annotated",
    support_markers = lapply(final$support_markers %||% list(), as.character),
    claim_warnings = as.list(warnings),
    ## What the agent says the retrieved menu was missing, when it says so. A column rather
    ## than a sentence buried in the rationale, because "how often did retrieval not offer
    ## the answer" is a rate somebody has to be able to count.
    unlisted_identity = as.character(final$unlisted_identity %||% "")[1],
    candidates = as.character(.candidate_names(pool, cluster)),
    candidate_entries = entries,
    identity_groups = final$identity_groups %||% list(),
    turns = as.integer(outcome$turns %||% 0L),
    tool_calls = lapply(outcome$transcript %||% list(), function(item) item$call),
    turn_budget_exhausted = isTRUE(outcome$turn_budget_exhausted),
    exclusion_audit = outcome$exclusion_audit %||% list(),
    annotation_qc = .qc_status(review, nzchar(demotion$rejected)),
    quality_control = review
  )
}

results <- results[clusters_sorted]
scma_flush_response_cache()

saveRDS(
  list(
    tag = tag, context = context, llm_status = llm_status,
    models = list(
      annotator = scma_resolve_model(ANNOTATOR_MODEL),
      annotator_reasoning_effort = ANNOTATOR_EFFORT
    ),
    results = results,
    transcripts = lapply(conversation, function(conv) conv$transcript)
  ),
  file.path(CACHE, sprintf("%s_annotations.rds", tag))
)
statuses <- vapply(results, function(r) r$resolution_status, "")
turn_counts <- vapply(results, function(r) as.integer(r$turns %||% 0L), 0L)
warned <- sum(vapply(results, function(r) length(r$claim_warnings %||% list()) > 0L, TRUE))
cat(sprintf(
  paste0(
    "[done] %s: %d resolved, %d mixed, %d unresolved of %d clusters; ",
    "turns median %d, max %d; %d clusters carry a claim warning; %d not model-decided\n"
  ),
  tag, sum(statuses == RESOLVED), sum(statuses == MIXED), sum(statuses == UNRESOLVED),
  length(results),
  as.integer(stats::median(turn_counts[turn_counts > 0L] %||% 0L)),
  max(turn_counts, 0L), warned,
  sum(vapply(results, function(r) !identical(r$annotation_source, "cluster_annotation"), TRUE))
))
