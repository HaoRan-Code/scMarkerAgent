#!/usr/bin/env Rscript
# =============================================================================
# test_two_arm_packet_parity.R -- the two arms must build the same evidence.
#
# The two arms run the same algorithm through separate implementations, so what
# can silently drift is not a number but a FIELD or an ordering: one arm sorting a panel
# differently, or dropping a column the other keeps, means the two agents reason from
# different evidence while every other test still passes.
#
# The expectations here are the same literals the other arm's pool test asserts against
# the same synthetic dataset, so a change to one arm's packet has to be made in both.
#
# The pool builders read caches and call a model, so they are not run. The pure functions
# over an already-built pool are evaluated into a private environment and exercised
# directly.
#
#   Usage: Rscript test_two_arm_packet_parity.R
# =============================================================================
suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
})
## Packaged adaptation (the one departure from software/tests/test_two_arm_packet_parity.R):
## the pipeline is the installed copy under inst/scmarkeragent/, handed in by the
## testthat driver, not a source-tree layout resolved from --file.
package <- Sys.getenv("SCMA_PACKAGED_PIPELINE")
stopifnot(dir.exists(package))
## The pool and the loop are two files on this arm as they are two modules on the other,
## so both are read: the packet comes from the first and the answer checks from the second.
## The gate joins them because a marker row now carries its verdict, and a row that
## reported significance from a private copy of the thresholds is the drift this file
## exists to catch. The client joins them because the reply parser lives with the
## transport that records what it parsed, exactly as it does on the other arm.
targets <- file.path(
  package, "rflow",
  c("evidence_gate.R", "llm_client.R", "annotator_pool.R", "cluster_annotation.R")
)

## ---- every helper this arm calls must exist on this arm --------------------------
## Checked first, and statically, because it is the one drift the rest of this file
## cannot see. A helper written on the other arm and only CALLED here still parses, and
## every packet assertion below still passes; the arm errors the first time a cluster
## reaches that branch, which is how `.subtype_exclusive_markers` shipped uncalled for a
## release. Function-call tokens come from R's own parser rather than a pattern, so a
## name inside a comment or a string is not mistaken for a call.
arm <- list.files(file.path(package, "rflow"), pattern = "[.]R$", full.names = TRUE)
private <- function(names) grep("^[.][A-Za-z]", names, value = TRUE)
called <- unique(unlist(lapply(targets, function(target) {
  tokens <- utils::getParseData(parse(target, keep.source = TRUE))
  private(tokens$text[tokens$token == "SYMBOL_FUNCTION_CALL"])
})))
defined <- unique(unlist(lapply(arm, function(target) {
  sub(" <- function", "", regmatches(
    readLines(target, warn = FALSE),
    regexpr("^\\.?[A-Za-z_.]+ <- function", readLines(target, warn = FALSE))
  ))
})))
undefined <- setdiff(called, defined)
if (length(undefined)) {
  stop(sprintf(
    "R arm calls helpers it does not define: %s", paste(undefined, collapse = ", ")
  ))
}

## The same drift in a constant, which the helper check above cannot see: a configuration
## name is looked up only when the branch reading it runs, so an arm that never assigns it
## parses, passes every assertion in this file, and fails on the first cluster that reaches
## the branch. `SUBTYPE_EXCLUSIVE_REQUIRED` reached a 38-cluster run that way -- every LLM
## call was paid for before the post-hoc audit asked for the name. Constants need their own
## check precisely because this file INJECTS the configuration-derived ones below: the
## injection is what lets a pure function be exercised at all, and also what would hide a
## name the arm itself never assigns.
bound <- function(expression) {
  if (!is.call(expression)) {
    return(character(0))
  }
  here <- if (as.character(expression[[1]])[1] %in% c("<-", "=", "<<-") &&
    is.name(expression[[2]])) {
    as.character(expression[[2]])
  } else {
    character(0)
  }
  c(here, unlist(lapply(as.list(expression)[-1], bound)))
}
## Read from the same files as the helper check and resolved against the whole arm, and
## restricted to the upper-case spelling this package gives configuration: a data.table
## expression names its columns as bare symbols too, and `.N` yields one spelled `N`.
read <- unique(unlist(lapply(targets, function(target) {
  tokens <- utils::getParseData(parse(target, keep.source = TRUE))
  grep("^[A-Z][A-Z0-9_]*$", tokens$text[tokens$token == "SYMBOL"], value = TRUE)
})))
assigned <- unique(unlist(lapply(arm, function(target) {
  unlist(lapply(as.list(parse(target)), bound))
})))
unassigned <- setdiff(read, assigned)
if (length(unassigned)) {
  stop(sprintf(
    "R arm reads constants it does not assign: %s", paste(unassigned, collapse = ", ")
  ))
}

## Evaluate only what a pure function needs: `name <- function(...)`, `name <- c(...)`,
## and plain literal constants. Everything else in the scripts reads the filesystem or
## calls a model.
env <- new.env(parent = globalenv())
for (target in targets) {
  for (expression in as.list(parse(target))) {
    if (!is.call(expression) || !identical(as.character(expression[[1]]), "<-")) next
    value <- expression[[3]]
    literal <- is.numeric(value) || is.character(value) || is.logical(value)
    if (literal || (is.call(value) && as.character(value[[1]])[1] %in% c("function", "c"))) {
      try(eval(expression, env), silent = TRUE)
    }
  }
}

## The configuration-derived constants the validators read. Injecting them from the shared
## defaults is itself the assertion that both arms take them from one place.
config <- fromJSON(file.path(package, "config", "defaults.json"), simplifyVector = FALSE)
A <- config$cluster_annotation
env$ANNOTATOR_SCHEMA <- as.character(A$annotator_schema_version)
env$UNKNOWN <- as.character(A$unknown_token)
env$WARNING_GENES_SHOWN <- as.integer(A$warning_genes_shown)
env$SOURCE_GENES_PER_QUERY <- as.integer(A$source_genes_per_query)
env$SOURCES_PER_GENE <- as.integer(A$sources_per_marker)
env$SOURCE_MAX_CHARS <- as.integer(A$source_sentence_max_chars)
env$NEGATIVE_SOURCE_MIN_PCT_IN <- as.numeric(A$negative_source_min_pct_in)
env$CLUSTER_MARKER_ROWS <- as.integer(A$cluster_marker_rows)
env$JUDGE_SOURCES_PUSHED <- as.integer(A$judge_sources_pushed)
env$JUDGE_SOURCES_PER_GENE <- as.integer(A$judge_sources_per_gene)
env$JUDGE_QUERY_TURNS <- as.integer(A$judge_query_turns)
env$SCHEMA_RETRIES <- as.integer(A$schema_retries)
env$POOL_SEARCH_ROWS <- as.integer(A$pool_search_rows)
env$QC_PASSED <- "passed"
env$QC_REVISED <- "passed_after_revision"
env$QC_DEMOTED <- "demoted_to_parent"
env$QC_FAILED <- "failed"
env$QC_UNCHECKED <- "not_checked"
env$NOT_AVAILABLE <- as.character(config$output$not_applicable_token)
env$SIG_LOG2FC <- as.numeric(config$evidence_gate$avg_log2fc_min_exclusive)
env$SIG_PCT <- as.numeric(config$evidence_gate$pct_in_min_exclusive)
env$SIG_PADJ <- as.numeric(config$evidence_gate$padj_max_exclusive)
env$CROSS_SPECIES <- isTRUE(config$features$cross_species_markers)
env$PROMPT_DIR <- file.path(package, "prompts")
env$PANEL_READING <- "panel_reading"
env$PANEL_READING_PLACEHOLDER <- "{{PANEL_READING}}"

## ---- the two things the render has to get right ---------------------------------
## Both judges receive the shared panel-reading block, so neither can drift into its own
## reading of what a raised marker is; and a curated sentence carrying a quote survives
## into the packet. R reads backslashes in a regex REPLACEMENT, so `\"` from `toJSON` was
## eaten and the packet reached the model as text that no longer parsed as JSON -- on the
## two clusters whose sources happened to quote, with nothing in the log to show for it.
local({
  shared <- trimws(paste(
    readLines(file.path(env$PROMPT_DIR, "panel_reading.txt"), warn = FALSE),
    collapse = "\n"
  ))
  for (name in c("cluster_annotator", "subtype_judge", "annotation_judge")) {
    rendered <- env$.prompt(name)
    stopifnot(
      !grepl(env$PANEL_READING_PLACEHOLDER, rendered, fixed = TRUE),
      grepl(shared, rendered, fixed = TRUE)
    )
  }
  quoted <- list(sources = list(GENE = "a marker of the \"inflammatory\" subset"))
  packet <- env$.render("head {{EVIDENCE_PACKET_JSON}}", quoted)
  stopifnot(identical(
    fromJSON(sub("^head ", "", packet), simplifyVector = FALSE),
    quoted
  ))
})

## ---- one reply, two objects ------------------------------------------------------
## A template asks for one object and a model occasionally sends two: a query, and then the
## answer it would have given without waiting for the observation. Reading the whole string
## parsed as neither, and on human_liver c24 that cost the cluster its check -- filed
## `not_checked`, which reads as "no judge was reachable" when a judge had plainly replied.
## The FIRST object wins, so the question is honoured and the self-answer is dropped.
local({
  both <- paste0(
    '{"action":"query","tool":"sources","args":{"genes":["SI"]}}\n',
    '{"verdict":"supported","reason":"answered without waiting"}'
  )
  first <- env$scma_parse_json(both)
  stopifnot(identical(first$action, "query"), is.null(first$verdict))
  ## And the cases the previous reader already handled have to keep working: a fenced block,
  ## a brace inside a string, an escaped quote.
  one <- env$scma_parse_json('{"verdict":"supported","reason":"a } brace in a string"}')
  stopifnot(identical(one$verdict, "supported"),
            identical(one$reason, "a } brace in a string"))
  fenced <- env$scma_parse_json('```json\n{"verdict":"contradicted"}\n```')
  stopifnot(identical(fenced$verdict, "contradicted"))
  stopifnot(identical(env$scma_parse_json('{"s":"say \\"hi\\" }"}')$s, 'say "hi" }'))
  stopifnot(is.null(env$scma_parse_json("not json at all")))
})

stopifnot(
  is.function(env$.cluster_packet),
  is.function(env$.run_tool),
  is.function(env$.claim_warnings),
  is.function(env$.valid_final),
  is.function(env$.is_raised),
  is.function(env$.run_judge_tool),
  identical(env$TOOLS, c(
    "sources", "gene_across_clusters", "candidates_with_gene", "pool_search"
  )),
  # A judge gets fewer tools than the annotator, and that is the gate: `candidates_with_gene`
  # and `pool_search` both name identities other than the one under test, and a judge holding
  # a menu it could prefer from is a second annotator.
  identical(env$JUDGE_TOOLS, c("sources", "gene_across_clusters")),
  identical(env$PCT_DECIMALS, 1L),
  identical(env$QUOTE_TOLERANCE, 0.15)
)

## ---- the same synthetic dataset the shared fixture builds ----------------------
## `significant` follows from the row: the synthetic DE table puts every padj at 1e-9, so
## the gate reduces here to a positive fold change over a tenth of the cluster, and the
## fixture states the outcome rather than recomputing it. `specificity` is a constant on
## this dataset -- its job in the packet is to be a separate axis from `n_pub`, and one
## value is enough to assert the arms carry it.
marker <- function(gene, polarity, n_pub, tier, pct_in, pct_out, lfc, significant) {
  list(
    gene = gene, polarity = polarity, n_pub = as.integer(n_pub), tier = tier,
    pct_in = pct_in, pct_out = pct_out, avg_log2FC = lfc,
    significant = significant, specificity = 0.5
  )
}
enterocyte <- list(
  cell_type = "enterocyte", retrieval_rank = 1L,
  markers = list(
    marker("VIL1", "positive", 53, "high", 70.8, 72.1, -0.22, FALSE),
    marker("SI", "positive", 34, "high", 86.5, 50.1, 0.01, TRUE),
    marker("FABP2", "positive", 24, "high", 92.4, 68.5, -0.10, FALSE),
    marker("RBP2", "positive", 3, "low", 95.8, 45.3, 1.43, TRUE),
    marker("PTPRC", "negative", 1, "low", 0.0, 1.2, -0.04, FALSE)
  ),
  exclusion_sources = list(),
  unmeasured_curated_genes = list(count = 0L, genes = list()),
  program = list(median_in = 0.6, median_out = 0.4)
)
endocrine <- list(
  cell_type = "enteroendocrine cell", retrieval_rank = 2L,
  markers = list(
    marker("CHGA", "positive", 160, "high", 2.3, 4.1, -2.41, FALSE),
    marker("NEUROD1", "positive", 15, "high", 0.0, 2.2, -0.22, FALSE),
    marker("PROX1", "positive", 4, "high", 77.7, 21.7, 1.19, TRUE)
  ),
  exclusion_sources = list(),
  unmeasured_curated_genes = list(count = 1L, genes = list("MISSING1")),
  program = list(median_in = 0.6, median_out = 0.4)
)
tuft <- list(
  cell_type = "tuft cell", retrieval_rank = 3L,
  markers = list(
    marker("POU2F3", "positive", 40, "high", 1.0, 40.0, -3.00, FALSE),
    marker("PROX1", "positive", 1, "low", 77.7, 21.7, 1.19, TRUE)
  ),
  exclusion_sources = list(),
  unmeasured_curated_genes = list(count = 0L, genes = list()),
  program = list(median_in = 0.6, median_out = 0.4)
)
pool <- list(
  context = list(
    species = "Human", tissue = "intestine", disease = "Normal",
    development_stage = env$NOT_AVAILABLE, clusters_in_dataset = 2L
  ),
  clusters = list(`1` = list(
    cluster_id = "1", status = "pool", n_cells = 288L,
    candidates = list(enterocyte, endocrine, tuft)
  )),
  gene_clusters = list(FABP2 = list(
    list(cluster_id = "1", pct_in = 92.4, pct_out = 68.5),
    list(cluster_id = "2", pct_in = 10.0, pct_out = 80.0)
  )),
  gene_carriers = list(PROX1 = list(
    list(cell_type = "enteroendocrine cell", n_pub = 4L, tier = "high", polarity = "positive"),
    list(cell_type = "tuft cell", n_pub = 1L, tier = "low", polarity = "positive")
  )),
  native_gene = list(FABP2 = "FABP2", PROX1 = "PROX1")
)

## ---- the panel builder ---------------------------------------------------------
## The fixture above hardcodes the panel the agent reads; this rebuilds it from the same
## curated rows and DE measurements the production caches hold, so the two arms are shown
## to derive that panel and not just to format it. The curated table carries columns named
## `candidate` and `gene` exactly as production does, because a local sharing one of those
## names is evaluated as the whole column inside `.(...)` and turns the lookup into a
## cartesian join whose first row then becomes every marker's polarity and n_pub.
curated <- data.table::data.table(
  candidate = c(rep("enterocyte", 5), rep("enteroendocrine cell", 4), rep("tuft cell", 2)),
  gene = c(
    "VIL1", "SI", "FABP2", "RBP2", "PTPRC",
    "CHGA", "NEUROD1", "PROX1", "MISSING1",
    "POU2F3", "PROX1"
  ),
  marker_polarity = c(
    "positive", "positive", "positive", "positive", "negative",
    "positive", "positive", "positive", "positive",
    "positive", "positive"
  ),
  n_pub = c(53L, 34L, 24L, 3L, 1L, 160L, 15L, 4L, 2L, 40L, 1L),
  tier = c(
    "high", "high", "high", "low", "low",
    "high", "high", "high", "medium",
    "high", "low"
  )
)
curated[, gene_key := toupper(gene)]
data.table::setkey(curated, candidate, gene_key)
de_genes <- c("FABP2", "VIL1", "SI", "RBP2", "CHGA", "NEUROD1", "PROX1", "POU2F3", "PTPRC")
de_pct_in <- c(92.4, 70.8, 86.5, 95.8, 2.3, 0.0, 77.7, 1.0, 0.0)
de_pct_out <- c(68.5, 72.1, 50.1, 45.3, 4.1, 2.2, 21.7, 40.0, 1.2)
de_lfc <- c(-0.10, -0.22, 0.01, 1.43, -2.41, -0.22, 1.19, -3.00, -0.04)
## One padj for every row, as in the shared fixture, so the gate reduces here to a
## positive fold change over a tenth of the cluster and the two arms are compared on
## building the field at all rather than on a threshold sweep.
de_padj <- rep(1e-9, length(de_genes))
de_keys <- paste("1", de_genes, sep = "\x1f")
measured <- list(
  pct_in = setNames(de_pct_in, de_keys),
  pct_out = setNames(de_pct_out, de_keys),
  avg_log2FC = setNames(de_lfc, de_keys),
  padj = setNames(de_padj, de_keys)
)
native_of <- setNames(de_genes, de_genes)
specificity <- as.list(setNames(rep(0.5, length(de_genes)), de_genes))
built <- env$.panel_rows(
  "enterocyte",
  c("VIL1", "SI", "FABP2", "RBP2"), "PTPRC",
  "1", measured, curated, native_of, specificity
)
## `significant` is compared too, because it is derived rather than copied: an arm that
## read the gate from its own thresholds would agree on every other field and disagree
## here, which is exactly the drift this file exists to catch.
compared <- c("gene", "polarity", "n_pub", "tier", "significant", "specificity")
stopifnot(identical(
  lapply(built, function(row) row[compared]),
  lapply(enterocyte$markers, function(row) row[compared])
))
## MISSING1 is curated for the identity and absent from the DE table, so it contributes no
## row: a gene with no measurement supports nothing and refutes nothing.
stopifnot(length(env$.panel_rows(
  "enteroendocrine cell", c("CHGA", "NEUROD1", "PROX1", "MISSING1"), character(0),
  "1", measured, curated, native_of, specificity
)) == 3L)

## ---- the opening packet --------------------------------------------------------
packet <- env$.cluster_packet(pool, "1")
stopifnot(identical(packet$query$cluster_id, "1"))
stopifnot(identical(packet$query$cells_in_cluster, 288L))
stopifnot(identical(
  vapply(packet$candidates, function(entry) entry$cell_type, ""),
  c("enterocyte", "enteroendocrine cell", "tuft cell")
))
## The same measurements keyed by gene: ordered by contrast rather than by publication
## support, and naming every candidate here that claims the gene. Both orderings are
## asserted because the two arms build them from different sorts -- one over a dict, one
## over a named list -- and a tie broken differently would put the arms on two packets.
cluster_markers <- packet$cluster_markers
stopifnot(identical(
  vapply(cluster_markers, function(row) row$gene, ""),
  c("PROX1", "RBP2", "SI")
))
stopifnot(identical(
  names(cluster_markers[[1]]),
  c("gene", "pct_in", "pct_out", "specificity", "claimed_by")
))
stopifnot(identical(
  lapply(cluster_markers[[1]]$claimed_by, function(claim) claim$cell_type),
  list("enteroendocrine cell", "tuft cell")
))
stopifnot(identical(cluster_markers[[1]]$claimed_by[[1]]$n_pub, 4L))
stopifnot(identical(cluster_markers[[2]]$claimed_by[[1]]$cell_type, "enterocyte"))
# A candidate shows its complete panel and nothing derived from it. No gate column:
# whether the gene is up-regulated is what avg_log2FC already says.
stopifnot(identical(
  names(packet$candidates[[1]]),
  c("cell_type", "markers", "exclusion_sources", "unmeasured_curated_genes", "program")
))
## `tier` is not among them: the resource derives it from `n_pub` by a fixed rule, so it
## is the same fact twice, and the packet's field list never described it. The row still
## carries it for `marker_evidence.csv`, which is why this asserts the packet projection
## rather than the row.
stopifnot(identical(
  names(packet$candidates[[1]]$markers[[1]]),
  c("gene", "polarity", "n_pub", "specificity",
    "pct_in", "pct_out", "avg_log2FC", "significant")
))
# Best-corroborated first, so an agent reading from the top meets the gene the resource
# says defines the identity before any incidental marker.
stopifnot(identical(
  vapply(packet$candidates[[1]]$markers, function(m) m$gene, ""),
  c("VIL1", "SI", "FABP2", "RBP2", "PTPRC")
))
# Percent, on the DE table's own scale.
stopifnot(identical(packet$candidates[[1]]$markers[[3]]$pct_in, 92.4))
# The fixture's one negative marker is undetected, so no provenance ships with it and the
# packet still carries no sentence. The mirror case is asserted against the builder below.
stopifnot(identical(packet$candidates[[1]]$exclusion_sources, list()))
stopifnot(!grepl("sentence", toJSON(packet, auto_unbox = TRUE), fixed = TRUE))

## ---- provenance behind a detected exclusion ------------------------------------
## Same threshold, same shape, same call as the other arm: the sentences follow the
## measurement, so raising the negative marker above the threshold is what ships them.
fake_sources <- function(candidate, gene, k) {
  list(list(pmid = "1", pmcid = "PMC1",
            sentence = sprintf("%s is absent from %s.", gene, candidate)))
}
raised_negative <- enterocyte$markers
raised_negative[[5]]$pct_in <- 66.1
stopifnot(identical(
  env$.exclusion_sources("enterocyte", raised_negative, fake_sources),
  list(PTPRC = list(list(pmid = "1", pmcid = "PMC1",
                         sentence = "PTPRC is absent from enterocyte.")))
))
## Under the threshold nothing ships, and a positive marker never does whatever its value.
stopifnot(identical(
  env$.exclusion_sources("enterocyte", enterocyte$markers, fake_sources), list()
))

## ---- the packet quality control is judged from -------------------------------
## The narrowing is the contract, so it is asserted as an exact key set rather than by
## spot-checking: a judge handed `candidates` or the retrieval order could prefer a
## different label, and the disagreement would say nothing about the label under test.
## The same assertion is made against the other arm's builder.
judge <- env$.judge_packet(pool, "1", "enterocyte", fake_sources)
stopifnot(setequal(
  names(judge),
  c("query", "label_under_test", "panel", "unmeasured_curated_genes", "sources")
))
stopifnot(identical(judge$label_under_test, "enterocyte"))
stopifnot(identical(judge$query$cluster_id, "1"))
## The whole definition ordered by publication support, raised or not: VIL1 and FABP2 are
## detected in most of the cluster AND in most cells outside it, and that they fail to
## separate is exactly what a judge has to see. The undetected negative marker is left
## out -- an exclusion nobody detects excludes nothing.
stopifnot(identical(
  vapply(judge$panel, function(row) row$gene, ""), c("VIL1", "SI", "FABP2", "RBP2")
))
stopifnot(setequal(names(judge$sources), c("VIL1", "SI", "FABP2", "RBP2")))
## With a parent the question changes from "is this plausible" to "is this separable",
## so the parent's panel travels with it and the genes absent from it are named.
finer <- env$.judge_packet(pool, "1", "tuft cell", fake_sources, parent = "enterocyte")
stopifnot(identical(finer$parent_label, "enterocyte"))
stopifnot(identical(finer$raised_and_absent_from_parent_panel, "PROX1"))
stopifnot(identical(
  vapply(finer$parent_panel, function(row) row$gene, ""),
  c("VIL1", "SI", "FABP2", "RBP2")
))
## POU2F3 defines tuft cells and sits at 1.0 here: the row that decides this is not one.
stopifnot(identical(
  vapply(finer$panel, function(row) row$gene, ""), c("POU2F3", "PROX1")
))

stopifnot(identical(
  env$.exclusion_sources("enterocyte", raised_negative, NULL),
  list(PTPRC = list())
))

## ---- raised is exactly a positive fold change ----------------------------------
stopifnot(
  isTRUE(env$.is_raised(list(avg_log2FC = 0.01))),
  !isTRUE(env$.is_raised(list(avg_log2FC = 0))),
  !isTRUE(env$.is_raised(list(avg_log2FC = -0.22))),
  !isTRUE(env$.is_raised(list()))
)

## ---- tools ---------------------------------------------------------------------
across <- env$.run_tool(pool, "1", "gene_across_clusters", list(gene = "fabp2"))
stopifnot(identical(across$gene, "FABP2"))
stopifnot(identical(
  across$clusters,
  list(
    list(cluster_id = "1", pct_in = 92.4, pct_out = 68.5),
    list(cluster_id = "2", pct_in = 10.0, pct_out = 80.0)
  )
))
carriers <- env$.run_tool(pool, "1", "candidates_with_gene", list(gene = "PROX1"))
stopifnot(identical(
  vapply(carriers$candidates, function(row) row$cell_type, ""),
  c("enteroendocrine cell", "tuft cell")
))
stub <- function(cell_type, gene, k) {
  list(list(pmid = "1", pmcid = "PMC1", sentence = sprintf("%s marks %s.", gene, cell_type)))
}
sources <- env$.run_tool(
  pool, "1", "sources", list(candidate = "enterocyte", genes = list("fabp2")), stub
)
stopifnot(identical(sources$sources$FABP2[[1]]$pmcid, "PMC1"))
stopifnot(isFALSE(sources$truncated))
many <- env$.run_tool(
  pool, "1", "sources",
  list(
    candidate = "enterocyte",
    genes = as.list(sprintf("G%d", seq_len(env$SOURCE_GENES_PER_QUERY + 4L)))
  ),
  stub
)
stopifnot(length(many$sources) == env$SOURCE_GENES_PER_QUERY, isTRUE(many$truncated))
## The budget is what the template promises, and 64% of the calls in the shipped run asked
## for exactly it -- so a cluster wanting more paid a second turn, which re-sends the whole
## 26k-token packet to buy a few sentences.
stopifnot(identical(env$SOURCE_GENES_PER_QUERY, 12L))
unknown <- env$.run_tool(pool, "1", "gene_stats", list(genes = list("ALB")))
stopifnot(grepl("unknown tool", unknown$error))

## ---- warnings ------------------------------------------------------------------
## Read over the whole measured positive panel, not over the gene the answer quoted:
## a claim can be written on one incidental marker while the genes the resource says
## define that cell sit flat. Each line arrives paired with its identity, because 212
## curated names carry a colon of their own and parsing the name back out of the line
## filed those warnings against identities that do not exist.
lines <- env$.claim_warnings(pool, "1", list(list("cooc", "enteroendocrine cell")))
stopifnot(identical(
  lines,
  list(list(
    "enteroendocrine cell",
    "cooc enteroendocrine cell: not_raised 2/3 | top_n_pub: CHGA(2.3/4.1),NEUROD1(0.0/2.2) | +0 more"
  ))
))
## A correct call on a dominant population is flagged too, and not corrected: three of
## these four positive markers really are unraised, and the high pct_out beside the high
## pct_in is in the line for a reader to weigh.
dominant <- env$.claim_warnings(pool, "1", list(list("selected", "enterocyte")))
stopifnot(identical(
  dominant,
  list(list(
    "enterocyte",
    "selected enterocyte: not_raised 2/4 | top_n_pub: VIL1(70.8/72.1),FABP2(92.4/68.5) | +0 more"
  ))
))
## The identity travels verbatim even when it holds the delimiter the old parse cut at.
colon_name <- "BCR::ABL1 positive primitive cell"
colon_pool <- list(clusters = list(`1` = list(candidates = list(list(
  cell_type = colon_name,
  markers = list(list(
    gene = "ABL1", polarity = "positive",
    avg_log2FC = -1.0, pct_in = 1.0, pct_out = 2.0
  ))
)))))
colon <- env$.claim_warnings(colon_pool, "1", list(list("selected", colon_name)))
stopifnot(identical(
  colon,
  list(list(
    colon_name,
    sprintf("selected %s: not_raised 1/1 | top_n_pub: ABL1(1.0/2.0) | +0 more", colon_name)
  ))
))

## ---- the answer -----------------------------------------------------------------
final <- list(
  action = "final", schema_version = env$ANNOTATOR_SCHEMA,
  selected = "enterocyte", subtype = "", state = "",
  co_occurring_identities = list(),
  claim_evidence = list(list(
    identity = "enterocyte", decisive_gene = "FABP2", pct_in = 92.4, pct_out = 68.5
  )),
  identity_groups = list(
    list(identity = "absorptive", candidates = list("enterocyte")),
    list(identity = "secretory", candidates = list("enteroendocrine cell", "tuft cell"))
  ),
  support_markers = list("FABP2"), confidence = "high",
  reason = "FABP2 is detected in 92.4% of these cells."
)
stopifnot(isTRUE(env$.valid_final(final, pool, "1")))

# A gene the resource never tied to that cell type cannot carry the claim.
borrowed <- final
borrowed$claim_evidence <- list(list(
  identity = "enterocyte", decisive_gene = "CHGA", pct_in = 2.3, pct_out = 4.1
))
stopifnot(!isTRUE(env$.valid_final(borrowed, pool, "1")))

# A number the agent typed rather than read is not an answer.
misquoted <- final
misquoted$claim_evidence <- list(list(
  identity = "enterocyte", decisive_gene = "FABP2", pct_in = 99.0, pct_out = 68.5
))
stopifnot(!isTRUE(env$.valid_final(misquoted, pool, "1")))
rounded <- final
rounded$claim_evidence <- list(list(
  identity = "enterocyte", decisive_gene = "FABP2", pct_in = 92.44, pct_out = 68.5
))
stopifnot(isTRUE(env$.valid_final(rounded, pool, "1")))

# Every claimed identity needs its own entry, co-occurring ones included.
mixture <- final
mixture$co_occurring_identities <- list("enteroendocrine cell")
mixture$identity_groups <- list(
  list(identity = "absorptive", candidates = list("enterocyte")),
  list(identity = "endocrine", candidates = list("enteroendocrine cell")),
  list(identity = "tuft", candidates = list("tuft cell"))
)
stopifnot(!isTRUE(env$.valid_final(mixture, pool, "1")))
mixture$claim_evidence <- c(mixture$claim_evidence, list(list(
  identity = "enteroendocrine cell", decisive_gene = "PROX1", pct_in = 77.7, pct_out = 21.7
)))
stopifnot(isTRUE(env$.valid_final(mixture, pool, "1")))

# The grouping has to cover every supplied candidate exactly once.
ungrouped <- final
ungrouped$identity_groups <- list(list(identity = "absorptive", candidates = list("enterocyte")))
stopifnot(!isTRUE(env$.valid_final(ungrouped, pool, "1")))

# Support markers are genes listed for the selected identity.
stolen <- final
stolen$support_markers <- list("CHGA")
stopifnot(!isTRUE(env$.valid_final(stolen, pool, "1")))

# A query turn and a final turn are told apart.
stopifnot(
  isTRUE(env$.valid_query(list(action = "query", tool = "sources", args = list()))),
  !isTRUE(env$.valid_query(final))
)

# An observation appends only itself: the packet must stay byte-identical from turn to turn.
block <- env$.observation_block(1L, across)
stopifnot(grepl("^\n\n# Observation 1\n", block))
stopifnot(!grepl("cells_in_cluster", block, fixed = TRUE))

## ---- the finer name has to earn itself -----------------------------------------
## Whether a marker is exclusive is a question about the PARENT, so the same refinement
## answers differently under two of them. Under enteroendocrine, the only tuft marker this
## cluster raises is PROX1, which the parent's own definition also contains: that is the
## parent read at a narrower name. Under enterocyte, whose definition does not contain it,
## the same PROX1 is the refinement's own.
stopifnot(identical(
  env$.subtype_exclusive_markers(pool, "1", "enteroendocrine cell", "tuft cell"),
  character(0)
))
stopifnot(identical(
  env$.subtype_exclusive_markers(pool, "1", "enterocyte", "tuft cell"), "PROX1"
))
## POU2F3 belongs to tuft alone and is measured at 1.0 against 40.0 here, so it counts for
## nothing: the comparison is over the markers the cluster RAISES.
stopifnot(!("POU2F3" %in% env$.subtype_exclusive_markers(pool, "1", "enterocyte", "tuft cell")))

## ---- the exclusion audit, both sides -------------------------------------------
## Facts only, and the demoted identity is reported even when it carries nothing: an
## empty list there is a measured absence, and absence is the half of the comparison the
## selected identity's own rows cannot supply.
detected <- pool
detected$clusters[["1"]]$candidates[[1]]$markers[[5]]$pct_in <- 66.1
audit <- env$.binding_exclusions(
  detected, "1", "enterocyte",
  list(list(identity = "absorptive", candidates = list("enterocyte"))),
  list("tuft cell")
)
stopifnot(identical(
  vapply(audit, function(block) block$role, ""), c("selected", "co_occurring")
))
stopifnot(identical(audit[[1]]$detected_exclusions[[1]]$gene, "PTPRC"))
stopifnot(identical(audit[[1]]$detected_exclusions[[1]]$curated_negative_for, "enterocyte"))
stopifnot(identical(audit[[2]]$identity, "tuft cell"))
stopifnot(identical(length(audit[[2]]$detected_exclusions), 0L))
## Undetected, nothing binds and no observation is raised at all.
stopifnot(identical(
  env$.binding_exclusions(pool, "1", "enterocyte", list(), list("tuft cell")), list()
))

## ---- how a label leaves quality control --------------------------------------
## The same five words on both arms, decided the same way. `not_checked` is the one that
## matters most: an offline replay reaching no judge must not report every cluster as
## having survived a check nobody ran. The other arm asserts this against the same
## cases.
stopifnot(identical(env$.qc_status(list(checked = FALSE), FALSE), "not_checked"))
stopifnot(identical(env$.qc_status(list(), FALSE), "not_checked"))
carried <- list(
  checked = TRUE, selected_round = 1L, judged_label = "enterocyte",
  annotation_verdict = list(verdict = "supported")
)
stopifnot(identical(env$.qc_status(carried, FALSE), "passed"))
stopifnot(identical(
  env$.qc_status(modifyList(carried, list(selected_round = 2L)), FALSE),
  "passed_after_revision"
))
stopifnot(identical(env$.qc_status(carried, TRUE), "demoted_to_parent"))
stopifnot(identical(
  env$.qc_status(
    modifyList(carried, list(annotation_verdict = list(verdict = "contradicted"))), FALSE
  ),
  "failed"
))
## A level wrong in a nameable direction still identified the lineage, so it outranks a
## contradiction when the rounds have to choose which answer ships.
stopifnot(env$.verdict_rank("supported") > env$.verdict_rank("too_specific"))
stopifnot(env$.verdict_rank("too_specific") > env$.verdict_rank("insufficient_evidence"))
stopifnot(env$.verdict_rank("insufficient_evidence") > env$.verdict_rank("contradicted"))
stopifnot(identical(env$.verdict_rank("contradicted"), env$.verdict_rank("nonsense")))

## What goes back to the annotator: the verdict and the genes, and nothing that would
## re-send evidence it already holds.
feedback <- env$.qc_feedback(list(
  demoted_subtype = "tuft cell", judged_label = "enterocyte",
  subtype_verdict = list(
    verdict = "not_separable", conflicting_markers = list("PROX1"),
    reason = "PROX1 is curated for the parent as well."
  ),
  annotation_verdict = list(
    verdict = "supported", conflicting_markers = list(), reason = "SI carries it."
  )
))
stopifnot(identical(feedback$subtype$verdict, "not_separable"))
stopifnot(identical(feedback$subtype$name, "tuft cell"))
stopifnot(identical(feedback$annotation$label, "enterocyte"))
stopifnot(grepl("candidates you already hold", feedback$note, fixed = TRUE))
stopifnot(!grepl(
  "unmeasured_curated_genes", toJSON(feedback, auto_unbox = TRUE), fixed = TRUE
))

## ---- the R review loop actually runs -----------------------------------------
## Parsing is not running. `SUBTYPE_EXCLUSIVE_REQUIRED` parsed fine for a release and
## then stopped a 38-cluster run on the first cluster that reached its branch, after
## every model call had been paid for. So the loop is driven here with a stubbed model:
## the assertions are modest, but every line of it executes.
env$UNKNOWN <- as.character(A$unknown_token)
env$scma_cached_call_llm_many <- function(prompts, api_url, api_key, ...) {
  lapply(prompts, function(prompt) {
    verdict <- if (grepl("Subtype Quality Control", prompt, fixed = TRUE)) {
      "not_separable"
    } else {
      "supported"
    }
    list(sprintf('{"verdict":"%s","reason":"stub","conflicting_markers":["PROX1"]}', verdict), NULL)
  })
}
env$.prompt <- function(name) {
  paste(
    readLines(file.path(package, "prompts", paste0(name, ".txt")), warn = FALSE),
    collapse = "\n"
  )
}
answer <- list(selected = "enterocyte", subtype = "tuft cell")
conversation <- list(`1` = list(trace_id = "t", turn = 1L))
## Keyed by cluster, because the R arm judges a batch of clusters at once and each one is
## served by its own stream -- sharing one would let a cluster's reading exhaust another's.
reviews <- env$.review_many(
  pool, "1", list(`1` = answer), list(`1` = fake_sources), conversation, "u", "k"
)
review <- reviews[["1"]]
## The refinement was refused, so it is dropped and the parent is what gets checked.
stopifnot(identical(review$demoted_subtype, "tuft cell"))
stopifnot(identical(review$judged_label, "enterocyte"))
stopifnot(isTRUE(review$checked))
## The parent it fell back to is supported, so this is a finished answer, not a failure
## to re-open: nothing is retried and no second turn is paid for. The column still records
## that a refinement was dropped.
stopifnot(isTRUE(review$passed))
stopifnot(isFALSE(review$retry))
stopifnot(identical(review$rank, 3L))
stopifnot(identical(env$.qc_status(
  c(review, list(selected_round = 1L)), nzchar(review$demoted_subtype)
), "demoted_to_parent"))

cat("two-arm packet parity contract OK\n")
