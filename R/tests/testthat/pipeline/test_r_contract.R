#!/usr/bin/env Rscript

## Packaged adaptation (the one departure from software/tests/test_r_contract.R):
## the pipeline is the installed copy under inst/scmarkeragent/, handed in by the
## testthat driver, not a source-tree layout resolved from --file.
rflow <- file.path(Sys.getenv("SCMA_PACKAGED_PIPELINE"), "rflow")
stopifnot(dir.exists(rflow))
for (path in list.files(rflow, pattern = "[.]R$", full.names = TRUE)) {
  parse(path)
}

source(file.path(rflow, "config.R"))
stopifnot(CFG$schema_version == "scmarkeragent-defaults-v1")
stopifnot(isTRUE(CFG$features$cluster_annotation))
# The retrieval admission controls are a contract between the arms, not constants of the
# R arm: assert that R reads exactly what the shared configuration declares, so a sweep
# of the config moves both arms together instead of leaving this test pinned to a number.
stopifnot(TOP_CANDIDATES == as.integer(CFG$retrieval$top_candidates))
stopifnot(RU_MIN_HITS == as.integer(CFG$retrieval$min_significant_hits))
stopifnot(RU_MIN_POOL_FLOOR == as.integer(CFG$retrieval$min_pool_floor))
stopifnot(TOP_CANDIDATES >= 1L, RU_MIN_HITS >= 1L, RU_MIN_POOL_FLOOR >= 1L)
stopifnot(NOT_AVAILABLE == "N/A")
stopifnot(na_display(NULL) == "N/A", na_display("") == "N/A", na_display("CL:1") == "CL:1")

# The R arm must carry no Cell Ontology resource, module or graph operation.
stopifnot(is.null(CFG$resources$cell_ontology))
stopifnot(!exists("OBO"))
stopifnot(!file.exists(file.path(rflow, "cl_ontology.R")))
stopifnot(!file.exists(file.path(rflow, "decisive_conflict.R")))
stopifnot(!exists(paste0("DB_", "CURATION")))

# The retrieval gate carries no AUC term and reads the genome-wide BH column. Its fold
# change bound is zero: "up-regulated one versus rest" is the whole of what it means, and
# both arms and the annotator's own wording read that one number from this configuration.
source(file.path(rflow, "evidence_gate.R"))
stopifnot(!exists("SIG_AUC"))
stopifnot(SIG_LOG2FC == 0, SIG_PCT == 0.1, SIG_PADJ == 0.05)
stopifnot(SIG_LOG2FC == as.numeric(CFG$evidence_gate$avg_log2fc_min_exclusive))
stopifnot(isTRUE(sig_pass(0.5, 0.5, 0.001)))
stopifnot(isTRUE(sig_pass(0.05, 0.5, 0.001)))
stopifnot(!isTRUE(sig_pass(0, 0.5, 0.001)))
stopifnot(!isTRUE(sig_pass(-0.2, 0.5, 0.001)))
stopifnot(!isTRUE(sig_pass(0.5, 0.05, 0.001)))
stopifnot(!isTRUE(sig_pass(0.5, 0.5, 0.5)))
stopifnot(!isTRUE(sig_pass(NA, 0.5, 0.001)))

Sys.unsetenv(c("OPENAI_API_KEY", "OPENAI_BASE_URL"))
Sys.setenv(
  SCMA_OFFLINE = "0",
  OPENAI_API_KEY = "unit-test-openai-key-not-secret"
)
credentials <- .resolve_credentials()
stopifnot(
  credentials$key == "unit-test-openai-key-not-secret",
  credentials$url == "https://api.openai.com/v1/responses",
  credentials$source == "environment",
  .provider_mode() == "official_openai"
)
Sys.setenv(
  OPENAI_BASE_URL = "https://openai-compatible.invalid/v1"
)
credentials <- .resolve_credentials()
stopifnot(
  credentials$url == "https://openai-compatible.invalid/v1/responses",
  .provider_mode() == "custom_openai_base"
)
Sys.setenv(SCMA_OFFLINE = "1")
credentials <- .resolve_credentials()
stopifnot(
  credentials$key == "",
  credentials$url == "",
  credentials$source == "offline"
)
Sys.setenv(SCMA_OFFLINE = "0")

# Both arms render the same templates from the same directory. One agent decides a
# cluster, so there is one annotator prompt; the other two check what it decided and are
# handed a single label with its own panel, never the candidates together.
stopifnot(dir.exists(PROMPT_DIR))
stopifnot(setequal(
  basename(list.files(PROMPT_DIR, pattern = "[.]txt$")),
  c(
    "cluster_annotator.txt", "subtype_judge.txt", "annotation_judge.txt",
    "panel_reading.txt", "cross_species_note.txt"
  )
))
# The evidence packet has to be the LAST thing in the prompt. Everything before it is
# identical across every call, which is the only part a provider can serve from a cached
# prefix; it is also why a tool result is appended rather than folded into the packet.
text <- paste(
  readLines(file.path(PROMPT_DIR, "cluster_annotator.txt"), warn = FALSE),
  collapse = "\n"
)
stopifnot(grepl("\\{\\{EVIDENCE_PACKET_JSON\\}\\}", text))
stopifnot(nchar(trimws(sub(".*\\{\\{EVIDENCE_PACKET_JSON\\}\\}", "", text))) == 0L)

# How a marker row is read is one question, so it is one file. Both judges include it and
# the annotator holds it verbatim: two readings of the same number is how a label's own
# best-published marker came to count against it in one template and be argued away in
# another.
shared <- trimws(paste(
  readLines(file.path(PROMPT_DIR, "panel_reading.txt"), warn = FALSE),
  collapse = "\n"
))
stopifnot(startsWith(shared, "[R1]"), grepl(shared, text, fixed = TRUE))
for (name in c("subtype_judge.txt", "annotation_judge.txt")) {
  judge <- paste(readLines(file.path(PROMPT_DIR, name), warn = FALSE), collapse = "\n")
  stopifnot(grepl(PANEL_READING_PLACEHOLDER, judge, fixed = TRUE))
}

# One agent decides a cluster, so the configuration names exactly one model and one
# schema for it. Asserting the whole key set, rather than the presence of a few keys,
# is what keeps a second model from being added here without both arms being taught it.
A <- CFG$cluster_annotation
stopifnot(nzchar(as.character(A$annotator_model)))
stopifnot(nzchar(as.character(A$annotator_schema_version)))
stopifnot(setequal(names(A), c(
  "annotator_schema_version", "annotator_model", "annotator_reasoning_effort",
  "sources_per_marker", "source_genes_per_query", "source_sentence_max_chars",
  # A marker's sentences are served in batches rather than as a fixed shortlist, so both
  # arms need the batch cap and the seed that fixes which arbitrary order they come in.
  # Different seeds on the two arms would mean different evidence for the same cluster.
  "source_batches_per_marker", "source_order_seed",
  "pool_search_rows",
  "unmeasured_genes_shown", "warning_genes_shown", "negative_source_min_pct_in",
  "cluster_marker_rows",
  "subtype_exclusive_markers_required", "max_turns", "schema_retries",
  "max_rationale_chars", "unknown_token",
  # Quality control reads its model and its bounds from the same block, so an arm cannot
  # judge with one model while the other judges with another.
  "judge_model", "judge_reasoning_effort", "judge_max_rounds",
  # `judge_sources_pushed` caps the SENTENCES that ship, not the panel. It was
  # `judge_panel_genes` and it cut the panel, which put 10% of the genes the answers rested
  # on outside the judge's view entirely; the name changed with the meaning so a stale
  # config cannot quietly restore the cut.
  "judge_sources_pushed", "judge_sources_per_gene", "judge_query_turns"
)))
# Whatever the numbers are, a judge has to be able to ask, or code is the last word on what
# evidence exists -- and that is the arrangement that hid the rows.
stopifnot(as.integer(A$judge_query_turns) >= 1L)
cat("R production contract OK\n")
