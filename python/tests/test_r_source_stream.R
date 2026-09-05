#!/usr/bin/env Rscript
# =============================================================================
# test_r_source_stream.R -- the R arm serves source sentences as a stream, not a shortlist.
#
# The shared stream contract, plus the wiring a test of the server alone cannot see. The
# server itself was written long before anything called it: `scma_source_server()` existed
# and matched the other arm's server method for method, while `cluster_annotation.R` still
# answered every request with the same first k sentences and filed the second ask as a
# duplicate.
# Both arms parsed, both arms' unit tests passed, and the two agents read different
# evidence -- so what is asserted here is the CONNECTION as much as the behaviour.
#
#   Usage: Rscript test_r_source_stream.R
# =============================================================================
suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
})
args <- commandArgs(trailingOnly = FALSE)
match <- grep("^--file=", args, value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", match[1])), ".."))
package <- file.path(root, "src", "scmarkeragent")

targets <- file.path(package, "rflow", c("marker_sources.R", "annotator_pool.R"))

## Evaluate only what a pure function needs, as the packet-parity test does: everything
## else in these scripts reads the filesystem or calls a model.
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

config <- fromJSON(file.path(package, "config", "defaults.json"), simplifyVector = FALSE)
A <- config$cluster_annotation
env$SOURCES_PER_GENE <- as.integer(A$sources_per_marker)
env$SOURCE_GENES_PER_QUERY <- as.integer(A$source_genes_per_query)
env$SOURCE_MAX_CHARS <- as.integer(A$source_sentence_max_chars)
env$NOT_AVAILABLE <- as.character(config$output$not_applicable_token)

## ---- a pair with more sentences than one batch, over two publication counts --------
## The same shape the shared fixture builds, so the two arms are asserted against
## one fixture rather than each against its own.
rows <- data.table(
  cell_type = "enterocyte",
  gene_upper = rep(c("SI", "VIL1"), each = 7L),
  pmid = as.character(rep(1:7, 2L)),
  pmcid = paste0("PMC", rep(1:7, 2L)),
  source = paste0("s", rep(1:7, 2L)),
  n_pub_support = rep(c(9, 9, 9, 9, 9, 4, 4), 2L)
)
setorder(rows, cell_type, gene_upper, -n_pub_support, source)
setkey(rows, cell_type, gene_upper)
context <- list(rows = rows, seed = "seed")

sentences <- function(records) vapply(records, function(r) as.character(r$sentence), "")

## ---- the server: asking again reaches what the first answer left behind ------------
server <- env$scma_source_server(context, batch = 3L, max_batches = 0L)
first <- server$take("enterocyte", "SI")
second <- server$take("enterocyte", "SI")
stopifnot(
  length(intersect(sentences(first$sources), sentences(second$sources))) == 0L,
  first$remaining == 4L, isFALSE(first$exhausted),
  second$remaining == 1L
)

## The stream ends and says so, and every sentence is reached exactly once.
server <- env$scma_source_server(context, batch = 3L, max_batches = 0L)
seen <- character(0)
for (i in seq_len(4L)) {
  answer <- server$take("enterocyte", "SI")
  seen <- c(seen, sentences(answer$sources))
}
stopifnot(
  identical(sort(seen), sort(paste0("s", 1:7))),
  isTRUE(answer$exhausted), answer$remaining == 0L
)

## Higher publication counts are served first: nothing from the lower tier reaches the
## opening batch while five better-corroborated sentences are still unread.
server <- env$scma_source_server(context, batch = 3L, max_batches = 0L)
opening_batch <- sentences(server$take("enterocyte", "SI")$sources)
stopifnot(!("s6" %in% opening_batch), !("s7" %in% opening_batch))

## Sentences already on the page are not served again.
server <- env$scma_source_server(context, batch = 3L, max_batches = 0L)
shipped <- server$opening("enterocyte", "SI")
stopifnot(
  length(intersect(sentences(shipped), sentences(server$take("enterocyte", "SI")$sources)))
  == 0L
)

## A pair can only be drawn from so many times, and the cap is on the draws rather than on
## what exists -- the reply still says how much is behind it.
server <- env$scma_source_server(context, batch = 1L, max_batches = 2L)
invisible(server$take("enterocyte", "SI"))
invisible(server$take("enterocyte", "SI"))
stopped <- server$take("enterocyte", "SI")
stopifnot(
  length(stopped$sources) == 0L, isTRUE(stopped$limit_reached), stopped$remaining == 5L
)

## One pair running out does not affect another.
server <- env$scma_source_server(context, batch = 3L, max_batches = 1L)
invisible(server$take("enterocyte", "SI"))
stopifnot(length(server$take("enterocyte", "VIL1")$sources) == 3L)

## ---- the wiring: the tool serves the stream, and says what is left -----------------
pool <- list(clusters = list(`1` = list(candidates = list(
  list(cell_type = "enterocyte", exclusion_sources = list())
))))
served <- env$scma_source_server(context, batch = 3L, max_batches = 0L)
call <- list(candidate = "enterocyte", genes = list("si"))
one <- env$.tool_sources(pool, "1", call, served)
two <- env$.tool_sources(pool, "1", call, served)
stopifnot(
  length(intersect(sentences(one$sources$SI), sentences(two$sources$SI))) == 0L,
  ## The note appears only while something is still behind the pair, and it counts what is
  ## left rather than what was sent.
  identical(one$not_yet_shown$SI, "4 more"),
  identical(two$not_yet_shown$SI, "1 more"),
  identical(one$sources$SI[[1]]$pmcid, "PMC1")
)
last <- env$.tool_sources(pool, "1", call, served)
stopifnot(is.null(last$not_yet_shown))

## A plain fetch callable still works and carries no note: the pool is built once for every
## cluster, before any per-cluster server exists.
plain <- function(candidate, gene, k) {
  list(list(pmid = "1", pmcid = "PMC1", sentence = "s1"))
}
answer <- env$.tool_sources(pool, "1", call, plain)
stopifnot(
  identical(sentences(answer$sources$SI), "s1"),
  is.null(answer$not_yet_shown)
)

## ---- the wiring: what the packet already shows is not offered again ----------------
## Regression guard. `exclusion_sources` are resolved once when the pool is built, so a
## per-cluster server has not served them and would hand back the page the agent is
## reading. Registering them is what makes the first request go further.
shipped <- lapply(1:3, function(i) {
  list(pmid = as.character(i), pmcid = paste0("PMC", i), sentence = paste0("s", i))
})
with_packet <- list(clusters = list(`1` = list(candidates = list(
  list(cell_type = "enterocyte", exclusion_sources = list(SI = shipped))
))))
registered <- env$scma_source_server(context, batch = 3L, max_batches = 0L)
env$scma_register_packet_sources(registered, with_packet, "1")
after <- env$.tool_sources(with_packet, "1", call, registered)
stopifnot(
  length(intersect(sentences(after$sources$SI), paste0("s", 1:3))) == 0L,
  identical(after$not_yet_shown$SI, "1 more")
)

## Negative control: without the registration the same request replies with the page the
## packet already carries. If this ever stops holding, the assertion above has stopped
## testing anything.
unregistered <- env$scma_source_server(context, batch = 3L, max_batches = 0L)
stopifnot(
  length(intersect(
    sentences(env$.tool_sources(with_packet, "1", call, unregistered)$sources$SI),
    paste0("s", 1:3)
  )) > 0L
)

## ---- the wiring: neither loop may short-circuit `sources` as a repeat --------------
## Behavioural coverage would have to run a conversation, which needs a model. This reads
## the two arms instead and requires the exemption in both: it is a one-line deletion that
## turns "these three did not settle it" back into an unanswerable request, and no other
## test in this suite would notice.
loop_r <- paste(
  readLines(file.path(package, "rflow", "cluster_annotation.R"), warn = FALSE),
  collapse = "\n"
)
loop_py <- paste(
  readLines(file.path(package, "cluster_annotation.py"), warn = FALSE),
  collapse = "\n"
)
stopifnot(
  ## The annotator loop and the judge loop, on each arm.
  length(gregexpr('identical(as.character(parsed$tool)[1], "sources")', loop_r,
    fixed = TRUE)[[1]]) == 2L,
  length(gregexpr('parsed["tool"] == "sources"', loop_py, fixed = TRUE)[[1]]) == 2L
)

cat("R source stream contract OK\n")
