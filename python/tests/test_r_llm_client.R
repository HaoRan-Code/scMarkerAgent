#!/usr/bin/env Rscript
# =============================================================================
# test_r_llm_client.R -- the .rds arm talks to a provider by itself, and agrees with
# the other arm about what it said.
#
# This exercises the transport rather than describing it. A fake provider is started on a
# loopback port and the real client is pointed at it, so the lines that build a libcurl
# handle and read a response actually execute.
#
# Three assertions carry the parity claim, and none of them calls the other arm:
#
#   - the cache key is a LITERAL hex digest, the same literal the other arm's transport
#     test asserts, so a change to how either arm keys a prompt has to be made in both;
#   - a cache file written by the other arm replays here without a call;
#   - the request body carries the model, the prompt and the reasoning effort, and
#     nothing else -- asserted against what the server RECEIVED, not against what this
#     arm believes it sent.
#
#   Usage: Rscript test_r_llm_client.R
# =============================================================================
suppressPackageStartupMessages({
  library(jsonlite)
})
args <- commandArgs(trailingOnly = FALSE)
match <- grep("^--file=", args, value = TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", match[1])), ".."))
rflow <- file.path(root, "src", "scmarkeragent", "rflow")

FAKE_KEY <- "unit-test-openai-key-not-secret"
MODEL <- "unit-test-model"

work <- tempfile("scma-llm-")
dir.create(work, recursive = TRUE, showWarnings = FALSE)
Sys.unsetenv(c(
  "OPENAI_API_KEY", "OPENAI_BASE_URL", "SCMA_LLM_MODEL",
  "SCMA_LLM_CACHE", "SCMA_LLM_RAW_LOG"
))
## The credential is passed as an argument AND present in the environment, exactly as in a
## run: the audit record names where the key came from, and it reads that from here.
Sys.setenv(SCMA_OFFLINE = "0", SCMA_WORK_DIR = work, OPENAI_API_KEY = FAKE_KEY)
source(file.path(rflow, "config.R"))
## The arm's `%||%` lives with the pool, which the annotating stage sources before any of
## this runs. Sourced here for the same reason and from the same file, rather than
## redefined, because a private copy is how two spellings of "missing" start to differ.
source(file.path(rflow, "annotator_pool.R"))
source(file.path(rflow, "llm_client.R"))

## ---- a provider on loopback -----------------------------------------------------
## Answers in the Responses shape and echoes what it received, so the test can assert the
## wire format rather than this arm's own record of it. `slow:` resolves through a promise
## so the single-threaded server keeps accepting while one request is outstanding, which
## is what makes the concurrency assertion below mean anything.
server_script <- file.path(work, "fake_provider.R")
writeLines(c(
  'suppressPackageStartupMessages({',
  '  library(httpuv); library(jsonlite); library(later); library(promises)',
  '})',
  'port_file <- commandArgs(trailingOnly = TRUE)[1]',
  'seen <- new.env(parent = emptyenv())',
  'json <- function(status, body) list(',
  '  status = status, headers = list("Content-Type" = "application/json"), body = body',
  ')',
  'responses_shape <- function(text) json(200L, as.character(toJSON(list(',
  '  id = "response_test", object = "response",',
  '  output = list(list(type = "message", role = "assistant",',
  '    content = list(list(type = "output_text", text = text)))),',
  '  usage = list(input_tokens = 7L, output_tokens = 3L)',
  '), auto_unbox = TRUE)))',
  'app <- list(call = function(req) {',
  '  body <- fromJSON(rawToChar(req$rook.input$read()), simplifyVector = FALSE)',
  '  prompt <- as.character(body$input)[1]',
  '  seen_fields <- list(',
  '    request = body,',
  '    content_type = as.character(req$CONTENT_TYPE)[1],',
  '    expect = if (is.null(req$HTTP_EXPECT)) "" else as.character(req$HTTP_EXPECT)[1]',
  '  )',
  '  ## Reflected back only when asked for. A provider that echoed the credential into',
  '  ## every reply would put it in the audit log through the response rather than',
  '  ## through the client, and the assertion that the client never writes it would',
  '  ## then be untestable.',
  '  if (startsWith(prompt, "echo-auth")) {',
  '    seen_fields$authorization <- as.character(req$HTTP_AUTHORIZATION)[1]',
  '  }',
  '  echo <- as.character(toJSON(seen_fields, auto_unbox = TRUE))',
  '  if (startsWith(prompt, "always-fail")) return(json(500L, "{\\"error\\":\\"boom\\"}"))',
  '  if (startsWith(prompt, "fail-once")) {',
  '    n <- if (exists(prompt, envir = seen, inherits = FALSE)) {',
  '      get(prompt, envir = seen, inherits = FALSE)',
  '    } else 0L',
  '    assign(prompt, n + 1L, envir = seen)',
  '    if (n == 0L) return(json(500L, "{\\"error\\":\\"boom\\"}"))',
  '  }',
  '  if (startsWith(prompt, "empty")) return(responses_shape("   "))',
  '  if (startsWith(prompt, "slow")) {',
  '    return(promise(function(resolve, reject) {',
  '      later(function() resolve(responses_shape(echo)), delay = 0.4)',
  '    }))',
  '  }',
  '  responses_shape(echo)',
  '})',
  'port <- NA_integer_',
  'for (attempt in seq_len(60)) {',
  '  candidate <- sample(20000:60000, 1)',
  '  handle <- tryCatch(startServer("127.0.0.1", candidate, app), error = function(e) NULL)',
  '  if (!is.null(handle)) { port <- candidate; break }',
  '}',
  'if (is.na(port)) stop("no free loopback port")',
  'writeLines(as.character(port), port_file)',
  '## A hard deadline rather than a loop forever: a failed assertion in the parent leaves',
  '## this process orphaned, and an orphan holding a port is a test that fails once and',
  '## then fails differently every time after.',
  'deadline <- Sys.time() + 300',
  'while (Sys.time() < deadline) service(200)'
), server_script)

port_file <- file.path(work, "port")
## Its own streams, so the parent's stdout is not held open by a server that outlives a
## failed assertion.
system2(
  "Rscript", c(server_script, port_file),
  wait = FALSE,
  stdout = file.path(work, "fake_provider.log"),
  stderr = file.path(work, "fake_provider.log")
)
for (attempt in seq_len(200)) {
  if (file.exists(port_file) && length(readLines(port_file, warn = FALSE))) break
  Sys.sleep(0.1)
}
if (!file.exists(port_file)) stop("fake provider did not start")
PORT <- as.integer(readLines(port_file, warn = FALSE)[1])
URL <- sprintf("http://127.0.0.1:%d/v1/responses", PORT)
on.exit(system2("pkill", c("-f", server_script), stdout = NULL, stderr = NULL), add = TRUE)

## Each scenario gets its own cache and its own cold-call log, so "how many calls did this
## cost" is countable rather than inferred.
scenario <- function(name) {
  Sys.setenv(
    SCMA_LLM_CACHE = file.path(work, paste0(name, "-cache.json")),
    SCMA_LLM_RAW_LOG = file.path(work, paste0(name, "-cold.jsonl"))
  )
}
cold_calls <- function(name) {
  path <- file.path(work, paste0(name, "-cold.jsonl"))
  if (!file.exists(path)) {
    return(list())
  }
  lapply(readLines(path, warn = FALSE), fromJSON, simplifyVector = FALSE)
}

## ---- the cache key is a shared literal --------------------------------------------
## The same hex string the other arm's transport test asserts. Neither arm computes the
## other's; they are pinned to one value, which is what lets a run finished on one replay
## on the other.
stopifnot(identical(
  .llm_cache_key("ping", "high", MODEL),
  "d738f945517ec7cb7c8ffe1e671d929b9637e3d59b22e08e8514114bb26c896a"
))
stopifnot(identical(
  .llm_cache_key("\u03b2 cell \u2014 \u4e2d\u6587", "high", MODEL),
  "fb18c245d3bf92d71d00e774a3de264bf3fd26092adacee5166d46aa5435375f"
))
## No effort is the empty string, not the word NULL and not a dropped field.
stopifnot(identical(
  .llm_cache_key("ping", NULL, MODEL),
  "16991eeca85954e48c4c92431a3ca4c082f0f6fe735fddc425a311ba479859a2"
))
## The environment override wins over what a stage names, and the key follows the model
## that will actually answer. So under an override every stage keys on the pinned model --
## a stage's own name cannot reach the key -- and the key differs from the unpinned one,
## which is what stops a pinned run replaying a reply another model produced.
Sys.setenv(SCMA_LLM_MODEL = "pinned-model")
stopifnot(identical(scma_resolve_model(MODEL), "pinned-model"))
stopifnot(identical(.llm_cache_key("ping", "high", MODEL), .llm_cache_key("ping", "high")))
stopifnot(!identical(
  .llm_cache_key("ping", "high", MODEL),
  "d738f945517ec7cb7c8ffe1e671d929b9637e3d59b22e08e8514114bb26c896a"
))
Sys.unsetenv("SCMA_LLM_MODEL")
stopifnot(identical(scma_resolve_model(MODEL), MODEL))
stopifnot(identical(scma_resolve_model(NULL), as.character(CFG$llm$model)))

## ---- an unset transport control reads the configuration ----------------------------
## It arrives as NULL from one caller and as 0 from another, and both mean the same thing
## on the other arm. Zero retries read literally is a call that never happens at all.
stopifnot(identical(.llm_setting(NULL, 7L), 7L))
stopifnot(identical(.llm_setting(0L, 7L), 7L))
stopifnot(identical(.llm_setting(3L, 7L), 3L))

## ---- nothing that identifies the caller reaches an error message ------------------
leaky <- sprintf(
  "POST %s failed: Authorization: Bearer %s rejected", URL, FAKE_KEY
)
safe <- .llm_safe_error(leaky, FAKE_KEY, URL)
stopifnot(!grepl(FAKE_KEY, safe, fixed = TRUE))
stopifnot(!grepl(URL, safe, fixed = TRUE))
stopifnot(grepl("<redacted", safe, fixed = TRUE))
stopifnot(nchar(.llm_safe_error(strrep("x", 5000))) == 2000L)

## ---- what actually goes over the wire ------------------------------------------------
## Asserted against what the SERVER received rather than against this arm's own record of
## what it meant to send.
scenario("wire")
wire <- scma_cached_call_llm(
  "echo-auth:1", URL, FAKE_KEY,
  retries = 1L, reasoning_effort = "high", model = MODEL
)
stopifnot(is.null(wire[[2]]), isFALSE(wire[[3]]))
seen <- fromJSON(wire[[1]], simplifyVector = FALSE)
## Three fields and no sampling control of any kind; a temperature or a seed here would
## make a replay a different question from the one that was asked.
stopifnot(setequal(names(seen$request), c("model", "input", "reasoning")))
stopifnot(identical(seen$request$model, MODEL))
stopifnot(identical(seen$request$input, "echo-auth:1"))
stopifnot(identical(seen$request$reasoning$effort, "high"))
stopifnot(identical(seen$authorization, paste("Bearer", FAKE_KEY)))
stopifnot(startsWith(seen$content_type, "application/json"))
## libcurl asks for `100 Continue` on any body past a kilobyte unless it is told not to,
## and an endpoint that ignores the probe costs a second per call before the body is sent.
stopifnot(identical(seen$expect, ""))

## No effort named means no `reasoning` key at all, rather than an empty one.
plain <- scma_cached_call_llm("echo-auth:2", URL, FAKE_KEY, retries = 1L, model = MODEL)
stopifnot(setequal(
  names(fromJSON(plain[[1]], simplifyVector = FALSE)$request), c("model", "input")
))

## ---- one real call, and the audit record it leaves -------------------------------------
scenario("one")
answer <- scma_cached_call_llm(
  "ping", URL, FAKE_KEY,
  retries = 1L, reasoning_effort = "high", model = MODEL,
  trace_id = "trace-1", turn_index = 3L
)
stopifnot(is.null(answer[[2]]), isFALSE(answer[[3]]))
records <- cold_calls("one")
stopifnot(length(records) == 1L)
record <- records[[1]]
stopifnot(identical(record$schema_version, "llm-cold-call-v1"))
stopifnot(identical(record$prompt, "ping"))
stopifnot(identical(
  record$prompt_sha256,
  "758d61f26a44448384e5c4468a0dcb7a2abe456067b0f7b505bc28b9411fe931"
))
stopifnot(identical(record$model, MODEL))
stopifnot(identical(record$reasoning_effort, "high"))
## The trajectory fields identify which conversation a prompt belongs to and in what
## order. They are in the record and in neither the request nor the key.
stopifnot(identical(record$trace_id, "trace-1"), identical(record$turn_index, 3L))
stopifnot(identical(record$attempt, 1L), identical(record$http_status, 200L))
stopifnot(identical(record$config_sha256, CONFIG_SHA256))
stopifnot(identical(record$credential_source, "environment"))
stopifnot(identical(record$provider_response$usage$input_tokens, 7L))
stopifnot(identical(record$cost$reason, "pricing_not_configured"))
stopifnot(is.null(record$error))
stopifnot(!is.null(record$parsed_json))
## The credential is not a field of the record and must not arrive inside one either.
stopifnot(!grepl(FAKE_KEY, paste(readLines(
  file.path(work, "one-cold.jsonl"), warn = FALSE
), collapse = "\n"), fixed = TRUE))

## ---- a repeat costs nothing ---------------------------------------------------------
replay <- scma_cached_call_llm(
  "ping", URL, FAKE_KEY,
  retries = 1L, reasoning_effort = "high", model = MODEL
)
stopifnot(isTRUE(replay[[3]]), identical(replay[[1]], answer[[1]]))
stopifnot(length(cold_calls("one")) == 1L)

## The cache on disk is an object keyed by that same digest, which is the whole of what
## the two arms have to agree on.
written <- fromJSON(file.path(work, "one-cache.json"), simplifyVector = FALSE)
stopifnot(identical(
  names(written), .llm_cache_key("ping", "high", MODEL)
))

## ---- a cache the other arm wrote ----------------------------------------------------
## Written here as bytes rather than by calling the other arm, and in that arm's own
## layout: spaces after the separators, non-ASCII left as UTF-8 rather than escaped.
scenario("interop")
interop_cache <- file.path(work, "interop-cache.json")
handle <- file(interop_cache, open = "wb")
writeBin(charToRaw(enc2utf8(paste0(
  '{"', .llm_cache_key("\u03b2 cell \u2014 \u4e2d\u6587", "high", MODEL), '": ',
  '"{\\"selected\\": \\"\u03b2 cell\\"}"}'
))), handle)
close(handle)
.LLM_CACHE$path <- NULL
.LLM_CACHE$data <- NULL
offline <- scma_cached_call_llm(
  "\u03b2 cell \u2014 \u4e2d\u6587", "", "",
  reasoning_effort = "high", model = MODEL
)
stopifnot(isTRUE(offline[[3]]), is.null(offline[[2]]))
stopifnot(identical(
  fromJSON(offline[[1]], simplifyVector = FALSE)$selected, "\u03b2 cell"
))
stopifnot(length(cold_calls("interop")) == 0L)

## A miss with no credentials is the answer, not an error to retry: the annotating loop
## reads this exact string to decide a cluster was never offered to a model.
missing <- scma_cached_call_llm("never seen", "", "", model = MODEL)
stopifnot(is.null(missing[[1]]), identical(missing[[2]], "cache_miss_no_credentials"))
stopifnot(length(cold_calls("interop")) == 0L)

## A cache file that parses but is not an object is refused rather than half-read.
scenario("broken")
writeLines('["not", "an", "object"]', file.path(work, "broken-cache.json"))
.LLM_CACHE$path <- NULL
.LLM_CACHE$data <- NULL
stopifnot(inherits(tryCatch(.llm_load_cache(), error = function(e) e), "error"))

## ---- retry ---------------------------------------------------------------------------
## Two attempts for one prompt, the first refused. Both are in the log, because the record
## of what a run cost is the record of every attempt it made.
scenario("retry")
recovered <- scma_cached_call_llm(
  "fail-once:a", URL, FAKE_KEY,
  retries = 2L, reasoning_effort = "high", model = MODEL
)
stopifnot(is.null(recovered[[2]]), !is.null(recovered[[1]]))
attempts <- cold_calls("retry")
stopifnot(length(attempts) == 2L)
stopifnot(identical(attempts[[1]]$attempt, 1L), identical(attempts[[1]]$error, "HTTP 500"))
stopifnot(identical(attempts[[2]]$attempt, 2L), is.null(attempts[[2]]$error))

## A prompt that never succeeds spends its whole budget and returns the last error. It is
## not cached: a failure must not replay as an answer.
scenario("exhausted")
refused <- scma_cached_call_llm(
  "always-fail", URL, FAKE_KEY,
  retries = 2L, reasoning_effort = "high", model = MODEL
)
stopifnot(is.null(refused[[1]]), identical(refused[[2]], "HTTP 500"))
stopifnot(length(cold_calls("exhausted")) == 2L)
stopifnot(!file.exists(file.path(work, "exhausted-cache.json")))

## A 200 carrying no text is a failure too, and says which kind.
scenario("empty")
blank <- scma_cached_call_llm(
  "empty:one", URL, FAKE_KEY,
  retries = 1L, reasoning_effort = "high", model = MODEL
)
stopifnot(is.null(blank[[1]]), identical(blank[[2]], "empty content"))

## ---- a batch ---------------------------------------------------------------------------
## The answers come back in the order the prompts went out, each with its own trajectory
## label. One pool serves them all, so a mixed-up slot would be invisible in any single
## call and wrong in every one.
scenario("batch")
prompts <- as.list(sprintf("batch-%02d", seq_len(12)))
batch <- scma_cached_call_llm_many(
  prompts, URL, FAKE_KEY,
  retries = 1L, reasoning_effort = "high", model = MODEL, max_workers = 12L,
  trace_ids = as.list(sprintf("t%02d", seq_along(prompts))),
  turn_indexes = as.list(seq_along(prompts))
)
stopifnot(length(batch) == length(prompts))
for (index in seq_along(prompts)) {
  stopifnot(is.null(batch[[index]][[2]]))
  stopifnot(identical(
    fromJSON(batch[[index]][[1]], simplifyVector = FALSE)$request$input, prompts[[index]]
  ))
}
labels <- cold_calls("batch")
stopifnot(length(labels) == length(prompts))
stopifnot(setequal(
  vapply(labels, function(r) paste(r$trace_id, r$turn_index), ""),
  sprintf("t%02d %d", seq_along(prompts), seq_along(prompts))
))

## Every prompt in the batch is now cached, so asking again costs no call at all.
again <- scma_cached_call_llm_many(
  prompts, URL, FAKE_KEY,
  retries = 1L, reasoning_effort = "high", model = MODEL
)
stopifnot(identical(again, batch))
stopifnot(length(cold_calls("batch")) == length(prompts))

## ---- the batch is actually concurrent ---------------------------------------------------
## Eight requests the provider holds for 0.4s each. In flight together they finish in about
## half a second; one at a time they cannot finish in under 3.2, so the bound separates
## them by a wide margin without pinning the test to a machine's speed.
scenario("parallel")
started <- Sys.time()
held <- scma_cached_call_llm_many(
  as.list(sprintf("slow:%d", seq_len(8))), URL, FAKE_KEY,
  retries = 1L, reasoning_effort = "high", model = MODEL, max_workers = 8L
)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
stopifnot(length(held) == 8L)
stopifnot(all(vapply(held, function(pair) is.null(pair[[2]]), TRUE)))
if (elapsed >= 2) {
  stop(sprintf("eight concurrent 0.4s calls took %.1fs; the pool is serialising", elapsed))
}

cat("R LLM client contract OK\n")
