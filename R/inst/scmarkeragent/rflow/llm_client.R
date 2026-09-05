# =============================================================================
# llm_client.R -- transport, replay cache and cold-call audit for the .rds arm.
#
# Native transport: this arm reaches the provider itself.
#
# Three things have to agree between the arms, and do because each COMPUTES them from the
# same inputs rather than borrowing the other's implementation:
#
#   - the cache key, sha256("<model>\x1f<effort>\x1f<prompt>") over UTF-8, so a run
#     finished on either arm replays on the other without a call;
#   - the cold-call record, `llm-cold-call-v1`, appended to one JSONL file;
#   - the request body, which carries the model, the prompt, and the reasoning effort
#     when one is asked for -- no temperature, no seed, no sampling controls at all.
#
# Concurrency is the one place the arms genuinely differ. R has no threads, so here the
# BATCH is the primitive and a single call is the batch of one, rather than one call per
# thread. libcurl drives every request from this one process with at most
# `max_workers` sockets in flight, and the retry proceeds in rounds. Per prompt that is
# the same attempt count against the same backoff; only the scheduling differs, in that
# a round's backoff waits for the round.
#
# Sourced by cluster_annotation.R after config.R; defines no side effects of its own.
# =============================================================================

suppressPackageStartupMessages({
  library(curl)
  library(digest)
  library(jsonlite)
})

LLM_MODEL <- as.character(CFG$llm$model)
LLM_THREADS <- as.integer(CFG$llm$threads)
LLM_TIMEOUT <- as.integer(CFG$llm$timeout_seconds)
LLM_RETRIES <- as.integer(CFG$llm$retries)
## Which configuration produced a reply, recorded beside it. Hashed once: the file cannot
## change under a running stage, and both arms hash the same bytes.
CONFIG_SHA256 <- digest(file = CONFIG_PATH, algo = "sha256")

.LLM_CACHE <- new.env(parent = emptyenv())
.LLM_CACHE$path <- NULL
.LLM_CACHE$data <- NULL
.LLM_CACHE$dirty <- FALSE

#' The model one call will actually use.
#'
#' A caller passes the model its own stage configuration names; the environment override
#' wins over that, and the packaged default applies when neither is given. One variable
#' pins every stage to one model, which is what a reproduction wants; a per-stage choice
#' belongs in the configuration, so there is deliberately no per-stage variable.
scma_resolve_model <- function(model = NULL) {
  override <- trimws(Sys.getenv("SCMA_LLM_MODEL"))
  if (nzchar(override)) {
    return(override)
  }
  named <- if (is.null(model) || !length(model)) "" else as.character(model)[1]
  if (is.na(named) || !nzchar(named)) LLM_MODEL else named
}

#' The FIRST complete JSON object in a reply, or NULL.
#'
#' The first and not the widest span. A template asks for one object and a model
#' occasionally sends two -- a query, and then the answer it would have given without
#' waiting for the observation. Matching from the first brace to the last swallowed both
#' and parsed as neither, so a reply that plainly asked a question was scored as no reply
#' at all. Taking the first honours the question and drops the answer that skipped it,
#' which is the right way round: the observation was requested, so the verdict that
#' ignored it is not the one to keep. R has no `raw_decode`, so the first object is found
#' by walking the braces that are not inside a string.
scma_parse_json <- function(text) {
  if (is.null(text) || !length(text)) {
    return(NULL)
  }
  text <- trimws(gsub("^```(json)?|```$", "", trimws(as.character(text)[1])))
  whole <- tryCatch(fromJSON(text, simplifyVector = FALSE), error = function(e) NULL)
  if (!is.null(whole)) {
    return(whole)
  }
  start <- regexpr("{", text, fixed = TRUE)
  if (start < 0) {
    return(NULL)
  }
  chars <- strsplit(substring(text, start), "")[[1]]
  depth <- 0L
  in_string <- FALSE
  escaped <- FALSE
  for (position in seq_along(chars)) {
    ch <- chars[position]
    if (in_string) {
      if (escaped) {
        escaped <- FALSE
      } else if (identical(ch, "\\")) {
        escaped <- TRUE
      } else if (identical(ch, "\"")) {
        in_string <- FALSE
      }
      next
    }
    if (identical(ch, "\"")) {
      in_string <- TRUE
    } else if (identical(ch, "{")) {
      depth <- depth + 1L
    } else if (identical(ch, "}")) {
      depth <- depth - 1L
      if (depth == 0L) {
        return(tryCatch(
          fromJSON(paste(chars[seq_len(position)], collapse = ""), simplifyVector = FALSE),
          error = function(e) NULL
        ))
      }
    }
  }
  NULL
}

#' An error message with the credential and the endpoint taken out of it.
#'
#' Applied to everything that reaches the audit record or the console, because a
#' transport error quotes the request it failed on.
.llm_safe_error <- function(value, api_key = "", api_url = "") {
  text <- paste(as.character(value), collapse = " ")
  if (!length(text) || is.na(text)) {
    text <- ""
  }
  for (secret in c(api_key, api_url)) {
    if (nzchar(secret)) {
      text <- gsub(secret, "<redacted>", text, fixed = TRUE)
    }
  }
  text <- gsub(
    "(?i)authorization\\s*[:=]\\s*bearer\\s+\\S+",
    "Authorization: Bearer <redacted>", text,
    perl = TRUE
  )
  text <- gsub("https?://[^[:space:]\"']+", "<redacted_url>", text)
  substr(text, 1L, 2000L)
}

.llm_now_utc <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6+00:00", tz = "UTC")
}

## A transport control, or the configured one: an unset control arrives as NULL from one
## caller and as 0 from another, and both mean "use the configuration" rather than
## "no timeout" or "no attempts".
.llm_setting <- function(value, default) {
  number <- suppressWarnings(as.integer(value %||% NA_integer_)[1])
  if (is.na(number) || number == 0L) as.integer(default) else number
}

.llm_raw_log_path <- function() {
  value <- Sys.getenv("SCMA_LLM_RAW_LOG")
  if (nzchar(value)) value else file.path(CACHE, "llm_cold_calls.jsonl")
}

## One line per attempt, written as bytes: a prompt carries curated sentences in any
## script, and a connection that re-encoded to the session locale would corrupt the
## record of what was actually sent.
.llm_write_raw_record <- function(record) {
  path <- .llm_raw_log_path()
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  line <- enc2utf8(as.character(toJSON(
    record,
    auto_unbox = TRUE, null = "null", na = "null", digits = NA
  )))
  handle <- file(path, open = "ab")
  on.exit(close(handle))
  writeLines(line, handle, useBytes = TRUE)
  invisible(NULL)
}

## A provider response body, parsed strictly. Strictly, and not with the lenient reader
## the MODEL's replies get: a reply is prose that contains an object, but a response body
## either is one or is something to record verbatim, and picking the first brace out of an
## error page would file a fragment of it as the provider's answer.
.llm_json_body <- function(text) {
  value <- tryCatch(fromJSON(text, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(value)) list(raw_text = text) else value
}

#' The assistant text in a Responses payload, or NULL when it carries none.
.llm_response_text <- function(payload) {
  if (!is.list(payload) || is.null(names(payload))) {
    return(NULL)
  }
  direct <- payload$output_text
  if (is.character(direct) && length(direct) == 1L) {
    return(direct)
  }
  for (item in payload$output %||% list()) {
    if (!is.list(item)) next
    for (content in item$content %||% list()) {
      if (!is.list(content)) next
      type <- as.character(content$type %||% "")[1]
      if (!(type %in% c("output_text", "text"))) next
      text <- content$text
      if (is.character(text) && length(text) == 1L && nzchar(trimws(text))) {
        return(text)
      }
    }
  }
  NULL
}

.llm_request_body <- function(prompt, model, reasoning_effort) {
  body <- list(model = model, input = prompt)
  effort <- as.character(reasoning_effort %||% "")[1]
  if (!is.na(effort) && nzchar(effort)) {
    body$reasoning <- list(effort = effort)
  }
  body
}

.llm_item <- function(prompt, model, reasoning_effort, trace_id, turn_index) {
  text <- enc2utf8(as.character(prompt)[1])
  list(
    prompt = text,
    prompt_sha256 = digest(text, algo = "sha256", serialize = FALSE),
    request = .llm_request_body(text, model, reasoning_effort),
    reasoning_effort = if (is.null(reasoning_effort) || !length(reasoning_effort)) {
      NULL
    } else {
      as.character(reasoning_effort)[1]
    },
    trace_id = trace_id,
    turn_index = if (is.null(turn_index) || !length(turn_index)) {
      NULL
    } else {
      as.integer(turn_index)[1]
    }
  )
}

## One concurrent pass over the requests still pending. Every request in a pass is in
## flight at once, bounded by `workers` sockets; libcurl returns when the last of them has
## answered or failed. Multiplexing is off so that one socket is one request and the bound
## means what it says.
.llm_round <- function(items, api_url, api_key, timeout, workers) {
  outcomes <- vector("list", length(items))
  pool <- new_pool(total_con = workers, host_con = workers, multiplex = FALSE)
  for (position in seq_along(items)) {
    local({
      slot <- position
      started <- .llm_now_utc()
      handle <- new_handle()
      handle_setheaders(
        handle,
        "Authorization" = paste0("Bearer ", api_key),
        "Content-Type" = "application/json",
        # libcurl asks for `100 Continue` on any body past a kilobyte and every prompt is
        # far past it. An endpoint that never answers the probe costs a second per call
        # before the body is sent at all, and no arm asks for it.
        "Expect" = ""
      )
      body <- charToRaw(enc2utf8(as.character(toJSON(
        items[[slot]]$request,
        auto_unbox = TRUE, null = "null", digits = NA
      ))))
      handle_setopt(
        handle,
        url = api_url,
        post = TRUE, postfields = body, postfieldsize = length(body),
        timeout = timeout, connecttimeout = min(timeout, 60L)
      )
      multi_add(
        handle,
        done = function(response) {
          text <- rawToChar(response$content)
          Encoding(text) <- "UTF-8"
          outcomes[[slot]] <<- list(
            started = started, status = as.integer(response$status_code),
            text = text, failure = NULL
          )
        },
        fail = function(message) {
          outcomes[[slot]] <<- list(
            started = started, status = NULL, text = NULL,
            failure = as.character(message)[1]
          )
        },
        pool = pool
      )
    })
  }
  multi_run(pool = pool)
  outcomes
}

## The retry loop, shared by the single call and the batch: `retries` attempts per request
## against the bounded exponential backoff every arm uses. A request that succeeds leaves
## the loop and the rest carry on without it.
.llm_transport <- function(items, api_url, api_key, timeout, retries, workers) {
  contents <- vector("list", length(items))
  errors <- as.list(rep("no attempt", length(items)))
  pending <- seq_along(items)
  for (attempt in seq_len(max(0L, retries))) {
    if (!length(pending)) break
    outcomes <- .llm_round(
      items[pending], api_url, api_key, timeout, min(workers, length(pending))
    )
    still <- integer(0)
    for (position in seq_along(pending)) {
      index <- pending[position]
      item <- items[[index]]
      outcome <- outcomes[[position]]
      payload <- NULL
      content <- NULL
      failure <- NULL
      if (is.null(outcome)) {
        # libcurl returned neither a response nor a failure for this handle. Recorded as a
        # failed attempt rather than quietly retried, so the log still shows the round.
        failure <- "no response from transport"
      } else if (!is.null(outcome$failure)) {
        failure <- outcome$failure
      } else {
        payload <- .llm_json_body(outcome$text)
        if (identical(outcome$status, 200L)) {
          candidate <- .llm_response_text(payload)
          if (!is.null(candidate) && nzchar(trimws(candidate))) {
            content <- candidate
          } else {
            failure <- "empty content"
          }
        } else {
          failure <- sprintf("HTTP %d", outcome$status)
        }
      }
      failure <- if (is.null(failure)) NULL else .llm_safe_error(failure, api_key, api_url)
      .llm_write_raw_record(list(
        schema_version = "llm-cold-call-v1",
        timestamp_utc = outcome$started %||% .llm_now_utc(),
        prompt_sha256 = item$prompt_sha256,
        prompt = item$prompt,
        model = item$request$model,
        reasoning_effort = item$reasoning_effort,
        trace_id = item$trace_id,
        turn_index = item$turn_index,
        config_sha256 = CONFIG_SHA256,
        credential_source = .resolve_credentials()$source,
        provider_mode = .provider_mode(),
        request = item$request,
        attempt = attempt,
        http_status = outcome$status,
        provider_response = payload,
        content = content,
        parsed_json = if (is.null(content)) NULL else scma_parse_json(content),
        usage = if (is.list(payload)) payload$usage else NULL,
        cost = list(value = NULL, currency = NULL, reason = "pricing_not_configured"),
        error = failure
      ))
      if (is.null(failure)) {
        contents[[index]] <- content
        errors[index] <- list(NULL)
      } else {
        errors[[index]] <- failure
        still <- c(still, index)
      }
    }
    pending <- still
    if (length(pending) && attempt < retries) {
      Sys.sleep(min(2^(attempt - 1L), 30))
    }
  }
  lapply(seq_along(items), function(index) list(contents[[index]], errors[[index]]))
}

#' Call the configured endpoint once, with bounded exponential retry.
#'
#' Returns `list(content, error)`; `error` is NULL when the call succeeded.
#'
#' `trace_id` and `turn_index` are written to the audit record and to nothing else. One
#' annotation decision is several calls -- the agent asks for evidence and is called
#' again -- so without them the cold-call log is a pile of prompts with no way to tell
#' which belong to the same cluster's conversation or in what order. They stay out of the
#' request body and out of the cache key: they identify a run's trajectory, not the
#' question asked, and a replay must still hit on the prompt alone.
scma_call_llm <- function(prompt, api_url, api_key, timeout = NULL, retries = NULL,
                          reasoning_effort = NULL, model = NULL, trace_id = NULL,
                          turn_index = NULL) {
  item <- .llm_item(
    prompt, scma_resolve_model(model), reasoning_effort, trace_id, turn_index
  )
  .llm_transport(
    list(item), api_url, api_key,
    .llm_setting(timeout, LLM_TIMEOUT), .llm_setting(retries, LLM_RETRIES),
    workers = 1L
  )[[1]]
}

## ---- replay cache ------------------------------------------------------------
## Keyed by (model, reasoning effort, full prompt). A hit replays the exact raw response
## without an API call, which is what lets a finished run be re-reported, re-plotted and
## re-audited for nothing.

scma_response_cache_path <- function() {
  value <- Sys.getenv("SCMA_LLM_CACHE")
  if (nzchar(value)) value else file.path(CACHE, "llm_response_cache.json")
}

.llm_cache_key <- function(prompt, reasoning_effort, model = NULL) {
  effort <- as.character(reasoning_effort %||% "")[1]
  digest(
    enc2utf8(paste(
      scma_resolve_model(model),
      if (is.na(effort)) "" else effort,
      as.character(prompt)[1],
      sep = "\u001f"
    )),
    algo = "sha256", serialize = FALSE
  )
}

## An environment rather than a named list: a run adds thousands of entries one at a time,
## and growing a named list copies the whole thing on every insert.
.llm_load_cache <- function() {
  path <- scma_response_cache_path()
  if (identical(.LLM_CACHE$path, path) && !is.null(.LLM_CACHE$data)) {
    return(.LLM_CACHE$data)
  }
  data <- new.env(parent = emptyenv())
  if (file.exists(path)) {
    text <- rawToChar(readBin(path, "raw", file.info(path)$size))
    Encoding(text) <- "UTF-8"
    loaded <- tryCatch(
      fromJSON(text, simplifyVector = FALSE),
      error = function(e) {
        stop(sprintf("invalid LLM response cache: %s", path), call. = FALSE)
      }
    )
    # `names()` is NULL for a JSON array and character(0) for an empty object, which is
    # what separates a cache file from a file that merely parses.
    if (!is.list(loaded) || is.null(names(loaded))) {
      stop(sprintf("LLM response cache is not an object: %s", path), call. = FALSE)
    }
    for (key in names(loaded)) {
      assign(key, as.character(loaded[[key]])[1], envir = data)
    }
  }
  .LLM_CACHE$path <- path
  .LLM_CACHE$data <- data
  .LLM_CACHE$dirty <- FALSE
  data
}

## Written whole through a temporary file and renamed, so a run interrupted mid-write
## leaves the previous cache intact rather than a truncated one that would refuse to parse
## on the next start.
scma_flush_response_cache <- function() {
  if (!isTRUE(.LLM_CACHE$dirty) || is.null(.LLM_CACHE$data)) {
    return(invisible(NULL))
  }
  path <- .LLM_CACHE$path %||% scma_response_cache_path()
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  keys <- sort(ls(.LLM_CACHE$data, all.names = TRUE))
  text <- if (!length(keys)) {
    "{}"
  } else {
    as.character(toJSON(
      setNames(lapply(keys, function(key) get(key, envir = .LLM_CACHE$data)), keys),
      auto_unbox = TRUE
    ))
  }
  temporary <- paste0(path, ".tmp")
  handle <- file(temporary, open = "wb")
  writeBin(charToRaw(enc2utf8(text)), handle)
  close(handle)
  file.rename(temporary, path)
  .LLM_CACHE$dirty <- FALSE
  invisible(NULL)
}

## The cache lookup, and the calls the misses need, for a whole batch at once. Both public
## entry points are this function; a single call is the batch of one.
.llm_cached_batch <- function(prompts, api_url, api_key, timeout, retries,
                              reasoning_effort, model, max_workers,
                              trace_ids, turn_indexes) {
  traces <- if (is.null(trace_ids)) vector("list", length(prompts)) else as.list(trace_ids)
  turns <- if (is.null(turn_indexes)) {
    vector("list", length(prompts))
  } else {
    as.list(turn_indexes)
  }
  resolved <- scma_resolve_model(model)
  cache <- .llm_load_cache()
  answers <- vector("list", length(prompts))
  replayed <- rep(FALSE, length(prompts))
  keys <- vapply(
    prompts, function(prompt) .llm_cache_key(prompt, reasoning_effort, model), ""
  )
  pending <- integer(0)
  for (index in seq_along(prompts)) {
    if (exists(keys[index], envir = cache, inherits = FALSE)) {
      answers[[index]] <- list(get(keys[index], envir = cache, inherits = FALSE), NULL)
      replayed[index] <- TRUE
    } else if (!nzchar(api_url %||% "") || !nzchar(api_key %||% "")) {
      # A miss with no credentials is not a failure to retry: it is the whole answer in
      # offline mode, and the annotating loop reads exactly this string to say so.
      answers[[index]] <- list(NULL, "cache_miss_no_credentials")
    } else {
      pending <- c(pending, index)
    }
  }
  if (!length(pending)) {
    return(list(answers = answers, replayed = replayed))
  }
  items <- lapply(pending, function(index) {
    .llm_item(
      prompts[[index]], resolved, reasoning_effort, traces[[index]], turns[[index]]
    )
  })
  results <- .llm_transport(
    items, api_url, api_key,
    .llm_setting(timeout, LLM_TIMEOUT), .llm_setting(retries, LLM_RETRIES),
    workers = max(1L, .llm_setting(max_workers, LLM_THREADS))
  )
  persisted <- FALSE
  for (position in seq_along(pending)) {
    index <- pending[position]
    answers[[index]] <- results[[position]]
    content <- results[[position]][[1]]
    if (!is.null(content)) {
      assign(keys[index], content, envir = cache)
      persisted <- TRUE
    }
  }
  if (persisted) {
    .LLM_CACHE$dirty <- TRUE
    # Once per batch rather than once per reply. A batch is one round of the annotating
    # loop, so an interrupted run loses at most the round it was in, and a 38-cluster
    # round does not rewrite the whole cache file 38 times.
    scma_flush_response_cache()
  }
  list(answers = answers, replayed = replayed)
}

#' Run one prompt against the cache, then the provider.
#'
#' Returns `list(content, error, from_cache)`.
scma_cached_call_llm <- function(prompt, api_url, api_key, timeout = NULL, retries = NULL,
                                 reasoning_effort = NULL, model = NULL, trace_id = NULL,
                                 turn_index = NULL) {
  batch <- .llm_cached_batch(
    list(prompt), api_url, api_key, timeout, retries, reasoning_effort, model,
    max_workers = 1L, trace_ids = list(trace_id), turn_indexes = list(turn_index)
  )
  c(batch$answers[[1]], list(batch$replayed[1]))
}

#' Run every prompt against the cache, then the provider, concurrently.
#'
#' Returns `list(content, error)` pairs in the order of `prompts`. `from_cache` is dropped
#' because a batched caller only needs the two.
#'
#' `trace_ids` and `turn_indexes`, when given, are per-prompt and go to the audit record
#' only. This is how a round is labelled: the loop advances one conversation per cluster
#' in lock step, so the prompts in one batch belong to DIFFERENT conversations at the same
#' turn, and a single scalar could not identify them.
scma_cached_call_llm_many <- function(prompts, api_url, api_key, timeout = NULL,
                                      retries = NULL, reasoning_effort = NULL,
                                      model = NULL, max_workers = NULL,
                                      trace_ids = NULL, turn_indexes = NULL) {
  prompts <- as.list(prompts)
  if (!length(prompts)) {
    return(list())
  }
  .llm_cached_batch(
    prompts, api_url, api_key, timeout, retries, reasoning_effort, model,
    max_workers, trace_ids, turn_indexes
  )$answers
}
