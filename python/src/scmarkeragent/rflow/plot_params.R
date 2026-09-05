# =============================================================================
# plot_params.R -- validated loader for production figure parameters, .rds arm.
#
# Reads config/plot_params.json, the one parameter file every arm draws from rather than
# each carrying its own copy of the numbers.
#
# Contract:
#   * A key that is ABSENT falls back to the caller's default (backward compatible).
#   * A key that is PRESENT but malformed (wrong type / out of range / bad choice)
#     fails LOUDLY naming the offending "group.key", never silently.
#
# The file is JSONC: `//` line-comments and trailing commas are allowed. Relocate it
# with SCMA_PLOT_PARAMS.
# =============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
})

.pp_script_dir <- function() {
  frame_files <- vapply(
    sys.frames(),
    function(frame) {
      value <- frame$ofile
      if (is.null(value)) "" else as.character(value)[1]
    },
    ""
  )
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files)) {
    return(dirname(normalizePath(frame_files[length(frame_files)])))
  }
  args <- commandArgs(trailingOnly = FALSE)
  match <- grep("^--file=", args, value = TRUE)
  if (length(match)) {
    return(dirname(normalizePath(sub("^--file=", "", match[1]))))
  }
  normalizePath(getwd())
}

# Resolved while this file is being sourced. Asking later would be too late: `source()`
# has left the call stack by then and the lookup falls back to the working directory.
.PP_SELF <- .pp_script_dir()

.PP_CACHE <- new.env(parent = emptyenv())
.PP_CACHE$path <- NULL
.PP_CACHE$data <- NULL

#' Resolved params file path.
scma_plot_params_path <- function() {
  override <- Sys.getenv("SCMA_PLOT_PARAMS", "")
  if (nzchar(override)) {
    return(override)
  }
  file.path(dirname(.PP_SELF), "config", "plot_params.json")
}

#' Remove `//` line-comments (not inside strings) and trailing commas -> strict JSON.
#'
#' Deterministic character-at-a-time string/escape tracking, so a `//` inside a quoted
#' value is kept and only a real comment is cut. Every arm's stripper cuts the same
#' characters.
.pp_strip_jsonc <- function(text) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  stripped <- vapply(lines, function(line) {
    chars <- strsplit(line, "", fixed = TRUE)[[1]]
    n <- length(chars)
    in_str <- FALSE
    esc <- FALSE
    buf <- character(0)
    i <- 1L
    while (i <= n) {
      ch <- chars[i]
      if (in_str) {
        buf <- c(buf, ch)
        if (esc) {
          esc <- FALSE
        } else if (ch == "\\") {
          esc <- TRUE
        } else if (ch == "\"") {
          in_str <- FALSE
        }
      } else {
        if (ch == "\"") {
          in_str <- TRUE
          buf <- c(buf, ch)
        } else if (ch == "/" && i < n && chars[i + 1L] == "/") {
          break # rest of line is a comment
        } else {
          buf <- c(buf, ch)
        }
      }
      i <- i + 1L
    }
    paste(buf, collapse = "")
  }, "", USE.NAMES = FALSE)
  gsub(",(\\s*[}\\]])", "\\1", paste(stripped, collapse = "\n"), perl = TRUE)
}

#' Parse the params file once (cached). An absent file means "all defaults" and returns
#' an empty list; a present but unparseable file is an error, never a silent default.
scma_plot_params_load <- function() {
  p <- scma_plot_params_path()
  if (!is.null(.PP_CACHE$path) && identical(.PP_CACHE$path, p) &&
    !is.null(.PP_CACHE$data)) {
    return(.PP_CACHE$data)
  }
  data <- list()
  if (file.exists(p)) {
    raw <- paste(readLines(p, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    data <- tryCatch(
      fromJSON(.pp_strip_jsonc(raw), simplifyVector = FALSE),
      error = function(e) {
        stop(
          "plot_params: ", p, " is not valid JSON/JSONC: ", conditionMessage(e),
          call. = FALSE
        )
      }
    )
    # An object parses to a named list; an array or scalar does not, and must not be
    # accepted as if it were a group table.
    if (!is.list(data) || (length(data) && is.null(names(data)))) {
      stop("plot_params: ", p, " top level must be an object", call. = FALSE)
    }
  }
  .PP_CACHE$path <- p
  .PP_CACHE$data <- data
  data
}

.PP_MISSING <- structure(list(), class = "scma_pp_missing")

.pp_raw <- function(group, key) {
  data <- scma_plot_params_load()
  bucket <- data[[group]]
  if (is.null(bucket) || !is.list(bucket)) {
    return(.PP_MISSING)
  }
  value <- bucket[[key]]
  if (is.null(value)) .PP_MISSING else value
}

.pp_is_missing <- function(value) inherits(value, "scma_pp_missing")

# JSON true/false parses to logical, and R treats logicals as numbers under arithmetic.
# Both arms therefore reject a boolean where a number is required, rather than reading
# `true` as 1.
.pp_is_number <- function(value) {
  length(value) == 1L && is.numeric(value) && !is.logical(value)
}

.pp_show <- function(value) {
  paste(format(unlist(value)), collapse = ", ")
}

#' A numeric parameter, optionally range-checked.
scma_pp_num <- function(group, key, default, lo = NULL, hi = NULL) {
  value <- .pp_raw(group, key)
  if (.pp_is_missing(value)) {
    return(default)
  }
  if (!.pp_is_number(value)) {
    stop(
      "plot_params: '", group, ".", key, "' must be a number, got ",
      .pp_show(value),
      call. = FALSE
    )
  }
  value <- as.numeric(value)
  if ((!is.null(lo) && value < lo) || (!is.null(hi) && value > hi)) {
    stop(
      "plot_params: '", group, ".", key, "'=", value, " out of range [",
      if (is.null(lo)) "None" else lo, ", ", if (is.null(hi)) "None" else hi, "]",
      call. = FALSE
    )
  }
  value
}

#' A true/false parameter.
scma_pp_boolean <- function(group, key, default) {
  value <- .pp_raw(group, key)
  if (.pp_is_missing(value)) {
    return(default)
  }
  if (!(length(value) == 1L && is.logical(value))) {
    stop(
      "plot_params: '", group, ".", key, "' must be true/false, got ",
      .pp_show(value),
      call. = FALSE
    )
  }
  as.logical(value)
}

#' A string parameter, optionally restricted to a choice set.
scma_pp_text <- function(group, key, default, choices = NULL) {
  value <- .pp_raw(group, key)
  if (.pp_is_missing(value)) {
    return(default)
  }
  if (!(length(value) == 1L && is.character(value))) {
    stop(
      "plot_params: '", group, ".", key, "' must be a string, got ",
      .pp_show(value),
      call. = FALSE
    )
  }
  value <- as.character(value)
  if (!is.null(choices) && !(value %in% choices)) {
    stop(
      "plot_params: '", group, ".", key, "'=", value, " not in ",
      paste(sort(choices), collapse = ", "),
      call. = FALSE
    )
  }
  value
}

#' A list of numbers, optionally of fixed length and range-checked element-wise.
scma_pp_num_list <- function(group, key, default, n = NULL, lo = NULL, hi = NULL) {
  value <- .pp_raw(group, key)
  if (.pp_is_missing(value)) {
    return(default)
  }
  items <- if (is.list(value)) value else as.list(value)
  if (!length(items) || !all(vapply(items, .pp_is_number, TRUE))) {
    stop(
      "plot_params: '", group, ".", key, "' must be a list of numbers, got ",
      .pp_show(value),
      call. = FALSE
    )
  }
  numbers <- vapply(items, as.numeric, 0)
  if (!is.null(n) && length(numbers) != n) {
    stop(
      "plot_params: '", group, ".", key, "' must have ", n, " numbers, got ",
      length(numbers),
      call. = FALSE
    )
  }
  for (x in numbers) {
    if ((!is.null(lo) && x < lo) || (!is.null(hi) && x > hi)) {
      stop(
        "plot_params: '", group, ".", key, "' element ", x, " out of range [",
        if (is.null(lo)) "None" else lo, ", ", if (is.null(hi)) "None" else hi, "]",
        call. = FALSE
      )
    }
  }
  numbers
}

#' A list of strings.
scma_pp_str_list <- function(group, key, default) {
  value <- .pp_raw(group, key)
  if (.pp_is_missing(value)) {
    return(default)
  }
  items <- if (is.list(value)) value else as.list(value)
  is_string <- function(one) length(one) == 1L && is.character(one)
  if (!length(items) || !all(vapply(items, is_string, TRUE))) {
    stop(
      "plot_params: '", group, ".", key, "' must be a list of strings, got ",
      .pp_show(value),
      call. = FALSE
    )
  }
  vapply(items, as.character, "")
}
