## Fetch and verify the curated marker resource.
##
## The resource is hundreds of megabytes and lives in a versioned archive rather than in
## the package, so getting it is a step a user has to take. It is also the thing an
## annotation is measured against, which is why every file is checked against a shipped
## SHA-256 index before the directory is declared usable: a truncated download is
## otherwise indistinguishable from a smaller database, and the run that follows would
## quietly annotate against half a resource.

.index_path <- function() {
  system.file("extdata", "resource_index.json", package = "scmarkeragent")
}

#' The shipped resource index
#'
#' @param path Optional path to an index file; defaults to the packaged one.
#' @return The index as a list: bundle name, total size and per-file size and SHA-256.
#' @export
resource_index <- function(path = .index_path()) {
  if (!nzchar(path) || !file.exists(path)) {
    stop(
      "resource index not found: this package was built without ",
      "required_files/resource_index.json",
      call. = FALSE
    )
  }
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

.sha256 <- function(path) {
  ## digest() would be one more dependency for one call, and openssl is not guaranteed
  ## either. sha256sum is in base R's tools on every platform R supports.
  unname(tools::sha256sum(path))
}

.check_files <- function(dest, index, quiet = FALSE) {
  problems <- character(0)
  for (name in names(index$files)) {
    expected <- index$files[[name]]
    path <- file.path(dest, name)
    if (!file.exists(path)) {
      problems <- c(problems, sprintf("missing: %s", name))
      next
    }
    size <- file.info(path)$size
    if (!isTRUE(size == expected$bytes)) {
      problems <- c(problems, sprintf(
        "wrong size: %s (%.0f bytes, expected %.0f)", name, size, expected$bytes
      ))
      next
    }
    if (!quiet) message("  checking ", name, " ...")
    if (!identical(.sha256(path), expected$sha256)) {
      problems <- c(problems, sprintf("checksum mismatch: %s", name))
    }
  }
  problems
}

#' Verify an existing resource directory
#'
#' @param dest Directory holding the unpacked resource.
#' @return `TRUE` invisibly when every indexed file is present and intact; otherwise it
#'   stops with the list of problems.
#' @export
verify_resources <- function(dest) {
  dest <- path.expand(dest)
  index <- resource_index()
  problems <- .check_files(dest, index)
  if (length(problems)) {
    stop(sprintf(
      "resource at %s is NOT usable (%d problem(s)):\n  %s",
      dest, length(problems), paste(problems, collapse = "\n  ")
    ), call. = FALSE)
  }
  message(sprintf(
    "resource at %s verified (%s, %d files)", dest, index$bundle, length(index$files)
  ))
  invisible(TRUE)
}

.unpack <- function(archive, dest, index) {
  staging <- file.path(dirname(dest), sprintf(".scmarkeragent-staging-%s", Sys.getpid()))
  dir.create(staging, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  utils::untar(archive, exdir = staging)
  ## The bundle may be archived with or without a top-level directory.
  first <- names(index$files)[1]
  roots <- c(staging, list.dirs(staging, recursive = FALSE))
  source <- roots[vapply(roots, function(r) file.exists(file.path(r, first)), logical(1))]
  if (!length(source)) {
    stop(sprintf(
      "archive does not contain the expected bundle (looked for %s)", first
    ), call. = FALSE)
  }
  source <- source[1]
  for (name in names(index$files)) {
    target <- file.path(dest, name)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    file.rename(file.path(source, name), target)
  }
}

#' Download the curated marker resource
#'
#' @param dest Directory to unpack into; pass it to a run as `resource_dir`.
#' @param url Archive URL. Defaults to the one recorded in the resource index; when no
#'   URL is configured yet, the function explains how to install the bundle by hand.
#' @return The destination path invisibly.
#' @export
download_resources <- function(dest, url = NULL) {
  dest <- path.expand(dest)
  index <- resource_index()
  human <- index$total_bytes / 1e9

  if (!length(.check_files(dest, index, quiet = TRUE))) {
    message(sprintf("resource already complete at %s (%s)", dest, index$bundle))
    return(invisible(dest))
  }

  if (is.null(url) || !nzchar(url)) url <- index$archive_url %||% ""
  if (!nzchar(url)) {
    stop(sprintf(
      paste0(
        "no archive URL is configured yet.\n",
        "  Bundle %s is ~%.1f GB across %d files.\n",
        "  Pass url = <archive>, or download the bundle by hand and unpack it into\n",
        "  %s\n",
        "  See required_files/README.md for the archive DOI."
      ), index$bundle, human, length(index$files), dest
    ), call. = FALSE)
  }

  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  archive <- tempfile(fileext = ".tar.gz")
  on.exit(unlink(archive), add = TRUE)
  message("downloading ", url)
  utils::download.file(url, archive, mode = "wb", quiet = FALSE)
  message("unpacking ...")
  .unpack(archive, dest, index)

  message("verifying ...")
  problems <- .check_files(dest, index)
  if (length(problems)) {
    stop(sprintf(
      "resource verification failed; the download is not usable:\n  %s",
      paste(problems, collapse = "\n  ")
    ), call. = FALSE)
  }
  message(sprintf("resource ready at %s (%s, %.1f GB)", dest, index$bundle, human))
  message(sprintf("  use it with: resource_dir = \"%s\"", dest))
  invisible(dest)
}

`%||%` <- function(x, y) if (is.null(x) || !length(x) || !nzchar(as.character(x)[1])) y else x
