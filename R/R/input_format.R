## Which execution arm reads a given input format.
##
## The two arms are not interchangeable back-ends for one reader: each reads its native
## object directly, so a Seurat `.rds` is annotated by this package and an `.h5ad` by the
## arm that owns that format, and neither format is converted to the other first.

ARM_BY_SUFFIX <- list(".h5ad" = "h5ad", ".rds" = "r")

#' Validate an input file before a run starts
#'
#' Checks that the file exists and that its format belongs to this arm, and explains the
#' way out when it does not. The check is at the entry point rather than inside the
#' reader because the two failures it prevents look nothing alike to a user: an `.h5ad`
#' handed to this arm dies inside `readRDS` with a message about an unknown format, and a
#' path typo dies even later, after preprocessing has already written a cache directory.
#' Both are the same mistake, and both are answerable here in one sentence.
#'
#' @param path Path to the input object.
#' @param expected_arm Which arm is asking; `"r"` for this package.
#' @return The normalized path, invisibly usable as the validated input.
#' @export
check_input <- function(path, expected_arm = "r") {
  path <- if (length(path)) as.character(path)[1] else ""
  if (is.na(path) || !nzchar(path)) {
    stop("no input given: pass the path to a Seurat .rds file", call. = FALSE)
  }
  suffix <- tolower(sub("^.*(\\.[^.]+)$", "\\1", basename(path)))
  supported <- paste(names(ARM_BY_SUFFIX), collapse = ", ")
  arm <- ARM_BY_SUFFIX[[suffix]]
  if (is.null(arm)) {
    stop(sprintf(
      paste0(
        "unsupported input format '%s': %s\n",
        "  scMarkerAgent reads %s only, and never converts between them.\n",
        "  Convert the object first, keeping the RAW COUNTS, because the pipeline",
        " starts from counts."
      ), suffix, path, supported
    ), call. = FALSE)
  }
  if (!identical(arm, expected_arm)) {
    stop(sprintf(
      paste0(
        "%s input belongs to the %s arm, not the %s arm: %s\n",
        "  This package reads Seurat .rds objects, through scmarkeragent::annotate().\n",
        "  The two arms read their own object directly; neither converts the",
        " other's format."
      ), suffix, toupper(arm), toupper(expected_arm), path
    ), call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(sprintf("input file does not exist: %s", path), call. = FALSE)
  }
  if (dir.exists(path)) {
    stop(sprintf("input path is not a file: %s", path), call. = FALSE)
  }
  normalizePath(path)
}
