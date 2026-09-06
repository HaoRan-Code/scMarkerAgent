## The public entry point of the R arm.
##
## The pipeline itself ships under inst/scmarkeragent/ -- the same rflow sources the
## benchmarks ran, assembled in by the release builder. This wrapper only validates
## what a user hands it, points the pipeline's environment at those inputs, runs the
## packaged run.R in a fresh Rscript process (the pipeline sequences its own stages
## and writes its own manifest), and returns that manifest parsed.

## Everything the packaged pipeline attaches or ::-calls unconditionally, answered
## before a run spawns the Rscript that would otherwise die mid-run with a subprocess
## traceback. These are declared Imports, so an install from a repository brings them
## in; the preflight still exists because the pipeline runs in a fresh Rscript whose
## library path may differ from the session's, and because a package can be removed
## after install. zip is listed because the delivered results archive needs it and the
## report stage stops without it; the guarded optionals (ggrastr, png) are not: the
## pipeline degrades without them and records what it skipped.
##
## presto is not on CRAN. It is served, at the commit the released pipeline was
## validated with, by the package's own CRAN-like repository on R-universe, which is
## what install.packages() resolves it from; the git line is the fallback when that
## repository cannot be reached. Neither uses the GitHub REST API (anonymous quota:
## 60 calls per hour per IP address), which is why presto is deliberately NOT declared
## in a Remotes: field -- remotes consults the API for that field on every install.
.presto_install_command <- paste0(
  "install.packages(\"presto\", repos = c(\"https://haoran-code.r-universe.dev\", ",
  "\"https://cloud.r-project.org\"))\n  or, without that repository: ",
  "remotes::install_git(\"https://github.com/immunogenomics/presto.git\", ",
  "ref = \"a24772a135c7895a8183b007376050556c60a05b\")"
)
.pipeline_dependencies <- function() {
  c(
    curl = requireNamespace("curl", quietly = TRUE),
    data.table = requireNamespace("data.table", quietly = TRUE),
    digest = requireNamespace("digest", quietly = TRUE),
    ggplot2 = requireNamespace("ggplot2", quietly = TRUE),
    Matrix = requireNamespace("Matrix", quietly = TRUE),
    scales = requireNamespace("scales", quietly = TRUE),
    Seurat = requireNamespace("Seurat", quietly = TRUE),
    presto = requireNamespace("presto", quietly = TRUE),
    zip = requireNamespace("zip", quietly = TRUE)
  )
}

.stop_if_pipeline_dependencies_missing <- function(available = .pipeline_dependencies()) {
  if (all(available)) {
    return(invisible(TRUE))
  }
  missing <- names(available)[!available]
  hint <- if ("presto" %in% missing) {
    paste0("\n  presto is not on CRAN; install it with:\n  ", .presto_install_command)
  } else {
    ""
  }
  stop(sprintf(
    "the pipeline needs packages that are not installed: %s%s",
    paste(missing, collapse = ", "), hint
  ), call. = FALSE)
}

#' Annotate a Seurat .rds dataset
#'
#' Runs the full scMarkerAgent pipeline on one Seurat object: QC, normalization,
#' Leiden clustering (always recomputed), genome-wide differential expression,
#' candidate retrieval from the curated marker resource, the annotating agent with
#' its quality-control judges, and the delivered result package (tables, audit
#' sidecar, figures, offline viewer and ZIP).
#'
#' @param input Path to a Seurat `.rds` file carrying raw counts.
#' @param tag Short name for the run; names the results directory.
#' @param species `"Human"`, `"Mouse"` or `"Rat"`.
#' @param tissue Free-text tissue, resolved against UBERON.
#' @param counts_source Where the raw counts live, as `"<assay>/<layer>"`
#'   (for example `"RNA/counts"`). Required and never guessed.
#' @param disease Disease context; several values may be given as a character
#'   vector. Defaults to `"Normal"`.
#' @param resource_dir The downloaded marker resource directory
#'   (see [download_resources()]). Verified before the run starts.
#' @param work_dir Cache and results root; defaults to `./.scmarkeragent`.
#' @param clustering_resolution Leiden resolution; `NULL` keeps the shared default.
#' @param offline Run the deterministic stages only, without contacting a model.
#' @param cluster_annotation Set `FALSE` to skip the annotating agent; the label is
#'   then the top of the retrieval order, and every field recording how the label
#'   was produced says so.
#' @return The parsed `run_manifest.json`, invisibly: run status, per-stage
#'   outcomes, and the paths of every delivered file.
#' @export
annotate <- function(input,
                     tag = "dataset",
                     species,
                     tissue,
                     counts_source,
                     disease = "Normal",
                     resource_dir = NULL,
                     work_dir = NULL,
                     clustering_resolution = NULL,
                     offline = FALSE,
                     cluster_annotation = TRUE) {
  input <- check_input(input, expected_arm = "r")
  if (missing(species)) {
    stop("species is required: \"Human\", \"Mouse\" or \"Rat\"", call. = FALSE)
  }
  species <- match.arg(as.character(species)[1], c("Human", "Mouse", "Rat"))
  if (missing(tissue) || !nzchar(trimws(as.character(tissue)[1]))) {
    stop("tissue is required; it is resolved against UBERON", call. = FALSE)
  }
  if (missing(counts_source) || !nzchar(trimws(as.character(counts_source)[1]))) {
    stop(
      paste0(
        "counts_source is required, as '<assay>/<layer>' (for example \"RNA/counts\").\n",
        "  The pipeline starts from raw counts, and silently annotating a normalized\n",
        "  matrix as if it were counts is the one failure that produces plausible\n",
        "  output from wrong input."
      ),
      call. = FALSE
    )
  }
  run_script <- system.file(
    "scmarkeragent", "rflow", "run.R",
    package = "scmarkeragent", mustWork = TRUE
  )
  if (!is.null(resource_dir)) {
    resource_dir <- normalizePath(path.expand(resource_dir), mustWork = TRUE)
    verify_resources(resource_dir)
  }
  .stop_if_pipeline_dependencies_missing()
  work_dir <- normalizePath(
    path.expand(work_dir %||% file.path(getwd(), ".scmarkeragent")),
    mustWork = FALSE
  )

  environment_variables <- c(
    sprintf("SCMA_INPUT_RDS=%s", input),
    sprintf("SCMA_TAG=%s", as.character(tag)[1]),
    sprintf("SCMA_SPECIES=%s", species),
    sprintf("SCMA_TISSUE=%s", as.character(tissue)[1]),
    sprintf("SCMA_DISEASE=%s", paste(as.character(disease), collapse = "|")),
    sprintf("SCMA_COUNTS_SOURCE=%s", as.character(counts_source)[1]),
    sprintf("SCMA_WORK_DIR=%s", work_dir),
    sprintf("SCMA_OFFLINE=%s", if (isTRUE(offline)) "1" else "0"),
    sprintf(
      "SCMA_CLUSTER_ANNOTATION=%s", if (isTRUE(cluster_annotation)) "1" else "0"
    )
  )
  if (!is.null(resource_dir)) {
    environment_variables <- c(
      environment_variables, sprintf("SCMA_RESOURCE_DIR=%s", resource_dir)
    )
  }
  if (!is.null(clustering_resolution)) {
    environment_variables <- c(
      environment_variables,
      sprintf("SCMA_CLUSTERING_RESOLUTION=%s", as.numeric(clustering_resolution))
    )
  }

  status <- system2(
    file.path(R.home("bin"), "Rscript"), shQuote(run_script),
    env = environment_variables
  )

  manifest_path <- file.path(
    work_dir, "results", as.character(tag)[1], "run_manifest.json"
  )
  manifest <- if (file.exists(manifest_path)) {
    jsonlite::fromJSON(manifest_path, simplifyVector = TRUE)
  } else {
    NULL
  }
  if (!identical(status, 0L)) {
    detail <- if (!is.null(manifest) && !is.null(manifest$error)) {
      paste0("\n  ", manifest$error)
    } else {
      ""
    }
    stop(
      sprintf("scMarkerAgent run failed (exit %s)%s", status, detail),
      call. = FALSE
    )
  }
  if (is.null(manifest)) {
    stop(sprintf("run finished but wrote no manifest: %s", manifest_path),
      call. = FALSE
    )
  }
  invisible(manifest)
}
