#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(data.table)
})
.sd <- local({
  args <- commandArgs(FALSE)
  match <- grep("^--file=", args, value = TRUE)
  dirname(normalizePath(sub("^--file=", "", match[1])))
})
source(file.path(.sd, "config.R"))
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

## Before any stage, and before a results directory exists: the wrong file, or the right
## file handed to the wrong arm, is answerable in one sentence here rather than inside a
## reader three stages later.
INPUT_RDS <- .check_input(INPUT_RDS, "r")

tag <- INPUT_TAG
outdir <- file.path(RESULTS, tag)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
manifest_path <- file.path(outdir, "run_manifest.json")
manifest <- list(
  schema_version = "run-manifest-v1",
  tag = tag,
  arm = "R",
  status = "running",
  stages = list()
)
save_manifest <- function() {
  write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
}
# Every stage runs. The old checkpoint reused a stage whenever its output files merely
# existed, without checking the input, the configuration or the marker database, so a rerun
# could silently annotate against a stale cache. Reuse of the preprocessing artifacts is a
# benchmark-only decision made explicitly by the benchmark rerun driver, never inferred
# here.
run_stage <- function(script, args = character()) {
  status <- system2(
    "Rscript", c(file.path(.sd, script), args),
    stdout = "", stderr = ""
  )
  manifest$stages[[script]] <<- list(
    status = if (status == 0L) "completed" else "failed",
    exit_code = status
  )
  save_manifest()
  if (status != 0L) stop(script, " failed with status ", status)
}

resolved <- CFG
resolved$resolved_paths <- list(
  package_root = PACKAGE_ROOT,
  resource_dir = RESOURCE_DIR,
  cache_dir = CACHE,
  results_dir = RESULTS,
  markers = DB_CSV,
  sources = DB_SOURCES,
  uberon_ontology = OBO_UBERON,
  ortholog_dir = ORTHO_DIR,
  prompt_dir = PROMPT_DIR
)
resolved$runtime <- list(
  input = normalizePath(INPUT_RDS),
  counts_source = COUNTS_SOURCE,
  species = INPUT_SPECIES,
  tissue = INPUT_TISSUE,
  disease = INPUT_DISEASE,
  clustering_resolution = LEIDEN_RESOLUTION,
  cross_species = CROSS_SPECIES,
  compute_umap = isTRUE(CFG$features$compute_umap),
  cluster_annotation = CLUSTER_ANNOTATION_ENABLED,
  credential_source = .resolve_credentials()$source,
  provider_mode = .provider_mode(),
  offline = tolower(Sys.getenv("SCMA_OFFLINE", "0")) %in%
    c("1", "true", "yes", "on")
)
manifest$runtime <- resolved$runtime
write_json(resolved, file.path(outdir, CFG$output$resolved_config_file),
  auto_unbox = TRUE, pretty = TRUE, null = "null"
)

tryCatch(
  {
    run_stage("preprocessing.R", tag)
    run_stage("candidate_scoring.R", tag)
    run_stage("cluster_annotation.R", tag)
    if (!(tolower(Sys.getenv("SCMA_NO_REPORT", "0")) %in%
      c("1", "true", "yes", "on"))) {
      # The report layer delivers the same files and the same fields every arm
      # delivers. Runs in-process because it is this arm's own code, not a stage of
      # its own.
      outputs <- tryCatch(
        {
          source(file.path(.sd, "reporting.R"))
          scma_generate_report(tag, outdir, cache = CACHE)
        },
        error = function(error) {
          manifest$stages$report <<- list(status = "failed", exit_code = 1L)
          save_manifest()
          stop("R report failed: ", conditionMessage(error), call. = FALSE)
        }
      )
      manifest$stages$report <- list(status = "completed", exit_code = 0L)
      manifest$outputs <- lapply(outputs, as.character)
      save_manifest()
    }
    result <- readRDS(file.path(CACHE, sprintf("%s_annotations.rds", tag)))$results
    manifest$clusters <- length(result)
    manifest$annotation_source <- as.list(sort(unique(vapply(
      result, function(record) as.character(record$annotation_source), ""
    ))))
    manifest$llm_status <- as.list(sort(unique(vapply(
      result, function(record) as.character(record$llm_status), ""
    ))))
    manifest$status <- "completed"
    manifest$llm_raw_log <- file.path(CACHE, "llm_cold_calls.jsonl")
    save_manifest()
  },
  error = function(error) {
    manifest$status <- "failed"
    manifest$error <- conditionMessage(error)
    save_manifest()
    stop(error)
  }
)

cat(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE, null = "null"), "\n")
