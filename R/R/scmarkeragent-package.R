#' @keywords internal
#' @details
#' The package is a thin, validated entry point over the two-arm scMarkerAgent
#' pipeline. The pipeline itself -- the same `rflow/` sources the benchmarks ran,
#' with the shared configuration, prompts, schemas and viewer assets -- ships
#' under `inst/scmarkeragent/` and runs in its own `Rscript` process, so this
#' package's namespace stays small: [annotate()] to run, [check_input()] to
#' route formats, and [download_resources()] / [verify_resources()] /
#' [resource_index()] to manage the curated marker resource.
"_PACKAGE"
