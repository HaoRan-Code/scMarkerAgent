suppressPackageStartupMessages(library(jsonlite))

.env <- function(key, default = "") {
  value <- Sys.getenv(key)
  if (nzchar(value)) value else default
}
.env_flag <- function(key, default) {
  value <- tolower(trimws(Sys.getenv(key)))
  if (!nzchar(value)) default else value %in% c("1", "true", "yes", "on")
}
.openai_endpoint <- function(base_url) {
  base <- sub("/+$", "", trimws(base_url))
  if (grepl("/responses$", base)) base else paste0(base, "/responses")
}
.resolve_credentials <- function() {
  if (tolower(.env("SCMA_OFFLINE", "0")) %in%
    c("1", "true", "yes", "on")) {
    return(list(key = "", url = "", source = "offline"))
  }
  key <- Sys.getenv("OPENAI_API_KEY")
  base <- .env("OPENAI_BASE_URL", "https://api.openai.com/v1")
  list(
    key = key,
    url = .openai_endpoint(base),
    source = if (nzchar(key)) "environment" else "none"
  )
}
.provider_mode <- function() {
  configured <- sub(
    "/+$", "", trimws(Sys.getenv("OPENAI_BASE_URL"))
  )
  if (!nzchar(configured) ||
    configured == "https://api.openai.com/v1") {
    "official_openai"
  } else {
    "custom_openai_base"
  }
}
.script_dir <- function() {
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
  if (length(match) != 1L) {
    stop("unable to resolve packaged R source directory")
  }
  dirname(normalizePath(sub("^--file=", "", match[1])))
}
.normalize_disease_name <- function(x) {
  value <- tolower(trimws(as.character(x)))
  value[is.na(x)] <- ""
  value
}
.normalize_disease_query <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  value <- unique(.normalize_disease_name(x))
  value <- value[nzchar(value)]
  if (!length(value) || "all" %in% value) NULL else value
}
.disease_exact_match <- function(db_labels, query) {
  target <- .normalize_disease_query(query)
  if (is.null(target)) {
    rep(TRUE, length(db_labels))
  } else {
    .normalize_disease_name(db_labels) %in% target
  }
}

SELF <- .script_dir()
PACKAGE_ROOT <- dirname(SELF)
CONFIG_PATH <- normalizePath(
  .env("SCMA_CONFIG", file.path(PACKAGE_ROOT, "config", "defaults.json"))
)
CFG <- fromJSON(CONFIG_PATH, simplifyVector = TRUE)
STATIC_RESOURCE_DIR <- normalizePath(
  file.path(PACKAGE_ROOT, "resources", "static"),
  mustWork = FALSE
)
RESOURCE_DIR <- normalizePath(
  .env("SCMA_RESOURCE_DIR", STATIC_RESOURCE_DIR),
  mustWork = FALSE
)
.resource <- function(relative) {
  candidate <- file.path(RESOURCE_DIR, relative)
  if (file.exists(candidate)) {
    candidate
  } else {
    file.path(STATIC_RESOURCE_DIR, relative)
  }
}

DB_CSV <- .resource(CFG$resources$markers)
DB_SOURCES <- .resource(CFG$resources$sources)
# There is deliberately no Cell Ontology resource. A candidate IS its curated free-text
# cell-type string; the CL id travels to the output as a display column only.
OBO_UBERON <- .resource(CFG$resources$uberon_ontology)
ORTHO_DIR <- if (dir.exists(file.path(RESOURCE_DIR, "ortholog"))) {
  file.path(RESOURCE_DIR, "ortholog")
} else {
  file.path(STATIC_RESOURCE_DIR, "ortholog")
}

WORK_DIR <- normalizePath(
  .env("SCMA_WORK_DIR", file.path(getwd(), ".scmarkeragent")),
  mustWork = FALSE
)
CACHE <- normalizePath(.env("SCMA_CACHE", file.path(WORK_DIR, "cache")),
  mustWork = FALSE
)
OUT <- CACHE
RESULTS <- normalizePath(.env("SCMA_RESULTS", file.path(WORK_DIR, "results")),
  mustWork = FALSE
)
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)

## Which arm reads which format, stated on this arm too. The R pipeline has to refuse an
## .h5ad by itself: routing lives in the entry points, and this arm has more than one
## (the packaged R entry point, and a user sourcing run.R directly). A check that exists
## only in the caller is a check the pipeline does not have.
ARM_BY_SUFFIX <- list(".h5ad" = "h5ad", ".rds" = "r")

.check_input <- function(path, expected_arm = "r") {
  path <- if (length(path)) as.character(path)[1] else ""
  if (is.na(path) || !nzchar(path)) {
    stop("no input given: set SCMA_INPUT_RDS to a Seurat .rds file", call. = FALSE)
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
        "  This pipeline reads Seurat .rds objects only.\n",
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

INPUT_RDS <- .env("SCMA_INPUT_RDS", "")
INPUT_TAG <- .env("SCMA_TAG", "dataset")
INPUT_SPECIES <- .env("SCMA_SPECIES", "")
INPUT_TISSUE <- .env("SCMA_TISSUE", "")
INPUT_DISEASE <- strsplit(.env("SCMA_DISEASE", "Normal"), "\\|")[[1]]
# Where the raw counts live, as '<assay>/<layer>'. Required: the counts matrix is never
# guessed, and the declaration is echoed into the run configuration.
COUNTS_SOURCE <- .env("SCMA_COUNTS_SOURCE", "")

QC_MIN_GENES <- as.integer(CFG$preprocessing$qc_min_genes)
QC_MIN_CELLS <- as.integer(CFG$preprocessing$qc_min_cells)
QC_MAX_MT <- as.numeric(CFG$preprocessing$qc_max_mt_percent)
# Clustering is always recomputed; only its resolution is settable, and it comes from the
# shared configuration unless the caller overrides it.
LEIDEN_RESOLUTION <- as.numeric(.env(
  "SCMA_CLUSTERING_RESOLUTION", as.character(CFG$preprocessing$leiden_resolution)
))
EXCLUDE_IN_VITRO <- isTRUE(CFG$features$exclude_in_vitro)
# Read from the shared configuration by BOTH arms; neither hard-codes them.
CROSS_SPECIES <- isTRUE(CFG$features$cross_species_markers)
CORROBORATION_ONLY <- isTRUE(CFG$retrieval$corroboration_only)
# The configuration carries the default and the environment variable carries the caller's
# override, which is how `--no-cluster-annotation` reaches this arm. Both CLIs spell the
# flag the same way, so turning the annotating agent off means the same thing in either.
CLUSTER_ANNOTATION_ENABLED <- .env_flag(
  "SCMA_CLUSTER_ANNOTATION", isTRUE(CFG$features$cluster_annotation)
)
TISSUE_ROOT <- as.character(CFG$eligibility$tissue_root)
RU_MIN_ELIG_GENES <- as.integer(CFG$retrieval$min_eligible_genes)
TOP_CANDIDATES <- as.integer(CFG$retrieval$top_candidates)
# Candidate admission: how many of a candidate's own positive markers must be
# significantly up-regulated in a cluster, and the pool size below which that
# requirement is relaxed one step at a time. Both arms read them from here.
RU_MIN_HITS <- as.integer(CFG$retrieval$min_significant_hits)
RU_MIN_POOL_FLOOR <- as.integer(CFG$retrieval$min_pool_floor)
if (RU_MIN_HITS < 1L || RU_MIN_POOL_FLOOR < 1L || TOP_CANDIDATES < 1L) {
  stop("retrieval admission controls must be positive integers")
}
MIN_CORROBORATING_PUBLICATIONS <- as.integer(
  CFG$eligibility$min_corroborating_publications
)
CORROBORATING_TIERS <- as.character(
  CFG$eligibility$corroborating_tiers
)

# Every table writes this token wherever a cell does not apply, so a blank never has to
# be read as either "absent" or "not applicable".
NOT_AVAILABLE <- CFG$output$not_applicable_token
na_display <- function(value) {
  if (is.null(value) || length(value) == 0L) {
    return(NOT_AVAILABLE)
  }
  text <- trimws(as.character(value)[1])
  if (is.na(text) || !nzchar(text) || tolower(text) %in% c("nan", "none", "na")) {
    NOT_AVAILABLE
  } else {
    text
  }
}
# How a label left quality control. Defined here rather than beside the judges because
# the interchange writes the field too, and a second spelling of `not_checked` on that
# side would be indistinguishable from a check that ran and passed.
QC_PASSED <- "passed"
QC_REVISED <- "passed_after_revision"
QC_DEMOTED <- "demoted_to_parent"
QC_FAILED <- "failed"
QC_UNCHECKED <- "not_checked"

FIGURES_DIR <- "figures"
FIGDATA_DIR <- "figure_data"
TABLE_A_FILE <- CFG$output$table_a_file
TABLE_B_FILE <- CFG$output$table_b_file
METRICS_FILE <- "metrics.csv"
STATS_FILE <- "stats_reproducibility.txt"

PROMPT_DIR <- file.path(PACKAGE_ROOT, "prompts")
## The panel-reading rules every template shares, and the marker that pulls them in.
## Serving sentences in batches needs the batch cap and the seed that fixes which
## arbitrary order they arrive in. A different seed on the two arms would mean different
## evidence for the same cluster, so both read them from this one file.
SOURCES_PER_MARKER <- as.integer(CFG$cluster_annotation$sources_per_marker)
SOURCE_BATCHES_PER_MARKER <- as.integer(
  CFG$cluster_annotation$source_batches_per_marker
)
SOURCE_ORDER_SEED <- as.character(CFG$cluster_annotation$source_order_seed)

PANEL_READING <- "panel_reading"
PANEL_READING_PLACEHOLDER <- "{{PANEL_READING}}"
