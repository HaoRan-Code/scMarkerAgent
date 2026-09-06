# The validation surface of annotate(): everything here must stop before a pipeline
# process is spawned, so no test needs the marker resource, a model, or more than an
# empty temp file. The pipeline itself is covered by test-pipeline-contract.R.

seurat_like_file <- function() {
  rds <- tempfile(fileext = ".rds")
  file.create(rds)
  rds
}

test_that("species is required and closed to the three supported names", {
  rds <- seurat_like_file()
  expect_error(
    annotate(rds, tissue = "blood", counts_source = "RNA/counts"),
    "species is required"
  )
  expect_error(
    annotate(rds,
      species = "Zebrafish", tissue = "blood", counts_source = "RNA/counts"
    ),
    "'arg' should be one of"
  )
})

test_that("tissue is required, and blank does not count", {
  rds <- seurat_like_file()
  expect_error(
    annotate(rds, species = "Human", counts_source = "RNA/counts"),
    "tissue is required"
  )
  expect_error(
    annotate(rds, species = "Human", tissue = "  ", counts_source = "RNA/counts"),
    "tissue is required"
  )
})

test_that("counts_source is required and the error says why it is never guessed", {
  rds <- seurat_like_file()
  expect_error(
    annotate(rds, species = "Human", tissue = "blood"),
    "counts_source is required"
  )
  expect_error(
    annotate(rds, species = "Human", tissue = "blood", counts_source = ""),
    "raw counts"
  )
})

test_that("the input is routed through check_input before anything else", {
  h5 <- tempfile(fileext = ".h5ad")
  file.create(h5)
  expect_error(
    annotate(h5, species = "Human", tissue = "blood", counts_source = "X"),
    "H5AD arm"
  )
})

test_that("a resource_dir that does not exist stops the run before it starts", {
  rds <- seurat_like_file()
  expect_error(
    annotate(rds,
      species = "Human", tissue = "blood", counts_source = "RNA/counts",
      resource_dir = file.path(tempdir(), "no-such-resource-dir")
    )
  )
})

test_that("a resource_dir that fails verification is refused", {
  rds <- seurat_like_file()
  empty <- tempfile("scma_res_")
  dir.create(empty)
  expect_error(
    annotate(rds,
      species = "Human", tissue = "blood", counts_source = "RNA/counts",
      resource_dir = empty
    ),
    "NOT usable"
  )
})

test_that("the dependency preflight names what is missing, and presto gets the hint", {
  expect_true(.stop_if_pipeline_dependencies_missing(c(Seurat = TRUE, presto = TRUE)))
  expect_error(
    .stop_if_pipeline_dependencies_missing(c(Seurat = FALSE, presto = TRUE)),
    "Seurat"
  )
  hint <- tryCatch(
    .stop_if_pipeline_dependencies_missing(c(Seurat = TRUE, presto = FALSE)),
    error = conditionMessage
  )
  expect_match(hint, "haoran-code.r-universe.dev", fixed = TRUE)
  expect_match(hint, "install_git.*immunogenomics/presto.git.*a24772a135c7895a8183b007376050556c60a05b")
})

test_that("presto is not resolved through the GitHub API", {
  # A Remotes: field makes every remotes::install_* call consult api.github.com, whose
  # anonymous quota is 60 calls per hour per IP address; the install then fails behind a
  # shared address even when presto is already installed. presto is served by the
  # package's R-universe repository instead, declared as an additional repository.
  description <- read.dcf(system.file("DESCRIPTION", package = "scmarkeragent"))
  expect_false("Remotes" %in% colnames(description))
  expect_match(description[, "Additional_repositories"], "haoran-code.r-universe.dev", fixed = TRUE)
  imports <- trimws(strsplit(description[, "Imports"], ",")[[1]])
  expect_true(any(grepl("^presto", imports)))
  expect_true(any(grepl("^Seurat", imports)))
})
