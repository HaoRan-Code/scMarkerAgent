# The packaged pipeline under inst/scmarkeragent/ is the product this package exists to
# run, so the software tree's R-arm contract tests run against the installed copy here.
# The four scripts under tests/testthat/pipeline/ are software/tests/test_r_*.R with the
# departures marked "Packaged adaptation" inline: the pipeline directory arrives in
# SCMA_PACKAGED_PIPELINE (set below from the installed package) instead of being resolved
# from a source-tree layout, and the source-stream test asserts only this arm's half of
# one cross-arm wiring check, because the other arm's sources do not ship in this
# package. Everything is offline: fake credentials, and a loopback httpuv stub for the
# transport test. Each script runs in its own Rscript in a temporary directory, because
# the pipeline configuration creates a .scmarkeragent/ workspace wherever it starts.

pipeline_dir <- system.file("scmarkeragent", package = "scmarkeragent")

run_pipeline_test <- function(script) {
  path <- normalizePath(testthat::test_path("pipeline", script), mustWork = TRUE)
  log <- tempfile(fileext = ".log")
  ## The script gets its own TMPDIR under ours: its scratch (and the logs of any
  ## helper process it spawns, like the transport test's fake provider) stays
  ## reachable for the failure message instead of dying with the child session.
  scratch <- tempfile(paste0("scratch-", sub("[.]R$", "", script), "-"))
  dir.create(scratch, recursive = TRUE, showWarnings = FALSE)
  old_dir <- setwd(scratch)
  old_env <- Sys.getenv(
    c("SCMA_PACKAGED_PIPELINE", "TMPDIR"),
    unset = NA, names = TRUE
  )
  on.exit({
    setwd(old_dir)
    for (name in names(old_env)) {
      if (is.na(old_env[[name]])) {
        Sys.unsetenv(name)
      } else {
        restored <- list(old_env[[name]])
        names(restored) <- name
        do.call(Sys.setenv, restored)
      }
    }
  }, add = TRUE)
  Sys.setenv(SCMA_PACKAGED_PIPELINE = pipeline_dir, TMPDIR = scratch)
  status <- system2(
    file.path(R.home("bin"), "Rscript"), shQuote(path),
    stdout = log, stderr = log
  )
  helper_logs <- list.files(
    scratch, pattern = "[.]log$", recursive = TRUE, full.names = TRUE
  )
  helper_text <- vapply(
    helper_logs,
    function(f) {
      sprintf(
        "--- %s ---\n%s", basename(f),
        paste(readLines(f, warn = FALSE), collapse = "\n")
      )
    },
    ""
  )
  list(
    status = status,
    log = paste(
      c(paste(readLines(log, warn = FALSE), collapse = "\n"), helper_text),
      collapse = "\n"
    )
  )
}

expect_pipeline_test_passes <- function(script) {
  result <- run_pipeline_test(script)
  expect(
    identical(result$status, 0L),
    sprintf("%s exited with status %s:\n%s", script, result$status, result$log)
  )
}

test_that("the packaged pipeline holds the R production contract", {
  expect_pipeline_test_passes("test_r_contract.R")
})

test_that("the packaged LLM client speaks the shared transport contract", {
  for (needed in c("curl", "digest", "httpuv", "later", "promises")) {
    skip_if_not_installed(needed)
  }
  expect_pipeline_test_passes("test_r_llm_client.R")
})

test_that("the packaged source stream serves sentences as a stream", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("digest")
  expect_pipeline_test_passes("test_r_source_stream.R")
})

test_that("the packaged evidence packet matches the two-arm parity fixture", {
  skip_if_not_installed("data.table")
  expect_pipeline_test_passes("test_two_arm_packet_parity.R")
})
