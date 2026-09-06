test_that("the shipped index describes the bundle the packages verify against", {
  index <- resource_index()
  expect_identical(index$schema_version, "scmarkeragent-resource-index-v1")
  expect_true(index$total_bytes > 0)
  expect_true("scmarkeragent_curated_markers.csv" %in% names(index$files))
  for (entry in index$files) {
    expect_match(entry$sha256, "^[0-9a-f]{64}$")
  }
})

test_that("the shipped index carries the published archive URL", {
  # download_resources() falls back on this field; without it a fresh install cannot
  # fetch the resource unless the user is handed a URL by hand.
  index <- resource_index()
  expect_match(index$archive_url, "^https://")
})

test_that("the shipped index is the same file as the repository's download contract", {
  canonical <- testthat::test_path("..", "..", "..", "required_files", "resource_index.json")
  skip_if_not(file.exists(canonical), "not running inside the repository checkout")
  expect_identical(
    unname(tools::sha256sum(canonical)),
    unname(tools::sha256sum(.index_path()))
  )
})

test_that("an incomplete directory is refused rather than half-used", {
  dest <- tempfile("scma_res_")
  dir.create(dest)
  expect_error(verify_resources(dest), "NOT usable")
})

test_that("without a configured URL the downloader explains the manual route", {
  index <- resource_index()
  index$archive_url <- NULL
  local_mocked_bindings(resource_index = function(...) index)
  dest <- tempfile("scma_res_")
  expect_error(download_resources(dest), "no archive URL is configured")
})
