test_that("the shipped index describes the bundle the packages verify against", {
  index <- resource_index()
  expect_identical(index$schema_version, "scmarkeragent-resource-index-v1")
  expect_true(index$total_bytes > 0)
  expect_true("scmarkeragent_curated_markers.csv" %in% names(index$files))
  for (entry in index$files) {
    expect_match(entry$sha256, "^[0-9a-f]{64}$")
  }
})

test_that("an incomplete directory is refused rather than half-used", {
  dest <- tempfile("scma_res_")
  dir.create(dest)
  expect_error(verify_resources(dest), "NOT usable")
})

test_that("without a configured URL the downloader explains the manual route", {
  dest <- tempfile("scma_res_")
  expect_error(download_resources(dest), "no archive URL is configured")
})
