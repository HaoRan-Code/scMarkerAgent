test_that("an .h5ad object is routed to the other arm, with the way out named", {
  h5 <- tempfile(fileext = ".h5ad")
  file.create(h5)
  expect_error(check_input(h5, "r"), "H5AD arm")
  expect_error(check_input(h5, "r"), "scmarkeragent::annotate")
})

test_that("an unsupported format is refused rather than guessed", {
  loom <- tempfile(fileext = ".loom")
  file.create(loom)
  expect_error(check_input(loom, "r"), "unsupported input format")
  expect_error(check_input(loom, "r"), "\\.rds")
})

test_that("a missing file is reported as missing", {
  expect_error(check_input(file.path(tempdir(), "absent.rds"), "r"), "does not exist")
})

test_that("no input at all says what to pass", {
  expect_error(check_input(""), "no input given")
})

test_that("a Seurat object passes and comes back normalized", {
  rds <- tempfile(fileext = ".rds")
  file.create(rds)
  expect_identical(check_input(rds, "r"), normalizePath(rds))
})
