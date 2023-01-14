library(blockCV)

context("Test the availablbe data in the package")

test_that("test that the data exist", {

  expect_true(file.exists(system.file("extdata/au", "bio_5.tif", package = "blockCV")))
  expect_true(file.exists(system.file("extdata", "species.csv", package = "blockCV")))

})


test_that("test that the data is corrrect", {
  skip_on_cran()

  pa <- read.csv(system.file("extdata", "species.csv", package = "blockCV"))
  r <- system.file("extdata/au/", package = "blockCV") |>
    list.files(full.names = TRUE) |>
    terra::rast()

  expect_true(exists("pa"))
  expect_true(nrow(pa) > 100)
  # expect_equal(dim(pa), c(254, 3))
  expect_is(pa, "data.frame")
  expect_true(exists("r"))
  expect_equal(terra::nlyr(r), 4)
  expect_is(r, "SpatRaster")

})
