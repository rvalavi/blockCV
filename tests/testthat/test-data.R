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
  expect_true(nrow(pa) > 300)
  expect_s3_class(pa, "data.frame")
  expect_true(exists("r"))
  expect_equal(terra::nlyr(r), 4)
  expect_s4_class(r, "SpatRaster")

})
