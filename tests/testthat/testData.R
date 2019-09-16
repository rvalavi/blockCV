library(blockCV)

context("Test the availablbe data in the package")

test_that("test that the data exist", {

  expect_true(file.exists(system.file("extdata", "awt.grd", package = "blockCV")))
  expect_true(file.exists(system.file("extdata", "PA.csv", package = "blockCV")))
  expect_true(file.exists(system.file("extdata", "PB.csv", package = "blockCV")))

})


test_that("test that the data is corrrect", {
  skip_on_cran()

  PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
  PB <- read.csv(system.file("extdata", "PB.csv", package = "blockCV"))
  awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))

  expect_true(exists("PA"))
  expect_equal(dim(PA), c(254, 3))
  expect_is(PA, "data.frame")
  expect_true(exists("PB"))
  expect_equal(dim(PB), c(10116, 3))
  expect_is(PB, "data.frame")
  expect_true(exists("awt"))
  expect_equal(raster::nlayers(awt), 10)
  expect_is(awt, "RasterBrick")

})
