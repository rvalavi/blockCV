library(blockCV)

context("envBlock function")

expect_names <- c("folds",
                  "biomodTable",
                  "k",
                  "species",
                  "records")

test_that("test that environmental blocking function with rasterBlock, standard and species column", {

  awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))
  PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
  pa_data <- sp::SpatialPointsDataFrame(PA[,c("x", "y")], PA, proj4string=crs(awt))

  # environmental clustering
  eb <- envBlock(rasterLayer = awt,
                 speciesData = pa_data,
                 species = "Species", # name of the column with species data
                 k = 5,
                 standardization = "standard",
                 rasterBlock = TRUE)

  expect_true(exists("eb"))
  expect_is(eb, "EnvironmentalBlock")
  expect_equal(names(eb), expect_names)
  expect_equal(length(eb$folds), 5)
  expect_is(eb$folds, "list")
  expect_equal(eb$k, 5)
  expect_is(eb$species, "character")
  expect_equal(dim(eb$records), c(5, 4))
  expect_true(
    !all(eb$records == 0)
  )

})

test_that("test that environmental blocking function with no rasterBlock, normalize and no species column", {

  awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))
  PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
  pa_data <- sp::SpatialPointsDataFrame(PA[,c("x", "y")], PA, proj4string=crs(awt))

  # environmental clustering
  eb2 <- envBlock(rasterLayer = awt,
                  speciesData = pa_data,
                  k = 5,
                  standardization = "normal",
                  rasterBlock = FALSE)

  expect_true(exists("eb2"))
  expect_is(eb2, "EnvironmentalBlock")
  expect_equal(names(eb2), expect_names)
  expect_equal(length(eb2$folds), 5)
  expect_is(eb2$folds, "list")
  expect_equal(eb2$k, 5)
  expect_null(eb2$species)
  expect_equal(dim(eb2$records), c(5, 2))
  expect_true(
    !all(eb2$records == 0)
  )

  expect_equal(print.EnvironmentalBlock(eb2), "EnvironmentalBlock")
  expect_output(summary.EnvironmentalBlock(eb2))

})
