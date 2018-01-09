library(blockCV)

context("envBlock function")

expect_names <- c("folds",
                  "biomodTable",
                  "k",
                  "species",
                  "records")

test_that("test that environmental blocking function works properly", {

  awt <- raster::brick(system.file("extdata", "awt.tif", package = "blockCV"))
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
