library(blockCV)

context("spatialBlock function")

expect_names <- c("folds",
                  "biomodTable",
                  "k",
                  "blocks",
                  "species",
                  "range",
                  "plots",
                  "records")

test_that("test spatiaBlock function with random assingment and raster file", {

  awt <- raster::brick(system.file("extdata", "awt.tif", package = "blockCV"))
  PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
  pa_data <- sp::SpatialPointsDataFrame(PA[,c("x", "y")], PA, proj4string=crs(awt))

  sb1 <- spatialBlock(speciesData = pa_data,
                      species = "Species",
                      rasterLayer = awt,
                      rows = 5,
                      cols = 8,
                      k = 5,
                      selection = 'random',
                      iteration = 25,
                      numLimit = NULL,
                      biomod2Format = TRUE,
                      xOffset = 0.3,
                      yOffset = 0,
                      showBlocks = FALSE,
                      progress = TRUE)

  expect_true(exists("sb1"))
  expect_is(sb1, "SpatialBlock")
  expect_equal(names(sb1), expect_names)
  expect_equal(length(sb1$folds), 5)
  expect_is(sb1$folds, "list")
  expect_is(sb1$biomodTable, "matrix")
  expect_equal(dim(sb1$biomodTable), c(nrow(pa_data), 5))
  expect_equal(sb1$k, 5)
  expect_is(sb1$blocks, "SpatialPolygonsDataFrame")
  expect_is(sb1$species, "character")
  expect_null(sb1$range)
  expect_is(sb1$plots, "ggplot")
  expect_equal(dim(sb1$records), c(5, 4))
  expect_true(
    !all(sb1$records == 0)
  )

})

test_that("test spatiaBlock function with systematic assingment and no raster file", {

  awt <- raster::brick(system.file("extdata", "awt.tif", package = "blockCV"))
  PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
  pa_data <- sp::SpatialPointsDataFrame(PA[,c("x", "y")], PA, proj4string=crs(awt))

  sb2 <- spatialBlock(speciesData = pa_data,
                      rasterLayer = awt,
                      theRange = 68000,
                      k = 5,
                      selection = 'systematic',
                      numLimit = 1,
                      biomod2Format = TRUE,
                      showBlocks = TRUE)

  expect_true(exists("sb2"))
  expect_is(sb2, "SpatialBlock")
  expect_equal(names(sb2), expect_names)
  expect_equal(length(sb2$folds), 5)
  expect_is(sb2$folds, "list")
  expect_is(sb2$biomodTable, "matrix")
  expect_equal(dim(sb2$biomodTable), c(nrow(pa_data), 5))
  expect_equal(sb2$k, 5)
  expect_is(sb2$blocks, "SpatialPolygonsDataFrame")
  expect_null(sb2$species)
  expect_is(sb2$range, "numeric")
  expect_is(sb2$plots, "ggplot")
  expect_equal(dim(sb2$records), c(5, 2))
  expect_true(
    !all(sb2$records == 0)
  )

  expect_equal(print.SpatialBlock(sb2), "SpatialBlock")
  expect_message(plot.SpatialBlock(sb2))
  expect_output(summary.SpatialBlock(sb2))

})
