library(blockCV)

context("Spatial Autocorrelation Finder")

expect_names <- c("range",
                  "rangeTable",
                  "plots",
                  "sampleNumber",
                  "variograms")

test_that("test spatialAutoRange function with multi-layer raster in parallel", {

  awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))
  nl <- raster::nlayers(awt)

  range1 <- spatialAutoRange(rasterLayer = awt,
                             sampleNumber = 5000,
                             doParallel = TRUE,
                             nCores = NULL,
                             plotVariograms = FALSE,
                             showPlots = FALSE)

  expect_true(exists("range1"))
  expect_is(range1, "SpatialAutoRange")
  expect_equal(names(range1), expect_names)
  expect_equal(nrow(range1$rangeTable), nl)
  expect_equal(length(range1$variograms), nl)
  expect_equal(dim(range1$rangeTable), c(nl, 3))
  expect_is(range1$plots[[1]], "ggplot")
  expect_is(range1$variograms[[1]], "autofitVariogram")
  expect_is(range1$variograms, "list")
  expect_is(range1$rangeTable, "data.frame")
  expect_is(range1$sampleNumber, "numeric")
  expect_is(range1$range, "numeric")
  expect_true(range1$range >= 0)
  expect_true(!is.null(range1$variograms))
  expect_true(
    all(names(awt) %in% range1$rangeTable$layers)
  )

  expect_equal(print.SpatialAutoRange(range1), c("SpatialAutoRange", "list"))
  expect_silent(plot.SpatialAutoRange(range1))
  expect_output(summary.SpatialAutoRange(range1))

})


test_that("test spatialAutoRange function with multi-layer raster without parallel processing", {

  awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))
  nl <- raster::nlayers(awt)
  raster::crs(awt) <- NA

  range2 <- spatialAutoRange(rasterLayer = awt,
                             sampleNumber = 5000,
                             doParallel = FALSE,
                             plotVariograms = TRUE,
                             progress = TRUE,
                             showPlots = TRUE)

  expect_true(exists("range2"))
  expect_is(range2, "SpatialAutoRange")
  expect_equal(names(range2), expect_names)
  expect_equal(nrow(range2$rangeTable), nl)
  expect_equal(length(range2$variograms), nl)
  expect_equal(dim(range2$rangeTable), c(nl, 3))
  expect_is(range2$plots[[1]], "ggplot")
  expect_is(range2$variograms[[1]], "autofitVariogram")
  expect_is(range2$variograms, "list")
  expect_is(range2$rangeTable, "data.frame")
  expect_is(range2$sampleNumber, "numeric")
  expect_is(range2$range, "numeric")
  expect_true(range2$range >= 0)
  expect_true(!is.null(range2$variograms))
  expect_true(
    all(names(awt) %in% range2$rangeTable$layers)
  )

  expect_equal(print.SpatialAutoRange(range2), c("SpatialAutoRange", "list"))
  expect_silent(plot.SpatialAutoRange(range2))
  expect_output(summary.SpatialAutoRange(range2))

})

test_that("test spatialAutoRange function with single-layer raster", {

  awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))

  range3 <- spatialAutoRange(rasterLayer = awt[[1]],
                             sampleNumber = 5000,
                             plotVariograms = TRUE,
                             showPlots = TRUE)

  expect_true(exists("range3"))
  expect_is(range3, "SpatialAutoRange")
  expect_is(range3$plots[[1]], "ggplot")
  expect_is(range3$variograms, "autofitVariogram")
  expect_is(range3$sampleNumber, "numeric")
  expect_is(range3$range, "numeric")
  expect_true(range3$range >= 0)
  expect_true(!is.null(range3$variograms))

  expect_equal(print.SpatialAutoRange(range3), c("SpatialAutoRange", "list"))
  expect_silent(plot.SpatialAutoRange(range3))

})
