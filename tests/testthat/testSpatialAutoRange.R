library(blockCV)

context("Spatial Autocorrelation Finder")

expect_names <- c("range",
                  "rangeTable",
                  "plots",
                  "sampleNumber",
                  "variograms")

awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))
nl <- raster::nlayers(awt)

test_that("test spatialAutoRange function with multi-layer raster in parallel", {

  range1 <- spatialAutoRange(rasterLayer = awt,
                             sampleNumber = 1000,
                             doParallel = TRUE,
                             nCores = NULL,
                             plotVariograms = TRUE,
                             showPlots = TRUE)

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

  raster::crs(awt) <- NA

  expect_warning(range2 <- spatialAutoRange(rasterLayer = awt,
                             sampleNumber = 1000,
                             progress = TRUE,
                             doParallel = FALSE,
                             plotVariograms = FALSE,
                             showPlots = FALSE)

  )
  expect_true(exists("range2"))

})

test_that("test spatialAutoRange function with single-layer raster", {

  range3 <- spatialAutoRange(rasterLayer = awt[[1]],
                             sampleNumber = 1000,
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
