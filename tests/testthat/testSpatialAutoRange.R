library(blockCV)

context("Spatial Autocorrelation Finder")

expect_names <- c("range",
                  "rangeTable",
                  "plots",
                  "sampleNumber",
                  "variograms")

test_that("test spatialAutoRange function", {

  awt <- raster::brick(system.file("extdata", "awt.tif", package = "blockCV"))
  names(awt) <- c("bc01",  "bc04",  "bc05",  "bc06",  "bc12",  "bc15",  "bc17",  "bc20",  "bc33", "slope", "topo")
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
})

