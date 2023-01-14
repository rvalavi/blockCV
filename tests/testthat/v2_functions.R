library(blockCV)
library(automap)

context("Functions of v2")

awt <- system.file("extdata/au/", package = "blockCV") |>
  list.files(full.names = TRUE) |>
  terra::rast()
nl <- terra::nlyr(awt)

pa_data <- read.csv(system.file("extdata/", "species.csv", package = "blockCV")) |>
  sf::st_as_sf(coords = c("x", "y"), crs = 7845)

# spatialAutoRange --------------------------------------------------------
expect_names <- c("range",
                  "rangeTable",
                  "plots",
                  "sampleNumber",
                  "variograms")


test_that("test spatialAutoRange function works", {
  skip_on_cran()

  range1 <- spatialAutoRange(rasterLayer = awt,
                             sampleNumber = 1000,
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

  expect_equal(print.SpatialAutoRange(range1), "SpatialAutoRange")
  expect_silent(plot.SpatialAutoRange(range1))
  expect_output(summary.SpatialAutoRange(range1))

})

test_that("test spatialAutoRange function with x", {
  skip_on_cran()

  range3 <- spatialAutoRange(rasterLayer = awt,
                             speciesData = pa_data,
                             progress = FALSE,
                             showPlots = TRUE)

  expect_true(exists("range3"))
  expect_is(range3, "SpatialAutoRange")
  expect_is(range3$plots[[1]], "ggplot")
  expect_is(range3$variograms[[1]], "autofitVariogram")
  expect_is(range3$range, "numeric")
  expect_true(range3$range >= 0)
  expect_true(!is.null(range3$variograms))

  expect_equal(print.SpatialAutoRange(range3), "SpatialAutoRange")
  expect_silent(plot.SpatialAutoRange(range3))

})


# spatialBlock ------------------------------------------------------------
expect_names <- c("folds",
                  "foldID",
                  "biomodTable",
                  "k",
                  "blocks",
                  "species",
                  "range",
                  "plots",
                  "records")


test_that("test spatiaBlock function works", {
  skip_on_cran()

  set.seed(1000)
  sb1 <- spatialBlock(speciesData = pa_data,
                      species = "occ",
                      rasterLayer = awt,
                      theRange = 450000,
                      k = 5,
                      selection = "random",
                      iteration = 5,
                      biomod2Format = TRUE,
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
  expect_is(sb1$range, "numeric")
  expect_equal(dim(sb1$records), c(5, 4))
  expect_true(
    !all(sb1$records == 0)
  )

})


# buffering ---------------------------------------------------------------
expect_names <- c("folds",
                  "k",
                  "species",
                  "range",
                  "dataType",
                  "records")

test_that("test buffering function works", {
  skip_on_cran()

  # buffering with presence-absence data
  bf1 <- buffering(speciesData= pa_data,
                   species= "occ",
                   theRange= 450000,
                   spDataType = "PA",
                   progress = TRUE)

  expect_true(exists("bf1"))
  expect_is(bf1, "BufferedBlock")
  expect_equal(names(bf1), expect_names)
  expect_equal(length(bf1$folds), nrow(pa_data))
  expect_is(bf1$folds, "list")
  expect_is(bf1$k, "integer")
  expect_is(bf1$species, "character")
  expect_is(bf1$range, "numeric")
  expect_equal(bf1$dataType, "PA")
  expect_equal(dim(bf1$records), c(nrow(pa_data), 4))
  expect_true(
    !all(bf1$records == 0)
  )

})


# envBlock ----------------------------------------------------------------
expect_names <- c("folds",
                  "foldID",
                  "biomodTable",
                  "k",
                  "species",
                  "records")

test_that("test that environmental blocking function works", {
  skip_on_cran()

  # environmental clustering
  eb <- envBlock(rasterLayer = awt,
                 speciesData = pa_data,
                 species = "occ", # name of the column with species data
                 k = 5,
                 standardization = "standard",
                 rasterBlock = FALSE)

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
