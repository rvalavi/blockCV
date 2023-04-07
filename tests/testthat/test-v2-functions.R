aus <- system.file("extdata/au/", package = "blockCV") |>
  list.files(full.names = TRUE) |>
  terra::rast()
nl <- terra::nlyr(aus)

pa_data <- read.csv(system.file("extdata/", "species.csv", package = "blockCV")) |>
  sf::st_as_sf(coords = c("x", "y"), crs = 7845)


# spatialAutoRange --------------------------------------------------------
expect_names <- c("range",
                  "rangeTable",
                  "plots",
                  "sampleNumber",
                  "variograms")


test_that("test spatialAutoRange function works", {
  # skip_on_cran()

  range1 <- spatialAutoRange(rasterLayer = aus,
                             sampleNumber = 1000,
                             showPlots = TRUE)

  expect_true(exists("range1"))
  expect_s3_class(range1, "SpatialAutoRange")
  expect_equal(names(range1), expect_names)
  expect_equal(nrow(range1$rangeTable), nl)
  expect_equal(length(range1$variograms), nl)
  expect_equal(dim(range1$rangeTable), c(nl, 3))
  expect_s3_class(range1$plots[[1]], "ggplot")
  expect_s3_class(range1$variograms[[1]], "autofitVariogram")
  expect_s3_class(range1$rangeTable, "data.frame")
  expect_type(range1$variograms, "list")
  expect_type(range1$sampleNumber, "double")
  expect_type(range1$range, "double")
  expect_true(range1$range >= 0)
  expect_true(!is.null(range1$variograms))
  expect_true(
    all(names(aus) %in% range1$rangeTable$layers)
  )

  expect_equal(print.SpatialAutoRange(range1), "SpatialAutoRange")
  expect_silent(plot.SpatialAutoRange(range1))
  expect_output(summary.SpatialAutoRange(range1))

})

test_that("test spatialAutoRange function with x", {
  # skip_on_cran()

  range3 <- spatialAutoRange(rasterLayer = aus,
                             speciesData = pa_data,
                             progress = FALSE,
                             showPlots = TRUE)

  expect_true(exists("range3"))
  expect_s3_class(range3, "SpatialAutoRange")
  expect_s3_class(range3$plots[[1]], "ggplot")
  expect_s3_class(range3$variograms[[1]], "autofitVariogram")
  expect_type(range3$range, "double")
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
  # skip_on_cran()

  set.seed(1000)
  sb1 <- spatialBlock(speciesData = pa_data,
                      species = "occ",
                      rasterLayer = aus,
                      theRange = 450000,
                      k = 5,
                      selection = "random",
                      iteration = 5,
                      biomod2Format = TRUE,
                      showBlocks = TRUE,
                      progress = TRUE)

  expect_true(exists("sb1"))
  expect_s3_class(sb1, "SpatialBlock")
  expect_equal(names(sb1), expect_names)
  expect_equal(length(sb1$folds), 5)
  expect_type(sb1$folds, "list")
  expect_type(sb1$biomodTable, "logical")
  expect_equal(dim(sb1$biomodTable), c(nrow(pa_data), 5))
  expect_equal(sb1$k, 5)
  # expect_s4_class(sb1$blocks, "SpatialPolygonsDataFrame")
  expect_s3_class(sb1$blocks, "sf")
  expect_type(sb1$species, "character")
  expect_type(sb1$range, "double")
  expect_equal(dim(sb1$records), c(5, 4))
  expect_true(
    !all(sb1$records == 0)
  )

  expect_output(print.SpatialBlock(sb1))
  expect_output(summary.SpatialBlock(sb1))
  expect_message(plot.SpatialBlock(sb1))

})


# buffering ---------------------------------------------------------------
expect_names <- c("folds",
                  "k",
                  "species",
                  "range",
                  "dataType",
                  "records")

test_that("test buffering function works", {
  # skip_on_cran()

  bf1 <- buffering(speciesData= pa_data,
                   species= "occ",
                   theRange= 450000,
                   spDataType = "PA",
                   progress = TRUE)

  expect_true(exists("bf1"))
  expect_s3_class(bf1, "BufferedBlock")
  expect_equal(names(bf1), expect_names)
  expect_equal(length(bf1$folds), nrow(pa_data))
  expect_type(bf1$folds, "list")
  expect_type(bf1$k, "integer")
  expect_type(bf1$species, "character")
  expect_type(bf1$range, "double")
  expect_equal(bf1$dataType, "PA")
  expect_equal(dim(bf1$records), c(nrow(pa_data), 4))
  expect_true(
    !all(bf1$records == 0)
  )

  expect_output(print.BufferedBlock(bf1))
  expect_output(summary.BufferedBlock(bf1))

})


# envBlock ----------------------------------------------------------------
expect_names <- c("folds",
                  "foldID",
                  "biomodTable",
                  "k",
                  "species",
                  "records")

test_that("test that environmental blocking function works", {
  # skip_on_cran()

  eb <- envBlock(rasterLayer = aus,
                 speciesData = pa_data,
                 species = "occ",
                 k = 5,
                 standardization = "standard",
                 rasterBlock = FALSE)

  expect_true(exists("eb"))
  expect_s3_class(eb, "EnvironmentalBlock")
  expect_equal(names(eb), expect_names)
  expect_equal(length(eb$folds), 5)
  expect_type(eb$folds, "list")
  expect_equal(eb$k, 5)
  expect_type(eb$species, "character")
  expect_equal(dim(eb$records), c(5, 4))
  expect_true(
    !all(eb$records == 0)
  )

  expect_output(print.EnvironmentalBlock(eb))
  expect_output(summary.EnvironmentalBlock(eb))

})
