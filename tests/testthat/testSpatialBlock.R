library(blockCV)

context("spatialBlock function")

expect_names <- c("folds",
                  "foldID",
                  "biomodTable",
                  "k",
                  "blocks",
                  "species",
                  "range",
                  "plots",
                  "records")

awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))
PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
pa_data <- sf::st_as_sf(PA, coords = c("x", "y"), crs = raster::projection(awt))

# r <- awt[[1]]
# r[r > 0] <- 1
# bound <- raster::rasterToPolygons(r, dissolve = TRUE)

test_that("test spatiaBlock function with random assingment and raster file", {

  set.seed(1000)
  # expect_warning(
    sb1 <- spatialBlock(speciesData = pa_data,
                        species = "Species",
                        rasterLayer = awt,
                        theRange = 70000,
                        k = 5,
                        selection = "random",
                        # border = bound,
                        iteration = 5,
                        numLimit = 0,
                        biomod2Format = TRUE,
                        xOffset = 0.3,
                        yOffset = 0.2,
                        showBlocks = FALSE,
                        progress = TRUE)
  # )

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
  expect_null(sb1$plots)
  expect_equal(dim(sb1$records), c(5, 4))
  expect_true(
    !all(sb1$records == 0)
  )

})


test_that("test spatiaBlock function with systematic assingment and no raster file", {

  sb2 <- spatialBlock(speciesData = sf::as_Spatial(pa_data),
                      rows = 5,
                      cols = 8,
                      k = 5,
                      selection = "systematic",
                      biomod2Format = FALSE,
                      showBlocks = TRUE)

  expect_true(exists("sb2"))
  expect_is(sb2, "SpatialBlock")
  expect_equal(names(sb2), expect_names)
  expect_equal(length(sb2$folds), 5)
  expect_is(sb2$folds, "list")
  expect_null(sb2$biomodTable, "matrix")
  expect_equal(sb2$k, 5)
  expect_is(sb2$blocks, "SpatialPolygonsDataFrame")
  expect_null(sb2$species)
  expect_null(sb2$range)
  expect_is(sb2$plots, "ggplot")
  expect_equal(dim(sb2$records), c(5, 2))
  expect_true(
    !all(sb2$records == 0)
  )

  expect_equal(print.SpatialBlock(sb2), "SpatialBlock")
  expect_message(plot.SpatialBlock(sb2))
  expect_output(summary.SpatialBlock(sb2))

})


test_that("test spatiaBlock function with non-numeric iteration", {

  sb2 <- spatialBlock(speciesData = pa_data,
                      rows = 5,
                      cols = 8,
                      k = 5,
                      selection = "random",
                      iteration = "12", # not numeric
                      # numLimit = 2,
                      biomod2Format = FALSE,
                      showBlocks = TRUE)

  expect_true(exists("sb2"))
  expect_is(sb2, "SpatialBlock")
  expect_equal(names(sb2), expect_names)
  expect_equal(length(sb2$folds), 5)
  expect_is(sb2$folds, "list")
  expect_null(sb2$biomodTable, "matrix")
  expect_equal(sb2$k, 5)
  expect_is(sb2$blocks, "SpatialPolygonsDataFrame")
  expect_null(sb2$species)
  expect_null(sb2$range)
  expect_is(sb2$plots, "ggplot")
  expect_equal(dim(sb2$records), c(5, 2))
  expect_true(
    !all(sb2$records == 0)
  )

  expect_equal(print.SpatialBlock(sb2), "SpatialBlock")
  expect_message(plot.SpatialBlock(sb2))
  expect_output(summary.SpatialBlock(sb2))

})

test_that("test spatiaBlock with checkerboard assingment and only row blocks", {

  sb3 <- spatialBlock(speciesData = sf::as_Spatial(pa_data),
                      rows = 5,
                      selection = "checkerboard",
                      biomod2Format = FALSE,
                      showBlocks = TRUE)

  expect_true(exists("sb3"))
  expect_is(sb3, "SpatialBlock")
  expect_equal(names(sb3), expect_names)
  expect_equal(length(sb3$folds), 2)
  expect_is(sb3$folds, "list")
  expect_null(sb3$biomodTable)
  expect_equal(sb3$k, 2)
  expect_is(sb3$blocks, "SpatialPolygonsDataFrame")
  expect_null(sb3$species)
  expect_null(sb3$range)
  expect_is(sb3$plots, "ggplot")
  expect_equal(dim(sb3$records), c(2, 2))
  expect_true(
    !all(sb3$records == 0)
  )

  expect_equal(print.SpatialBlock(sb3), "SpatialBlock")
  expect_message(plot.SpatialBlock(sb3))
  expect_output(summary.SpatialBlock(sb3))

})


test_that("test spatiaBlock with user-defined blocks", {

  sb <- spatialBlock(speciesData = pa_data,
                     theRange = 70000,
                     selection = "random",
                     iteration = 1,
                     biomod2Format = FALSE,
                     showBlocks = FALSE)

  sb4 <- spatialBlock(speciesData = pa_data,
                      blocks = sb$blocks[,-2],
                      selection = "checkerboard",
                      biomod2Format = FALSE,
                      showBlocks = FALSE)

  expect_true(exists("sb4"))
  expect_is(sb4, "SpatialBlock")
  expect_equal(names(sb4), expect_names)
  expect_equal(length(sb4$folds), 2)
  expect_is(sb4$folds, "list")
  expect_null(sb4$biomodTable)
  expect_equal(sb4$k, 2)
  expect_is(sb4$blocks, "SpatialPolygonsDataFrame")
  expect_null(sb4$species)
  expect_null(sb4$range)
  expect_null(sb4$plots)
  expect_equal(dim(sb4$records), c(2, 2))
  expect_true(
    !all(sb4$records == 0)
  )

  expect_equal(print.SpatialBlock(sb4), "SpatialBlock")
  expect_message(plot.SpatialBlock(sb4))
  expect_output(summary.SpatialBlock(sb4))

})

test_that("test spatiaBlock function with NULL numLimit", {

  sb2 <- spatialBlock(speciesData = pa_data,
                      species = "Species",
                      rows = 5,
                      cols = 8,
                      k = 5,
                      selection = "random",
                      iteration = 5,
                      numLimit = NULL,
                      biomod2Format = FALSE,
                      showBlocks = FALSE)

  expect_true(exists("sb2"))
  expect_is(sb2, "SpatialBlock")
  expect_equal(names(sb2), expect_names)
  expect_equal(dim(sb2$records), c(5, 4))
  expect_true(
    !all(sb2$records == 0)
  )

})

test_that("test spatialBlock failur: number of blocks, wrong selection", {

  expect_error(spatialBlock(speciesData = pa_data,
                            cols = 5,
                            k = 15, # very high k
                            selection = "random"))

  expect_error(spatialBlock(speciesData = pa_data,
                            cols = 5,
                            rows = 8,
                            k = 1, # very low k
                            selection = "random"))

  expect_error(spatialBlock(speciesData = pa_data,
                            cols = 5,
                            k = 5,
                            selection = "rand")) # wrong selection

})

test_that("test spatialBlock failur: wrong species data", {

  expect_error(
    spatialBlock(speciesData = awt, # wrong speceis data
                 species = "Species",
                 rasterLayer = awt,
                 theRange = 70000,
                 k = 5)
  )

})

test_that("test spatialBlock failur: wrong user-defined blocks", {

  expect_error(
    spatialBlock(speciesData = pa_data, # wrong user-defined blocks
                 species = "Species",
                 rasterLayer = awt,
                 blocks = r,
                 theRange = 70000,
                 k = 5)
  )

})


test_that("test spatialBlock with no speceis column match", {

  expect_warning(
  sb4 <- spatialBlock(speciesData = pa_data,
                      species = "response", # wrong response
                      rasterLayer = awt,
                      maskBySpecies = FALSE,
                      theRange = 70000,
                      k = 5)
  )

  expect_true(exists("sb4"))
  expect_is(sb4, "SpatialBlock")
  expect_equal(names(sb4), expect_names)
  expect_equal(length(sb4$folds), 5)
  expect_equal(sb4$k, 5)
  expect_equal(dim(sb4$records), c(5, 2))
  expect_true(
    !all(sb4$records == 0)
  )

})

# skip these on CRAN ------------------------------------------------------

test_that("test spatialBlock failur: wrong user-defined border", {
  skip_on_cran()

  expect_error(
    spatialBlock(speciesData = pa_data, # wrong speceis data
                 species = "Species",
                 rasterLayer = awt,
                 border = r,
                 theRange = 70000,
                 k = 5)
  )

})

test_that("test spatialBlock with no smaller mask raster", {
  skip_on_cran()

  ext <- c(xmin = 206961.5,
           xmax = 481078.8,
           ymin = 8038112,
           ymax = 8308669 )
  awt2 <- raster::crop(awt, ext)
  expect_warning(
    sb4 <- spatialBlock(speciesData = pa_data,
                        species = "response", # wrong response
                        rasterLayer = awt2,
                        theRange = 70000,
                        k = 5)
  )
  expect_true(exists("sb4"))

})
