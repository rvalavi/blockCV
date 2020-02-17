library(blockCV)

context("envBlock function")

expect_names <- c("folds",
                  "foldID",
                  "biomodTable",
                  "k",
                  "species",
                  "records")

awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))
awt <- awt[[1:3]]
PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
pa_data <- sf::st_as_sf(PA, coords = c("x", "y"), crs = raster::projection(awt))

test_that("test that environmental blocking function with rasterBlock, standard and species column", {

  # environmental clustering
  expect_warning( # there sould be not enough records in the folds
  eb <- envBlock(rasterLayer = awt,
                 speciesData = pa_data,
                 species = "Species", # name of the column with species data
                 k = 5,
                 standardization = "standard",
                 rasterBlock = TRUE)
  )

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

  # environmental clustering
  eb2 <- envBlock(rasterLayer = awt,
                  speciesData = sf::as_Spatial(pa_data),
                  k = 5,
                  standardization = "normal",
                  biomod2Format = FALSE,
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


test_that("test that environmental blocking with no standardisation and species column", {

  expect_warning(
  eb <- envBlock(rasterLayer = awt,
                 speciesData = pa_data,
                 species = "response", # wrong column name
                 k = 5,
                 standardization = "none",
                 rasterBlock = TRUE)
  )

  expect_true(exists("eb"))
  expect_is(eb, "EnvironmentalBlock")
  expect_equal(names(eb), expect_names)
  expect_equal(length(eb$folds), 5)
  expect_is(eb$folds, "list")
  expect_equal(eb$k, 5)
  expect_equal(dim(eb$records), c(5, 2))
  expect_true(
    !all(eb$records == 0)
  )

})


test_that("test environmental blocking with no spatial or sf object", {

  expect_error(
    envBlock(rasterLayer = awt,
             speciesData = awt, # no spatial or sf object
             k = 5,
             standardization = "normal",
             rasterBlock = FALSE)

  )

})

test_that("test environmental blocking with no raster data", {
  expect_error(
    envBlock(rasterLayer = pa_data, # no raster data
             speciesData = pa_data,
             k = 5,
             standardization = "normal",
             rasterBlock = FALSE)

  )

})

test_that("test environmental blocking with in-complete raster - NA extraction", {

  ext <- c(xmin = 206961.5,
           xmax = 481078.8,
           ymin = 8038112,
           ymax = 8308669 )
  awt2 <- raster::crop(awt, ext)
  expect_error(
    envBlock(rasterLayer = awt2, # in-complete raster - NA extraction
             speciesData = pa_data,
             k = 5,
             standardization = "normal",
             rasterBlock = FALSE)

  )

})

test_that("test environmental blocking with too many folds", {

  expect_warning(
    envBlock(rasterLayer = awt,
             speciesData = pa_data[sample(nrow(pa_data), 50),],
             k = 15,
             standardization = "normal",
             numLimit = 10,
             rasterBlock = FALSE)
  )

})

