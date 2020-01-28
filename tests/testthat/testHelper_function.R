library(blockCV)

context("external helper functions fully")

awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))
PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
pa_data <- sf::st_as_sf(PA, coords = c("x", "y"), crs = raster::projection(awt))

# test_that("helper function with no CRS", {
#   skip_on_cran()
#
#   suppressWarnings(awt_wgs <- raster::projectRaster(from = awt, crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#   raster::projection(awt_wgs) <- NA
#
#   expect_warning(
#     net <- blockCV:::rasterNet(x = awt_wgs, resolution = 70000, mask = TRUE)
#   )
#
#   expect_true(exists("net"))
#   expect_is(net, "sf")
#
# })

test_that("helper function with no species data", {

  net <- blockCV:::rasterNet(x = pa_data, resolution = 70000, mask = TRUE)

  expect_true(exists("net"))
  expect_is(net, "sf")

})

test_that("helper function error with no block size", {

  expect_error(
    blockCV:::rasterNet(x = awt, mask = TRUE)
  )

})

test_that("helper function error with wrong ofset", {

  expect_error(
    blockCV:::rasterNet(x = pa_data, xOffset = 3)
  )

  expect_error(
    blockCV:::rasterNet(x = pa_data, yOffset = 3)
  )

})
