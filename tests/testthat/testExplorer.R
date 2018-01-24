library(blockCV)

context("Explorer function")

test_that("test that explorer function with presence-absence data", {

  awt <- raster::brick(system.file("extdata", "awt.tif", package = "blockCV"))
  PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
  pa_data <- sp::SpatialPointsDataFrame(PA[,c("x", "y")], PA, proj4string=raster::crs(awt))

  sb <- spatialBlock(speciesData = pa_data,
                     species = "Species",
                     rasterLayer = awt,
                     theRange = 66000,
                     k = 5,
                     selection = 'random',
                     iteration = 250,
                     numLimit = NULL,
                     biomod2Format = TRUE)

  expect_silent(foldExplorer(sb, awt, pa_data))

  # buffering with presence-absence data
  bf <- buffering(speciesData= pa_data,
                  species= "Species", # to count the number of presences and absences
                  theRange= 66500,
                  spDataType = "PA",
                  progress = T)

  expect_silent(foldExplorer(bf, awt, pa_data))

  # environmental clustering
  eb <- envBlock(rasterLayer = awt,
                 speciesData = pa_data,
                 species = "Species",
                 k = 5)

  expect_silent(foldExplorer(eb, awt, pa_data))

  expect_silent(rangeExplorer(rasterLayer = awt,
                speciesData = pa_data,
                species = "Species",
                rangeTable = NULL,
                minRange = 30000, # limit the search domain
                maxRange = 100000))


})
