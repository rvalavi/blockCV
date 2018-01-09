library(blockCV)

context("Buffering function")

expect_names <- c("folds",
                  "k",
                  "species",
                  "range",
                  "dataType",
                  "records")

test_that("test that buffering function works properly with presence-absence data", {

  PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
  Zone55s <- "+proj=utm +zone=55 +south +ellps=GRS80 +units=m +no_defs"
  pa_data <- sp::SpatialPointsDataFrame(PA[,c("x", "y")], PA, proj4string=sp::CRS(Zone55s))

  # buffering with presence-absence data
  bf1 <- buffering(speciesData= pa_data,
                   species= "Species", # to count the number of presences and absences
                   theRange= 68000,
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


test_that("test that buffering function works properly with presence-background data", {

  PB <- read.csv(system.file("extdata", "PB.csv", package = "blockCV"))
  Zone55s <- "+proj=utm +zone=55 +south +ellps=GRS80 +units=m +no_defs"
  pb_data <- sp::SpatialPointsDataFrame(PB[,c("x", "y")], PB, proj4string=sp::CRS(Zone55s))

  # buffering with presence-background data
  bf2 <- buffering(speciesData= pb_data,
                   species= "Species",
                   theRange= 68000,
                   spDataType = "PB",
                   progress = FALSE)

  expect_true(exists("bf2"))
  expect_is(bf2, "BufferedBlock")
  expect_equal(names(bf2), expect_names)
  expect_equal(length(bf2$folds), 116)
  expect_is(bf2$folds, "list")
  expect_is(bf2$k, "integer")
  expect_is(bf2$species, "character")
  expect_is(bf2$range, "numeric")
  expect_equal(bf2$dataType, "PB")
  expect_equal(dim(bf2$records), c(116, 4))
  expect_true(
    !all(bf2$records == 0)
  )

})

test_that("test that buffering function works properly with no species specified", {

  PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
  Zone55s <- "+proj=utm +zone=55 +south +ellps=GRS80 +units=m +no_defs"
  pa_data <- sp::SpatialPointsDataFrame(PA[,c("x", "y")], PA, proj4string=sp::CRS(Zone55s))

  # buffering with presence-absence data
  bf3 <- buffering(speciesData= pa_data,
                   theRange= 68000,
                   spDataType = "PA",
                   progress = TRUE)

  expect_true(exists("bf3"))
  expect_is(bf3, "BufferedBlock")
  expect_equal(names(bf3), expect_names)
  expect_equal(length(bf3$folds), nrow(pa_data))
  expect_is(bf3$folds, "list")
  expect_is(bf3$k, "integer")
  expect_null(bf3$species)
  expect_is(bf3$range, "numeric")
  expect_equal(bf3$dataType, "PA")
  expect_equal(dim(bf3$records), c(nrow(pa_data), 2))
  expect_true(
    !all(bf3$records == 0)
  )

  expect_equal(print.BufferedBlock(bf3), "BufferedBlock")
  expect_output(summary.BufferedBlock(bf3))

})
