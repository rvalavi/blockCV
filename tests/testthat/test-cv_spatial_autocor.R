aus <- system.file("extdata/au/", package = "blockCV") |>
  list.files(full.names = TRUE) |>
  terra::rast()
nl <- terra::nlyr(aus)

pa_data <- read.csv(system.file("extdata/", "species.csv", package = "blockCV")) |>
  sf::st_as_sf(coords = c("x", "y"), crs = 7845)



expect_names <- c("range",
                  "range_table",
                  "plots",
                  "num_sample",
                  "variograms")


test_that("test cv_spatial_autocor function works", {

  sac <- cv_spatial_autocor(
    r = aus,
    num_sample = 1000,
    plot = TRUE
  )

  expect_true(exists("sac"))
  expect_s3_class(sac, "cv_spatial_autocor")
  expect_equal(names(sac), expect_names)
  expect_equal(nrow(sac$range_table), nl)
  expect_equal(length(sac$variograms), nl)
  expect_equal(dim(sac$range_table), c(nl, 3))
  expect_s3_class(sac$plots[[1]], "ggplot")
  expect_s3_class(sac$variograms[[1]], "autofitVariogram")
  expect_s3_class(sac$range_table, "data.frame")
  expect_type(sac$variograms, "list")
  expect_type(sac$num_sample, "double")
  expect_type(sac$range, "double")
  expect_true(sac$range >= 0)
  expect_true(!is.null(sac$variograms))
  expect_true(all(names(aus) %in% sac$range_table$layers))

  expect_equal(print.cv_spatial_autocor(sac), "cv_spatial_autocor")
  expect_silent(plot.cv_spatial_autocor(sac))
  expect_output(summary.cv_spatial_autocor(sac))

})

test_that("test cv_spatial_autocor function with x", {

  sac <- cv_spatial_autocor(
    x = pa_data,
    column = "occ",
    progress = FALSE,
    plot = TRUE
  )

  expect_true(exists("sac"))
  expect_s3_class(sac, "cv_spatial_autocor")
  expect_s3_class(sac$variograms[[1]], "autofitVariogram")
  expect_type(sac$range, "double")
  expect_true(sac$range >= 0)
  expect_true(!is.null(sac$variograms))

})

test_that("test cv_spatial_autocor function works with wgs crs", {


  sac <- cv_spatial_autocor(
    r = terra::project(aus, "epsg:4326"),
    num_sample = 1000,
    plot = TRUE
  )

  expect_true(exists("sac"))
  expect_s3_class(sac, "cv_spatial_autocor")
  expect_equal(names(sac), expect_names)
  expect_equal(nrow(sac$range_table), nl)
  expect_equal(length(sac$variograms), nl)
  expect_equal(dim(sac$range_table), c(nl, 3))
  expect_s3_class(sac$plots[[1]], "ggplot")
  expect_s3_class(sac$variograms[[1]], "autofitVariogram")
  expect_s3_class(sac$range_table, "data.frame")
  expect_type(sac$variograms, "list")
  expect_type(sac$num_sample, "double")
  expect_type(sac$range, "double")
  expect_true(sac$range >= 0)
  expect_true(!is.null(sac$variograms))
  expect_true(all(names(aus) %in% sac$range_table$layers))

  expect_equal(print.cv_spatial_autocor(sac), "cv_spatial_autocor")
  expect_silent(plot.cv_spatial_autocor(sac))
  expect_output(summary.cv_spatial_autocor(sac))

})


test_that("test cv_spatial_autocor function gets error with x and without column", {

  expect_error(
    cv_spatial_autocor(
      x = pa_data,
      progress = FALSE,
      plot = TRUE
    )
  )

})
