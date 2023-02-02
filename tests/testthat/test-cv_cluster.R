expect_names <- c("folds_list",
                  "folds_ids",
                  "biomod_table",
                  "k",
                  "column",
                  "type",
                  "records")


aus <- system.file("extdata/au/", package = "blockCV") |>
  list.files(full.names = TRUE) |>
  terra::rast()

pa_data <- read.csv(system.file("extdata/", "species.csv", package = "blockCV")) |>
  sf::st_as_sf(coords = c("x", "y"), crs = 7845)


test_that("test that environmental cluster function with raster_cluster", {

  # environmental clustering
  set.seed(42)
  eb <- cv_cluster(r = aus,
                   x = pa_data,
                   column = "occ",
                   k = 3,
                   scale = TRUE,
                   raster_cluster = TRUE)

  expect_true(exists("eb"))
  expect_s3_class(eb, "cv_cluster")
  expect_equal(names(eb), expect_names)
  expect_equal(length(eb$folds_list), 3)
  expect_type(eb$folds_list, "list")
  expect_equal(eb$k, 3)
  expect_type(eb$column, "character")
  expect_equal(dim(eb$records), c(3, 4))
  expect_true(
    !all(eb$records == 0)
  )

})

test_that("test that spacial cluster function with no column", {

  # spatial clustering
  set.seed(42)
  eb <- cv_cluster(x = sf::as_Spatial(pa_data),
                   k = 5,
                   biomod2 = FALSE)

  expect_true(exists("eb"))
  expect_s3_class(eb, "cv_cluster")
  expect_equal(names(eb), expect_names)
  expect_equal(length(eb$folds_list), 5)
  expect_type(eb$folds_list, "list")
  expect_equal(eb$k, 5)
  expect_null(eb$column)
  expect_equal(dim(eb$records), c(5, 2))
  expect_true(
    !all(eb$records == 0)
  )

  expect_equal(print.cv_cluster(eb), "cv_cluster")
  expect_output(summary.cv_cluster(eb))

})


test_that("test that environmental cluster with no scale and wrong column", {

  set.seed(42)
  expect_warning(
    eb <- cv_cluster(r = aus,
                     x = pa_data,
                     column = "response", # wrong column name
                     k = 5,
                     scale = FALSE,
                     raster_cluster = TRUE,
                     algorithm = "MacQueen")
  )

  expect_true(exists("eb"))
  expect_s3_class(eb, "cv_cluster")
  expect_equal(names(eb), expect_names)
  expect_equal(length(eb$folds_list), 5)
  expect_type(eb$folds_list, "list")
  expect_equal(eb$k, 5)
  expect_equal(dim(eb$records), c(5, 2))
  expect_true(
    !all(eb$records == 0)
  )

})


test_that("test environmental cluster with no spatial or sf object", {

  expect_error(
    cv_cluster(r = aus,
               x = aus, # no spatial or sf object
               k = 5,
               scale = TRUE,
               raster_cluster = FALSE)

  )

})

test_that("test environmental cluster with no raster data", {
  expect_error(
    cv_cluster(r = data.frame(x = 1), # no raster data
               x = pa_data)

  )

})

