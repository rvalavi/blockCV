expect_names <- c("folds_list",
                  "folds_ids",
                  "biomod_table",
                  "k",
                  "size",
                  "column",
                  "blocks",
                  "records")

aus <- system.file("extdata/au/", package = "blockCV") |>
  list.files(full.names = TRUE) |>
  terra::rast()

pa_data <- read.csv(system.file("extdata/", "species.csv", package = "blockCV")) |>
  sf::st_as_sf(coords = c("x", "y"), crs = 7845)
pa_data <- pa_data[1:200, ]

test_that("test cv_spatial function with random assingment and raster file", {

  set.seed(1000)
  scv <- cv_spatial(x = pa_data,
                    column = "occ",
                    r = aus,
                    size = 450000,
                    k = 5,
                    selection = "random",
                    iteration = 5,
                    biomod2 = TRUE,
                    offset = c(0.3, 0.2),
                    plot = FALSE,
                    progress = TRUE)

  expect_true(exists("scv"))
  expect_s3_class(scv, "cv_spatial")
  expect_equal(names(scv), expect_names)
  expect_equal(length(scv$folds_list), 5)
  expect_type(scv$folds_list, "list")
  expect_type(scv$biomod_table, "logical")
  expect_equal(dim(scv$biomod_table), c(nrow(pa_data), 5))
  expect_equal(scv$k, 5)
  expect_s3_class(scv$blocks, "sf")
  expect_type(scv$column, "character")
  expect_type(scv$size, "double")
  expect_equal(dim(scv$records), c(5, 4))
  expect_true(
    !all(scv$records == 0)
  )

})


test_that("test cv_spatial function with systematic assingment and no raster file", {

  scv <- cv_spatial(x = sf::as_Spatial(pa_data),
                    rows_cols = c(10, 10),
                    k = 5,
                    selection = "systematic",
                    biomod2 = FALSE,
                    plot = TRUE)

  expect_true(exists("scv"))
  expect_s3_class(scv, "cv_spatial")
  expect_equal(names(scv), expect_names)
  expect_equal(length(scv$folds_list), 5)
  expect_type(scv$folds_list, "list")
  expect_null(scv$biomod_table, "matrix")
  expect_equal(scv$k, 5)
  expect_s3_class(scv$blocks, "sf")
  expect_null(scv$column)
  expect_null(scv$size)
  expect_equal(dim(scv$records), c(5, 2))
  expect_true(
    !all(scv$records == 0)
  )

  expect_equal(print.cv_spatial(scv), "cv_spatial")
  expect_message(plot.cv_spatial(scv))
  expect_output(summary.cv_spatial(scv))

})


test_that("test cv_spatial function with non-numeric iteration", {

  scv <- cv_spatial(x = pa_data,
                    hexagon = FALSE,
                    k = 5,
                    selection = "random",
                    iteration = "2", # not numeric
                    biomod2 = FALSE,
                    plot = FALSE)

  expect_true(exists("scv"))
  expect_s3_class(scv, "cv_spatial")
  expect_equal(names(scv), expect_names)
  expect_equal(length(scv$folds_list), 5)
  expect_type(scv$folds_list, "list")
  expect_null(scv$biomod_table, "matrix")
  expect_equal(scv$k, 5)
  expect_s3_class(scv$blocks, "sf")
  expect_null(scv$column)
  expect_null(scv$size)
  expect_equal(dim(scv$records), c(5, 2))
  expect_true(
    !all(scv$records == 0)
  )

  expect_equal(print.cv_spatial(scv), "cv_spatial")
  expect_message(plot.cv_spatial(scv))
  expect_output(summary.cv_spatial(scv))

})

test_that("test cv_spatial with checkerboard assingment and only row blocks", {

  scv <- cv_spatial(x = sf::as_Spatial(pa_data),
                    hexagon = FALSE,
                    rows_cols = c(5, 1),
                    selection = "checkerboard",
                    biomod2 = FALSE,
                    plot = TRUE)

  expect_true(exists("scv"))
  expect_s3_class(scv, "cv_spatial")
  expect_equal(names(scv), expect_names)
  expect_equal(length(scv$folds_list), 2)
  expect_type(scv$folds_list, "list")
  expect_null(scv$biomod_table)
  expect_equal(scv$k, 2)
  expect_s3_class(scv$blocks, "sf")
  expect_null(scv$column)
  expect_null(scv$size)
  expect_equal(dim(scv$records), c(2, 2))
  expect_true(
    !all(scv$records == 0)
  )

  expect_equal(print.cv_spatial(scv), "cv_spatial")
  expect_message(plot.cv_spatial(scv))
  expect_output(summary.cv_spatial(scv))

})


test_that("test cv_spatial with user-defined blocks", {

  sb <- cv_spatial(x = pa_data,
                   size = 450000,
                   selection = "random",
                   iteration = 1,
                   biomod2 = FALSE,
                   plot = FALSE)

  scv <- cv_spatial(x = pa_data,
                    user_blocks = sb$blocks,
                    selection = "random",
                    iteration = 5,
                    biomod2 = FALSE,
                    progress = FALSE,
                    plot = FALSE)

  expect_true(exists("scv"))
  expect_s3_class(scv, "cv_spatial")
  expect_equal(names(scv), expect_names)
  expect_equal(length(scv$folds_list), 5)
  expect_type(scv$folds_list, "list")
  expect_null(scv$biomod_table)
  expect_equal(scv$k, 5)
  expect_s3_class(scv$blocks, "sf")
  expect_null(scv$column)
  expect_null(scv$size)
  expect_equal(dim(scv$records), c(5, 2))
  expect_true(
    !all(scv$records == 0)
  )

  expect_equal(print.cv_spatial(scv), "cv_spatial")
  expect_message(plot.cv_spatial(scv))
  expect_output(summary.cv_spatial(scv))

})


test_that("test cv_spatial failur: number of blocks, wrong selection", {

  expect_error(cv_spatial(x = pa_data,
                          rows_cols = c(2, 2),
                          hexagon = FALSE,
                          k = 15, # very high k
                          selection = "random"))

  expect_error(cv_spatial(x = pa_data,
                          k = 1, # very low k
                          selection = "random"))

  expect_error(cv_spatial(x = pa_data,
                          k = 5,
                          selection = "chance")) # wrong selection

})


test_that("test cv_spatial failur: wrong user-defined blocks", {

  expect_error(
    cv_spatial(x = pa_data, # wrong user-defined blocks
               column = "occ",
               r = aus,
               user_blocks = data.frame(a = 1),
               size = 450000,
               k = 5)
  )

})


test_that("test cv_spatial with no speceis column match", {

  expect_warning(
    scv <- cv_spatial(x = pa_data,
                      column = "response", # wrong response
                      r = aus,
                      size = 450000,
                      k = 5)
  )

  expect_true(exists("scv"))
  expect_s3_class(scv, "cv_spatial")
  expect_equal(names(scv), expect_names)
  expect_equal(length(scv$folds_list), 5)
  expect_equal(scv$k, 5)
  expect_equal(dim(scv$records), c(5, 2))
  expect_true(
    !all(scv$records == 0)
  )

})
