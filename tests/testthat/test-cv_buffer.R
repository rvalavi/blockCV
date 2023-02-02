expect_names <- c("folds_list",
                  "k",
                  "column",
                  "size",
                  "presence_bg",
                  "records")


pa_data <- read.csv(system.file("extdata/", "species.csv", package = "blockCV")) |>
  sf::st_as_sf(coords = c("x", "y"), crs = 7845)
pa_data <- pa_data[1:200, ]

test_that("test that cv_buffer function works properly with presence-absence data",
          {
            bloo <- cv_buffer(
              x = pa_data,
              column = "occ",
              size = 250000,
              presence_bg = FALSE,
              progress = TRUE
            )

            expect_true(exists("bloo"))
            expect_s3_class(bloo, "cv_buffer")
            expect_equal(names(bloo), expect_names)
            expect_equal(length(bloo$folds), nrow(pa_data))
            expect_type(bloo$folds, "list")
            expect_type(bloo$k, "integer")
            expect_type(bloo$column, "character")
            expect_type(bloo$size, "double")
            expect_equal(bloo$presence_bg, FALSE)
            expect_equal(dim(bloo$records), c(nrow(pa_data), 4))
            expect_true(!all(bloo$records == 0))

})


test_that("test that cv_buffer function works properly with presence-background data",
          {
            bloo <- cv_buffer(
              x = pa_data,
              column = "occ",
              size = 250000,
              presence_bg = TRUE,
              progress = FALSE
            )

            expect_true(exists("bloo"))
            expect_s3_class(bloo, "cv_buffer")
            expect_equal(names(bloo), expect_names)
            expect_equal(length(bloo$folds), sum(pa_data$occ))
            expect_type(bloo$folds, "list")
            expect_type(bloo$k, "integer")
            expect_type(bloo$column, "character")
            expect_type(bloo$size, "double")
            expect_equal(bloo$presence_bg, TRUE)
            expect_equal(dim(bloo$records), c(sum(pa_data$occ), 4))
            expect_true(!all(bloo$records == 0))

})

test_that("test that cv_buffer function works properly with no species specified",
          {
            bloo <- cv_buffer(
              x = sf::as_Spatial(pa_data),
              size = 250000
            )

            expect_true(exists("bloo"))
            expect_s3_class(bloo, "cv_buffer")
            expect_equal(names(bloo), expect_names)
            expect_equal(length(bloo$folds), nrow(pa_data))
            expect_type(bloo$folds, "list")
            expect_type(bloo$k, "integer")
            expect_null(bloo$species)
            expect_type(bloo$size, "double")
            expect_equal(bloo$presence_bg, FALSE)
            expect_equal(dim(bloo$records), c(nrow(pa_data), 2))
            expect_true(!all(bloo$records == 0))

            expect_equal(print.cv_buffer(bloo), "cv_buffer")
            expect_output(summary.cv_buffer(bloo))

})


test_that("test cv_buffer function with no matching species column", {
  expect_warning(
    bloo <- cv_buffer(
      x = pa_data,
      column = "response",
      size = 250000,
      presence_bg = FALSE,
      progress = TRUE
    )
  )

  expect_true(exists("bloo"))
  expect_equal(dim(bloo$records), c(nrow(pa_data), 2))
  expect_true(!all(bloo$records == 0))

})

test_that("test cv_buffer function with no spatial column data", {
  expect_error(
    cv_buffer(
      x = "pa_data",
      size = 250000
    )
  )

})

test_that("test cv_buffer function to have sptial points with no CRS", {
  sf::st_crs(pa_data) <- NA

  expect_error(
    cv_buffer(
      x = pa_data,
      column = "occ",
      size = 250000
    )
  )

})
