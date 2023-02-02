aus <- system.file("extdata/au/", package = "blockCV") |>
  list.files(full.names = TRUE) |>
  terra::rast()

pa_data <- read.csv(system.file("extdata/", "species.csv", package = "blockCV")) |>
  sf::st_as_sf(coords = c("x", "y"), crs = 7845)
pa_data <- pa_data[1:200, ]

test_that("test that cv_similarity function works with cv_spatil",
          {
            scv <- cv_spatial(x = pa_data,
                              selection = "random",
                              iteration = 1,
                              biomod2 = FALSE,
                              plot = FALSE,
                              report = FALSE,
                              progress = FALSE)

            plt <- cv_similarity(cv = scv, x = pa_data, r = aus)

            expect_true(exists("plt"))
            expect_s3_class(plt, "ggplot")
            expect_type(plt, "list")

})


test_that("test that cv_similarity function works with cv_buffer",
          {
            bloo <- cv_buffer(x = pa_data,
                              size = 250000,
                              progress = FALSE,
                              report = FALSE)

            plt <- cv_similarity(cv = bloo, x = pa_data, r = aus)

            expect_true(exists("plt"))
            expect_s3_class(plt, "ggplot")
            expect_type(plt, "list")

})
