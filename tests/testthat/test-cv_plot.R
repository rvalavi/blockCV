pa_data <- read.csv(system.file("extdata/", "species.csv", package = "blockCV")) |>
  sf::st_as_sf(coords = c("x", "y"), crs = 7845)


test_that("test that cv_plot function works",
          {
            scv <- cv_spatial(x = pa_data,
                              size = 450000,
                              k = 5,
                              selection = "random",
                              iteration = 1,
                              biomod2 = FALSE,
                              plot = FALSE,
                              progress = FALSE)

            plt <- cv_plot(cv = scv, x = pa_data)

            expect_true(exists("plt"))
            expect_s3_class(plt, "ggplot")
            expect_type(plt, "list")

})
