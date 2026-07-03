pa_data <- read.csv(system.file("extdata/", "species.csv", package = "blockCV")) |>
    sf::st_as_sf(coords = c("x", "y"), crs = 7845)


test_that("test that cv_plot function works", {
    scv <- cv_spatial(
        x = pa_data,
        size = 450000,
        k = 5,
        selection = "random",
        iteration = 1,
        biomod2 = FALSE,
        plot = FALSE,
        progress = FALSE
    )

    plt <- cv_plot(cv = scv, x = pa_data)

    expect_true(exists("plt"))
    expect_true(ggplot2::is_ggplot(plt))

})


test_that("cv_plot combine_folds shows a single fold map for k-fold objects", {
    scv <- cv_spatial(
        x = pa_data,
        size = 450000,
        k = 5,
        selection = "random",
        iteration = 1,
        biomod2 = FALSE,
        plot = FALSE,
        progress = FALSE
    )

    # default auto-generated fold palette
    plt <- cv_plot(cv = scv, x = pa_data, combine_folds = TRUE)
    expect_true(ggplot2::is_ggplot(plt))

    # user-supplied fold colours
    plt2 <- cv_plot(cv = scv, x = pa_data, combine_folds = TRUE,
                    fold_colors = grDevices::rainbow(5))
    expect_true(ggplot2::is_ggplot(plt2))

    # leave-one-out objects are not supported
    bloo <- cv_buffer(
        x = pa_data[1:150, ],
        column = "occ",
        size = 250000,
        presence_bg = FALSE,
        progress = FALSE
    )
    expect_error(cv_plot(cv = bloo, x = pa_data[1:150, ], combine_folds = TRUE))

})
