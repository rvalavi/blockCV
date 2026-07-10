pa_data <- sf::st_as_sf(
    read.csv(system.file("extdata/", "species.csv", package = "blockCV")),
    coords = c("x", "y"),
    crs = 7845
)


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

    # user-supplied fold colours must cover all folds
    expect_error(
        cv_plot(cv = scv, x = pa_data, combine_folds = TRUE,
                fold_colors = grDevices::rainbow(4)),
        "'fold_colors' must provide at least 5 colours"
    )

    # x must match the observations used to create the folds
    expect_error(
        cv_plot(cv = scv, x = pa_data[1:150, ], combine_folds = TRUE),
        "Number of rows in 'x' does not match the folds in 'cv'!"
    )

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


test_that("point-based plots explain and use the original sample data", {
    bloo <- cv_buffer(
        x = pa_data[1:50, ],
        column = "occ",
        size = 250000,
        presence_bg = FALSE,
        progress = FALSE
    )

    expect_error(
        cv_plot(cv = bloo),
        "The original sample data are required"
    )
    expect_error(
        plot(bloo),
        "The original sample data are required"
    )

    plt <- plot(bloo, pa_data[1:50, ], num_plots = 1:2)
    expect_true(ggplot2::is_ggplot(plt))

    plt_data <- plot(bloo, data = pa_data[1:50, ], num_plots = 1)
    expect_true(ggplot2::is_ggplot(plt_data))

    cl <- cv_cluster(
        x = pa_data[1:50, ],
        column = "occ",
        k = 3,
        report = FALSE
    )
    expect_error(
        plot(cl),
        "The original sample data are required"
    )
    plt_cluster <- plot(cl, pa_data[1:50, ], combine_folds = TRUE)
    expect_true(ggplot2::is_ggplot(plt_cluster))
})


test_that("cv_plot only notes the background fade when bg_alpha < points_alpha", {
    clpb <- cv_cluster(
        x = pa_data[1:100, ],
        column = "occ",
        k = 3,
        presence_bg = TRUE,
        report = FALSE
    )

    # background more transparent than presences: the fade note is drawn
    p_fade <- cv_plot(cv = clpb, x = pa_data[1:100, ], combine_folds = TRUE,
                      points_alpha = 0.7, bg_alpha = 0.1)
    expect_match(p_fade$labels$caption, "faded")

    # bg_alpha == points_alpha disables the fade, so no note
    p_equal <- cv_plot(cv = clpb, x = pa_data[1:100, ], combine_folds = TRUE,
                       points_alpha = 0.7, bg_alpha = 0.7)
    expect_null(p_equal$labels$caption)

    # bg_alpha above points_alpha is capped (with a warning) and not called faded
    expect_warning(
        p_cap <- cv_plot(cv = clpb, x = pa_data[1:100, ], combine_folds = TRUE,
                         points_alpha = 0.2, bg_alpha = 0.9),
        "cannot exceed"
    )
    expect_null(p_cap$labels$caption)
})
