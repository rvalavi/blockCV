aus <- system.file("extdata/au/", package = "blockCV") |>
    list.files(full.names = TRUE) |>
    terra::rast()

pa_data <- read.csv(system.file("extdata/", "species.csv", package = "blockCV")) |>
    sf::st_as_sf(coords = c("x", "y"), crs = 7845)
pa_data <- pa_data[1:200, ]

scv <- cv_spatial(
    x = pa_data,
    k = 5,
    selection = "random",
    iteration = 1,
    balance = FALSE,
    biomod2 = FALSE,
    plot = FALSE,
    report = FALSE,
    progress = FALSE
)


test_that("cv_distance works on a k-fold object and attaches W", {

    plt <- cv_distance(cv = scv, x = pa_data, r = aus, num_sample = 2000,
                       add_random = TRUE, plot = FALSE)

    expect_true(ggplot2::is_ggplot(plt))
    W <- attr(plt, "W")
    expect_true(is.numeric(W))
    expect_true(all(c("CV", "LOO", "Random") %in% names(W)))
    expect_true(all(W >= 0))
})


test_that("cv_distance works in feature space", {

    plt <- cv_distance(cv = scv, x = pa_data, r = aus, space = "feature",
                       num_sample = 2000, plot = FALSE)

    expect_true(ggplot2::is_ggplot(plt))
    expect_true("CV" %in% names(attr(plt, "W")))
})


test_that("cv_distance works on a leave-one-out object (no Random curve)", {
    bloo <- cv_buffer(
        x = pa_data,
        size = 250000,
        progress = FALSE,
        report = FALSE
    )

    expect_message(
        plt <- cv_distance(cv = bloo, x = pa_data, r = aus, num_sample = 2000, plot = FALSE),
        "leave-one-out"
    )
    expect_true(ggplot2::is_ggplot(plt))
    W <- attr(plt, "W")
    expect_false("Random" %in% names(W))
    expect_true(all(c("CV", "LOO") %in% names(W)))
})


test_that("cv_distance errors on a non-cv object", {
    expect_error(cv_distance(cv = list(a = 1), x = pa_data, r = aus, plot = FALSE))
})
