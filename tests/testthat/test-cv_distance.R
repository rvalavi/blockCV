aus <- terra::rast(
    list.files(system.file("extdata/au/", package = "blockCV"), full.names = TRUE)
)

pa_data <- sf::st_as_sf(
    read.csv(system.file("extdata/", "species.csv", package = "blockCV")),
    coords = c("x", "y"),
    crs = 7845
)
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


test_that("cv_distance returns a cv_distance object with W", {

    res <- cv_distance(cv = scv, x = pa_data, r = aus, num_sample = 2000,
                       add_random = TRUE, plot = FALSE)

    expect_s3_class(res, "cv_distance")
    expect_named(res, c("distances", "W", "plot"))
    W <- res$W
    expect_true(is.numeric(W))
    expect_true(all(c("CV", "LOO", "Random") %in% names(W)))
    expect_true(all(W >= 0))
})


test_that("cv_distance returns a per-fold distance summary", {

    res <- cv_distance(cv = scv, x = pa_data, r = aus, num_sample = 2000,
                       add_random = FALSE, plot = FALSE)

    d <- res$distances
    expect_s3_class(d, "data.frame")
    # one row per fold
    expect_equal(nrow(d), length(scv$folds_list))
    expect_setequal(d$fold, seq_along(scv$folds_list))
    expect_named(d, c("fold", "n_test", "min", "q1", "median", "q3", "max", "pct_below_pred"))
    # quartiles are ordered and distances non-negative
    expect_true(all(d$min <= d$median & d$median <= d$max))
    expect_true(all(d$min >= 0))
    # pct_below_pred is a valid percentage
    expect_true(all(d$pct_below_pred >= 0 & d$pct_below_pred <= 100))
})


test_that("plot toggle controls drawing only, not the returned object", {

    # the ggplot is always built and returned, whether or not it is drawn
    res0 <- cv_distance(cv = scv, x = pa_data, r = aus, num_sample = 2000,
                        seed = 1, plot = FALSE)
    res1 <- cv_distance(cv = scv, x = pa_data, r = aus, num_sample = 2000,
                        seed = 1, plot = TRUE)
    expect_true(ggplot2::is_ggplot(res0$plot))
    expect_true(ggplot2::is_ggplot(res1$plot))
    # same data regardless of the plot toggle
    expect_equal(res0$distances, res1$distances)
    expect_equal(res0$W, res1$W)
})


test_that("cv_distance works in feature space", {

    res <- cv_distance(cv = scv, x = pa_data, r = aus, space = "feature",
                       num_sample = 2000, plot = FALSE)

    expect_s3_class(res, "cv_distance")
    expect_true("CV" %in% names(res$W))
})


test_that("cv_distance works on a leave-one-out object (no Random curve)", {
    bloo <- cv_buffer(
        x = pa_data,
        size = 250000,
        progress = FALSE,
        report = FALSE
    )

    expect_message(
        res <- cv_distance(cv = bloo, x = pa_data, r = aus, num_sample = 2000, plot = FALSE),
        "leave-one-out"
    )
    expect_s3_class(res, "cv_distance")
    W <- res$W
    expect_false("Random" %in% names(W))
    expect_true(all(c("CV", "LOO") %in% names(W)))
})


test_that("cv_distance errors on a non-cv object", {
    expect_error(cv_distance(cv = list(a = 1), x = pa_data, r = aus, plot = FALSE))
})
