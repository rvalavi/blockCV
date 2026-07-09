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
    column = "occ",
    k = 5,
    selection = "random",
    iteration = 1,
    biomod2 = FALSE,
    plot = FALSE,
    report = FALSE,
    progress = FALSE
)


test_that("cv_summary returns structural info without a raster", {
    s <- cv_summary(scv)

    expect_s3_class(s, "cv_summary")
    expect_named(s, c("n_folds", "is_loo", "records", "distances", "novelty", "warnings"))
    expect_equal(s$n_folds, length(scv$folds_list))
    expect_false(s$is_loo)
    expect_s3_class(s$records, "data.frame")
    expect_null(s$distances)
    expect_null(s$novelty)
    expect_s3_class(s$warnings, "data.frame")
    expect_named(s$warnings, c("fold", "type", "message"))
})


test_that("cv_summary adds distance and novelty diagnostics when a raster is given", {
    s <- cv_summary(scv, x = pa_data, r = aus, num_sample = 2000, seed = 1, progress = FALSE)

    expect_s3_class(s$distances, "data.frame")
    expect_s3_class(s$novelty, "data.frame")
    expect_equal(nrow(s$distances), length(scv$folds_list))
    expect_equal(nrow(s$novelty), length(scv$folds_list))
})


test_that("cv_summary needs x for the raster diagnostics", {
    expect_error(cv_summary(scv, r = aus), "x")
})


test_that("structural warnings flag degenerate folds (as data, not conditions)", {
    records <- data.frame(
        train_0 = c(100, 100, 100, 100),
        train_1 = c(  0,  60,  60,  60),  # fold1: class 1 absent from training
        test_0  = c( 10,  10,   2,  10),  # fold3: only 2 test points (tiny)
        test_1  = c(  0,  30,   1,  29)   # fold1: single-class test; fold3: imbalance vs ~15
    )
    fake <- list(folds_list = vector("list", 4), records = records, column = "occ")

    w <- expect_no_warning(
        blockCV:::.cv_warnings(fake, distances = NULL, min_test = 5L, is_loo = FALSE)
    )
    expect_s3_class(w, "data.frame")
    types <- w$type
    expect_true("class_missing_train" %in% types)
    expect_true("single_class_test" %in% types)
    expect_true("tiny_test" %in% types)
    expect_true("imbalance" %in% types)

    # leave-one-out objects skip the per-fold structural checks
    expect_equal(nrow(blockCV:::.cv_warnings(fake, is_loo = TRUE)), 0L)

    # high_leakage only fires when a distance summary is supplied
    dist <- data.frame(fold = 1:4, pct_below_pred = c(95, 40, 92, 10))
    w2 <- blockCV:::.cv_warnings(fake, distances = dist, min_test = 5L, is_loo = FALSE)
    expect_setequal(w2$fold[w2$type == "high_leakage"], c(1, 3))
})


test_that("cv_summary works on a leave-one-out object", {
    bloo <- cv_buffer(x = pa_data, size = 250000, progress = FALSE, report = FALSE)

    s <- cv_summary(bloo)
    expect_true(s$is_loo)
    # structural warnings are skipped for leave-one-out designs
    expect_equal(nrow(s$warnings), 0L)
})
