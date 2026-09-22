# internal input-validation guards in R/checks.R

test_that(".check_x_matches_cv errors when x and cv sizes disagree", {
    cv <- structure(list(folds_list = list(list(1:8, 9:10))), class = "cv_spatial")

    expect_error(
        blockCV:::.check_x_matches_cv(data.frame(a = 1:5), cv),
        "does not match"
    )
    expect_true(blockCV:::.check_x_matches_cv(data.frame(a = 1:10), cv))
})


test_that(".check_num_plots keeps in-range folds and errors when none remain", {
    expect_equal(blockCV:::.check_num_plots(c(0L, 2L, 99L), k = 5L), 2L)
    expect_error(
        blockCV:::.check_num_plots(c(0L, 99L), k = 5L),
        "None of the requested"
    )
})


test_that(".check_within errors when points fall outside the raster extent", {
    r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    pts <- sf::st_as_sf(
        data.frame(x = c(10, 20), y = c(10, 20)),
        coords = c("x", "y")
    )
    expect_error(blockCV:::.check_within(pts, r), "outside the raster extent")
})


test_that(".check_group_col validates the grouping-column argument", {
    x <- data.frame(site = c("a", "b"))

    expect_error(blockCV:::.check_group_col(NULL, x), "must be provided")
    expect_error(blockCV:::.check_group_col(c("site", "other"), x), "single column name")
})


test_that(".presence_index requires a column for presence-background data", {
    x <- data.frame(occ = c(0, 1, 1))
    expect_error(blockCV:::.presence_index(x, NULL, TRUE), "must be provided")
})


test_that(".check_ext rejects a degenerate raster extent", {
    r <- terra::rast(nrows = 2, ncols = 2)
    terra::ext(r) <- c(0, 0, 0, 1)  # xmin == xmax
    expect_error(blockCV:::.check_ext(r), "Invalid raster extent")
})


test_that(".check_kmeans_k caps clusters at the number of data points", {
    expect_error(blockCV:::.check_kmeans_k(10, 5), "Cannot create")
    expect_true(blockCV:::.check_kmeans_k(3, 5))
})
