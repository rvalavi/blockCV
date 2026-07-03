expect_names <- c(
    "folds_list",
    "folds_ids",
    "biomod_table",
    "k",
    "column",
    "type",
    "space",
    "q",
    "W",
    "Gij",
    "Gj",
    "Gjstar",
    "plot",
    "records"
)


aus <- system.file("extdata/au/", package = "blockCV") |>
    list.files(full.names = TRUE) |>
    terra::rast()

pa_data <- read.csv(system.file("extdata/", "species.csv", package = "blockCV")) |>
    sf::st_as_sf(coords = c("x", "y"), crs = 7845)
pa_data <- pa_data[1:100, ]


test_that("test that cv_knndm works with the default block clustering", {
    set.seed(1)
    knn <- cv_knndm(
        x = pa_data,
        column = "occ",
        r = aus,
        k = 5,
        num_sample = 3000,
        plot = FALSE,
        report = FALSE
    )

    expect_true(exists("knn"))
    expect_s3_class(knn, "cv_knndm")
    expect_equal(names(knn), expect_names)
    expect_equal(knn$k, 5L)
    expect_type(knn$k, "integer")
    expect_equal(length(knn$folds_list), 5L)
    expect_equal(length(knn$folds_ids), nrow(pa_data))
    expect_false(anyNA(knn$folds_ids))
    expect_equal(sort(unique(knn$folds_ids)), 1:5)
    expect_type(knn$W, "double")
    expect_equal(knn$space, "geographical")
    # regular sampling of a non-rectangular raster returns up to num_sample points
    expect_true(length(knn$Gij) > 0 && length(knn$Gij) <= 3000)
    expect_length(knn$Gj, nrow(pa_data))
    expect_length(knn$Gjstar, nrow(pa_data))
    expect_equal(dim(knn$records), c(5L, 4L))
    expect_equal(dim(knn$biomod_table), c(nrow(pa_data), 5L))
    expect_true(!all(knn$records == 0))
    # no fold exceeds maxp
    expect_true(max(table(knn$folds_ids)) <= 0.5 * nrow(pa_data))
})


test_that("test that cv_knndm works with hierarchical and kmeans clustering", {
    set.seed(1)
    kh <- cv_knndm(
        x = pa_data, r = aus, k = 4, clustering = "hierarchical",
        num_sample = 3000, nk_len = 40, plot = FALSE, report = FALSE
    )
    expect_s3_class(kh, "cv_knndm")
    expect_equal(names(kh), expect_names)
    expect_equal(length(kh$folds_list), 4L)
    expect_equal(dim(kh$records), c(4L, 2L))

    set.seed(1)
    km <- cv_knndm(
        x = pa_data, r = aus, k = 4, clustering = "kmeans",
        num_sample = 3000, nk_len = 40, plot = FALSE, report = FALSE
    )
    expect_s3_class(km, "cv_knndm")
    expect_equal(sort(unique(km$folds_ids)), 1:4)

    expect_equal(print(kh), "cv_knndm")
    expect_output(summary(kh))
})


test_that("test that cv_knndm works in feature space", {
    set.seed(1)
    kf <- cv_knndm(
        x = pa_data, r = aus, k = 4, space = "feature",
        clustering = "hierarchical", num_sample = 2000, nk_len = 40,
        plot = FALSE, report = FALSE
    )
    expect_s3_class(kf, "cv_knndm")
    expect_equal(kf$space, "feature")
    expect_equal(sort(unique(kf$folds_ids)), 1:4)

    # blocks is coerced to hierarchical in feature space
    expect_message(
        cv_knndm(x = pa_data, r = aus, k = 4, space = "feature",
                 clustering = "blocks", num_sample = 2000, nk_len = 20,
                 plot = FALSE, report = FALSE)
    )
})


test_that("test that cv_knndm accepts explicit prediction points", {
    set.seed(1)
    pp <- sf::st_as_sf(terra::spatSample(aus[[1]], size = 1000, method = "regular",
                                         na.rm = TRUE, as.points = TRUE))
    knn <- cv_knndm(
        x = pa_data, predpoints = pp, k = 5, clustering = "hierarchical",
        nk_len = 30, plot = FALSE, report = FALSE
    )
    expect_s3_class(knn, "cv_knndm")
    expect_length(knn$Gij, nrow(pp))
})


test_that("test that cv_knndm is reproducible with a seed", {
    a <- cv_knndm(x = pa_data, r = aus, k = 5, num_sample = 2000, seed = 42, plot = FALSE, report = FALSE)
    b <- cv_knndm(x = pa_data, r = aus, k = 5, num_sample = 2000, seed = 42, plot = FALSE, report = FALSE)
    expect_identical(a$folds_ids, b$folds_ids)
})


test_that("test cv_knndm errors with invalid arguments", {
    # invalid x
    expect_error(cv_knndm(x = "pa_data", r = aus, k = 5))
    # k too small
    expect_error(cv_knndm(x = pa_data, r = aus, k = 1, num_sample = 2000))
    # maxp outside (1/k, 1)
    expect_error(cv_knndm(x = pa_data, r = aus, k = 5, maxp = 2, num_sample = 2000))
    expect_error(cv_knndm(x = pa_data, r = aus, k = 5, maxp = 0.1, num_sample = 2000))
    # no prediction source
    expect_error(cv_knndm(x = pa_data, k = 5))
})


test_that("test cv_knndm errors with mismatched or missing CRS", {
    nocrs <- pa_data
    sf::st_crs(nocrs) <- NA
    expect_error(cv_knndm(x = nocrs, column = "occ", r = aus, k = 5))

    reproj <- sf::st_transform(pa_data, 4326)
    expect_error(cv_knndm(x = reproj, r = aus, k = 5, num_sample = 2000))
})
