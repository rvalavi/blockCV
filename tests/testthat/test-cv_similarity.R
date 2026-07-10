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
    selection = "random",
    iteration = 1,
    biomod2 = FALSE,
    plot = FALSE,
    report = FALSE,
    progress = FALSE
)


test_that("cv_similarity returns a cv_similarity object with MESS", {

    res <- cv_similarity(cv = scv, x = pa_data, r = aus, method = "MESS",
                         plot = FALSE, progress = FALSE)

    expect_s3_class(res, "cv_similarity")
    expect_named(res, c("extrapolation", "overall", "plot", "pbg"))
    expect_true(ggplot2::is_ggplot(res$plot))
})

test_that("cv_similarity is silent by default in non-interactive sessions", {
    skip_if(interactive())

    expect_silent(res <- cv_similarity(cv = scv, x = pa_data, r = aus,
                                       method = "MESS", plot = FALSE))
    expect_s3_class(res, "cv_similarity")
})

test_that("cv_similarity works with L1", {

    res <- cv_similarity(cv = scv, x = pa_data, r = aus, method = "L1",
                         plot = FALSE, progress = FALSE)

    expect_s3_class(res, "cv_similarity")
    expect_true(ggplot2::is_ggplot(res$plot))
})

test_that("cv_similarity works with L2", {

    res <- cv_similarity(cv = scv, x = pa_data, r = aus, method = "L2",
                         plot = FALSE, progress = FALSE)

    expect_s3_class(res, "cv_similarity")
    expect_true(ggplot2::is_ggplot(res$plot))
})

test_that("cv_similarity returns an extrapolation table and a spatial map", {

    res <- cv_similarity(cv = scv, x = pa_data, r = aus, method = "MESS",
                         plot = FALSE, progress = FALSE)
    extr <- res$extrapolation
    expect_s3_class(extr, "data.frame")
    expect_true(all(c("fold", "n_test", "pct_novel", "limiting_var") %in% names(extr)))
    expect_true(is.numeric(res$overall))

    mp <- cv_similarity(cv = scv, x = pa_data, r = aus, method = "MESS", type = "map",
                        plot = FALSE, progress = FALSE)
    expect_s3_class(mp, "cv_similarity")
    expect_true(ggplot2::is_ggplot(mp$plot))
    expect_true(any(vapply(mp$plot$layers, function(l) inherits(l$geom, "GeomSf"), logical(1))))
})


test_that("cv_similarity works with cv_buffer", {
    bloo <- cv_buffer(
        x = pa_data,
        size = 250000,
        progress = FALSE,
        report = FALSE
    )

    res <- cv_similarity(cv = bloo, x = pa_data, r = aus, plot = FALSE, progress = FALSE)

    expect_s3_class(res, "cv_similarity")
    expect_true(ggplot2::is_ggplot(res$plot))
})
