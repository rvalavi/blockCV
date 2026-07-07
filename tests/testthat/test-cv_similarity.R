aus <- system.file("extdata/au/", package = "blockCV") |>
    list.files(full.names = TRUE) |>
    terra::rast()

pa_data <- read.csv(system.file("extdata/", "species.csv", package = "blockCV")) |>
    sf::st_as_sf(coords = c("x", "y"), crs = 7845)
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


test_that("test that cv_similarity function works with cv_spatil and MESS", {

    plt <- cv_similarity(cv = scv, x = pa_data, r = aus, method = "MESS")

    expect_true(exists("plt"))
    expect_true(ggplot2::is_ggplot(plt))

})

test_that("test that cv_similarity function works with cv_spatil and L1", {

    plt <- cv_similarity(cv = scv, x = pa_data, r = aus, method = "L1")

    expect_true(exists("plt"))
    expect_true(ggplot2::is_ggplot(plt))

})

test_that("test that cv_similarity function works with cv_spatil and L2", {

    plt <- cv_similarity(cv = scv, x = pa_data, r = aus, , method = "L2")

    expect_true(exists("plt"))
    expect_true(ggplot2::is_ggplot(plt))

})

test_that("test that cv_similarity returns an extrapolation table and a spatial map", {

    plt <- cv_similarity(cv = scv, x = pa_data, r = aus, method = "MESS")
    extr <- attr(plt, "extrapolation")
    expect_s3_class(extr, "data.frame")
    expect_true(all(c("fold", "n_test", "pct_novel", "limiting_var") %in% names(extr)))

    mp <- cv_similarity(cv = scv, x = pa_data, r = aus, method = "MESS", type = "map")
    expect_true(ggplot2::is_ggplot(mp))
    expect_true(any(vapply(mp$layers, function(l) inherits(l$geom, "GeomSf"), logical(1))))

})


test_that("test that cv_similarity function works with cv_buffer", {
    bloo <- cv_buffer(
        x = pa_data,
        size = 250000,
        progress = FALSE,
        report = FALSE
    )

    plt <- cv_similarity(cv = bloo, x = pa_data, r = aus)

    expect_true(exists("plt"))
    expect_true(ggplot2::is_ggplot(plt))

})

