test_that("cv functions are silent by default in non-interactive sessions", {
    skip_if(interactive())

    pts <- sf::st_as_sf(
        data.frame(
            x = rep(seq_len(6), each = 6),
            y = rep(seq_len(6), times = 6),
            occ = rep(0:1, length.out = 36)
        ),
        coords = c("x", "y"),
        crs = 3857
    )

    r <- terra::rast(
        xmin = 0,
        xmax = 7,
        ymin = 0,
        ymax = 7,
        nrows = 7,
        ncols = 7,
        crs = "EPSG:3857"
    )
    terra::values(r) <- seq_len(terra::ncell(r))

    expect_silent(sp <- cv_spatial(
        x = pts,
        rows_cols = c(3, 3),
        k = 3,
        hexagon = FALSE,
        selection = "systematic"
    ))
    expect_s3_class(sp, "cv_spatial")

    expect_silent(cl <- cv_cluster(x = pts, k = 3, biomod2 = FALSE))
    expect_s3_class(cl, "cv_cluster")

    expect_silent(buf <- cv_buffer(x = pts, size = 1.5))
    expect_s3_class(buf, "cv_buffer")

    expect_silent(nn <- cv_nndm(x = pts, r = r, size = 2, num_sample = 20))
    expect_s3_class(nn, "cv_nndm")
})


test_that("report and plot flags suppress output without dropping records", {
    pts <- sf::st_as_sf(
        data.frame(
            x = rep(seq_len(6), each = 6),
            y = rep(seq_len(6), times = 6),
            occ = rep(0:1, length.out = 36)
        ),
        coords = c("x", "y"),
        crs = 3857
    )

    r <- terra::rast(
        xmin = 0,
        xmax = 7,
        ymin = 0,
        ymax = 7,
        nrows = 7,
        ncols = 7,
        crs = "EPSG:3857"
    )
    terra::values(r) <- seq_len(terra::ncell(r))

    expect_silent(sp <- cv_spatial(
        x = pts,
        rows_cols = c(3, 3),
        k = 3,
        hexagon = FALSE,
        selection = "systematic",
        progress = FALSE,
        report = FALSE,
        plot = FALSE
    ))
    expect_false(is.null(sp$records))

    expect_silent(cl <- cv_cluster(
        x = pts,
        k = 3,
        biomod2 = FALSE,
        report = FALSE
    ))
    expect_false(is.null(cl$records))

    expect_silent(buf <- cv_buffer(
        x = pts,
        size = 1.5,
        progress = FALSE,
        report = FALSE
    ))
    expect_false(is.null(buf$records))

    expect_silent(nn <- cv_nndm(
        x = pts,
        r = r,
        size = 2,
        num_sample = 20,
        report = FALSE,
        plot = FALSE
    ))
    expect_false(is.null(nn$records))
})
