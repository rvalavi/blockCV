test_that("presence-background validation requires numeric 0/1 values", {
    x <- data.frame(occ = c(0, 1, 1, 0))

    expect_equal(.presence_index(x, "occ", TRUE), c(2L, 3L))
    expect_equal(.presence_index(x, "occ", FALSE), seq_len(nrow(x)))

    fractional <- x
    fractional$occ[2] <- 0.5
    expect_error(
        .presence_index(fractional, "occ", TRUE),
        "0s.*1s"
    )

    missing <- x
    missing$occ[2] <- NA_real_
    expect_error(
        .presence_index(missing, "occ", TRUE),
        "missing values"
    )

    text <- x
    text$occ <- as.character(text$occ)
    expect_error(
        .presence_index(text, "occ", TRUE),
        "numeric"
    )
})


test_that("public presence-background constructors reject invalid responses", {
    x <- data.frame(
        x = c(0, 1000, 2000, 0, 1000, 2000),
        y = c(0, 0, 0, 1000, 1000, 1000),
        occ = c(1, 0, 1, 0, 1, 0)
    )
    x <- sf::st_as_sf(x, coords = c("x", "y"), crs = 3857)
    model_domain <- sf::st_as_sf(sf::st_as_sfc(sf::st_bbox(x)))

    fractional <- x
    fractional$occ[2] <- 0.5
    expect_error(
        cv_buffer(
            x = fractional,
            column = "occ",
            size = 1000,
            presence_bg = TRUE,
            progress = FALSE
        ),
        "0s.*1s"
    )
    expect_error(
        cv_nndm(
            x = fractional,
            column = "occ",
            model_domain = model_domain,
            size = 1000,
            num_sample = 20,
            presence_bg = TRUE,
            plot = FALSE,
            report = FALSE
        ),
        "0s.*1s"
    )

    missing <- x
    missing$occ[2] <- NA_real_
    expect_error(
        cv_buffer(
            x = missing,
            column = "occ",
            size = 1000,
            presence_bg = TRUE,
            progress = FALSE
        ),
        "missing values"
    )
    expect_error(
        cv_nndm(
            x = missing,
            column = "occ",
            model_domain = model_domain,
            size = 1000,
            num_sample = 20,
            presence_bg = TRUE,
            plot = FALSE,
            report = FALSE
        ),
        "missing values"
    )
})
