expect_names <- c(
    "folds_list",
    "k",
    "column",
    "size",
    "exclusion",
    "plot",
    "presence_bg",
    "records"
)


aus <- terra::rast(
    list.files(system.file("extdata/au/", package = "blockCV"), full.names = TRUE)
)

pa_data <- sf::st_as_sf(
    read.csv(system.file("extdata/", "species.csv", package = "blockCV")),
    coords = c("x", "y"),
    crs = 7845
)
pa_data <- pa_data[1:100, ]


test_that("test that cv_nndm function works properly with presence-absence data", {
    bloo <- cv_nndm(
        x = pa_data,
        column = "occ",
        r = aus,
        size = 250000,
        num_sample = 3000,
        presence_bg = FALSE
    )

    expect_true(exists("bloo"))
    expect_s3_class(bloo, "cv_nndm")
    expect_equal(names(bloo), expect_names)
    expect_equal(length(bloo$folds_list), nrow(pa_data))
    expect_type(bloo$folds_list, "list")
    expect_type(bloo$k, "integer")
    expect_type(bloo$column, "character")
    expect_type(bloo$size, "double")
    expect_equal(bloo$presence_bg, FALSE)
    expect_equal(dim(bloo$records), c(nrow(pa_data), 4))
    expect_true(!all(bloo$records == 0))

    # the exclusion radius reported for each fold must reproduce that fold's split
    expect_s3_class(bloo$exclusion, "data.frame")
    expect_equal(names(bloo$exclusion), c("fold", "test_id", "exclusion_distance"))
    expect_equal(nrow(bloo$exclusion), bloo$k)
    expect_true(all(bloo$exclusion$exclusion_distance > 0))
    dmat <- sf::st_distance(pa_data)
    units(dmat) <- NULL
    for(i in seq_len(bloo$k)){
        expect_equal(
            as.numeric(which(dmat[bloo$exclusion$test_id[i], ] > bloo$exclusion$exclusion_distance[i])),
            sort(bloo$folds_list[[i]][[1]])
        )
        expect_true(bloo$exclusion$test_id[i] %in% bloo$folds_list[[i]][[2]])
    }

})


test_that("test that cv_nndm function works properly with presence-background data", {
    bloo <- cv_nndm(
        x = pa_data,
        column = "occ",
        r = aus,
        size = 250000,
        sampling = "regular",
        presence_bg = TRUE
    )

    expect_true(exists("bloo"))
    expect_s3_class(bloo, "cv_nndm")
    expect_equal(names(bloo), expect_names)
    expect_equal(length(bloo$folds_list), sum(pa_data$occ))
    expect_type(bloo$folds_list, "list")
    expect_type(bloo$k, "integer")
    expect_type(bloo$column, "character")
    expect_type(bloo$size, "double")
    expect_equal(bloo$presence_bg, TRUE)
    expect_equal(dim(bloo$records), c(sum(pa_data$occ), 4))
    expect_true(!all(bloo$records == 0))

    # in presence-background, every fold's exclusion is keyed to a presence record
    expect_equal(nrow(bloo$exclusion), bloo$k)
    expect_true(all(pa_data$occ[bloo$exclusion$test_id] == 1))

})


test_that("test that cv_nndm function works properly with no species specified", {
    bloo <- cv_nndm(
        x = sf::as_Spatial(pa_data),
        r = aus,
        size = 250000,
        num_sample = 3000,
        plot = FALSE
    )

    expect_true(exists("bloo"))
    expect_s3_class(bloo, "cv_nndm")
    expect_equal(names(bloo), expect_names)
    expect_equal(length(bloo$folds_list), nrow(pa_data))
    expect_type(bloo$folds_list, "list")
    expect_type(bloo$k, "integer")
    expect_null(bloo$species)
    expect_type(bloo$size, "double")
    expect_equal(bloo$presence_bg, FALSE)
    expect_equal(dim(bloo$records), c(nrow(pa_data), 2))
    expect_true(!all(bloo$records == 0))

    expect_output(print(bloo), "blockCV cv_nndm")
    expect_output(summary(bloo))

})


test_that("test cv_nndm function with no matching species column", {
    expect_warning(
        bloo <- cv_nndm(
            x = pa_data,
            column = "response",
            r = aus,
            size = 250000,
            num_sample = 3000,
            presence_bg = FALSE
        )
    )

    expect_true(exists("bloo"))
    expect_equal(dim(bloo$records), c(nrow(pa_data), 2))
    expect_true(!all(bloo$records == 0))

})


test_that("test cv_nndm function with no spatial column data", {
    expect_error(
        cv_nndm(
            x = "pa_data",
            r = aus,
            size = 250000
        )
    )

})

test_that("cv_nndm falls back to Euclidean distance when x has no CRS", {
    nocrs <- pa_data
    sf::st_crs(nocrs) <- NA
    dom <- sf::st_as_sfc(sf::st_bbox(nocrs))

    expect_warning(
        bloo <- cv_nndm(
            x = nocrs,
            column = "occ",
            model_domain = dom,
            size = 250000,
            num_sample = 500,
            plot = FALSE,
            report = FALSE
        ),
        "CRS undefined"
    )

    expect_s3_class(bloo, "cv_nndm")
    expect_equal(length(bloo$folds_list), nrow(nocrs))

})

test_that("cv_nndm errors when x and r have mismatched CRS", {
    nocrs <- pa_data
    sf::st_crs(nocrs) <- NA

    expect_error(
        suppressWarnings(
            cv_nndm(
                x = nocrs,
                column = "occ",
                r = aus,
                size = 250000
            )
        )
    )

})
