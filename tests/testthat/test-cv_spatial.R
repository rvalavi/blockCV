expect_names <- c("folds_list",
                  "folds_ids",
                  "biomod_table",
                  "k",
                  "size",
                  "column",
                  "presence_bg",
                  "blocks",
                  "records")

aus <- system.file("extdata/au/", package = "blockCV") |>
    list.files(full.names = TRUE) |>
    terra::rast()

pa_data <- read.csv(system.file("extdata/", "species.csv", package = "blockCV")) |>
    sf::st_as_sf(coords = c("x", "y"), crs = 7845)
pa_data <- pa_data[1:200, ]

test_that("test cv_spatial function with random assingment and raster file", {

    set.seed(1000)
    scv <- cv_spatial(
        x = pa_data,
        column = "occ",
        r = aus,
        size = 450000,
        k = 5,
        selection = "random",
        iteration = 5,
        biomod2 = TRUE,
        offset = c(0.3, 0.2),
        plot = FALSE,
        progress = TRUE
    )

    expect_true(exists("scv"))
    expect_s3_class(scv, "cv_spatial")
    expect_equal(names(scv), expect_names)
    expect_equal(length(scv$folds_list), 5)
    expect_type(scv$folds_list, "list")
    expect_type(scv$biomod_table, "logical")
    expect_equal(dim(scv$biomod_table), c(nrow(pa_data), 5))
    expect_equal(scv$k, 5)
    expect_s3_class(scv$blocks, "sf")
    expect_type(scv$column, "character")
    expect_type(scv$size, "double")
    expect_equal(dim(scv$records), c(5, 4))
    expect_true(
        !all(scv$records == 0)
    )

})

test_that("continuous column is balanced with quantile bins", {

    cont_data <- pa_data
    set.seed(101)
    cont_data$biomass <- stats::rnorm(nrow(cont_data))

    set.seed(102)
    scv <- cv_spatial(
        x = cont_data,
        column = "biomass",
        rows_cols = c(20, 20),
        hexagon = FALSE,
        k = 4,
        selection = "random",
        iteration = 20,
        biomod2 = FALSE,
        progress = FALSE,
        plot = FALSE
    )

    expect_equal(dim(scv$records), c(4, 8))
    expect_equal(
        names(scv$records),
        c(paste0("train_Q", 1:4), paste0("test_Q", 1:4))
    )
    bins <- attr(scv$records, "column_bins")
    expect_equal(nrow(bins), 4)
    expect_equal(attr(bins, "requested_bins"), 4L)
    expect_true(all(bins$lower <= bins$upper))

})

test_that("numeric multi-class column is not quantile-binned", {

    class_data <- pa_data
    class_data$class5 <- rep(1:5, length.out = nrow(class_data))

    scv <- cv_spatial(
        x = class_data,
        column = "class5",
        rows_cols = c(20, 20),
        hexagon = FALSE,
        k = 5,
        selection = "random",
        iteration = 5,
        biomod2 = FALSE,
        progress = FALSE,
        plot = FALSE
    )

    expect_equal(dim(scv$records), c(5, 10))
    expect_null(attr(scv$records, "column_bins"))

})

test_that("num_bins = NULL disables quantile binning for a continuous column", {

    # deterministic, non-integer, > 15 unique values (would be binned by default)
    cont_data <- pa_data
    cont_data$biomass <- rep(seq(0.5, 20.5, by = 1), length.out = nrow(cont_data))
    n_unique <- length(unique(cont_data$biomass))

    set.seed(104)
    scv <- suppressWarnings(
        cv_spatial(
            x = cont_data,
            column = "biomass",
            rows_cols = c(20, 20),
            hexagon = FALSE,
            k = 4,
            selection = "random",
            iteration = 5,
            biomod2 = FALSE,
            num_bins = NULL,
            progress = FALSE,
            plot = FALSE
        )
    )

    # every unique value becomes its own class, so no quantile bins are stored
    expect_null(attr(scv$records, "column_bins"))
    expect_equal(ncol(scv$records), n_unique * 2)

})


test_that("test cv_spatial function with systematic assingment and no raster file", {

    scv <- cv_spatial(
        x = sf::as_Spatial(pa_data),
        rows_cols = c(10, 10),
        k = 5,
        selection = "systematic",
        biomod2 = FALSE,
        plot = TRUE
    )

    expect_true(exists("scv"))
    expect_s3_class(scv, "cv_spatial")
    expect_equal(names(scv), expect_names)
    expect_equal(length(scv$folds_list), 5)
    expect_type(scv$folds_list, "list")
    expect_null(scv$biomod_table, "matrix")
    expect_equal(scv$k, 5)
    expect_s3_class(scv$blocks, "sf")
    expect_null(scv$column)
    expect_null(scv$size)
    expect_equal(dim(scv$records), c(5, 2))
    expect_true(
        !all(scv$records == 0)
    )

    expect_equal(print(scv), "cv_spatial")
    expect_message(plot(scv))
    expect_output(summary(scv))

})


test_that("test cv_spatial function with non-numeric iteration", {

    scv <- cv_spatial(
        x = pa_data,
        hexagon = FALSE,
        k = 5,
        selection = "random",
        iteration = "2", # not numeric
        biomod2 = FALSE,
        plot = FALSE
    )

    expect_true(exists("scv"))
    expect_s3_class(scv, "cv_spatial")
    expect_equal(names(scv), expect_names)
    expect_equal(length(scv$folds_list), 5)
    expect_type(scv$folds_list, "list")
    expect_null(scv$biomod_table, "matrix")
    expect_equal(scv$k, 5)
    expect_s3_class(scv$blocks, "sf")
    expect_null(scv$column)
    expect_null(scv$size)
    expect_equal(dim(scv$records), c(5, 2))
    expect_true(
        !all(scv$records == 0)
    )

    expect_equal(print(scv), "cv_spatial")
    expect_message(plot(scv))
    expect_output(summary(scv))

})

test_that("test cv_spatial with checkerboard assingment and only row blocks", {

    scv <- cv_spatial(
        x = sf::as_Spatial(pa_data),
        hexagon = FALSE,
        rows_cols = c(5, 1),
        selection = "checkerboard",
        biomod2 = FALSE,
        plot = TRUE
    )

    expect_true(exists("scv"))
    expect_s3_class(scv, "cv_spatial")
    expect_equal(names(scv), expect_names)
    expect_equal(length(scv$folds_list), 2)
    expect_type(scv$folds_list, "list")
    expect_null(scv$biomod_table)
    expect_equal(scv$k, 2)
    expect_s3_class(scv$blocks, "sf")
    expect_null(scv$column)
    expect_null(scv$size)
    expect_equal(dim(scv$records), c(2, 2))
    expect_true(
        !all(scv$records == 0)
    )

    expect_equal(print(scv), "cv_spatial")
    expect_message(plot(scv))
    expect_output(summary(scv))

})


test_that("test cv_spatial with user-defined blocks", {

    # make a spatial block polygon
    user_poly <- .make_blocks(x_obj = pa_data, blocksize = 450000)

    scv <- cv_spatial(
        x = pa_data,
        user_blocks = user_poly,
        selection = "random",
        iteration = 5,
        biomod2 = FALSE,
        progress = FALSE,
        plot = FALSE
    )

    expect_true(exists("scv"))
    expect_s3_class(scv, "cv_spatial")
    expect_equal(names(scv), expect_names)
    expect_equal(length(scv$folds_list), 5)
    expect_type(scv$folds_list, "list")
    expect_null(scv$biomod_table)
    expect_equal(scv$k, 5)
    expect_s3_class(scv$blocks, "sf")
    expect_null(scv$column)
    expect_null(scv$size)
    expect_equal(dim(scv$records), c(5, 2))
    expect_true(
        !all(scv$records == 0)
    )

    expect_equal(print(scv), "cv_spatial")
    expect_message(plot(scv))
    expect_output(summary(scv))

})


test_that("test cv_spatial failur: number of blocks, wrong selection", {

    expect_error(cv_spatial(
        x = pa_data,
        rows_cols = c(2, 2),
        hexagon = FALSE,
        k = 15,
        # very high k
        selection = "random"
    ))

    expect_error(cv_spatial(
        x = pa_data,
        k = 1,
        # very low k
        selection = "random"
    ))

    expect_error(cv_spatial(x = pa_data, k = 5, selection = "chance")) # wrong selection

})


test_that("test cv_spatial failur: wrong user-defined blocks", {

    expect_error(
        cv_spatial(
            x = pa_data,
            # wrong user-defined blocks
            column = "occ",
            r = aus,
            user_blocks = data.frame(a = 1),
            size = 450000,
            k = 5
        )
    )

})


test_that("test cv_spatial with no speceis column match", {


    expect_warning(scv <- cv_spatial(
        x = pa_data,
        column = "response",
        # wrong response
        r = aus,
        size = 450000,
        k = 5
    ))

    expect_true(exists("scv"))
    expect_s3_class(scv, "cv_spatial")
    expect_equal(names(scv), expect_names)
    expect_equal(length(scv$folds_list), 5)
    expect_equal(scv$k, 5)
    expect_equal(dim(scv$records), c(5, 2))
    expect_true(
        !all(scv$records == 0)
    )

})
