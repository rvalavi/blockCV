expect_names <- c(
    "folds_list",
    "folds_ids",
    "biomod_table",
    "k",
    "column",
    "presence_bg",
    "type",
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


test_that("test that environmental cluster function with raster_cluster", {
    # environmental clustering
    set.seed(42)
    eb <- cv_cluster(
        x = pa_data,
        column = "occ",
        r = aus,
        k = 3,
        scale = TRUE,
        raster_cluster = TRUE
    )

    expect_true(exists("eb"))
    expect_s3_class(eb, "cv_cluster")
    expect_equal(names(eb), expect_names)
    expect_equal(length(eb$folds_list), 3)
    expect_type(eb$folds_list, "list")
    expect_equal(eb$k, 3)
    expect_type(eb$column, "character")
    expect_equal(dim(eb$records), c(3, 4))
    expect_true(!all(eb$records == 0))

})

test_that("test that spacial cluster function with no column", {
    # spatial clustering
    set.seed(42)
    eb <- cv_cluster(x = pa_data, k = 5, biomod2 = FALSE)

    expect_true(exists("eb"))
    expect_s3_class(eb, "cv_cluster")
    expect_equal(names(eb), expect_names)
    expect_equal(length(eb$folds_list), 5)
    expect_type(eb$folds_list, "list")
    expect_equal(eb$k, 5)
    expect_null(eb$column)
    expect_equal(dim(eb$records), c(5, 2))
    expect_true(!all(eb$records == 0))

    expect_equal(print(eb), "cv_cluster")
    expect_output(summary(eb))

})

test_that("continuous column is reported with quantile bins", {

    cont_data <- pa_data
    set.seed(201)
    cont_data$biomass <- stats::rnorm(nrow(cont_data))

    set.seed(202)
    eb <- cv_cluster(
        x = cont_data,
        column = "biomass",
        k = 5,
        biomod2 = FALSE,
        num_bins = 4
    )

    expect_equal(dim(eb$records), c(5, 8))
    expect_equal(
        names(eb$records),
        c(paste0("train_Q", 1:4), paste0("test_Q", 1:4))
    )
    bins <- attr(eb$records, "column_bins")
    expect_equal(nrow(bins), 4)
    expect_equal(attr(bins, "requested_bins"), 4L)

})

test_that("continuous quantile bins tolerate tied breaks", {

    cont_data <- pa_data
    n <- nrow(cont_data)
    n0 <- floor(n * 0.4)
    n1 <- floor(n * 0.4)
    cont_data$biomass <- c(
        rep(0, n0),
        rep(1, n1),
        seq(2, 3, length.out = n - n0 - n1)
    )

    set.seed(203)
    eb <- cv_cluster(
        x = cont_data,
        column = "biomass",
        k = 5,
        biomod2 = FALSE,
        num_bins = 4
    )

    bins <- attr(eb$records, "column_bins")
    expect_equal(attr(bins, "requested_bins"), 4L)
    expect_lt(nrow(bins), 4)
    expect_equal(ncol(eb$records), nrow(bins) * 2)

})


test_that("test that environmental cluster with no scale and wrong column", {
    set.seed(42)
    expect_warning(
        eb <- cv_cluster(
            x = pa_data,
            column = "response",
            # wrong column name
            r = aus,
            k = 5,
            scale = FALSE,
            raster_cluster = TRUE,
            algorithm = "MacQueen"
        )
    )

    expect_true(exists("eb"))
    expect_s3_class(eb, "cv_cluster")
    expect_equal(names(eb), expect_names)
    expect_equal(length(eb$folds_list), 5)
    expect_type(eb$folds_list, "list")
    expect_equal(eb$k, 5)
    expect_equal(dim(eb$records), c(5, 2))
    expect_true(!all(eb$records == 0))

})


test_that("test environmental cluster with no spatial or sf object", {
    expect_error(cv_cluster(
        x = aus,
        # no spatial or sf object
        r = aus,
        k = 5,
        scale = TRUE,
        raster_cluster = FALSE
    ))

})

test_that("test environmental cluster with no raster data", {
    expect_error(cv_cluster(
        r = data.frame(x = 1),
        # no raster data
        x = pa_data
    ))

})


test_that("categorical (factor) raster layers are rejected", {
    r_cat <- aus[[1]]
    terra::values(r_cat) <- sample(1:3, terra::ncell(r_cat), replace = TRUE)
    levels(r_cat) <- data.frame(id = 1:3, cover = c("a", "b", "c"))
    names(r_cat) <- "landcover"

    expect_error(
        cv_cluster(x = pa_data, r = c(aus, r_cat), k = 3),
        "numeric covariates"
    )
})


test_that("spatial_weight blends geography into environmental clustering", {
    set.seed(42)
    ec <- cv_cluster(x = pa_data, column = "occ", r = aus, k = 5,
                     scale = TRUE, biomod2 = FALSE, spatial_weight = 0.4)

    expect_s3_class(ec, "cv_cluster")
    expect_equal(ec$type, "Environmental Cluster")
    expect_equal(length(ec$folds_list), 5)

    # spatial_weight = 0 reproduces the default (pure environmental) folds exactly
    set.seed(42)
    w0  <- cv_cluster(x = pa_data, column = "occ", r = aus, k = 5,
                      scale = TRUE, biomod2 = FALSE, spatial_weight = 0)
    set.seed(42)
    def <- cv_cluster(x = pa_data, column = "occ", r = aus, k = 5,
                      scale = TRUE, biomod2 = FALSE)
    expect_identical(w0$folds_ids, def$folds_ids)

    # spatial_weight = 1 runs and clusters purely on the (standardised) coordinates
    set.seed(42)
    w1 <- cv_cluster(x = pa_data, column = "occ", r = aus, k = 5,
                     scale = TRUE, biomod2 = FALSE, spatial_weight = 1)
    expect_equal(length(w1$folds_list), 5)

    # higher spatial_weight yields more spatially compact folds (0 > 0.4 >= 1)
    xy <- scale(sf::st_coordinates(pa_data))
    spread <- function(cv){
        mean(sapply(split(as.data.frame(xy), cv$folds_ids), function(d)
            if (nrow(d) < 2) 0 else mean(dist(as.matrix(d)))))
    }
    expect_lt(spread(ec), spread(w0))
    expect_lte(spread(w1), spread(ec))
})


test_that("spatial_weight works with raster_cluster = TRUE", {
    set.seed(42)
    w0 <- cv_cluster(x = pa_data, column = "occ", r = aus, k = 5, scale = TRUE,
                     raster_cluster = TRUE, biomod2 = FALSE, spatial_weight = 0)
    set.seed(42)
    def <- cv_cluster(x = pa_data, column = "occ", r = aus, k = 5, scale = TRUE,
                      raster_cluster = TRUE, biomod2 = FALSE)
    # spatial_weight = 0 preserves the original raster_cluster path exactly
    expect_identical(w0$folds_ids, def$folds_ids)

    set.seed(42)
    ec <- cv_cluster(x = pa_data, column = "occ", r = aus, k = 5, scale = TRUE,
                     raster_cluster = TRUE, biomod2 = FALSE, spatial_weight = 0.5)
    expect_s3_class(ec, "cv_cluster")
    expect_equal(length(ec$folds_list), 5)
    expect_equal(length(ec$folds_ids), nrow(pa_data))
})


test_that("spatial_weight with scale = FALSE does not scale the covariates", {
    # scale = FALSE must keep covariates on native units; only coordinates are scaled
    set.seed(42)
    ec <- cv_cluster(x = pa_data, column = "occ", r = aus, k = 5,
                     scale = FALSE, biomod2 = FALSE, spatial_weight = 0.4)
    expect_s3_class(ec, "cv_cluster")
    expect_equal(length(ec$folds_list), 5)
})


test_that("spatial_weight is validated only when r is supplied, else ignored", {
    # values outside [0, 1] error when a raster is provided (both ends of the range)
    expect_error(
        cv_cluster(x = pa_data, r = aus, k = 3, spatial_weight = 1.5),
        "between 0 and 1"
    )
    expect_error(
        cv_cluster(x = pa_data, r = aus, k = 3, spatial_weight = -0.5),
        "between 0 and 1"
    )
    # non-scalar / non-numeric weights are also rejected
    expect_error(
        cv_cluster(x = pa_data, r = aus, k = 3, spatial_weight = c(0.2, 0.5)),
        "between 0 and 1"
    )
    # without a raster the weight is ignored (not validated) but warns when non-zero
    expect_warning(
        cv_cluster(x = pa_data, k = 5, biomod2 = FALSE, spatial_weight = 0.4),
        "only used for environmental clustering"
    )
    # an out-of-range weight is not validated (so does not error) when r is absent
    set.seed(42)
    expect_warning(
        eb <- cv_cluster(x = pa_data, k = 5, biomod2 = FALSE, spatial_weight = 1.5),
        "only used for environmental clustering"
    )
    expect_equal(length(eb$folds_list), 5)
})
