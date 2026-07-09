# cv_nndm must not leak points between train and test, and in the
# presence-background case the target test point must be a presence.
# Regression guard for the x_1s[i] vs bare-i index bug in the fold loop.

pa <- read.csv(system.file("extdata", "species.csv", package = "blockCV"))
pa_sf <- sf::st_as_sf(pa, coords = c("x", "y"), crs = 7845)
covar <- terra::rast(system.file("extdata/au/bio_5.tif", package = "blockCV"))

# every fold's train and test indices must be disjoint
expect_no_leakage <- function(cv){
    for(f in cv$folds_list){
        expect_length(intersect(f[[1]], f[[2]]), 0)
    }
}


test_that("cv_nndm presence-background folds have no train/test leakage", {
    skip_on_cran()
    nndm <- cv_nndm(
        x = pa_sf, column = "occ", r = covar,
        size = 350000, num_sample = 2000, sampling = "regular",
        min_train = 0.1, presence_bg = TRUE,
        plot = FALSE, report = FALSE
    )

    expect_no_leakage(nndm)

    # the target test point (first test index) is the presence being held out
    x_1s <- which(pa$occ == 1)
    for(p in seq_along(nndm$folds_list)){
        target <- nndm$folds_list[[p]][[2]][1]
        expect_equal(target, x_1s[p])
        expect_equal(pa$occ[target], 1)
    }
})


test_that("cv_nndm presence-background with add_bg has no leakage", {
    skip_on_cran()
    nndm <- cv_nndm(
        x = pa_sf, column = "occ", r = covar,
        size = 350000, num_sample = 2000, sampling = "regular",
        min_train = 0.1, presence_bg = TRUE, add_bg = TRUE,
        plot = FALSE, report = FALSE
    )

    expect_no_leakage(nndm)

    # first test index is still the target presence; any extras are backgrounds
    x_1s <- which(pa$occ == 1)
    for(p in seq_along(nndm$folds_list)){
        test_set <- nndm$folds_list[[p]][[2]]
        expect_equal(test_set[1], x_1s[p])
        if(length(test_set) > 1){
            expect_true(all(pa$occ[test_set[-1]] == 0))
        }
    }
})


test_that("cv_nndm non-presence-background folds have no leakage", {
    skip_on_cran()
    nndm <- cv_nndm(
        x = pa_sf, column = "occ", r = covar,
        size = 350000, num_sample = 2000, sampling = "regular",
        min_train = 0.1, presence_bg = FALSE,
        plot = FALSE, report = FALSE
    )

    expect_no_leakage(nndm)
})
