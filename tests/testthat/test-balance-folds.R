# the fold-balancing objective shared by cv_spatial (random selection) and
# cv_cluster(balance = TRUE); see .balance_score / .balance_folds in utils.R

pa <- read.csv(system.file("extdata", "species.csv", package = "blockCV"))
pa_sf <- sf::st_as_sf(pa, coords = c("x", "y"), crs = 7845)

# build a records table for a given fold vector (helper for the unit tests)
records_for <- function(folds, response, k){
    tt <- .records_table(k, response)
    for(p in seq_len(k)){
        test_set <- which(folds == p)
        train_set <- which(folds != p)
        tt <- .records_table_row(tt, p, train_set, test_set, response)
    }
    tt
}


test_that(".balance_score works for the single-column presence-background case", {
    # a length-1 test_cols must not trip the diag()-of-a-scalar identity gotcha
    resp <- .column_response(pa, "occ")
    tt <- records_for((seq_len(nrow(pa)) %% 5L) + 1L, resp, k = 5)

    score <- .balance_score(tt, "test_1", k = 5)
    expect_length(score, 2)
    expect_true(all(is.finite(score)))

    # a perfectly even split scores zero on both components
    even <- .records_table(4, .column_response(data.frame(occ = rep(0:1, each = 8)), "occ"))
    even[, c("test_0", "test_1")] <- 2   # 2 of each class in every one of 4 folds
    zero <- .balance_score(even, c("test_0", "test_1"), k = 4)
    expect_equal(zero, c(0, 0))
})


test_that("the objective beats the old min/sd rule on the packaged data", {
    resp <- .column_response(pa, "occ")
    k <- 5
    n <- nrow(pa)
    test_cols <- grep("^test", names(.records_table(k, resp)), value = TRUE)

    # a shared pool of candidate assignments scored by both rules
    set.seed(123)
    tables <- lapply(seq_len(200), function(i){
        records_for(sample(rep_len(seq_len(k), n)), resp, k)
    })

    # old rule: ratchet up the minimum cell, ratchet down the sd of *all* cells,
    # requiring both to improve at once (the behaviour being replaced)
    old_pick <- NA_integer_; min_num <- 0; max_sd <- Inf
    for(i in seq_along(tables)){
        obj <- tables[[i]]
        if(min(obj) >= min_num && stats::sd(unlist(obj)) < max_sd){
            old_pick <- i; min_num <- min(obj); max_sd <- stats::sd(unlist(obj))
        }
    }

    # new rule: lexicographic argmin on (eligible empty cells, imbalance)
    scores <- vapply(tables, .balance_score, numeric(2), test_cols = test_cols, k = k)
    new_pick <- order(scores[1, ], scores[2, ])[1]

    imbalance <- function(tt) .balance_score(tt, test_cols, k)[2]
    tvar <- function(tt) sum(apply(as.matrix(tt[, test_cols, drop = FALSE]), 2, stats::var))

    # (b) the chosen split is at least as balanced as the old rule's pick
    expect_lte(imbalance(tables[[new_pick]]), imbalance(tables[[old_pick]]))
    expect_lte(tvar(tables[[new_pick]]), tvar(tables[[old_pick]]))
})


test_that("cv_cluster(balance = TRUE) minimises eligible empty test-class cells", {
    set.seed(42)
    cl <- cv_cluster(
        x = pa_sf,
        column = "occ",
        k = 5,
        balance = TRUE,
        report = FALSE,
        progress = FALSE
    )

    rec <- cl$records
    test_cols <- grep("^test", names(rec), value = TRUE)
    tt <- as.matrix(rec[, test_cols, drop = FALSE])
    n_c <- colSums(tt)

    # every class with at least k records appears in every test fold for this data
    fillable <- tt[, n_c >= 5, drop = FALSE]
    expect_true(all(fillable > 0))
})
