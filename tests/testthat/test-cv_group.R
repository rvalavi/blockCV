group_test_data <- function(){
    sizes <- c(2L, 3L, 4L, 2L, 5L, 3L)
    groups <- rep(sprintf("site_%02d", seq_along(sizes)), sizes)
    dat <- data.frame(
        site = groups,
        occ = rep(c(0L, 1L), length.out = length(groups)),
        x = seq_along(groups),
        y = (seq_along(groups) %% 5L) + rep(seq_along(sizes), sizes),
        stringsAsFactors = FALSE
    )
    sf::st_as_sf(dat, coords = c("x", "y"), crs = 7845)
}

balanced_group_test_data <- function(){
    group_ids <- sprintf("site_%02d", 1:9)
    group_occ <- c(0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L)
    groups <- rep(group_ids, each = 3L)
    dat <- data.frame(
        site = groups,
        occ = rep(group_occ, each = 3L),
        x = seq_along(groups),
        y = rep(seq_along(group_ids), each = 3L),
        stringsAsFactors = FALSE
    )
    sf::st_as_sf(dat, coords = c("x", "y"), crs = 7845)
}

expect_group_integrity <- function(cv, x, group_col = "site"){
    groups <- as.character(x[[group_col]])

    row_tests <- unlist(lapply(cv$folds_list, `[[`, 2L), use.names = FALSE)
    expect_length(row_tests, nrow(x))
    expect_setequal(row_tests, seq_len(nrow(x)))

    group_folds <- tapply(cv$folds_ids, groups, function(z) unique(stats::na.omit(z)))
    expect_true(all(lengths(group_folds) == 1L))

    for(i in seq_along(cv$folds_list)){
        train_groups <- unique(groups[cv$folds_list[[i]][[1L]]])
        test_groups <- unique(groups[cv$folds_list[[i]][[2L]]])
        expect_false(any(test_groups %in% train_groups))
    }
}


test_that("cv_group keeps each group intact", {
    pts <- group_test_data()

    res <- cv_group(
        x = pts,
        group_col = "site",
        column = "occ",
        k = 3,
        biomod2 = FALSE,
        report = FALSE
    )

    expect_s3_class(res, "cv_group")
    expect_equal(res$group_col, "site")
    expect_equal(res$presence_bg, FALSE)
    expect_equal(res$k, 3)
    expect_equal(res$type, "Grouped (merged to k folds)")
    expect_group_integrity(res, pts)
})


test_that("cv_group leave-group-out mode leaves one group out per fold", {
    pts <- group_test_data()
    group_sizes <- table(pts$site)

    res <- cv_group(
        x = pts,
        group_col = "site",
        column = "occ",
        biomod2 = FALSE,
        report = FALSE
    )

    expect_s3_class(res, "cv_group")
    expect_equal(res$type, "Leave-group-out")
    expect_equal(res$k, length(group_sizes))
    expect_equal(res$n_groups, length(group_sizes))
    expect_equal(length(res$folds_list), length(group_sizes))
    expect_equal(dim(res$records), c(length(group_sizes), 4L))
    expect_null(res$biomod_table)
    expect_group_integrity(res, pts)

    test_groups <- vapply(res$folds_list, function(fold){
        unique(as.character(pts$site[fold[[2L]]]))
    }, character(1))
    expect_setequal(test_groups, names(group_sizes))
    expect_equal(unname(vapply(res$folds_list, function(fold) length(fold[[2L]]), integer(1))),
                 unname(as.integer(group_sizes[test_groups])))
})


test_that("cv_group merges groups into k folds without splitting groups", {
    pts <- group_test_data()

    res <- cv_group(
        x = pts,
        group_col = "site",
        column = "occ",
        k = 3,
        biomod2 = TRUE,
        report = FALSE
    )

    expect_equal(res$type, "Grouped (merged to k folds)")
    expect_equal(res$k, 3)
    expect_equal(length(res$folds_list), 3)
    expect_equal(dim(res$biomod_table), c(nrow(pts), 3L))
    expect_group_integrity(res, pts)

    group_folds <- tapply(res$folds_ids, pts$site, unique)
    expect_equal(sort(as.integer(table(unlist(group_folds)))), c(2L, 2L, 2L))

    for(i in seq_len(res$k)){
        test_set <- res$folds_list[[i]][[2L]]
        train_set <- res$folds_list[[i]][[1L]]
        expect_true(all(!res$biomod_table[test_set, i]))
        expect_true(all(res$biomod_table[train_set, i]))
    }
})


test_that("cv_group balance mode improves class coverage while preserving groups", {
    pts <- balanced_group_test_data()

    expect_warning(
        unbalanced <- cv_group(
            x = pts,
            group_col = "site",
            column = "occ",
            k = 3,
            balance = FALSE,
            biomod2 = FALSE,
            report = FALSE
        ),
        "zero records"
    )

    expect_silent(
        balanced <- cv_group(
            x = pts,
            group_col = "site",
            column = "occ",
            k = 3,
            balance = TRUE,
            iteration = 200,
            seed = 42,
            biomod2 = FALSE,
            report = FALSE,
            progress = FALSE
        )
    )

    expect_group_integrity(unbalanced, pts)
    expect_group_integrity(balanced, pts)
    expect_equal(balanced$type, "Grouped (merged to k folds)")

    test_cols <- grep("^test", names(balanced$records), value = TRUE)
    expect_gt(sum(unbalanced$records[, test_cols] == 0), 0)
    expect_true(all(balanced$records[, test_cols] > 0))

    balanced_again <- cv_group(
        x = pts,
        group_col = "site",
        column = "occ",
        k = 3,
        balance = TRUE,
        iteration = 200,
        seed = 42,
        biomod2 = FALSE,
        report = FALSE,
        progress = FALSE
    )
    expect_identical(balanced_again$folds_ids, balanced$folds_ids)
})


test_that("cv_group rejects missing or invalid groups", {
    pts <- group_test_data()

    expect_error(
        cv_group(x = pts, group_col = "missing_site", report = FALSE),
        "There is no column named 'missing_site'"
    )

    pts_na <- pts
    pts_na$site[1] <- NA_character_
    expect_error(
        cv_group(x = pts_na, group_col = "site", report = FALSE),
        "contains missing values"
    )

    pts_one <- pts
    pts_one$site <- "only_site"
    expect_error(
        cv_group(x = pts_one, group_col = "site", report = FALSE),
        "at least two distinct groups"
    )
})


test_that("cv_group plotting requires data and returns ggplot objects", {
    pts <- group_test_data()
    res <- cv_group(
        x = pts,
        group_col = "site",
        column = "occ",
        k = 3,
        biomod2 = FALSE,
        report = FALSE
    )

    expect_error(
        plot(res),
        "The original sample data are required"
    )

    plt <- cv_plot(cv = res, x = pts, num_plots = 1:2)
    expect_true(ggplot2::is_ggplot(plt))

    combined <- cv_plot(cv = res, x = pts, combine_folds = TRUE)
    expect_true(ggplot2::is_ggplot(combined))

    from_plot <- plot(res, data = pts, combine_folds = TRUE)
    expect_true(ggplot2::is_ggplot(from_plot))
})
