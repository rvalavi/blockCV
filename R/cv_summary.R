#' Summarise the quality of a set of cross-validation folds
#'
#' A one-call diagnostic that gathers the fold-quality signals scattered across \code{blockCV} into a
#' single object: fold sizes and class prevalence (always), the test-to-nearest-train distance and
#' environmental-novelty diagnostics (when a raster or prediction domain is supplied), and a list of
#' automatically detected warnings about degenerate folds.
#'
#' The structural part (sizes, class prevalence and the structural warnings) is read from the \code{cv}
#' object itself and needs no raster. When \code{r} (or \code{pred_points}/\code{model_domain}) is supplied,
#' the per-fold distance diagnostic of \code{\link{cv_distance}} is added; when the covariate raster \code{r}
#' is supplied, the per-fold environmental-novelty diagnostic of \code{\link{cv_similarity}} is added as well.
#' These reuse the same computations as those functions.
#'
#' The warnings are returned as data (a data.frame), not raised as \code{warning()} conditions, so they can
#' be inspected programmatically. The structural warnings are a k-fold concept and are skipped for
#' leave-one-out objects (\code{\link{cv_buffer}}, \code{\link{cv_nndm}}), whose folds hold a single test
#' point by design. The flagged issues are:
#' \itemize{
#'   \item{\code{empty_test} - a fold with no test points.}
#'   \item{\code{single_class_test} - a fold whose test set contains a single class (breaks AUC and other
#'   class-wise metrics).}
#'   \item{\code{class_missing_train} - a class entirely absent from a fold's training set.}
#'   \item{\code{tiny_test} - a fold with fewer than \code{min_test} test points (unstable per-fold estimate).}
#'   \item{\code{imbalance} - a class severely under-represented in a fold's test set relative to an even split.}
#'   \item{\code{high_leakage} - a fold where at least 90\% of test points sit nearer to training than the
#'   median prediction distance (only available when the distance diagnostic is computed).}
#' }
#'
#' @param cv a \code{blockCV} cross-validation object, i.e. the output of \code{\link{cv_spatial}},
#' \code{\link{cv_cluster}}, \code{\link{cv_group}}, \code{\link{cv_buffer}}, \code{\link{cv_nndm}} or
#' \code{\link{cv_knndm}}.
#' @param x a simple features (sf) object of the sample points used to create \code{cv}. Required only for
#' the distance/novelty diagnostics (when \code{r}, \code{pred_points} or \code{model_domain} is supplied).
#' @param r a terra SpatRaster of the environmental covariates. When supplied, both the distance and the
#' environmental-novelty diagnostics are computed. Optional.
#' @param pred_points,model_domain optional prediction points or a model-domain polygon used by the distance
#' diagnostic; see \code{\link{cv_distance}}. Supplying either (without \code{r}) adds the distance diagnostic
#' but not the novelty diagnostic.
#' @param method the similarity method passed to \code{\link{cv_similarity}}: \code{"MESS"}, \code{"L1"} or
#' \code{"L2"}.
#' @param num_sample number of random raster samples used by the distance and novelty diagnostics.
#' @param min_test integer; folds with fewer than this many test points are flagged as \code{tiny_test}.
#' @param seed integer; an optional random seed for the raster-sampling baselines.
#' @param progress logical; whether to show a progress bar for the novelty diagnostic.
#'
#' @seealso \code{\link{cv_distance}} and \code{\link{cv_similarity}} for the individual diagnostics, and
#' \code{\link{cv_plot}} to visualise the folds
#'
#' @return an object of class \code{cv_summary}: a list with
#' \itemize{
#'   \item{\code{n_folds} - the number of folds.}
#'   \item{\code{is_loo} - whether \code{cv} is a leave-one-out object.}
#'   \item{\code{records} - the per-fold train/test counts per class (fold sizes and class prevalence).}
#'   \item{\code{distances} - the per-fold test-to-nearest-train distance summary, or \code{NULL}.}
#'   \item{\code{novelty} - the per-fold environmental-novelty (extrapolation) summary, or \code{NULL}.}
#'   \item{\code{warnings} - a data.frame of flagged degenerate folds (\code{fold}, \code{type},
#'   \code{message}); zero rows when nothing is flagged.}
#' }
#' @export
#'
#' @examples
#' \donttest{
#' library(blockCV)
#'
#' points <- read.csv(system.file("extdata/", "species.csv", package = "blockCV"))
#' pa_data <- sf::st_as_sf(points, coords = c("x", "y"), crs = 7845)
#' covar <- terra::rast(system.file("extdata/au/bio_5.tif", package = "blockCV"))
#'
#' sb <- cv_spatial(x = pa_data, column = "occ", size = 450000, k = 5, iteration = 1)
#'
#' # structural summary + warnings only (no raster)
#' cv_summary(sb)
#'
#' # add the distance and novelty diagnostics
#' cv_summary(sb, x = pa_data, r = covar, num_sample = 5000)
#' }
cv_summary <- function(
        cv,
        x = NULL,
        r = NULL,
        pred_points = NULL,
        model_domain = NULL,
        method = "MESS",
        num_sample = 10000L,
        min_test = 5L,
        seed = NULL,
        progress = FALSE
){
    .check_cv(cv)
    is_loo <- .is_loo(cv)

    # the distance diagnostic needs a prediction domain (r / pred_points / model_domain);
    # the novelty diagnostic needs the covariate raster r
    want_distance <- !is.null(r) || !is.null(pred_points) || !is.null(model_domain)
    want_novelty  <- !is.null(r)

    distances <- NULL
    novelty <- NULL
    pbg <- FALSE
    if(want_distance){
        if(is.null(x)){
            stop("'x' (the sample points) is required to compute the distance/novelty diagnostics.")
        }
        # only the per-fold summary is needed, so skip the random baseline (add_random = FALSE)
        dist_data <- .cv_distance_data(
            cv = cv, x = x, r = r, pred_points = pred_points, model_domain = model_domain,
            add_random = FALSE, num_sample = num_sample, seed = seed
        )
        distances <- dist_data$distances
        pbg <- isTRUE(dist_data$pbg)
        if(want_novelty){
            novelty <- .cv_similarity_data(
                cv = cv, x = x, r = r, method = method,
                num_sample = num_sample, seed = seed, progress = progress
            )$extrapolation
        }
    }

    out <- list(
        n_folds = length(cv$folds_list),
        is_loo = is_loo,
        records = cv$records,
        distances = distances,
        novelty = novelty,
        pbg = pbg,
        warnings = .cv_warnings(cv, distances = distances, min_test = min_test, is_loo = is_loo)
    )
    class(out) <- "cv_summary"
    out
}


# ---- structural warnings ---------------------------------------------------
# Detect degenerate folds from the cv object's train/test records (and, when
# available, the distance diagnostic). Returns the flags as data, never as
# warning() conditions. Structural per-fold checks are a k-fold concept, so they
# are skipped for leave-one-out objects whose folds hold a single test point by
# design.
.cv_warnings <- function(cv, distances = NULL, min_test = 5L, is_loo = FALSE){
    warns <- data.frame(fold = integer(), type = character(), message = character(),
                        stringsAsFactors = FALSE)
    add <- function(fold, type, message){
        warns <<- rbind(warns, data.frame(fold = fold, type = type, message = message,
                                          stringsAsFactors = FALSE))
    }

    records <- cv$records
    if(is_loo || is.null(records) || !nrow(records)) return(warns)

    nm <- names(records)
    test_cols <- grep("^test", nm, value = TRUE)
    train_cols <- grep("^train", nm, value = TRUE)
    k <- nrow(records)
    multiclass <- length(test_cols) >= 2

    for(i in seq_len(k)){
        test_counts <- unlist(records[i, test_cols, drop = FALSE], use.names = FALSE)
        train_counts <- unlist(records[i, train_cols, drop = FALSE], use.names = FALSE)
        n_test <- sum(test_counts)

        if(n_test == 0){
            add(i, "empty_test", "fold has no test points")
            next
        }
        if(multiclass){
            if(sum(test_counts > 0) <= 1){
                add(i, "single_class_test",
                    "test set has only one class (breaks AUC and class-wise metrics)")
            }
            if(any(train_counts == 0)){
                miss <- sub("^train_", "", train_cols[train_counts == 0])
                add(i, "class_missing_train",
                    sprintf("class(es) absent from training: %s", paste(miss, collapse = ", ")))
            }
        }
        if(n_test < min_test){
            add(i, "tiny_test",
                sprintf("only %d test point(s) (< %d); per-fold estimate is unstable", n_test, min_test))
        }
    }

    # severe imbalance: a class present but well below its even share of a fold's
    # test set. Only checked for classes common enough that an even split gives at
    # least ~2 per fold (n_class >= 2k); rarer classes cannot be split evenly and
    # are covered by the single_class_test / class_missing_train flags instead.
    if(multiclass && k > 1){
        for(col in test_cols){
            counts <- as.numeric(records[[col]])
            n_class <- sum(counts)
            if(n_class < 2 * k) next
            expected <- n_class / k
            cls <- sub("^test_", "", col)
            for(i in seq_len(k)){
                if(counts[i] > 0 && counts[i] < 0.25 * expected){
                    add(i, "imbalance",
                        sprintf("class '%s' under-represented in test (%d vs ~%.0f expected)",
                                cls, counts[i], round(expected)))
                }
            }
        }
    }

    # high leakage: only meaningful once the distance diagnostic has been computed
    if(!is.null(distances) && "pct_below_pred" %in% names(distances)){
        hi <- distances[!is.na(distances$pct_below_pred) & distances$pct_below_pred >= 90, , drop = FALSE]
        for(i in seq_len(nrow(hi))){
            add(hi$fold[i], "high_leakage",
                sprintf("%.0f%% of test points nearer to train than the median prediction distance (optimistic/leaky fold)",
                        hi$pct_below_pred[i]))
        }
    }

    if(nrow(warns)) warns[order(warns$fold, warns$type), , drop = FALSE] else warns
}


# print at most 'n' rows of a per-fold table, noting any that are hidden
.print_capped <- function(df, n = 12){
    if(nrow(df) > n){
        print(utils::head(df, n), row.names = FALSE)
        cat(sprintf("  ... %d more row(s)\n", nrow(df) - n))
    } else {
        print(df, row.names = FALSE)
    }
}


#' @export
#' @method print cv_summary
print.cv_summary <- function(x, ...){
    cat("blockCV fold-quality summary\n")
    cat(sprintf("Folds: %d%s\n", x$n_folds, if(isTRUE(x$is_loo)) " (leave-one-out)" else ""))

    if(!is.null(x$records)){
        if(isTRUE(x$is_loo)){
            cat(sprintf("\nLeave-one-out design: %d held-out point(s); per-fold size table omitted.\n",
                        nrow(x$records)))
        } else {
            cat("\nFold sizes and class prevalence (train/test counts):\n")
            print(x$records)
        }
    }
    if(!is.null(x$distances)){
        cat("\nTest-to-nearest-train distances (leakage):\n")
        .print_capped(x$distances)
    }
    if(!is.null(x$novelty)){
        cat("\nEnvironmental novelty (extrapolation):\n")
        .print_capped(x$novelty)
    }
    if(isTRUE(x$pbg) && (!is.null(x$distances) || !is.null(x$novelty))){
        cat("\nNote: presence-background object; distance and novelty diagnostics computed on presence points only (background samples excluded).\n")
    }

    cat("\nWarnings:\n")
    if(is.null(x$warnings) || !nrow(x$warnings)){
        cat("  none\n")
    } else {
        for(i in seq_len(nrow(x$warnings))){
            cat(sprintf("  [fold %s] %s: %s\n",
                        x$warnings$fold[i], x$warnings$type[i], x$warnings$message[i]))
        }
    }
    invisible(x)
}
