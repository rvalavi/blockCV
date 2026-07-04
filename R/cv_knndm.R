#' Use the k-fold Nearest Neighbour Distance Matching (kNNDM) to separate train and test folds
#'
#' A k-fold version of the Nearest Neighbour Distance Matching algorithm (Linnenbrink et al., 2024).
#' Unlike \code{\link{cv_nndm}}, which is a leave-one-out (LOO) strategy producing as many folds as
#' there are points, kNNDM produces a small number of \code{k} folds (cheap to fit models on) while
#' still matching the nearest neighbour distance distribution function between the test and training
#' data to the one between the target prediction points and the training data.
#'
#' The method groups the sample points and assigns whole groups to folds, tuning the degree of
#' aggregation so that the resulting test-to-train nearest neighbour distance distribution
#' (\eqn{\hat{G}_{j}^{*}}) best matches the prediction-to-sample distribution (\eqn{\hat{G}_{ij}}). The
#' match is scored with the Wasserstein-1 distance (\code{W} in the output; the area between the two
#' empirical distribution functions), and the grouping with the smallest value is returned. Following
#' Linnenbrink et al. (2024), a one-sided Kolmogorov-Smirnov test is first used to check whether the
#' points are clustered relative to the prediction area; if they are not, a random cross-validation is
#' returned. When a clustering is needed, the intermediate groups are ordered along the first principal
#' component of the coordinates and distributed across the \code{k} folds, and only groupings in which
#' no fold holds more than \code{maxp} of the data are considered.
#'
#' Three ways of generating the groups are available via \code{clustering}: \code{"hierarchical"} and
#' \code{"kmeans"} follow Linnenbrink et al. (2024), while \code{"blocks"} (the default) is a
#' \code{blockCV}-specific variant that groups the points with hexagonal (or square) spatial blocks, in
#' the spirit of \code{\link{cv_spatial}}. The \code{"blocks"} option is a \code{blockCV} extension and
#' is not part of the original kNNDM algorithm. When \code{space = "feature"} the distances are computed
#' in the (scaled) covariate space of \code{r} instead of the geographical space, and \code{"blocks"} is
#' not available.
#'
#' When \code{column} is supplied and \code{balance_classes = TRUE}, class balance is used as a
#' validity gate during the candidate scan: among groupings with every class present in every test
#' fold, the one with the smallest \code{W} is selected. If no class-complete grouping is available,
#' the overall smallest-\code{W} grouping is returned and the final records table reports the missing
#' class(es).
#'
#' @inheritParams cv_spatial
#' @param x a simple features (sf) or SpatialPoints object of spatial sample data (e.g., species
#' data or ground truth sample for image classification).
#' @param column character (optional). Indicating the name of the column in which response variable
#' (e.g. species data as a binary response i.e. 0s and 1s) is stored. This is used to report whether
#' all folds contain all classes, and by default to prefer class-complete candidate groupings.
#' @param r a terra SpatRaster object. This defines the area that the model is going to predict; when
#' neither \code{predpoints} nor \code{modeldomain} is supplied, prediction points are sampled from it.
#' It is also required (for the covariates) when \code{space = "feature"}.
#' @param predpoints a simple features (sf) object of prediction points (optional). If provided, these are
#' used directly as the prediction locations instead of sampling from \code{r} or \code{modeldomain}.
#' @param modeldomain an sf polygon of the prediction area (optional). If provided (and \code{predpoints}
#' is not), prediction points are sampled from it.
#' @param k integer value. The number of desired folds for cross-validation. The default is \code{k = 5}.
#' @param maxp numeric; strictly between \code{1/k} and 1. The maximum proportion of the points allowed
#' in any one fold; groupings that exceed it are discarded.
#' @param clustering character; the method used to group the points before allocating them to folds.
#' One of \code{"blocks"} (default; a \code{blockCV} spatial-blocks variant), \code{"hierarchical"}
#' (\code{\link[stats]{hclust}}), or \code{"kmeans"} (\code{\link[stats]{kmeans}}).
#' @param space character; \code{"geographical"} (default) matches distances in geographical space, or
#' \code{"feature"} matches distances in the (scaled) covariate space of \code{r}.
#' @param hexagon logical. Creates hexagonal (default) spatial blocks when \code{clustering = "blocks"}.
#' If \code{FALSE}, square blocks are used.
#' @param num_sample integer; the number of prediction points to sample from \code{r} or \code{modeldomain}
#' when \code{predpoints} is not supplied.
#' @param sampling either \code{"regular"} (default) or \code{"random"} for sampling prediction points.
#' @param nk_len integer; the number of candidate groupings (numbers of clusters / block sizes) to explore.
#' @param linkage character; the agglomeration method passed to \code{\link[stats]{hclust}} when
#' \code{clustering = "hierarchical"}.
#' @param scale logical; whether to scale the covariates when \code{space = "feature"}.
#' @param balance_classes logical (optional). If \code{TRUE}, and \code{column} is supplied, candidate
#' groupings with every class represented in every test fold are preferred. If no such grouping exists,
#' the best distance-matching grouping is returned and the usual zero-record warning is reported.
#' Defaults to \code{TRUE} when \code{column} is supplied and \code{FALSE} otherwise.
#' @param seed integer; a random seed for reproducibility.
#' @param plot logical; whether to plot the distance distribution functions. Defaults to \code{interactive()}.
#' @param report logical; whether to print summary of records in each fold. Defaults to \code{interactive()}.
#'
#' @seealso \code{\link{cv_nndm}}, \code{\link{cv_spatial}}, \code{\link{cv_cluster}},
#' \code{\link{cv_spatial_autocor}}, and \code{CAST::knndm()} for the authors'
#' CAST implementation of kNNDM.
#'
#' @references Linnenbrink, J., Milà, C., Ludwig, M., & Meyer, H. (2024). kNNDM CV: k-fold nearest
#' neighbour distance matching cross-validation for map accuracy estimation. Geoscientific Model
#' Development, 17(15), 5897-5912. https://doi.org/10.5194/gmd-17-5897-2024
#'
#' Milà, C., Mateu, J., Pebesma, E., & Meyer, H. (2022). Nearest neighbour distance matching
#' leave-one-out cross-validation for map validation. Methods in Ecology and Evolution, 13(6), 1304-1316.
#'
#' @return An object of class S3. A list of objects including:
#'    \itemize{
#'     \item{folds_list - a list containing the folds. Each fold has two vectors with the training (first) and testing (second) indices}
#'     \item{folds_ids - a vector of values indicating the number of the fold for each observation (each number corresponds to the same point in x)}
#'     \item{biomod_table - a matrix with the folds to be used in \pkg{biomod2} package}
#'     \item{k - number of the folds}
#'     \item{column - the name of the column if provided}
#'     \item{type - the clustering method used (or "random" if a random CV was returned)}
#'     \item{space - whether distances were matched in geographical or feature space}
#'     \item{q - the selected number of intermediate groups}
#'     \item{W - the Wasserstein statistic of the selected folds (lower means a better match)}
#'     \item{Gij - the prediction-to-sample nearest neighbour distances}
#'     \item{Gj - the sample-to-sample (leave-one-out) nearest neighbour distances}
#'     \item{Gjstar - the test-to-train nearest neighbour distances of the selected folds}
#'     \item{records - a table with the number of points in each category of training and testing}
#'     }
#' @export
#'
#' @examples
#' \donttest{
#' library(blockCV)
#'
#' # import presence-absence species data
#' points <- read.csv(system.file("extdata/", "species.csv", package = "blockCV"))
#' # make an sf object from data.frame
#' pa_data <- sf::st_as_sf(points, coords = c("x", "y"), crs = 7845)
#'
#' # load raster data
#' path <- system.file("extdata/au/bio_5.tif", package = "blockCV")
#' covar <- terra::rast(path)
#'
#' knndm <- cv_knndm(x = pa_data,
#'                   column = "occ", # optional
#'                   r = covar,
#'                   k = 5,
#'                   num_sample = 5000)
#'
#' }
cv_knndm <- function(
        x,
        column = NULL,
        r = NULL,
        predpoints = NULL,
        modeldomain = NULL,
        k = 5L,
        maxp = 0.5,
        clustering = "blocks",
        space = "geographical",
        hexagon = TRUE,
        num_sample = 1e4,
        sampling = "regular",
        nk_len = 100L,
        linkage = "ward.D2",
        scale = TRUE,
        biomod2 = TRUE,
        deg_to_metre = 111325,
        balance_classes = NULL,
        seed = NULL,
        plot = interactive(),
        report = interactive()
){

    # check the arguments
    clustering <- match.arg(clustering, choices = c("blocks", "hierarchical", "kmeans"))
    space <- match.arg(space, choices = c("geographical", "feature"))
    sampling <- match.arg(sampling, choices = c("regular", "random"))

    # check x is an sf object
    x <- .check_x(x)
    # is column in x?
    column <- .check_column(column, x)
    if(is.null(balance_classes)){
        balance_classes <- !is.null(column)
    }
    if(!is.logical(balance_classes) || length(balance_classes) != 1L || is.na(balance_classes)){
        stop("'balance_classes' must be TRUE or FALSE.")
    }
    if(balance_classes && is.null(column)){
        stop("'balance_classes = TRUE' requires 'column'.")
    }
    # x's CRS must be defined
    if(is.na(sf::st_crs(x))){
        stop("The coordinate reference system of 'x' must be defined.")
    }

    # spatial blocks are only meaningful in geographical space
    if(space == "feature" && clustering == "blocks"){
        clustering <- "hierarchical"
        message("Spatial 'blocks' clustering is not available in feature space; using 'hierarchical'.")
    }

    # check k and maxp
    n <- nrow(x)
    k <- as.integer(abs(k))
    if(k < 2 || k >= n) stop("'k' must be a natural number between 2 and the number of points.")
    if(!(maxp > 1 / k && maxp < 1)) stop("'maxp' must be strictly between 1/k and 1.")

    if(!is.null(seed)) set.seed(seed)

    # feature space requires the covariate raster
    if(space == "feature" && is.null(r)){
        stop("'r' (covariate raster) is required when space = 'feature'.")
    }
    # make sure r and x share the same CRS
    if(!is.null(r)){
        r <- .check_r(r)
        if(!.same_crs(x, terra::crs(r))){
            stop("The coordinate reference systems of 'x' and 'r' must match.")
        }
    }

    # resolve the prediction points -------------------------------------------
    predpts <- .knndm_predpoints(x, r, predpoints, modeldomain, num_sample, sampling)

    # distances and clustering coordinates ------------------------------------
    if(space == "geographical"){
        tdist <- sf::st_distance(x)
        units(tdist) <- NULL
        # prediction-to-sample nearest neighbour distances
        Gij <- sf::st_distance(predpts, x)
        units(Gij) <- NULL
        Gij <- apply(Gij, 1, min)
        # coordinates used for PC1 ordering and k-means
        clust_coords <- sf::st_coordinates(x)
        if(clustering == "blocks" && sf::st_is_longlat(x)){
            warning("Block sizes for geographic (lon/lat) coordinates are approximate (converted via 'deg_to_metre').")
        }
    } else{
        # feature space: distances in the (scaled) covariate space
        ft <- .knndm_features(x, r, predpts, scale)
        clust_coords <- ft$train
        tdist <- as.matrix(stats::dist(ft$train))
        Gij <- .nn_cross(ft$pred, ft$train)
    }

    # sample-to-nearest-other-sample (LOO) distances
    Gj <- vapply(seq_len(n), function(i) min(tdist[i, -i]), numeric(1))

    # random-CV gate (Linnenbrink et al., 2024) -------------------------------
    ks <- suppressWarnings(stats::ks.test(Gj, Gij, alternative = "greater"))
    if(ks$p.value >= 0.05){
        message("Gij <= Gj; a random cross-validation assignment is returned (clustering not required).")
        fold_vect <- sample(rep(seq_len(k), ceiling(n / k)), size = n)
        Gjstar <- .nn_diff_fold(tdist, fold_vect)
        W <- .wasserstein(Gjstar, Gij)
        method <- "random"
        q <- k

    } else{
        # search over candidate groupings -------------------------------------
        cand <- .knndm_candidates(clustering, k, n, nk_len, tdist, clust_coords,
                                  linkage, x, hexagon, deg_to_metre)

        best <- NULL
        best_complete <- NULL
        for(clust in cand){
            folds <- .merge_clusters(clust, clust_coords, k, maxp)
            if(is.null(folds)) next # violates maxp or cannot form k folds
            gjs <- .nn_diff_fold(tdist, folds)
            w <- .wasserstein(gjs, Gij)
            candidate <- list(
                W = w,
                fold_vect = folds,
                Gjstar = gjs,
                q = length(unique(clust))
            )
            if(is.null(best) || w < best$W){
                best <- candidate
            }
            if(balance_classes &&
               .knndm_class_complete(folds, x, column, k) &&
               (is.null(best_complete) || w < best_complete$W)){
                best_complete <- candidate
            }
        }
        # no valid kNNDM configuration was found after clustering was required
        if(is.null(best)){
            stop("No grouping satisfied 'maxp'. Increase 'maxp', reduce 'k', increase 'nk_len', or use a different clustering method.")
        } else{
            method <- clustering
        }
        if(balance_classes && !is.null(best_complete)){
            best <- best_complete
        }
        fold_vect <- best$fold_vect
        Gjstar <- best$Gjstar
        q <- best$q
        W <- best$W
    }

    # build folds, biomod table and records (k-fold output) -------------------
    fold_list <- list()
    biomod_table <- data.frame(RUN1 = rep(TRUE, n))
    for(i in seq_len(k)){
        test_set <- which(fold_vect == i)
        train_set <- which(fold_vect != i)
        fold_list[[i]] <- list(train_set, test_set)
        if(biomod2){ # creating a biomod2 CV.user.table for validation
            colm <- paste0("RUN", i)
            biomod_table[, colm] <- FALSE
            biomod_table[train_set, colm] <- TRUE
        }
    }

    # calculate train test table summary
    train_test_table <- .table_summary(fold_list, x, column, k)
    # give a warning if any fold is empty
    zerofolds <- which(apply(train_test_table, 1, function(x) any(x < 1)))
    if(length(zerofolds) > 0){
        if(length(zerofolds) > 1){
            warning("Folds ", paste(zerofolds, collapse = ", "), " have class(es) with zero records")
        } else{
            warning("Fold ", zerofolds, " has class(es) with zero records")
        }
    }
    if(report) print(train_test_table)

    # distance distribution functions for plotting
    r_range <- seq(0, max(c(Gij, Gjstar)), length.out = 200)
    plot_data <- data.frame(pr = stats::ecdf(Gij)(r_range),
                            lo = stats::ecdf(Gj)(r_range),
                            nn = stats::ecdf(Gjstar)(r_range),
                            r = r_range)

    plt <- ggplot2::ggplot(plot_data, ggplot2::aes(x = get("r"))) +
        ggplot2::geom_step(alpha = 0.7, linewidth = 1.2, ggplot2::aes(y = get("pr"), color = "Prediction")) +
        ggplot2::geom_step(alpha = 0.7, linewidth = 0.6, ggplot2::aes(y = get("lo"), color = "LOO")) +
        ggplot2::geom_step(alpha = 0.7, linewidth = 1.2, ggplot2::aes(y = get("nn"), color = "kNNDM")) +
        ggplot2::scale_color_manual(values = c("Prediction" = "#000000",
                                               "LOO" = "#56B4E9",
                                               "kNNDM" = "#E69F00")) +
        ggplot2::labs(color = "", x = "r", y = expression(G[r])) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.text = ggplot2::element_text(size = 12))

    if(plot) plot(plt)

    final_objs <- list(
        folds_list = fold_list,
        folds_ids = fold_vect,
        biomod_table = if(biomod2) as.matrix(biomod_table) else NULL,
        k = k,
        column = column,
        type = paste0("kNNDM (", method, ")"),
        space = space,
        q = q,
        W = W,
        Gij = Gij,
        Gj = Gj,
        Gjstar = Gjstar,
        plot = plt,
        records = train_test_table
    )

    class(final_objs) <- c("cv_knndm")
    return(final_objs)
}


# resolve prediction points from predpoints, modeldomain or a sampled raster
.knndm_predpoints <- function(x, r, predpoints, modeldomain, num_sample, sampling){
    if(!is.null(predpoints)){
        predpoints <- .check_x(predpoints, name = "predpoints")
        if(!.same_crs(x, predpoints)) stop("The coordinate reference systems of 'x' and 'predpoints' must match.")
        return(sf::st_geometry(predpoints))
    }
    if(!is.null(modeldomain)){
        modeldomain <- .check_x(modeldomain, name = "modeldomain")
        if(!.same_crs(x, modeldomain)) stop("The coordinate reference systems of 'x' and 'modeldomain' must match.")
        return(sf::st_sample(sf::st_geometry(modeldomain), size = num_sample, type = sampling))
    }
    if(is.null(r)) stop("Provide one of 'r', 'predpoints', or 'modeldomain' to define the prediction area.")
    rx <- tryCatch(
        terra::spatSample(r[[1]], size = num_sample, method = sampling,
                          values = FALSE, na.rm = TRUE, as.points = TRUE),
        error = function(cond) stop("Sampling the raster 'r' failed: ", conditionMessage(cond))
    )
    if(is.null(rx) || nrow(rx) == 0) stop("Sampling the raster 'r' returned no points.")
    return(sf::st_geometry(sf::st_as_sf(rx)))
}


# extract and scale covariates at the sample and prediction points (feature space)
.knndm_features <- function(x, r, predpts, scale){
    train <- terra::extract(r, x, ID = FALSE)
    if(anyNA(train)) stop("Some 'x' points fall outside 'r' (missing covariate values).")
    pred <- terra::extract(r, sf::st_sf(geometry = predpts), ID = FALSE)
    pred <- pred[stats::complete.cases(pred), , drop = FALSE]
    train <- as.matrix(train)
    pred <- as.matrix(pred)
    if(scale){
        train <- base::scale(train)
        ctr <- attr(train, "scaled:center")
        scl <- attr(train, "scaled:scale")
        pred <- sweep(sweep(pred, 2, ctr, "-"), 2, scl, "/")
        attr(train, "scaled:center") <- NULL
        attr(train, "scaled:scale") <- NULL
    }
    return(list(train = train, pred = pred))
}


# generate candidate group assignments for the requested clustering method
.knndm_candidates <- function(clustering, k, n, nk_len, tdist, clust_coords,
                              linkage, x, hexagon, deg_to_metre){
    cand <- list()
    if(clustering == "hierarchical"){
        hc <- stats::hclust(stats::as.dist(tdist), method = linkage)
        grid <- unique(round(exp(seq(log(k), log(max(k, n - 2)), length.out = nk_len))))
        for(q in grid) cand[[length(cand) + 1]] <- stats::cutree(hc, k = q)

    } else if(clustering == "kmeans"){
        if(sf::st_is_longlat(x)){
            warning("k-means clustering on geographic (lon/lat) coordinates is approximate; consider projecting 'x'.")
        }
        grid <- unique(round(exp(seq(log(k), log(max(k, n - 2)), length.out = nk_len))))
        for(q in grid){
            cand[[length(cand) + 1]] <- stats::kmeans(clust_coords, centers = q, iter.max = 100)$cluster
        }

    } else{ # blocks (blockCV variant)
        bb <- sf::st_bbox(x)
        if(sf::st_is_longlat(x)){
            dx <- (bb[["xmax"]] - bb[["xmin"]]) * deg_to_metre
            dy <- (bb[["ymax"]] - bb[["ymin"]]) * deg_to_metre
        } else{
            dx <- bb[["xmax"]] - bb[["xmin"]]
            dy <- bb[["ymax"]] - bb[["ymin"]]
        }
        diag_m <- sqrt(dx^2 + dy^2)
        nb <- unique(round(exp(seq(log(2), log(max(3, ceiling(sqrt(n)))), length.out = nk_len))))
        for(s in diag_m / nb){
            bid <- .points_to_blocks(x, s, hexagon, deg_to_metre)
            if(length(unique(bid)) < k) next
            cand[[length(cand) + 1]] <- bid
        }
    }
    return(cand)
}


# assign points to spatial blocks of a given size (reuses .make_blocks)
.points_to_blocks <- function(x, size, hexagon, degree){
    blocks <- .make_blocks(
        x_obj = x,
        blocksize = size,
        hexagonal = hexagon,
        extend_perc = 0,
        degree = degree
    )
    inter <- sf::st_intersects(sf::st_geometry(x), sf::st_geometry(blocks))
    bid <- vapply(inter, function(z) if(length(z)) z[1] else NA_integer_, integer(1))
    # nearest block for points that fall outside the tiling
    if(anyNA(bid)){
        na_idx <- which(is.na(bid))
        bid[na_idx] <- sf::st_nearest_feature(sf::st_geometry(x)[na_idx], sf::st_geometry(blocks))
    }
    return(as.integer(bid))
}


# check whether every class is represented in every test fold
.knndm_class_complete <- function(fold_vect, x, column, k){
    fold_list <- lapply(seq_len(k), function(i){
        list(which(fold_vect != i), which(fold_vect == i))
    })
    train_test_table <- .table_summary(fold_list, x, column, k)
    test_cols <- startsWith(names(train_test_table), "test_")
    all(as.matrix(train_test_table[, test_cols, drop = FALSE]) > 0)
}


# merge intermediate groups into k folds along the first principal component
# following CAST's kNNDM assignment. Returns NULL if a fold would exceed maxp.
.merge_clusters <- function(clust, coords, k, maxp){
    ids <- unique(clust)
    if(length(ids) < k) return(NULL)
    n <- length(clust)
    cap <- maxp * n
    target <- n / k

    # order groups along the first principal component of the coordinates
    pc1 <- stats::prcomp(coords, center = TRUE, scale. = FALSE, rank. = 1)$x[, 1]
    score <- vapply(ids, function(g) mean(pc1[clust == g]), numeric(1))
    size <- vapply(ids, function(g) sum(clust == g), numeric(1))
    ord <- order(score)
    ids <- ids[ord]
    size <- size[ord]

    fold_of <- integer(length(ids))
    next_fold <- 1L
    for(i in seq_along(ids)){
        if(size[i] >= target){
            if(next_fold > k) return(NULL)
            fold_of[i] <- next_fold
            next_fold <- next_fold + 1L
        }
    }

    remaining_folds <- setdiff(seq_len(k), fold_of[fold_of > 0L])
    missing <- fold_of == 0L
    if(any(missing)){
        if(!length(remaining_folds)) return(NULL)
        fold_of[missing] <- rep(remaining_folds, length.out = sum(missing))
    }
    out <- integer(n)
    for(i in seq_along(ids)) out[clust == ids[i]] <- fold_of[i]

    load <- tabulate(out, nbins = k)
    if(length(unique(out)) < k) return(NULL) # must form k non-empty folds
    if(any(load > cap)) return(NULL) # cannot satisfy maxp
    return(out)
}


# nearest neighbour distance to a point in a different fold (Gjstar)
.nn_diff_fold <- function(distmat, folds){
    vapply(seq_along(folds), function(i){
        oth <- folds != folds[i]
        if(any(oth)) min(distmat[i, oth]) else NA_real_
    }, numeric(1))
}


# minimum (Euclidean) distance from each row of 'a' to the rows of 'b'
.nn_cross <- function(a, b){
    tb <- t(b)
    vapply(seq_len(nrow(a)), function(i){
        sqrt(min(colSums((tb - a[i, ])^2)))
    }, numeric(1))
}


# Wasserstein-1 distance (area between the empirical distribution functions)
.wasserstein <- function(a, b){
    grid <- sort(unique(c(a, b)))
    fa <- stats::ecdf(a)(grid)
    fb <- stats::ecdf(b)(grid)
    w <- diff(grid)
    sum(abs(fa[-length(fa)] - fb[-length(fb)]) * w)
}


# compare the CRS of two spatial objects (accepts sf objects or a CRS/WKT)
.same_crs <- function(a, b){
    ca <- sf::st_crs(a)
    cb <- sf::st_crs(b)
    if(is.na(ca) || is.na(cb)) return(FALSE)
    isTRUE(ca == cb)
}


#' @export
#' @method print cv_knndm
print.cv_knndm <- function(x, ...){
    print(class(x))
}


#' @export
#' @method plot cv_knndm
plot.cv_knndm <- function(x, y, ...){
    plot(x$plot)
    message("Please use cv_plot function to plot each fold.")
}


#' @export
#' @method summary cv_knndm
summary.cv_knndm <- function(object, ...){
    cat("Number of recoreds in each training and testing fold:\n")
    print(object$records)
}
