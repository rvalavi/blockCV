#' Use the Nearest Neighbour Distance Matching (NNDM) to separate train and test folds
#'
#' A fast implementation of the Nearest Neighbour Distance Matching (NNDM) algorithm (Milà et al., 2022) in C++. Similar
#' to \code{\link{cv_buffer}}, this is a variation of leave-one-out (LOO) cross-validation. It tries to match the
#' nearest neighbour distance distribution function between the test and training data to the nearest neighbour
#' distance distribution function between the target prediction and training points (Milà et al., 2022).
#'
#' When working with presence-background data (presences and \emph{background} points sampled across the study
#' area to represent the available conditions rather than confirmed absences, specified by
#' \code{presence_bg = TRUE}), only presence records are used
#' for specifying the folds (recommended). The testing fold comprises only the target \emph{presence} point (optionally,
#' all background points within the distance are also included when \code{add_bg = TRUE}; this is the
#' distance that matches the nearest neighbour distance distribution function of training-testing presences and
#' training-presences and prediction points; often lower than \code{size}).
#' Any non-target presence points inside the distance are excluded.
#' All points (presence and background) outside of distance are used for the training set.
#' The methods cycles through all the presence data, so the number of folds is equal to
#' the number of presence points in the dataset.
#'
#' For all other types of data (including presence-absence, count, continuous, and multi-class)
#' set \code{presence_bg = FALE}, and the function behaves similar to the methods
#' explained by Milà and colleagues (2022).
#'
#' @param x a simple features (sf) or SpatialPoints object of spatial sample data (e.g., species
#' data or ground truth sample for image classification).
#' @inheritParams cv_buffer
#' @param r a terra SpatRaster object of a predictor variable (optional). This defines the area that
#' the model is going to predict; prediction points are sampled from it when neither \code{pred_points}
#' nor \code{model_domain} is supplied. One of \code{r}, \code{pred_points}, or \code{model_domain} is required.
#' @param pred_points a simple features (sf) object of prediction points (optional). If provided, these are
#' used directly as the prediction locations instead of sampling from \code{r} or \code{model_domain}.
#' @param model_domain an sf polygon of the prediction area (optional). If provided (and \code{pred_points}
#' is not), prediction points are sampled from it.
#' @param size numeric value of the range of spatial autocorrelation (the \code{phi} parameter).
#' This distance should be in \strong{metres} (or in the coordinate units of \code{x} when its CRS is
#' undefined). The range could be explored by \code{\link{cv_spatial_autocor}}.
#' @param num_sample integer; the number of sample points from predictor (\code{r}) to be used for calculating
#' the G function of prediction points.
#' @param sampling either \code{"random"} or \code{"regular"} for sampling prediction points.
#' When  \code{sampling = "regular"}, the actual number of samples might be less than \code{num_sample}
#' for non-rectangular rasters (points falling on no-value areas are removed).
#' @param min_train numeric; between 0 and 1. A constraint on the minimum proportion of train points in each fold.
#' @inheritParams cv_buffer
#' @inheritParams cv_buffer
#' @param plot logical; whether to plot the G functions. Defaults to \code{interactive()}.
#' @param report logical; whether to print summary of records in each fold.
#' Defaults to \code{interactive()}.
#'
#' @seealso \code{\link{cv_buffer}}, \code{\link{cv_knndm}} and \code{\link{cv_spatial_autocor}};
#' \code{\link{cv_plot}} to visualise, and \code{\link{cv_distance}} and \code{\link{cv_similarity}} to evaluate, the folds
#'
#' @references C. Milà, J. Mateu, E. Pebesma, and H. Meyer, Nearest Neighbour Distance Matching
#' Leave-One-Out Cross-Validation for map validation, Methods in Ecology and Evolution (2022).
#'
#' @return An object of class S3. A list of objects including:
#'     \itemize{
#'     \item{folds_list - a list containing the folds. Each fold has two vectors with the training (first) and testing (second) indices}
#'     \item{k - number of the folds}
#'     \item{size - the distance band to separated trainig and testing folds)}
#'     \item{column - the name of the column if provided}
#'     \item{presence_bg - whether this was treated as presence-background data}
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
#' nndm <- cv_nndm(x = pa_data,
#'                 column = "occ", # optional
#'                 r = covar,
#'                 size = 350000, # size in metres no matter the CRS
#'                 num_sample = 10000,
#'                 sampling = "regular",
#'                 min_train = 0.1)
#'
#' }
cv_nndm <- function(
        x,
        column = NULL,
        r = NULL,
        size,
        pred_points = NULL,
        model_domain = NULL,
        num_sample = 1e4,
        sampling = "random",
        min_train = 0.05,
        presence_bg = FALSE,
        add_bg = FALSE,
        n_bins = 4L,
        plot = interactive(),
        report = interactive()
){

    # check for sampling arg
    sampling <- match.arg(sampling, choices = c("random", "regular"))

    # check x is an sf object
    x <- .check_x(x)
    # is column in x?
    column <- .check_column(column, x)
    # CRS: fall back to planar (Euclidean) distances when undefined
    if(is.na(sf::st_crs(x))){
        warning("CRS undefined; assuming planar coordinates, distances in coordinate units")
    }

    # if a raster is provided, it must share x's CRS
    if(!is.null(r)){
        r <- .check_r(r)
        if(!.same_crs(x, terra::crs(r))){
            stop("The coordinate reference systems of 'x' and 'r' must match.")
        }
    }

    # resolve the prediction points from pred_points, model_domain, or a sampled raster
    rx <- .knndm_predpoints(x, r, pred_points, model_domain, num_sample, sampling)
    # convert to sf object
    rx <- sf::st_as_sf(rx)

    if(presence_bg){
        unqsp <- unique(x[, column, drop = TRUE])
        if(!is.numeric(unqsp) || any(unqsp < 0) || any(unqsp > 1)){
            stop("Presence-background option is only for species data with 0s (backgrounds/pseudo-absences) and 1s (presences).\n", "The data should be numeric.\n")
        }
        # indices of presences
        x_1s <- which(x[, column, drop = TRUE] == 1)
    } else{
        x_1s <- 1:nrow(x)
    }

    # the k
    n <- length(x_1s)

    Gij <- sf::st_distance(rx, x[x_1s, ])
    units(Gij) <- NULL
    Gij <- apply(Gij, 1, min)

    tdist <- sf::st_distance(x[x_1s, ])
    units(tdist) <- NULL
    full_distmat <- tdist # for calculating indices
    diag(tdist) <- NA
    Gj <- apply(tdist, 1, function(i) min(i, na.rm=TRUE))
    Gjstar <- as.numeric(Gj)
    # training points CDF
    tp_cdf <- stats::ecdf(Gjstar)

    # run nndm in cpp
    distmat <- nndm_cpp(
        X = tdist,
        Gij = Gij,
        phi = size,
        min_train = min_train
    )

    msize <- apply(distmat, 1, function(x) min(x, na.rm=TRUE))

    # get a complete dist matrix for PBG
    if(presence_bg){
        full_distmat <- sf::st_distance(x)
        units(full_distmat) <- NULL
    }

    # add background only if both true
    add_bg <- (presence_bg && add_bg)
    # Note: in presence-background, length of full-matrix is longer than msize
    fold_list <- lapply(1:n, function(i, pbag = add_bg){
        if(pbag){
            test_ids <- which(full_distmat[x_1s[i], ] <= msize[i])
            inside <- x[test_ids, column, drop = TRUE]
            test_set <- test_ids[which(inside == 0)]
            test_set <- c(i, test_set)
        } else{
            test_set <- i
        }
        list(as.numeric(which(full_distmat[x_1s[i], ] > msize[i])),
             as.numeric(test_set))
    }
    )
    # calculate train test table summary
    train_test_table <- .table_summary(fold_list, x, column, n, n_bins = n_bins)
    if(report){
        print(summary(train_test_table)[c(1,4,6), ])
    }

    # range for plotting
    r_range <- seq(0, size, length.out = 200)
    # prediction points CDF
    pp_cdf <- stats::ecdf(Gij)
    # NNDM CDF
    fp_cdf <- stats::ecdf(msize)
    # plotting datta
    plot_data <- data.frame(pr = pp_cdf(r_range),
                            lo = tp_cdf(r_range),
                            nn = fp_cdf(r_range),
                            r = r_range)

    plt <- ggplot2::ggplot(plot_data, ggplot2::aes(x = get("r"))) +
        ggplot2::geom_step(alpha = 0.7, linewidth = 1.2, ggplot2::aes(y = get("pr"), color = "Prediciton")) +
        ggplot2::geom_step(alpha = 0.7, linewidth = 0.6, ggplot2::aes(y = get("lo"), color = "LOO")) +
        ggplot2::geom_step(alpha = 0.7, linewidth = 1.2, ggplot2::aes(y = get("nn"), color = "NNDM")) +
        ggplot2::scale_color_manual(values = c("Prediciton" = "#000000",
                                               "LOO" = "#56B4E9",
                                               "NNDM" = "#E69F00")) +
        ggplot2::labs(color = "", x = "r", y = expression(G[r])) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.text = ggplot2::element_text(size = 12))

    if(plot) plot(plt)

    final_objs <- list(
        folds_list = fold_list,
        k = n,
        column = column,
        size = size,
        plot = plt,
        presence_bg = presence_bg,
        records = train_test_table
    )

    class(final_objs) <- c("cv_nndm")
    return(final_objs)
}

#' @export
#' @method print cv_nndm
print.cv_nndm <- function(x, ...){
    print(class(x))
}


#' @export
#' @method plot cv_nndm
plot.cv_nndm <- function(x, y, data = NULL, ...){
    if(!missing(y) || !is.null(data)){
        return(.plot_cv_fold_map(
            cv = x,
            y = if(!missing(y)) y else NULL,
            data = data,
            has_y = !missing(y),
            ...
        ))
    }

    plot(x$plot)
    invisible(x$plot)
}


#' @export
#' @method summary cv_nndm
summary.cv_nndm <- function(object, ...){
    cat("Summary of number of recoreds in each training and testing fold:\n")
    if(!is.null(object$records)){
        print(summary(object$records)[c(1,4,6), ])
    }
}
