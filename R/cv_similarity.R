#' Compute similarity measures to evaluate possible extrapolation in testing folds
#'
#' This function supports three similarity measurement methods:
#' Multivariate Environmental Similarity Surface (MESS), Manhattan distance (L1),
#' and Euclidean distance (L2). It is primarily designed for evaluating the similarity
#' between training and test datasets in environmental or spatial modeling applications.
#'
#'
#' The function support three method of similarity measures, the Multivariate Environmental
#' Similarity Surface (MESS), Manhattan distance (L1), and Euclidean distance (L2).
#'
#' The MESS is calculated as described in Elith et al. (2010). MESS represents
#' how similar a point in a testing fold is to a training fold (as a reference
#' set of points), with respect to a set of predictor variables in \code{r}.
#' The negative values are the sites where at least one variable has a value that is outside
#' the range of environments over the reference set, so these are novel environments.
#'
#' When using the L1 (Manhattan) or L2 (Euclidean) distance options (experimental), the
#' function performs the following steps for each test sample:
#'
#' \itemize{
#' \item{1. Calculates the minimum distance between each test sample and all training samples
#'    in the same fold using the selected metric (L1 or L2).}
#' \item{2. Calculates a baseline distance: the average of the minimum distances between a set
#'    of random background samples (defined by \code{num_sample}) and all training/test
#'    samples combined.}
#' \item{3. Computes a similarity score by subtracting the test sample’s minimum distance from
#'    the baseline average. A higher score indicates the test sample is more similar to
#'    the training data, while lower or negative scores indicate novelty.}
#' }
#'
#' This provides a simple, distance-based novelty metric, useful for assessing
#' extrapolation or dissimilarity in prediction scenarios.
#'
#' @inheritParams cv_plot
#' @param x a simple features (sf) or SpatialPoints object of the spatial sample data used for creating
#' the \code{cv} object.
#' @param r a terra SpatRaster object of environmental predictor that are going to be used for modelling. This
#' is used to calculate similarity between the training and testing points.
#' @param method the similarity method including: MESS, L1 and L2. Read the details section.
#' @param num_sample number of random samples from raster to calculate similarity distances (only for L1 and L2).
#' @param num_plot a vector of indices of folds.
#' @param jitter_width numeric; the width of jitter points.
#' @param points_size numeric; the size of points.
#' @param points_alpha numeric; the opacity of points
#' @param points_colors character; a character vector of colours for points
#' @inheritParams cv_spatial
#'
#' @seealso \code{\link{cv_spatial}}, \code{\link{cv_cluster}}, \code{\link{cv_buffer}}, and \code{\link{cv_nndm}}
#'
#' @references Elith, J., Kearney, M., & Phillips, S. (2010). The art of modelling range-shifting species: The art of modelling range-shifting species. Methods in Ecology and Evolution, 1(4), 330–342.
#'
#' @return a ggplot object
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
#' path <- system.file("extdata/au/", package = "blockCV")
#' files <- list.files(path, full.names = TRUE)
#' covars <- terra::rast(files)
#'
#' # hexagonal spatial blocking by specified size and random assignment
#' sb <- cv_spatial(x = pa_data,
#'                  column = "occ",
#'                  size = 450000,
#'                  k = 5,
#'                  iteration = 1)
#'
#' # compute extrapolation
#' cv_similarity(cv = sb, r = covars, x = pa_data)
#'
#' }
cv_similarity <- function(
        cv,
        x,
        r,
        num_plot = seq_along(cv$folds_list),
        method = "MESS",
        num_sample = 10000L,
        jitter_width = 0.1,
        points_size = 2,
        points_alpha = 0.7,
        points_colors = NULL,
        progress = TRUE
){

    # check required packages
    pkg <- c("ggplot2")
    .check_pkgs(pkg)
    # check x is an sf object
    x <- .check_x(x)
    # change the r to terra object
    r <- .check_r(r)
    # check points fall within the raster extent
    .check_within(x, r)

    # check cv obj
    .check_cv(cv)

    method <- match.arg(tolower(method), choices = c("mess", "l1", "l2"))
    MESS <- method == "mess"
    L1 <- method == "l1"

    # The iteration must be a natural number
    tryCatch(
        {
            num_plot <- abs(as.integer(num_plot))
            num_plot <- sort(num_plot)
        },
        error = function(cond) {
            message("'num_plot' must be a natural numbers.")
        }
    )

    # get the folds list
    folds_list <- cv$folds_list
    # length of the folds
    k <- length(folds_list)
    if(max(num_plot) > k){
        num_plot <- num_plot[num_plot <= k]
    }
    # extract the raster values
    points <- terra::extract(r, x, ID = FALSE)
    # to set as nrow for df; cv_buffer has only one target points unless P-BG
    n <- nrow(points)
    if(.is_loo(cv)){
        n <- ifelse(cv$presence_bg, nrow(points), 1)
    }
    # number of predictors
    m <- ncol(points)
    df <- data.frame(id = seq_len(n))
    # add progress bar
    if(progress) pb <- utils::txtProgressBar(min = 0, max = length(num_plot), style = 3)

    if (!MESS) {
        # scale a dataset based on params of another scaled dataset
        f <- function(x, s) {
            mns <- attr(s,"scaled:center")
            sds <- attr(s,"scaled:scale")
            for (i in names(mns)) {
                x[[i]] <- (x[[i]] - mns[i]) / sds[i]
            }
            attr(x, "scaled:center") <- NULL
            attr(x, "scaled:scale") <- NULL
            return(x)
        }
        rand_points <- terra::spatSample(r, size = num_sample, na.rm = TRUE)
        rand_points <- rand_points[stats::complete.cases(rand_points), ] # just make sure..
        # scale and centre the random points
        rand_points <- as.matrix(scale(rand_points))
        # scale and centre the samples points based on the same stats
        points <- as.matrix(f(x = points, s = rand_points))
        # remove the attributes
        attr(rand_points, "scaled:center") <- NULL
        attr(rand_points, "scaled:scale") <- NULL
    }

    # calculate MESS for testing data
    for(i in num_plot){
        df[, paste("Fold", i, sep = "")] <- NA
        train <- folds_list[[i]][[1]]
        test <- folds_list[[i]][[2]]
        if (MESS) {
            mes <- sapply(1:m, function(j) .messi3(points[test, j], points[train, j]))
            if(.is_loo(cv)){
                mmes <- min(mes)
            } else{
                mmes <- apply(mes, 1, min, na.rm = TRUE)
            }
        } else {
            mmes <- similarity_cpp(
                train_mat = points[train, ],
                test_mat = points[test,, drop=FALSE],
                rand_mat = rand_points,
                L1 = L1
            )
        }
        df[1:length(mmes), paste("Fold", i, sep = "")] <- mmes
        if(progress) utils::setTxtProgressBar(pb, i)
    }
    fold_names <- paste("Fold", num_plot, sep = "")
    # reshape for plotting
    mes_reshp <- stats::reshape(
        df,
        direction = "long",
        idvar = "id",
        varying = fold_names,
        times = fold_names,
        v.names =  "value",
        timevar = "folds"
    )

    # remove NAs
    mes_reshp <- mes_reshp[stats::complete.cases(mes_reshp), ]
    if(.is_loo(cv)) mes_reshp$folds <- as.numeric(substr(mes_reshp$folds, 5, 25))
    # get the max value for color legend
    maxabs <- max(abs(mes_reshp$value))
    # define point colors
    cols <- c("#D53E4F", "#FC8D59", "#FEE08B", "#FFFFBF", "#E6F598", "#99D594", "#3288BD")
    points_colors <- if(is.null(points_colors)) cols else points_colors

    # provide alternatives for class(cv)
    goem_buffer <- ggplot2::geom_point(size = points_size, alpha = points_alpha)
    geom_other <- ggplot2::geom_jitter(width = jitter_width, size = points_size, alpha = points_alpha)
    geom_vio <- ggplot2::geom_violin(
        ggplot2::aes(x = get("folds"), y = get("value"), group = get("folds")),
        inherit.aes = FALSE,
        fill = NA
    )
    # which geom to choose
    geom_exta <- if(.is_loo(cv)) goem_buffer else geom_other
    mes_reshp$folds <- factor(mes_reshp$folds)

    col_name <- switch(
        method,
        mess = "MESS",
        l1 = "L1 distance",
        l2 = "L2 distance"
    )
    y_name <- switch(
        method,
        mess = "MESS Values",
        l1 = "Distance differnce with random samples",
        l2 = "Distance differnce with random samples"
    )

    p1 <- ggplot2::ggplot(
        data = mes_reshp,
        ggplot2::aes(x = get("folds"), y = get("value"), colour = get("value"))) +
        ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
        geom_exta +
        switch(!.is_loo(cv), geom_vio) +
        ggplot2::scale_color_gradientn(colours = points_colors,
                                       limits = c(-maxabs, maxabs),
                                       na.value = "#BEBEBE03") +
        ggplot2::labs(x = "Folds", y = y_name, col = col_name) +
        ggplot2::theme_bw()

    return(p1)
}


# calculating mess for cv_similarity
# function borrowed from dismo:::.messi3
.messi3 <- function(p, v) {
    # seems 2-3 times faster than messi2
    v <- stats::na.omit(v)
    f <- 100*findInterval(p, sort(v)) / length(v)
    minv <- min(v)
    maxv <- max(v)
    res <- 2*f
    f[is.na(f)] <- -99
    i <- f>50 & f<100
    res[i] <- 200-res[i]

    i <- f==0
    res[i] <- 100*(p[i]-minv)/(maxv-minv)
    i <- f==100
    res[i] <- 100*(maxv-p[i])/(maxv-minv)
    res
}
