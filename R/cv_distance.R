#' Compare a cross-validation design to the prediction domain via nearest-neighbour distances
#'
#' A diagnostic that works on \emph{any} \code{blockCV} cross-validation object (the output of
#' \code{\link{cv_spatial}}, \code{\link{cv_cluster}}, \code{\link{cv_buffer}}, \code{\link{cv_nndm}}
#' or \code{\link{cv_knndm}}). It draws the same nearest-neighbour distance distribution plot produced
#' internally by \code{\link{cv_nndm}} and \code{\link{cv_knndm}}, but for a fold configuration that has
#' \emph{already} been generated. This lets you check whether a given fold-generation strategy actually
#' pushes the test-to-train distances towards the distances the model will face when predicting over the
#' target area. It does not compute pairwise distances for general use; it is a fold-quality diagnostic.
#'
#' Three (optionally four) empirical cumulative distribution functions of nearest-neighbour distances
#' are compared over a distance \code{r}:
#' \itemize{
#'   \item{\strong{Prediction} (\eqn{\hat{G}_{ij}}) - from each prediction point to its nearest sample
#'   point. This is the distance regime the model meets at prediction time and is the target to match.}
#'   \item{\strong{LOO} (\eqn{\hat{G}_{j}}) - from each sample point to its nearest \emph{other} sample
#'   point. This is the nearest-neighbour (leave-one-out) limit, i.e. the most optimistic, random-like case.}
#'   \item{\strong{CV} (\eqn{\hat{G}_{j}^{*}}) - from each test point to its nearest \emph{training} point
#'   \emph{within the same fold} of the supplied \code{cv} object. The closer this curve sits to the
#'   Prediction curve, the more the fold design mimics the prediction task.}
#'   \item{\strong{Random} - the CV curve of a random \code{k}-fold split with the same number of folds
#'   as \code{cv}, drawn only when \code{add_random = TRUE}. It shows where a naive random assignment lands,
#'   so the improvement of the supplied design is visible. For leave-one-out objects (\code{cv_buffer},
#'   \code{cv_nndm}) a random split coincides with the LOO curve, so it is omitted.}
#' }
#'
#' The gap between each curve and the Prediction curve is summarised by the Wasserstein-1 distance (the
#' area between the two empirical distribution functions; lower means a closer match). These values are
#' shown in the subtitle and attached to the returned plot as \code{attr(p, "W")}.
#'
#' When \code{space = "feature"} the distances are computed in the (optionally scaled) covariate space of
#' \code{r} instead of the geographical space, mirroring \code{\link{cv_knndm}}.
#'
#' When the supplied \code{cv} object was built with \code{presence_bg = TRUE} (in \code{\link{cv_spatial}},
#' \code{\link{cv_cluster}}, \code{\link{cv_knndm}}, \code{\link{cv_buffer}} or \code{\link{cv_nndm}}), all
#' three distance distributions are computed on the \emph{presences} only: the prediction, LOO and CV curves
#' are expressed relative to the presence-only training data rather than the (often random) \emph{background}
#' points (locations sampled across the study area to represent the available conditions rather than confirmed
#' absences), matching the nearest neighbour distance matching framing. This is read automatically from \code{cv}.
#'
#' @inheritParams cv_knndm
#' @param cv a \code{blockCV} cross-validation object, i.e. the output of \code{\link{cv_spatial}},
#' \code{\link{cv_cluster}}, \code{\link{cv_buffer}}, \code{\link{cv_nndm}} or \code{\link{cv_knndm}}.
#' @param x a simple features (sf) object of the spatial sample points used to create the \code{cv} object.
#' @param r a terra SpatRaster object. It defines the area the model predicts over; when neither
#' \code{pred_points} nor \code{model_domain} is supplied, prediction points are sampled from it. It is
#' also required (for the covariates) when \code{space = "feature"}. One of \code{r}, \code{pred_points},
#' or \code{model_domain} is required.
#' @param add_random logical; when \code{TRUE} (default), overlay the CV curve of a random \code{k}-fold
#' split with the same number of folds as \code{cv} (ignored for leave-one-out objects, see details).
#' @param num_random integer; the number of random \code{k}-fold splits used to estimate the random baseline
#' when \code{add_random = TRUE}. The plotted random curve is the mean across these splits, with a 10-90%
#' band.
#' @param plot logical; whether to draw the plot. The plot is always returned invisibly.
#'
#' @seealso \code{\link{cv_similarity}}, \code{\link{cv_knndm}}, \code{\link{cv_nndm}},
#' \code{\link{cv_spatial}}, \code{\link{cv_cluster}}, \code{\link{cv_buffer}}, and \code{\link{cv_plot}} to visualise the folds
#'
#' @references Milà, C., Mateu, J., Pebesma, E., & Meyer, H. (2022). Nearest neighbour distance matching
#' leave-one-out cross-validation for map validation. Methods in Ecology and Evolution, 13(6), 1304-1316.
#'
#' Linnenbrink, J., Milà, C., Ludwig, M., & Meyer, H. (2024). kNNDM CV: k-fold nearest neighbour distance
#' matching cross-validation for map accuracy estimation. Geoscientific Model Development, 17(15), 5897-5912.
#'
#' @return a ggplot object. The Wasserstein-1 distances of each curve to the Prediction curve are attached
#' as \code{attr(p, "W")} (a named numeric vector; lower means a closer match to the prediction domain).
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
#' # generate spatial folds
#' sb <- cv_spatial(x = pa_data, column = "occ", size = 450000, k = 5, iteration = 1)
#'
#' # assess how close the folds are to the prediction domain
#' cv_distance(cv = sb, x = pa_data, r = covar, num_sample = 5000)
#' }
cv_distance <- function(
        cv,
        x,
        r = NULL,
        pred_points = NULL,
        model_domain = NULL,
        space = "geographical",
        add_random = TRUE,
        num_random = 10L,
        num_sample = 10000L,
        sampling = "regular",
        scale = TRUE,
        seed = NULL,
        plot = TRUE
){
    # check required packages
    .check_pkgs("ggplot2")
    space <- match.arg(space, choices = c("geographical", "feature"))
    sampling <- match.arg(sampling, choices = c("regular", "random"))

    # check the cv and sample objects
    .check_cv(cv)
    x <- .check_x(x)
    # the fold indices in 'cv' must line up with the supplied sample points
    .check_x_matches_cv(x, cv)
    if(is.na(sf::st_crs(x))){
        stop("The coordinate reference system of 'x' must be defined.")
    }

    # feature space requires the covariate raster
    if(space == "feature" && is.null(r)){
        stop("'r' (covariate raster) is required when space = 'feature'.")
    }
    if(!is.null(r)){
        r <- .check_r(r)
        if(!.same_crs(x, terra::crs(r))){
            stop("The coordinate reference systems of 'x' and 'r' must match.")
        }
    }

    if(!is.null(seed)) set.seed(seed)

    n <- nrow(x)

    # resolve the prediction points -------------------------------------------
    predpts <- .knndm_predpoints(x, r, pred_points, model_domain, num_sample, sampling)

    # presence-background objects: express every distance relative to the presences (1s)
    pbg <- isTRUE(cv$presence_bg) && !is.null(cv$column) && cv$column %in% colnames(x)
    x_1s <- if(pbg) which(x[, cv$column, drop = TRUE] == 1) else seq_len(n)
    xp <- if(pbg) x[x_1s, ] else x
    np <- length(x_1s)
    # map full-data indices (as stored in folds_list) to the presence-only rows
    full_to_p <- match(seq_len(n), x_1s)

    # sample-to-sample distances and prediction-to-sample distances -----------
    if(space == "geographical"){
        tdist <- sf::st_distance(xp)
        units(tdist) <- NULL
        Gij <- sf::st_distance(predpts, xp)
        units(Gij) <- NULL
        Gij <- apply(Gij, 1, min)
    } else{
        ft <- .knndm_features(xp, r, predpts, scale)
        tdist <- as.matrix(stats::dist(ft$train))
        Gij <- .nn_cross(ft$pred, ft$train)
    }

    # sample-to-nearest-other-sample (LOO) distances
    Gj <- vapply(seq_len(np), function(i) min(tdist[i, -i]), numeric(1))

    # test-to-nearest-train distances of the supplied folds (Gjstar) ----------
    is_loo <- .is_loo(cv)
    Gjstar <- unlist(lapply(cv$folds_list, function(f){
        train <- f[[1]]
        test <- f[[2]]
        # presence-background: keep only the presences, mapped to the presence-only rows
        if(pbg){
            train <- full_to_p[train]; train <- train[!is.na(train)]
            test  <- full_to_p[test];  test  <- test[!is.na(test)]
        }
        if(!length(train) || !length(test)) return(numeric(0))
        apply(tdist[test, train, drop = FALSE], 1, min)
    }))

    # optional random k-fold curve(s) with matched number of folds -----------
    show_random <- isTRUE(add_random) && !is_loo
    num_random <- max(1L, as.integer(num_random))
    rand_draws <- NULL
    if(show_random){
        k <- length(cv$folds_list)
        rand_draws <- lapply(seq_len(num_random), function(j){
            rv <- sample(rep(seq_len(k), ceiling(np / k)), size = np)
            .nn_diff_fold(tdist, rv)
        })
    } else if(isTRUE(add_random) && is_loo){
        message("A random split of a leave-one-out object coincides with the LOO curve; 'Random' curve omitted.")
    }

    # Wasserstein-1 distances to the prediction distribution (lower is closer)
    W <- c(CV = .wasserstein(Gjstar, Gij), LOO = .wasserstein(Gj, Gij))
    if(show_random){
        W <- c(W, Random = mean(vapply(rand_draws, .wasserstein, numeric(1), b = Gij)))
    }

    # distance distribution functions for plotting ----------------------------
    max_r <- max(c(Gij, Gjstar, Gj, unlist(rand_draws)))
    r_range <- seq(0, max_r, length.out = 200)
    plot_data <- data.frame(r = r_range,
                            Prediction = stats::ecdf(Gij)(r_range),
                            LOO = stats::ecdf(Gj)(r_range),
                            CV = stats::ecdf(Gjstar)(r_range))
    if(show_random){
        # ecdf of each random draw on the shared grid: mean line + 10-90% band
        M <- vapply(rand_draws, function(g) stats::ecdf(g)(r_range), numeric(length(r_range)))
        plot_data$Random <- rowMeans(M)
        plot_data$Random_lo <- apply(M, 1, stats::quantile, probs = 0.1)
        plot_data$Random_hi <- apply(M, 1, stats::quantile, probs = 0.9)
    }

    cols <- c(Prediction = "#000000", LOO = "#56B4E9", Random = "#009E73", CV = "#E69F00")
    present <- c("Prediction", "LOO", if(show_random) "Random", "CV")
    w_lab <- paste(sprintf("%s: %.3g", names(W), W), collapse = "   ")

    plt <- ggplot2::ggplot(plot_data, ggplot2::aes(x = get("r"))) +
        ggplot2::geom_step(alpha = 0.7, linewidth = 1.2, ggplot2::aes(y = get("Prediction"), color = "Prediction")) +
        ggplot2::geom_step(alpha = 0.7, linewidth = 0.6, ggplot2::aes(y = get("LOO"), color = "LOO"))
    if(show_random){
        plt <- plt +
            ggplot2::geom_ribbon(ggplot2::aes(ymin = get("Random_lo"), ymax = get("Random_hi")),
                                 fill = "#009E73", alpha = 0.2) +
            ggplot2::geom_step(alpha = 0.7, linewidth = 0.8, linetype = "dashed",
                               ggplot2::aes(y = get("Random"), color = "Random"))
    }
    plt <- plt +
        ggplot2::geom_step(alpha = 0.7, linewidth = 1.2, ggplot2::aes(y = get("CV"), color = "CV")) +
        ggplot2::scale_color_manual(values = cols, breaks = present) +
        ggplot2::labs(color = "", x = "r", y = expression(G[r]),
                      subtitle = paste0("Wasserstein-1 to prediction (lower = closer)   ", w_lab),
                      caption = if(pbg) "Presence-background object: distances computed on presences only" else NULL) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.text = ggplot2::element_text(size = 12))

    attr(plt, "W") <- W
    if(plot) plot(plt)
    invisible(plt)
}

