#' Use the Nearest Neighbour Distance Matching (NNDM) to separate train and test folds
#'
#' A fast implementation of the Nearest Neighbour Distance Matching (NNDM) algorithm (Milà et al., 2022) in C++. Similar
#' to \code{\link{cv_buffer}}, this is a variation of leave-one-out (LOO) cross-validation. It tries to match the
#' nearest neighbour distance distribution function between the test and training data to the nearest neighbour
#' distance distribution function between the target prediction and training points (Milà et al., 2022).
#'
#' When working with presence-background (presence and pseudo-absence) species distribution
#' data (should be specified by \code{presence_background = TRUE} argument), only presence records are used
#' for specifying the folds (recommended). The testing fold comprises only the target \emph{presence} point (optionally,
#' all background points within the distance are also included when \code{add_background = TRUE}; this is the
#' distance that matches the nearest neighbour distance distribution function of training-testing presences and
#' training-presences and prediction points; often lower than \code{size}).
#' Any non-target presence points inside the distance are excluded.
#' All points (presence and background) outside of distance are used for the training set.
#' The methods cycles through all the presence data, so the number of folds is equal to
#' the number of presence points in the dataset.
#'
#' For all other types of data (including presence-absence, count, continuous, and multi-class)
#' set \code{presence_background = FALE}, and the function behaves similar to the methods
#' explained by Milà and colleagues (2022).
#'
#' @param x a simple features (sf) or SpatialPoints object of spatial sample data (e.g., species
#' data or ground truth sample for image classification).
#' @inheritParams cv_buffer
#' @param r a terra SpatRaster object of a predictor variable. This defines the area that model is going to predict.
#' @param size numeric value of the range of spatial autocorrelation (the \code{phi} parameter).
#' This distance should be in \strong{metres}. The range could be explored by \code{\link{cv_spatial_autocor}}.
#' @param num_sample integer; the number of sample points from predictor (\code{r}) to be used for calculating
#' the G function of prediction points.
#' @param sampling either \code{"random"} or \code{"regular"} for sampling prediction points.
#' When  \code{sampling = "regular"}, the actual number of samples might be less than \code{num_sample}
#' for non-rectangular rasters (points falling on no value areas are removed).
#' @param min_train numeric; between 0 and 1. A constraint on the minimum proportion of train points in each fold.
#' @inheritParams cv_buffer
#' @inheritParams cv_buffer
#' @param plot logical; whether to plot the G functions in ggplot.
#' @inheritParams cv_buffer
#'
#' @seealso \code{\link{cv_buffer}} and \code{\link{cv_spatial_autocor}}
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
#'     \item{presence_background - whether this was treated as presence-background data}
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
#' path <- system.file("extdata/au/", package = "blockCV")
#' files <- list.files(path, full.names = TRUE)
#' rasters <- terra::rast(files)
#'
#' nndm <- cv_nndm(x = pa_data,
#'                 column = "occ", # optional
#'                 r = rasters[[1]],
#'                 size = 350000, # size in metres no matter the CRS
#'                 num_sample = 10000,
#'                 sampling = "regular",
#'                 min_train = 0.1,)
#'
#' }
cv_nndm <- function(
    x,
    column = NULL,
    r,
    size,
    num_sample = 1e4,
    sampling = "random",
    min_train = 0.05,
    presence_background = FALSE,
    add_background = FALSE,
    plot = TRUE,
    print = TRUE
){

  # check for sampling arg
  sampling <- match.arg(sampling, choices = c("random", "regular"))

  # check x is an sf object
  x <- .x_check(x)
  # is column in x?
  column <- .column_check(column, x)
  # if(!is.null(column)){
  #   if(!column %in% colnames(x)){
  #     warning(sprintf("There is no column named '%s' in 'x'.\n", column))
  #     column <- NULL
  #   }
  # }

  # check r
  r <- .r_check(r)
  r <- r[[1]]

  # get the sample points form raster
  tryCatch(
    {
      rx <- terra::spatSample(x = r,
                              size = num_sample,
                              method = sampling,
                              values = FALSE,
                              na.rm = TRUE,
                              as.points = TRUE)
    },
    error = function(cond) {
      message("Sampling raster failed!")
    }
  )
  # convert to sf object
  rx <- sf::st_as_sf(rx)

  if(presence_background){
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
  diag(tdist) <- NA
  Gj <- apply(tdist, 1, function(i) min(i, na.rm=TRUE))
  Gjstar <- as.numeric(Gj)

  tp_cdf <- ecdf(Gjstar)

  # Start algorithm
  rmin <- min(Gjstar)
  imin <- which.min(Gjstar)[1]
  jmin <- which(tdist[imin, ] == rmin)

  # run nndm in cpp
  distmat <- nndm_cpp(
    X = tdist,
    Gjstar = Gjstar,
    Gij = Gij,
    rmin = rmin,
    imin = imin,
    jmin = as.numeric(jmin),
    phi = size,
    min_train = min_train
  )

  msize <- apply(distmat, 1, function(x) min(x, na.rm=TRUE))

  # distance matrix by sf
  if(presence_background){
    full_distmat <- sf::st_distance(x)
    units(full_distmat) <- NULL
  } else{
    full_distmat <- tdist
  }

  # if(progress) pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
  fold_list <- lapply(1:n, function(i, pbag = presence_background){
    if(pbag){
      test_ids <- which(full_distmat[i, ] <= msize[x_1s[i]])
      inside <- x[test_ids, column, drop = TRUE]
      test_set <- test_ids[which(inside == 0)]
      test_set <- c(i, test_set)
    } else{
      test_set <- i
    }
    # if(progress) utils::setTxtProgressBar(pb, i)
    list(which(full_distmat[i, ] > msize), test_set)
  }
  )

  # calculate train test table summary
  if(print){
    train_test_table <- .ttt(fold_list, x, column, n)
    print(summary(train_test_table)[c(1,4,6), ])
  }


  r_range <- seq(0, size, length.out = 200)
  pp_cdf <- ecdf(Gij)
  # tp_cdf <- ecdf(Gjstar)
  fp_cdf <- ecdf(msize)

  plot_data <- rbind(
    data.frame(value = pp_cdf(r_range), r = r_range, name = "Prediction"),
    data.frame(value = tp_cdf(r_range), r = r_range, name = "LOO"),
    data.frame(value = fp_cdf(r_range), r = r_range, name = "NNDM")
  )

  plt <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes_string(x="r", y="value", colour="name", size="name")) +
    ggplot2::geom_step(alpha = 0.8) +
    ggplot2::scale_size_manual(values = c(0.5, 1.1, 1.1)) +
    ggplot2::scale_colour_manual(values = c("#56B4E9", "#E69F00", "#000000")) +
    ggplot2::ylab(expression(paste(G[r]))) +
    ggplot2::labs(colour = "", size = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.text.align = 0,
                   legend.text = ggplot2::element_text(size = 12))


  plot(plt)

  # final_objs <- list(
  #   folds_list = fold_list,
  #   k = n,
  #   column = column,
  #   size = size,
  #   presence_background = presence_background,
  #   records = if(print) train_test_table else NULL
  # )
  #
  # class(final_objs) <- c("cv_nndm")
  # return(final_objs)
  return(distmat)
}

