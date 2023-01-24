#' Use the Nearest Neighbour Distance Matching (NNDM) to separate train and test folds
#'
#' @param x
#' @param column
#' @param r
#' @param size
#' @param num_sample
#' @param sampling
#' @param min_train
#' @param presence_background
#' @param add_background
#' @param plot
#' @param progress
#'
#' @return
#' @export
#'
#' @examples
#'
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
    progress = TRUE,
    print = TRUE
){

  sampling <- c("random", "regular")[1]

  # check x is an sf object
  x <- .x_check(x)
  # is column in x?
  if(!is.null(column)){
    if(!column %in% colnames(x)){
      warning(sprintf("There is no column named '%s' in 'x'.\n", column))
      column <- NULL
    }
  }

  # check r
  r <- .r_check(r)
  r <- r[[1]]

  # get the sample points form raster
  rx <- terra::spatSample(x = r,
                          size = num_sample,
                          method = sampling,
                          values = FALSE,
                          na.rm = TRUE,
                          as.points = TRUE)

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
  jmin <- which.min(Gjstar)[1]
  kmin <- which(tdist[jmin, ] == rmin)

  # run nndm in cpp
  distmat <- nndm_cpp(
    X = tdist,
    Gjstar = Gjstar,
    Gij = Gij,
    rmin = rmin,
    jmin = jmin,
    kmin = as.numeric(kmin),
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

  plot(r_range, pp_cdf(r_range), type = "l")
  points(r_range, tp_cdf(r_range), type = "l", col = "red")
  points(r_range, fp_cdf(r_range), type = "l", col = "blue")

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

  final_objs <- list(
    folds_list = fold_list,
    k = n,
    column = column,
    size = size,
    presence_background = presence_background,
    records = if(print) train_test_table else NULL
  )

  class(final_objs) <- c("cv_nndm")
  return(final_objs)
}

