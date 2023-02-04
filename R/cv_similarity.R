#' Compute similarity measures to evaluate possible extrapolation in testing folds
#'
#' This function computes multivariate environmental similarity surface (MESS) as described
#' in Elith et al. (2010). MESS represents how similar a point in a testing fold is to a training
#' fold (as a reference set of points), with respect to a set of predictor variables in \code{r}.
#' The negative values are the sites where at least one variable has a value that is outside
#' the range of environments over the reference set, so these are novel environments.
#'
#' @inheritParams cv_plot
#' @param x a simple features (sf) or SpatialPoints object of the spatial sample data used for creating
#' the \code{cv} object.
#' @param r a terra SpatRaster object of environmental predictor that are going to be used for modelling. This
#' is used to calculate similarity between the training and testing points.
#' @param num_plot a vector of indices of folds.
#' @param jitter_width numeric; the width of jitter points.
#' @param points_size numeric; the size of points.
#' @param points_alpha numeric; the opacity of points
#' @param points_colors character; a character vector of colours for points
#' @inheritParams cv_spatial
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
cv_similarity <- function(cv,
                          x,
                          r,
                          num_plot = seq_along(cv$folds_list),
                          jitter_width = 0.1,
                          points_size = 2,
                          points_alpha = 0.7,
                          points_colors = NULL,
                          progress = TRUE){

  # check required packages
  pkg <- c("ggplot2")
  .check_pkgs(pkg)
  # check x is an sf object
  x <- .check_x(x)
  # change the r to terra object
  r <- .check_r(r)

  # check cv obj
  if(!class(cv) %in% c("cv_spatial", "cv_cluster", "cv_buffer", "cv_nndm")){
    stop("'cv' must be a blockCV cv_* object.")
  }

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
  # calculate MESS for testing data
  for(i in num_plot){
    df[, paste("Fold", i, sep = "")] <- NA
    train <- folds_list[[i]][[1]]
    test <- folds_list[[i]][[2]]
    mes <- sapply(1:m, function(j) .messi3(points[test, j], points[train, j]))
    if(.is_loo(cv)){
      mmes <- min(mes)
    } else{
      mmes <- apply(mes, 1, min, na.rm = TRUE)
    }
    df[1:length(mmes), paste("Fold", i, sep = "")] <- mmes
    if(progress) utils::setTxtProgressBar(pb, i)
  }
  fold_names <- paste("Fold", num_plot, sep = "")
  # reshape for plotting
  mes_reshp <- stats::reshape(df,
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
  geom_vio <- ggplot2::geom_violin(fill = NA)
  # which geom to choose
  geom_exta <- if(.is_loo(cv)) goem_buffer else geom_other

  p1 <- ggplot2::ggplot(
    data = mes_reshp,
    # ggplot2::aes(x = rlang::.data$folds, y = rlang::.data$value, col = rlang::.data$value)) +
    ggplot2::aes(x = get("folds"), y = get("value"), col = get("value"))) +
    ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
    geom_exta +
    switch(!.is_loo(cv), geom_vio, NULL) +
    ggplot2::scale_color_gradientn(colours = points_colors,
                                   limits = c(-maxabs, maxabs),
                                   na.value = "#BEBEBE03") +
    ggplot2::labs(x = "Folds", y = "MESS Values", col = "MESS") +
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
