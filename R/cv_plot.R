#' Visualising folds created by blockCV in ggplot
#'
#' This function visualises the folds create by blockCV. It also accepts a raster
#' layer to be used as background in the output plot.
#'
#' @param cv a blockCV cv_* object; a \code{cv_spatial}, \code{cv_cluster}, \code{cv_buffer}
#' or \code{cv_nndm}
#' @param x a simple features (sf) or SpatialPoints object of the spatial sample data used for creating
#' the \code{cv} object. This could be empty when \code{cv} is a \code{cv_spatial} object.
#' @param r a terra SpatRaster object (optional). If provided, it will be used as background of the plots.
#' It also supports \emph{stars}, \emph{raster}, or path to a raster file on disk.
#' @param nrow integer; number of rows for facet plot
#' @param ncol integer; number of columns for facet plot
#' @param num_plots a vector of indices of folds; by default the first 10 are shown (if available).
#' You can choose any of the folds to be shown e.g. \code{1:3} or \code{c(2, 7, 16, 22)}
#' @param max_pixels integer; maximum number of pixels used for plotting \code{r}
#' @param raster_colors character; a character vector of colours for raster background e.g. \code{terrain.colors(20)}
#' @param points_colors character; two colours to be used for train and test points
#' @param points_alpha numeric; the opacity of points
#' @param label_size integer; size of fold labels when a \code{cv_spatial} object is used.
#' @param remove_na logical; whether to remove excluded points in \code{cv_buffer} from the plot
#'
#' @importFrom grDevices gray.colors
#' @return a ggplot object
#' @export
#'
#' @examples
#' \donttest{
#' library(blockCV)
#'
#' # import presence-absence species data
#' points <- read.csv(system.file("extdata/", "species.csv", package = "blockCV"))
#' pa_data <- sf::st_as_sf(points, coords = c("x", "y"), crs = 7845)
#'
#' # spatial clustering
#' sc <- cv_cluster(x = pa_data, k = 5)
#'
#' # now plot the create folds
#' cv_plot(cv = sc,
#'         x = pa_data, # sample points
#'         nrow = 2,
#'         points_alpha = 0.5)
#'
#' }
cv_plot <- function(
    cv,
    x,
    r = NULL,
    nrow = NULL,
    ncol = NULL,
    num_plots = 1:10,
    max_pixels = 3e5,
    remove_na = TRUE,
    raster_colors = gray.colors(10, alpha = 1),
    points_colors = c("#E69F00", "#56B4E9"),
    points_alpha = 0.7,
    label_size = 4
){
  # check for availability of ggplot2
  pkg <- c("ggplot2")
  .check_pkgs(pkg)

  if(!class(cv) %in% c("cv_spatial", "cv_cluster", "cv_buffer", "cv_nndm")){
    stop("'cv' must be a blockCV cv_* object.")
  }

  # check x is an sf object
  if(!missing(x)){
    x <- .check_x(x)
  }

  # change the r to terra object
  if(!is.null(r)){
    r <- .check_r(r)
    r <- r[[1]]
  }

  # make geom_tile for raster plots
  if(!is.null(r)){
    map_df <- terra::spatSample(r,
                                size = max_pixels,
                                method = "regular",
                                xy = TRUE,
                                na.rm = TRUE)
    colnames(map_df) <- c("x", "y", "value")

    geom_rast <- ggplot2::geom_tile(
      data = map_df,
      # ggplot2::aes(x = rlang::.data$x, y = rlang::.data$y, fill = rlang::.data$value))
      ggplot2::aes(x = get("x"), y = get("y"), fill = get("value")))
    geom_rast_col <- ggplot2::scale_fill_gradientn(colours = raster_colors)
  }
  # make geom_sf for spatial blocks
  if(methods::is(cv, "cv_spatial")){
    blocks <- cv$blocks
    geom_poly <- ggplot2::geom_sf(data = sf::st_geometry(blocks),
                                  inherit.aes = FALSE,
                                  colour = "red",
                                  fill = "orangered4",
                                  alpha = 0.04,
                                  linewidth = 0.2)
  }

  if(!missing(x)){
    x_long <- .x_to_long(x, cv, num_plot = num_plots)
    # exclude NAs from cv_buffer
    if(.is_loo(cv) && remove_na){
      x_long <- x_long[which(stats::complete.cases(x_long$value)), ]
    }
  } else{
    # stop if x is missing for buffer and cluster
    if(!methods::is(cv, "cv_spatial")) stop("'x' is required for plotting cv_cluster and cv_buffer.")
  }

  if(missing(x)){
    if(methods::is(cv, "cv_spatial")){

      p1 <- ggplot2::ggplot(data = blocks) +
        switch(!is.null(r), geom_rast, NULL) +
        switch(!is.null(r), geom_rast_col, NULL) +
        ggplot2::geom_sf(colour = "red",
                         fill = "orangered4",
                         alpha = 0.04,
                         linewidth = 0.2) +
        ggplot2::geom_sf_text(
          ggplot2::aes(label = get("folds")),
          size = label_size, fun.geometry = sf::st_centroid) +
        ggplot2::labs(x = "", y = "") + # or set the axes labes to NULL
        ggplot2::scale_x_continuous(guide = ggplot2::guide_axis(check.overlap = TRUE)) +
        ggplot2::theme_minimal() +
        ggplot2::guides(fill = "none")

    }

  } else{

    p1 <- ggplot2::ggplot(data = x_long) +
      switch(!is.null(r), geom_rast, NULL) +
      switch(!is.null(r), geom_rast_col, NULL) +
      switch(methods::is(cv, "cv_spatial"), geom_poly, NULL) +
      # ggplot2::geom_sf(ggplot2::aes(col = rlang::.data$value),
      ggplot2::geom_sf(ggplot2::aes(col = get("value")),
                       alpha = points_alpha) +
      ggplot2::scale_color_manual(values = points_colors, na.value = "#BEBEBE03") +
      # ggplot2::facet_wrap(~ rlang::.data$folds, nrow = nrow, ncol = ncol) +
      ggplot2::facet_wrap(~get("folds"), nrow = nrow, ncol = ncol) +
      ggplot2::labs(x = "", y = "", col = "") + # set the axes labes to NULL
      ggplot2::theme_bw() +
      ggplot2::guides(fill = "none")

  }

  return(p1)
}


# is it a LOO CV object?
.is_loo <- function(x){
  methods::is(x, "cv_buffer") || methods::is(x, "cv_nndm")
}

# transform x and fold numbers for plotting
.x_to_long <- function(x, cv, num_plot=1:10){
  # get the folds list
  folds_list <- cv$folds_list
  # The iteration must be a natural number
  tryCatch(
    {
      num_plot <- abs(as.integer(num_plot))
      num_plot <- sort(num_plot)
    },
    error = function(cond) {
      message("'num_plot' must be natural numbers.")
    }
  )
  # length of the folds
  k <- length(folds_list)
  if(max(num_plot) > k){
    num_plot <- num_plot[num_plot <= k]
  }
  # get the length of unique ids
  # if(methods::is(cv, "cv_buffer")){
  if(.is_loo(cv)){
    len <- length(unique(unlist(cv$folds_list)))
  } else{
    len <- length(unlist(folds_list[[1]]))
  }
  if(len != nrow(x)){
    stop("Number of rows in 'x' does not match the folds in 'cv'!")
  }
  # create a dataframe temp
  df <- data.frame(id = seq_len(len))
  # make the indices in x
  for(i in num_plot){
    df[, paste("Fold", i, sep = "")] <- NA
    test <- folds_list[[i]][[2]]
    train <- folds_list[[i]][[1]]
    df[test, paste("Fold", i, sep = "")] <- 0
    df[train, paste("Fold", i, sep = "")] <- 1
  }
  # get the geometry column name
  sf_colname <- attr(x, "sf_column")
  # cbind x and the df with fold ids
  xf <- cbind(x, df)
  # convert to dataframe for reshaping
  x_df <- as.data.frame(xf)
  # name of columns to rehspae long
  fold_names <- paste("Fold", num_plot, sep = "")
  # reshape x-df to long
  x_reshape <- stats::reshape(x_df,
                              direction = "long",
                              idvar = "id",
                              varying = fold_names,
                              times = fold_names,
                              v.names = "value",
                              timevar = "folds"
  )
  # convert back to sf
  x_long <- sf::st_as_sf(x_reshape, sf_column_name = sf_colname)
  # convert to factor for plotting
  x_long$value <- as.factor(x_long$value)
  levels(x_long$value) <- c("Test", "Train")

  return(x_long)
}
