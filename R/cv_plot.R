#' Visualising folds created by blockCV in ggplot
#'
#' This function visualises the folds create by blockCV. It also accepts a raster
#' layer to be used as background.
#'
#' @param cv a blockCV cv_* object; a \code{cv_spatial}, \code{cv_cluster} or \code{cv_buffer}
#' @param x a simple features (sf) or SpatialPoints object of the spatial sample data used for creating
#' the \code{cv} object. This could be empty when \code{cv} is a \code{cv_spatial} object.
#' @param r a terra SpatRaster object (optional). If provided, it will be used as background of the plots.
#' It also supports \emph{stars}, \emph{raster}, or path to a raster file on disk.
#' @param nrow integer; number of rows for facet plot
#' @param ncol integer; number of columns for facet plot
#' @param num_plots a vector of indices; for a cv_* with <= 10 folds it shows all the folds (e.g. \code{1:10}),
#' if the number of folds is higher, 10 folds are shown randomly. You can choose any of folds to be shown
#' e.g. \code{1:3} or \code{c(2, 7, 16, 22)}
#' @param max_pixels integer; maximum number of pixels used for plotting \code{r}
#' @param raster_colors a character vector of colours for raster background e.g. \code{terrain.colors(20)}
#' @param points_colors two colours to be used for train and test points
#' @param points_alpha the opacity of points
#' @param label_size integer; size of fold labels when a \code{cv_spatial} object is used.
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' library(blockCV)
#'
#' # import presence-absence species data
#' points <- read.csv(system.file("inst/extdata/", "species.csv", package = "blockCV"))
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
cv_plot <- function(
    cv,
    x,
    r = NULL,
    nrow = NULL,
    ncol = NULL,
    num_plots = sample(length(cv$folds_list), min(length(cv$folds_list), 10)),
    max_pixels = 3e5,
    raster_colors = gray.colors(10, alpha = 1),
    points_colors = c("#E69F00", "#56B4E9"),
    points_alpha = 0.7,
    label_size = 4
){
  # check for availability of ggplot2
  pkg <- c("ggplot2")
  .pkg_check(pkg)

  if(!class(cv) %in% c("cv_spatial", "cv_cluster", "cv_buffer")){
    stop("'cv' must be a blockCV cv_* object.")
  }

  # check x is an sf object
  if(!missing(x)){
    x <- .x_check(x)
  }

  # change the r to terra object
  if(!is.null(r)){
    r <- .r_check(r)
    r <- r[[1]]
  }

  # make geom_tile for raster plots
  if(!is.null(r)){
    map_df <- terra::spatSample(r,
                                size = max_pixels,
                                method="regular",
                                xy=TRUE,
                                na.rm=TRUE)
    colnames(map_df) <- c("x", "y", "value")

    geom_rast <- ggplot2::geom_tile(data = map_df,
                                    ggplot2::aes_string(x="x", y="y", fill="value"))
    geom_rast_col <- ggplot2::scale_fill_gradientn(colours = raster_colors)
  }
  # make geom_sf for spatial blocks
  if(class(cv) == "cv_spatial"){
    blocks <- cv$blocks
    geom_poly <- ggplot2::geom_sf(data = sf::st_geometry(blocks),
                                  inherit.aes = FALSE,
                                  colour = "red",
                                  fill = "orangered4",
                                  alpha = 0.04,
                                  size = 0.2)
  }

  if(!missing(x)){
    x_long <- .x_to_long(x, cv, num_plot = num_plots)
  } else{
    # stop if x is missing for buffer and cluster
    if(class(cv) != "cv_spatial") stop("'x' is required for plotting cv_cluster and cv_buffer.")
  }

  if(missing(x)){
    if(class(cv) == "cv_spatial"){

      p1 <- ggplot2::ggplot(data = blocks) +
        switch(!is.null(r), geom_rast, NULL) +
        switch(!is.null(r), geom_rast_col, NULL) +
        ggplot2::geom_sf(colour = "red",
                         fill = "orangered4",
                         alpha = 0.04,
                         size = 0.2) +
        ggplot2::geom_sf_text(ggplot2::aes_string(label = "folds"),
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
      switch(class(cv) == "cv_spatial", geom_poly, NULL) +
      ggplot2::geom_sf(ggplot2::aes_string(col = "value"),
                       alpha = points_alpha) +
      ggplot2::scale_color_manual(values = points_colors, na.value = "#BEBEBE03") +
      ggplot2::facet_wrap(~folds, nrow = nrow, ncol = ncol) +
      ggplot2::labs(x = "", y = "", col = "") + # set the axes labes to NULL
      ggplot2::theme_bw() +
      ggplot2::guides(fill = "none")

  }

  return(p1)
}
