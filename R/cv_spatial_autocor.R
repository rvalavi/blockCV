#' Measure spatial autocorrelation in spatial response data or predictor raster files
#'
#' This function provides a quantitative basis for choosing block size. The spatial autocorrelation in either the
#' spatial sample points or all continuous predictor variables available as raster layers is assessed and reported.
#' The response (as defined be \code{column}) in spatial sample points can be binary such as species distribution data,
#' or continuous response like soil organic carbon. The function estimates spatial autocorrelation \emph{ranges} of all input
#' raster layers or the response data. This is the range over which observations are independent and is determined by
#' constructing the empirical variogram, a fundamental geostatistical tool for measuring spatial autocorrelation.
#' The empirical variogram models the structure of spatial autocorrelation by measuring variability between all possible
#' pairs of points (O'Sullivan and Unwin, 2010). Results are plotted. See the details section for further information.
#'
#' The input raster layers should be continuous for computing the variograms and estimating the range of spatial
#' autocorrelation. The input rasters should also have a specified coordinate reference system. However, if the reference
#' system is not specified, the function attempts to guess it based on the extent of the map. It assumes an un-projected
#' reference system for layers with extent lying between -180 and 180.
#'
#' Variograms are calculated based on the distances between pairs of points, so un-projected rasters (in degrees) will
#' not give an accurate result (especially over large latitudinal extents). For un-projected rasters, \emph{the great circle distance}
#' (rather than Euclidean distance) is used to calculate the spatial distances between pairs of points. To
#' enable more accurate estimate, it is recommended to transform un-projected maps (geographic coordinate
#' system / latitude-longitude) to a projected metric reference system (e.g. UTM or Lambert) where it is possible.
#' See \code{\link[automap]{autofitVariogram}} from \pkg{automap} and \code{\link[gstat]{variogram}} from \pkg{gstat} packages
#' for further information.
#'
#' @param r a terra SpatRaster object. If provided (and \code{x} is missing), it will be used for to calculate range.
#' @param x a simple features (sf) or SpatialPoints object of spatial sample data (e.g., species binary or continuous date).
#' @param column character; indicating the name of the column in which response variable (e.g. species data as a binary
#'  response i.e. 0s and 1s) is stored for calculating spatial autocorrelation range. This supports multiple column names.
#' @param num_sample integer; the number of sample points of each raster layer to fit variogram models. It is 5000 by default,
#' however it can be increased by user to represent their region well (relevant to the extent and resolution of rasters).
#' @param deg_to_metre integer. The conversion rate of degrees to metres.
#' @param plot logical; whether to plot the results.
#' @param progress logical; whether to shows a progress bar.
#' @param ... additional option for \code{\link{cv_plot}}
#'
#' @seealso \code{\link{cv_block_size}}
#'
#' @references O'Sullivan, D., Unwin, D.J., (2010). Geographic Information Analysis, 2nd ed. John Wiley & Sons.
#'
#' Roberts et al., (2017). Cross-validation strategies for data with temporal, spatial, hierarchical,
#' or phylogenetic structure. Ecography. 40: 913-929.
#'
#' @return An object of class S3. A list object including:
#'     \itemize{
#'     \item{range - the suggested range (i.e. size), which is the median of all calculated ranges in case of 'r'.}
#'     \item{range_table - a table of input covariates names and their autocorrelation range}
#'     \item{plots - the output plot (the plot is shown by default)}
#'     \item{num_sample - number sample of 'r' used for analysis}
#'     \item{variograms - fitted variograms for all layers}
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
#' covars <- terra::rast(files)
#'
#' # spatial autocorrelation of a binary/continuous response
#' sac1 <- cv_spatial_autocor(x = pa_data,
#'                            column = "occ", # binary or continuous data
#'                            plot = TRUE)
#'
#'
#' # spatial autocorrelation of continuous raster files
#' sac2 <- cv_spatial_autocor(r = covars,
#'                            num_sample = 5000,
#'                            plot = TRUE)
#'
#' # show the result
#' summary(sac2)
#' }
cv_spatial_autocor <- function(r,
                               x,
                               column = NULL,
                               num_sample = 5000L,
                               deg_to_metre = 111325,
                               plot = TRUE,
                               progress = TRUE,
                               ... # extra arguments for cv_plot
){

  # check for availability of required packages
  pkg <- c("ggplot2", "cowplot", "automap", "terra")
  .check_pkgs(pkg)

  # check x is an sf object
  if(!missing(x)){
    x <- .check_x(x)
  }

  # x and column must be provided
  if(!missing(x) && is.null(column)){
    stop("When 'x' is provided, 'column' must also be provided. Otherwise, provide only 'r'.")
  }

  # is column in x?
  if(!missing(x) && !is.null(column)){
    if(!all(column %in% colnames(x))){
      wc <- which(!column %in% colnames(x))
      stop(sprintf("There is no column named '%s' in 'x'.\n", column[wc]))
    }
  }

  # change the r to terra object
  if(!missing(r)){
    r <- .check_r(r)
  }

  if(!missing(r) && missing(x)){
    if(any(terra::is.factor(r))){
      stop("'r' should not include any factor layers!")
    }

    nlayer <- terra::nlyr(r)

    # assume wgs84 if crs is na
    if(is.na(sf::st_crs(r))){
      mapext <- terra::ext(r)[1:4]
      if(all(mapext >= -180) && all(mapext <= 180)){
        terra::crs(r) <- "epsg:4326"
        warning("The input layer CRS is undefined! Based on the extent of the input map it is assumed to have geographic coordinate system (epsg:4326)")
      } else{
        stop("Please provide a raster object with a defined CRS.")
      }
    }

    # reduce the num_sample if the raster does not have enough cells
    if(terra::ncell(r) < 5 * num_sample){
      rp <- length(terra::cells(r))
      if(rp < num_sample){
        num_sample <- rp
        message("The num_sample reduced to ", num_sample, "; the total number of available cells.\n")
      }
    }

    # to be used for filtering blocks
    samp_point <- terra::spatSample(r[[1]],
                                    size = 5e4,
                                    method = "regular",
                                    as.points = TRUE,
                                    na.rm = TRUE)
  }

  if(!missing(x)){
    if(is.na(sf::st_crs(x))){
      stop("The coordinate reference system of 'x' must be defined.")
    }
    # no raster is provided
    nlayer <- length(column)
  }

  if(nlayer < 2) progress <- FALSE
  if(progress) pb <- utils::txtProgressBar(min = 0, max = nlayer, style = 3)

  # fitting wariogram models
  vario_list <- lapply(
    seq_len(nlayer),
    .fit_variogram,
    rr = switch(missing(x), r, NULL),
    xx = switch(!missing(x), x, NULL),
    column = column,
    num_sample = num_sample,
    progress = progress,
    pb = if(progress) pb else NULL
  )

  # make a dataframe from variograms data
  vario_data <- data.frame(layers = seq_len(nlayer), range = 1, sill = 1)
  for(v in seq_along(vario_list)){
    vario_data$layers[v] <- if(missing(x)) names(r)[v] else column[v]
    vario_data$range[v] <- vario_list[[v]]$var_model[2,3]
    vario_data$sill[v] <- vario_list[[v]]$var_model[2,2]
  }

  # order them for plotting
  vario_data <- vario_data[order(vario_data$range), ]

  x_obj <- if(missing(x)) sf::st_as_sf(samp_point) else x

  size <- the_range <- stats::median(vario_data$range)

  if(sf::st_is_longlat(x_obj)){
    vario_data$range <- vario_data$range * 1000
    the_range <- the_range * 1000
    size <- the_range / deg_to_metre
  }

  # make an object for plotting
  vis_block <- sf::st_make_grid(x_obj, cellsize = round(size), what = "polygons")
  vis_block <- sf::st_sf(vis_block[x_obj])
  vis_block$folds <- 1:nrow(vis_block)
  plot_data <- list(blocks = vis_block)
  class(plot_data) <- "cv_spatial"

  # update num_sample
  num_sample <- ifelse(missing(x), num_sample, nrow(x))

  if(nlayer > 1){
    p1 <- .make_bar_plot(vario_data, the_range, num_sample)
  }

  # plot the spatial blocks
  p2 <- cv_plot(
    cv = plot_data,
    r = switch(!missing(r), r, NULL),
    label_size = 0,
    ...
  )
  p2 <- p2 + ggplot2::ggtitle(label = "Spatial blocks",
                              subtitle = paste("Using",
                                               round(the_range),
                                               "metres block size"))


  if(plot){
    if(nlayer > 1) plot(cowplot::plot_grid(p1, p2)) else plot(p2)
  }

  final_list <- list(range = the_range,
                     range_table = vario_data,
                     plots = if(nlayer > 1) list(barchart = p1, map_plot = p2) else p2,
                     num_sample = num_sample,
                     variograms = vario_list)

  # specify the output class
  class(final_list) <- c("cv_spatial_autocor")
  return(final_list)
}


#' @export
#' @method print cv_spatial_autocor
print.cv_spatial_autocor <- function(x, ...){
  print(class(x))
}


#' @export
#' @method plot cv_spatial_autocor
plot.cv_spatial_autocor <- function(x, y, ...){
  if(length(x$plots) == 2){
    plot(cowplot::plot_grid(x$plots$barchart, x$plots$mapplot))
  } else{
    plot(x$plots$mapplot)
  }
}

#' @export
#' @method summary cv_spatial_autocor
summary.cv_spatial_autocor <- function(object, ...){
  cat("Summary statistics of spatial autocorrelation ranges of all input layers:\n")
  print(summary(object$rangeTable$range))
  print(object$rangeTable[,1:2])
}



# make a bar plot for cv_spatial_autocor
.make_bar_plot <- function(vario_data, the_range, ptnum){
  # change the scale to km
  vario_data$range <- vario_data$range / 1000
  the_range <- the_range / 1000

  p <- ggplot2::ggplot(
    data = vario_data,
    ggplot2::aes(y = get("range"),
                 x = stats::reorder(factor(get("layers")),
                                    get("range"),
                                    decreasing = FALSE),
                 color = get("range"))
  ) +
    # ggplot2::geom_bar(
    #   ggplot2::aes(x = stats::reorder(factor(get("layers")), get("range")),
    #                y = get("range"), fill = get("range")),
    #   stat = "identity", data = vario_data,) +
    ggplot2::geom_point(size = 4) +
    ggplot2::geom_segment(
      ggplot2::aes(x = get("layers"), xend = get("layers"), y = 0, yend = get("range")),
      linewidth = 1.5) +
    ggplot2::labs(x = "Variables", y = "Range (km)") +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Autocorrelation range", subtitle = paste("Based on", ptnum, "sample points"))+
    ggplot2::guides(color = "none") +
    ggplot2::geom_hline(yintercept = the_range, color = 'red', linewidth = 0.5, linetype = 2) +
    ggplot2::annotate("text", x = floor(nrow(vario_data) / 3),
                      y =  (the_range + (max(vario_data$range) / 20)),
                      angle = 270,
                      label = "Block size",
                      color = 'red') +
    ggplot2::coord_flip()

  return(p)
}


# auto-fit variogram models
.fit_variogram <- function(i,
                           xx = NULL,
                           rr = NULL,
                           column = NULL,
                           num_sample = 1e4,
                           progress = FALSE,
                           pb = NULL){
  if(is.null(xx)){
    points <- terra::spatSample(
      x = rr[[i]],
      size = num_sample,
      method = "random",
      as.points = TRUE,
      na.rm = TRUE
    )
    points <- sf::as_Spatial(sf::st_as_sf(points))
    names(points) <- "target"
  } else{
    points <- xx[column[i]] # [i] in case there are more column
    points <- sf::as_Spatial(points)
    names(points) <- "target"
  }
  fit_vario <- automap::autofitVariogram(target~1, points)
  if(progress) utils::setTxtProgressBar(pb, i)

  return(fit_vario)
}
