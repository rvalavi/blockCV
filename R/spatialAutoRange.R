#' Measure spatial autocorrelation in the predictor raster files
#'
#' This function is deprecated and will be removed in future updates! Please use \code{\link{cv_spatial_autocor}} instead!
#'
#'
#' @param rasterLayer A raster object of covariates to find spatial autocorrelation range.
#' @param sampleNumber Integer. The number of sample points of each raster layer to fit variogram models. It is 5000 by default,
#' however it can be increased by user to represent their region well (relevant to the extent and resolution of rasters).
#' @param border deprecated option!
#' @param speciesData A spatial or sf object (optional). If provided, the \code{sampleNumber} is ignored and
#' variograms are created based on species locations. This option is not recommended if the species data is not
#' evenly distributed across the whole study area and/or the number of records is low.
#' @param showPlots Logical. Show final plot of spatial blocks and autocorrelation ranges.
#' @param maxpixels Number of random pixels to select the blocks over the study area.
#' @param plotVariograms deprecated option!
#' @param doParallel deprecated option!
#' @param nCores deprecated option!
#' @param progress Logical. Shows progress bar. It works only when \code{doParallel = FALSE}.
#' @param degMetre Numeric. The conversion rate of metres to degree. This is for constructing spatial
#' blocks for visualisation. When the input map is in geographic coordinate system (decimal degrees), the block size is
#' calculated based on deviding the calculated \emph{range} by this value to convert to the input map's unit
#' (by default 111325; the standard distance of a degree in metres, on the Equator).
#'
#' @seealso \code{\link{cv_spatial_autocor}}
#'
#' @export
spatialAutoRange <- function(rasterLayer,
                             sampleNumber = 5000L,
                             border = NULL,
                             speciesData = NULL,
                             doParallel = NULL,
                             nCores = NULL,
                             showPlots = TRUE,
                             degMetre = 111325,
                             maxpixels = 1e+05,
                             plotVariograms = FALSE,
                             progress = TRUE){


  message("This function is deprecated! Please use 'cv_spatial_autocor' instead.")

  if(missing(rasterLayer)) stop("'rasterLayer' must be provided!")

  .check_pkgs(c("ggplot2", "terra"))
  # check r
  rasterLayer <- .check_r(rasterLayer, name = "rasterLayer")
  # check r layers
  if(terra::nlyr(rasterLayer) < 1){
    stop("'rasterLayer' is not a valid raster.")
  }

  # check x is an sf object
  if(!is.null(speciesData)){
    speciesData <- .check_x(speciesData, name = "speciesData")

    df <- terra::extract(rasterLayer, speciesData, ID = FALSE)
    speciesData <- cbind(speciesData, df)

    out <- cv_spatial_autocor(x = speciesData,
                              column = names(df),
                              deg_to_metre = degMetre,
                              plot = showPlots,
                              progress = progress)
  } else{

    out <- cv_spatial_autocor(r = rasterLayer,
                              num_sample = sampleNumber,
                              deg_to_metre = degMetre,
                              plot = showPlots,
                              progress = progress)
  }


  finalList <- list(range = out$range,
                    rangeTable = out$range_table,
                    plots = out$plots,
                    sampleNumber = out$num_sample,
                    variograms = out$variograms)

 
  # specify the output class
  class(finalList) <- c("SpatialAutoRange")
  return(finalList)
}


#' @export
#' @method print SpatialAutoRange
print.SpatialAutoRange <- function(x, ...){
  print(class(x))
}


#' @export
#' @method plot SpatialAutoRange
plot.SpatialAutoRange <- function(x, y, ...){
  if(length(x$plots) == 2){
    plot(cowplot::plot_grid(x$plots[[1]], x$plots[[2]]))
  } else{
    plot(x$plots[[1]])
  }
}

#' @export
#' @method summary SpatialAutoRange
summary.SpatialAutoRange <- function(object, ...){
  cat("Summary statistics of spatial autocorrelation ranges of all input layers:\n")
  print(summary(object$rangeTable$range))
  print(object$rangeTable[,1:2])
}
