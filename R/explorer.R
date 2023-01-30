#' Explore the generated folds
#'
#' This function is deprecated! Please use \code{\link{cv_plot}} function for plotting the folds.
#'
#' @param blocks deprecated!
#' @param rasterLayer deprecated!
#' @param speciesData deprecated!
foldExplorer <- function(blocks, rasterLayer, speciesData){
  stop(
    "This function is deprecated! Please use `cv_plot` function for plotting the folds."
  )
}



#' Explore spatial block size
#'
#' This function is deprecated and will be removed in future updates! Please use \code{\link{cv_block_size}} instead!
#'
#' @param rasterLayer raster layer for make plot
#' @param speciesData a simple features (sf) or SpatialPoints object containing species data (response variable). If provided, the species data will be shown on the map.
#' @param species character value indicating the name of the field in which the species data (response variable e.g. 0s and 1s) are stored.
#' If provided, species presence and absence data will be shown in different colours.
#' @param rangeTable deprecated option!
#' @param minRange a numeric value to set the minimum possible range for creating spatial blocks. It is used to limit the searching domain of
#' spatial block size.
#' @param maxRange a numeric value to set the maximum possible range for creating spatial blocks. It is used to limit the searching
#' domain of spatial block size.
#'
#' @seealso \code{\link{cv_block_size}}
#'
#' @export
rangeExplorer <- function(rasterLayer,
                          speciesData=NULL,
                          species=NULL,
                          rangeTable=NULL,
                          minRange=NULL,
                          maxRange=NULL){

  message("This function is deprecated! Please use 'cv_block_size' instead.")
  # check x is an sf object
  if(!is.null(speciesData)){
    speciesData <- .check_x(speciesData, name = "speciesData")
  }

  # is column in x?
  if(!is.null(species)){
    if(!species %in% colnames(speciesData)){
      warning(sprintf("There is no column named '%s' in 'speciesData'.\n", species))
      species <- NULL
    }
  }

  # check r
  rasterLayer <- .check_r(rasterLayer, name = "rasterLayer")

  cv_block_size(r = rasterLayer, # priority
                x = speciesData,
                column = species,
                min_size = minRange,
                max_size = maxRange)


}
