#' Use environmental clustering to separate train and test folds
#'
#' This function is deprecated and will be removed in future updates! Please use \code{\link{cv_cluster}} instead!
#'
#' @param rasterLayer A raster object of covariates to identify environmental groups.
#' @inheritParams buffering
#' @inheritParams spatialBlock
#' @param rasterBlock Logical. If TRUE, the clustering is done in the raster layer rather than species data. See details for
#' more information.
#' @param sampleNumber Integer. The number of samples from raster layers to build the clusters.
#' @param standardization Standardize input raster layers. Three possible inputs are "normal" (the default), "standard" and "none".
#' See details for more information.
#' @param numLimit Integer value. The minimum number of points in each category of data (\emph{train_0},
#' \emph{train_1}, \emph{test_0} and \emph{test_1}). Shows a message if the number of points
#' in any of the folds happens to be less than this number.
#' @param verbose Logical. To print the report of the recods per fold.
#'
#' @seealso \code{\link{cv_cluster}}
#'
#' @export
envBlock <- function(rasterLayer,
                     speciesData,
                     species = NULL,
                     k = 5,
                     standardization = "normal",
                     rasterBlock = TRUE,
                     sampleNumber = 10000,
                     biomod2Format = TRUE,
                     numLimit = 0,
                     verbose = TRUE){

  message("This function is deprecated! Please use 'cv_cluster' instead.")

  if(missing(rasterLayer)) stop("'rasterLayer' must br provided!")
  if(missing(speciesData)) stop("'speciesData' must br provided!")

  # check x is an sf object
  speciesData <- .check_x(speciesData, name = "speciesData")

  # is column in x?
  if(!is.null(species)){
    if(!species %in% colnames(speciesData)){
      warning(sprintf("There is no column named '%s' in 'speciesData'.\n", species))
      species <- NULL
    }
  }

  # check r
  rasterLayer <- .check_r(rasterLayer, name = "rasterLayer")
  # check r layers
  if(terra::nlyr(rasterLayer) < 1){
    stop("'rasterLayer' is not a valid raster.")
  }

  scale <- ifelse(standardization == "none", FALSE, TRUE)

  # scale?
  if(scale){
    tryCatch(
      {
        rasterLayer <- terra::scale(rasterLayer)
      },
      error = function(cond) {
        message("Normalising the raster failed!")
      }
    )
  }


  out <- cv_cluster(x = speciesData,
                    column = species,
                    r = rasterLayer,
                    k = k,
                    scale = scale,
                    raster_cluster = rasterBlock,
                    num_sample = sampleNumber,
                    biomod2 = biomod2Format,
                    report = verbose)

  theList <- list(
    folds = out$folds_list,
    foldID = out$folds_ids,
    biomodTable = out$biomod_table,
    k = k,
    species = out$column,
    records = out$records
  )


  class(theList) <- c("EnvironmentalBlock")
  return(theList)
}


#' @export
#' @method print EnvironmentalBlock
print.EnvironmentalBlock <- function(x, ...){
  print(class(x))
}

#' @export
#' @method summary EnvironmentalBlock
summary.EnvironmentalBlock <- function(object, ...){
  print("Number of recoreds in each category")
  print(object$records)
}
