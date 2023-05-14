#' Use spatial blocks to separate train and test folds
#'
#' This function is deprecated and will be removed in future updates! Please use \code{\link{cv_spatial}} instead!
#'
#' @inheritParams buffering
#' @param species Character (optional). Indicating the name of the column in which species data (response variable e.g. 0s and 1s) is stored.
#' This argument is used \emph{to make folds with evenly distributed records}. \strong{This option only works by random fold selection and with binary or
#' multi-class responses} e.g. species presence-absence/background or land cover classes for remote sensing image classification.
#' If \code{speceis = NULL} the response classes will be treated the same and only training and testing records
#' will be counted and balanced.
#' @param theRange Numeric value of the specified range by which blocks are created and training/testing data are separated.
#' This distance should be in \strong{metres}. The range could be explored by \code{spatialAutoRange()} and \code{rangeExplorer()} functions.
#' @param k Integer value. The number of desired folds for cross-validation. The default is \code{k = 5}.
#' @param selection Type of assignment of blocks into folds. Can be \strong{random} (default), \strong{systematic}, \strong{checkerboard}, or \strong{predefined}.
#' The checkerboard does not work with user-defined spatial blocks. If the selection = 'predefined', user-defined blocks and foldsCol must be supplied.
#' @param iteration Integer value. The number of attempts to create folds that fulfil the set requirement for minimum number
#' of points in each training and testing fold (for each response class e.g. \emph{train_0}, \emph{train_1}, \emph{test_0}
#' and \emph{test_1}), as specified by \code{species} and \code{numLimit} arguments.
#' @param blocks A sf or SpatialPolygons object to be used as the blocks (optional). This can be a user defined polygon and it must cover all
#' the species (response) points. If the selection = 'predefined', this argument (and foldsCol) must be supplied.
#' @param foldsCol Character. Indicating the name of the column (in user-defined blocks) in which the associated folds are stored.
#' This argument is necessary if you choose the 'predefined' selection.
#' @param numLimit deprecated option!
#' @param maskBySpecies Since version 1.1, this option is always set to \code{TRUE}.
#' @param degMetre Integer. The conversion rate of metres to degree. See the details section for more information.
#' @param rasterLayer A raster object for visualisation (optional). If provided, this will be used to specify the blocks covering the area.
#' @param border deprecated option!
#' @param showBlocks Logical. If TRUE the final blocks with fold numbers will be created with ggplot and plotted. A raster layer could be specified
#' in \code{rasterlayer} argument to be as background.
#' @param biomod2Format Logical. Creates a matrix of folds that can be directly used in the \pkg{biomod2} package as
#' a \emph{DataSplitTable} for cross-validation.
#' @param progress Logical. If TRUE shows a progress bar when \code{numLimit = NULL} in random fold selection.
#' @param rows Integer value by which the area is divided into latitudinal bins.
#' @param cols Integer value by which the area is divided into longitudinal bins.
#' @param xOffset Numeric value between \strong{0} and \strong{1} for shifting the blocks horizontally.
#' The value is the proportion of block size.
#' @param yOffset Numeric value between \strong{0} and \strong{1} for shifting the blocks vertically. The value is the proportion of block size.
#' @param extend numeric; This parameter specifies the percentage by which the map's extent is
#' expanded to increase the size of the square spatial blocks, ensuring that all points fall
#' within a block. The value should be a numeric between 0 and 5.
#' @param seed Integer. A random seed generator for reproducibility.
#' @param verbose Logical. To print the report of the recods per fold.
#'
#' @seealso \code{\link{cv_spatial}}
#'
#' @export
spatialBlock <- function(speciesData,
                         species = NULL,
                         rasterLayer = NULL,
                         theRange = NULL,
                         rows = NULL,
                         cols = NULL,
                         k = 5L,
                         selection = "random",
                         iteration = 100L,
                         blocks = NULL,
                         foldsCol = NULL,
                         numLimit = 0L,
                         maskBySpecies = TRUE, # remove in the next version
                         degMetre = 111325,
                         border = NULL,
                         showBlocks = TRUE,
                         biomod2Format = TRUE,
                         xOffset = 0,
                         yOffset = 0,
                         extend = 0,
                         seed = 42,
                         progress = TRUE,
                         verbose = TRUE){

  message("This function is deprecated! Please use 'cv_spatial' instead.")


  speciesData <- .check_x(speciesData, name = "speciesData")

  .check_pkgs("sf")

  if(!is.null(species)){
    if(!species %in% colnames(speciesData)){
      warning(sprintf("There is no column named '%s' in 'speciesData'.\n", species))
      species <- NULL
    }
  }

  # checks for pre-defined folds
  if(selection == "predefined"){
    if(is.null(foldsCol) || is.null(blocks)){
      stop("The 'blocks' and 'foldsCol' should be specified for 'predefined' selection")
    }
    if(!foldsCol %in% colnames(blocks)){
      stop(sprintf("There is no column named '%s' in 'blocks'.\n", foldsCol))
    }
    if(!is.numeric(blocks[,foldsCol, drop = TRUE])){
      stop("The fold numbers in 'foldsCol' must be integer numbers.")
    }
  }

  # change the r to terra object
  if(!is.null(rasterLayer)){
    rasterLayer <- .check_r(rasterLayer, name = "rasterLayer")
  }

  out <- cv_spatial(
    x = speciesData,
    column = species,
    r = rasterLayer,
    k = k,
    hexagon = FALSE,
    flat_top = FALSE,
    size = theRange,
    rows_cols = c(rows, cols),
    selection = selection,
    iteration = iteration,
    user_blocks = blocks,
    folds_column = foldsCol,
    deg_to_metre = degMetre,
    biomod2 = biomod2Format,
    offset = c(xOffset, yOffset),
    extend = extend,
    seed = seed,
    progress = progress,
    report = verbose,
    plot = showBlocks
  )

  theList <- list(folds = out$folds_list,
                  foldID = out$folds_ids,
                  biomodTable = out$biomod_table,
                  k = k,
                  # blocks = sf::as_Spatial(out$blocks),
                  blocks = out$blocks,
                  species = out$column,
                  range = out$size,
                  plots = if(showBlocks) cv_plot(out) else NULL,
                  records = out$records)

  class(theList) <- c("SpatialBlock")
  return(theList)
}


#' @export
#' @method print SpatialBlock
print.SpatialBlock <- function(x, ...){
  print(class(x))
}


#' @export
#' @method plot SpatialBlock
plot.SpatialBlock <- function(x, y, ...){
  if(is.null(x$plots)){
    plot(x$blocks)
  } else{
    plot(x$plots)
  }
  message("Please use foldExplorer function to plot each fold interactively.")
}


#' @export
#' @method summary SpatialBlock
summary.SpatialBlock <- function(object, ...){
  cat("Number of recoreds in each training and testing fold:\n")
  print(object$records)
}
