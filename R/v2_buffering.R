#' Use distance (buffer) around records to separate train and test folds
#'
#' This function is deprecated and will be removed in future updates! Please use \code{\link{cv_buffer}} instead!
#'
#' @param speciesData A simple features (sf) or SpatialPoints object containing species data (response variable).
#' @param species Character. Indicating the name of the field in which species data (binary response i.e. 0 and 1) is stored. If \code{speceis = NULL}
#' the presence and absence data (response variable) will be treated the same and only training and testing records will be counted. This can be used for multi-class responses
#' such as land cover classes for remote sensing image classification, but it is not necessary. \emph{Do not use this argument when the response variable is
#' continuous or count data}.
#' @param theRange Numeric value of the specified range by which the training and testing datasets are separated.
#' This distance should be in \strong{\emph{metres}} no matter what the coordinate system is. The range can  be explored by \code{\link{spatialAutoRange}}.
#' @param spDataType Character input indicating the type of species data. It can take two values, \strong{PA} for \emph{presence-absence} data and \strong{PB} for
#' \emph{presence-background} data, when \code{species} argument is not \code{NULL}. See the details section for more information on these two approaches.
#' @param addBG Logical. Add background points to the test set when \code{spDataType = "PB"}.
#' @param progress Logical. If TRUE a progress bar will be shown.
#'
#' @seealso \code{\link{cv_buffer}}
#'
#' @export
buffering <- function(speciesData,
                      species = NULL,
                      theRange,
                      spDataType = "PA",
                      addBG = TRUE,
                      progress = TRUE){


  message("This function is deprecated! Please use 'cv_buffer' instead.")

  # check x is an sf object
  speciesData <- .check_x(speciesData, name = "speciesData")

  # x's CRS must be defined
  if(is.na(sf::st_crs(speciesData))){
    stop("The coordinate reference system of 'speciesData' must be defined.")
  }
  # is column in x?
  if(!is.null(species)){
    if(!species %in% colnames(speciesData)){
      warning(sprintf("There is no column named '%s' in 'speciesData'.\n", species))
      species <- NULL
    }
  }

  if(is.null(species) && (spDataType != "PA")) stop("'species' must be provided for presence-background data.")


  out <- cv_buffer(x = speciesData,
                   column = species,
                   size = theRange,
                   presence_bg = ifelse(spDataType == "PA", FALSE, TRUE),
                   add_bg = addBG,
                   progress = progress)


  theList <- list(
    folds = out$folds_list,
    k = out$k,
    species = species,
    range = theRange,
    dataType = spDataType,
    records = out$records
  )


  class(theList) <- c("BufferedBlock")
  return(theList)
}


#' @export
#' @method print BufferedBlock
print.BufferedBlock <- function(x, ...){
  print(class(x))
}

#' @export
#' @method summary BufferedBlock
summary.BufferedBlock <- function(object, ...){
  print("Number of recoreds in each category")
  print(object$records)
}

