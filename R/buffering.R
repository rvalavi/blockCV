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
#' @return An object of class S3. A list of objects including:
#'     \itemize{
#'     \item{folds - a list containing the folds. Each fold has two vectors with the training (first) and testing (second) indices}
#'     \item{k - number of the folds}
#'     \item{range - the distance band to separated trainig and testing folds)}
#'     \item{species - the name of the species (column), if provided}
#'     \item{dataType - species data type}
#'     \item{records - a table with the number of points in each category of training and testing}
#'     }
#'
#' @export
buffering <- function(speciesData,
                      species = NULL,
                      theRange,
                      spDataType = "PA",
                      addBG = TRUE,
                      progress = TRUE){


  message("This function is deprecated! Please use 'cv_buffer' instead.")

  out <- cv_buffer(x = speciesData,
                   column = species,
                   size = theRange,
                   presence_background = ifelse(spDataType == "PA", FALSE, TRUE),
                   add_background = addBG,
                   progress = progress)

  theList <- list(
    folds = out$folds_list,
    k = out$k,
    species = species,
    range = theRange,
    dataType = spDataType,
    records = out$records
  )

  # ## change the sp objects to sf
  # if(methods::is(speciesData, "SpatialPoints")){
  #   speciesData <- sf::st_as_sf(speciesData)
  # } else if(!methods::is(speciesData, "sf")){
  #   stop("speciesData should be a sf or SpatialPoints object")
  # }
  # ## the crs must be defined
  # if(is.na(sf::st_crs(speciesData))){
  #   stop("The coordinate reference system of speciesData must be defined.")
  # }
  # ## check if species is a col in speciesData
  # if(!is.null(species)){
  #   if(species %in% colnames(speciesData) == FALSE){
  #     warning("There is no match between the columns name in 'speciesData' and 'species' argument (response variable).\n")
  #     species <- NULL
  #   }
  # }
  # dmatrix <- sf::st_distance(speciesData)
  # distuni <- dmatrix[1,1] # take the unit to avoid using units package
  # distuni[1] <- theRange
  # foldList <- list()
  # if(!is.null(species)){
  #   if(spDataType=="PB"){
  #     unqsp <- unique(speciesData[, species, drop = TRUE])
  #     if(!is.numeric(unqsp) || any(unqsp < 0) || any(unqsp > 1)){
  #       stop("Presence-background (PB) type is only for species data with 0s (backgrounds) and 1s (presences).\n", "The data should be numeric.\n")
  #     }
  #     prI <- which(speciesData[, species, drop = TRUE] == 1) # presence indices to loop through
  #     n <- length(prI)
  #     cl <- sort(unique(speciesData[, species, drop = TRUE]))
  #     clen <- length(cl)
  #     trainTestTable <- as.data.frame(matrix(0, nrow = n, ncol = clen * 2))
  #     names(trainTestTable) <- c(paste("train", cl, sep = "_"), paste("test", cl, sep = "_"))
  #     if(progress==TRUE){
  #       pb <- progress::progress_bar$new(format = " Progress [:bar] :percent in :elapsed",
  #                                        total=n, clear=FALSE, width=75) # add progress bar
  #     }
  #     j <- 0
  #     for(i in prI){ # loop through presences
  #       j <- j + 1
  #       trainSet <- which(dmatrix[i, ] > distuni)
  #       if(addBG==TRUE){
  #         test <- which(dmatrix[i, ] <= distuni)
  #         inside <- speciesData[test, ]
  #         testSet <- test[which(inside[, species, drop = TRUE] == 0)]
  #         testSet[length(testSet) + 1] <- i
  #       } else{
  #         testSet <- i
  #       }
  #       foldList[[j]] <- assign(paste0("fold", j), list(trainSet, testSet))
  #       countrain <- table(speciesData[trainSet ,species, drop = TRUE])
  #       countest <- table(speciesData[testSet ,species, drop = TRUE])
  #       trainTestTable[j, which(cl %in% names(countrain))] <- countrain
  #       trainTestTable[j, clen + which(cl %in% names(countest))] <- countest
  #       if(progress==TRUE){
  #         pb$tick() # update progress bar
  #       }
  #     }
  #   } else if(spDataType=="PA"){
  #     n <- nrow(speciesData)
  #     cl <- sort(unique(speciesData[, species, drop = TRUE]))
  #     clen <- length(cl)
  #     trainTestTable <- as.data.frame(matrix(0, nrow = n, ncol = clen * 2))
  #     names(trainTestTable) <- c(paste("train", cl, sep = "_"), paste("test", cl, sep = "_"))
  #     if(progress==TRUE){
  #       pb <- progress::progress_bar$new(format = " Progress [:bar] :percent in :elapsed",
  #                                        total=n, clear=FALSE, width=75) # add progress bar
  #     }
  #     for(i in seq_len(n)){
  #       trainSet <- which(dmatrix[i, ] > distuni)
  #       testSet <- i
  #       foldList[[i]] <- assign(paste0("fold", i), list(trainSet, testSet))
  #       countrain <- table(speciesData[trainSet ,species, drop = TRUE])
  #       countest <- table(speciesData[testSet ,species, drop = TRUE])
  #       trainTestTable[i, which(cl %in% names(countrain))] <- countrain
  #       trainTestTable[i, clen + which(cl %in% names(countest))] <- countest
  #       if(progress==TRUE){
  #         pb$tick() # update progress bar
  #       }
  #     }
  #   }
  # } else{ # data with no species column
  #   n <- nrow(speciesData)
  #   trainTestTable <- base::data.frame(train=rep(0, n), test=0)
  #   if(progress==TRUE){
  #     pb <- progress::progress_bar$new(format = " Progress [:bar] :percent in :elapsed",
  #                                      total=n, clear=FALSE, width=75) # add progress bar
  #   }
  #   for(i in seq_len(n)){
  #     trainSet <- which(dmatrix[i, ] > distuni)
  #     testSet <- i
  #     foldList[[i]] <- assign(paste0("fold", i), list(trainSet, testSet))
  #     trainTestTable$train[i] <- length(trainSet)
  #     trainTestTable$test[i] <- length(testSet)
  #     if(progress==TRUE){
  #       pb$tick() # update progress bar
  #     }
  #   }
  # }
  # theList <- list(folds=foldList, k=n, species=species, range=theRange, dataType=spDataType, records=trainTestTable)
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

