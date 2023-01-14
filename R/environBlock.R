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
#' @return An object of class S3. A list of objects including:
#'    \itemize{
#'     \item{folds - a list containing the folds. Each fold has two vectors with the training (first) and testing (second) indices}
#'     \item{foldID - a vector of values indicating the number of the fold for each observation (each number corresponds to the same point in species data)}
#'     \item{biomodTable - a matrix with the folds to be used in \pkg{biomod2} package}
#'     \item{k - number of the folds}
#'     \item{species - the name of the species (column), if provided}
#'     \item{records - a table with the number of points in each category of training and testing}
#'     }
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
  speciesData <- .x_check(speciesData, name = "speciesData")

  # is column in x?
  if(!is.null(species)){
    if(!species %in% colnames(speciesData)){
      warning(sprintf("There is no column named '%s' in 'speciesData'.\n", species))
      species <- NULL
    }
  }

  # check r
  rasterLayer <- .r_check(rasterLayer, name = "rasterLayer")
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
        message("Scaling the raster failed!")
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
                    print = verbose)

  theList <- list(
    folds = out$folds_list,
    foldID = out$folds_ids,
    biomodTable = out$biomod_table,
    k = k,
    species = out$column,
    records = out$records
  )


  # ## change the sp objects to sf
  # if(methods::is(speciesData, "SpatialPoints")){
  #   speciesData <- sf::st_as_sf(speciesData)
  # } else if(!methods::is(speciesData, "sf")){
  #   stop("speciesData should be a sf or SpatialPoints object")
  # }
  # ## check if species is a col in speciesData
  # if(!is.null(species)){
  #   if(species %in% colnames(speciesData) == FALSE){
  #     warning("There is no match between the columns name in 'speciesData' and 'species' argument (the response variable).\n")
  #     species <- NULL
  #   }
  # }
  # if(methods::is(rasterLayer, 'Raster')){
  #   if(raster::nlayers(rasterLayer) >= 1){
  #     foldList <- list()
  #     foldNum <- rep(NA, nrow(speciesData))
  #     biomodTable <- data.frame(RUN1=rep(TRUE, nrow(speciesData)))
  #     if(standardization=="standard"){
  #       rasterLayer <- standardize(rasterLayer)
  #     } else if(standardization=="normal"){
  #       rasterLayer <- normalize(rasterLayer)
  #     } else if(standardization=="none" || is.null(standardization)){
  #       rasterLayer <- rasterLayer
  #     }
  #     if(rasterBlock==TRUE){
  #       ext <- raster::extract(rasterLayer, speciesData)
  #       if(anyNA(ext)){
  #         stop("The input raster layer does not cover all the species points.")
  #       }
  #       # check number of raster cells
  #       if(raster::ncell(rasterLayer) < 10 * sampleNumber){
  #         rp <- length(which(!is.na(raster::values(rasterLayer[[1]]))))
  #         if(rp < sampleNumber){
  #           sampleNumber <- rp
  #           message("The sampleNumber reduced to ", sampleNumber, "; the total number of available cells.\n")
  #         }
  #       }
  #       sampr <- raster::sampleRandom(rasterLayer, size = sampleNumber)
  #       sampr <- sampr[stats::complete.cases(sampr), ]
  #       sampr <- rbind(ext, sampr)
  #       kms <- stats::kmeans(sampr, centers = k, iter.max = 500, nstart = 25)
  #       speciesData$fold <- kms$cluster[seq_len(nrow(speciesData))]
  #     } else{
  #       ext <- raster::extract(rasterLayer, speciesData)
  #       if(anyNA(ext)){
  #         stop("The input raster layer does not cover all the species points.")
  #       }
  #       kms <- stats::kmeans(ext, centers = k, iter.max = 500, nstart = 25)
  #       speciesData$fold <- kms$cluster
  #     }
  #     if(is.null(species)){
  #       trainTestTable <- data.frame(train = rep(0, k), test = 0)
  #     } else{
  #       cl <- sort(unique(speciesData[, species, drop = TRUE]))
  #       clen <- length(cl)
  #       trainTestTable <- as.data.frame(matrix(0, nrow = k, ncol = clen * 2))
  #       names(trainTestTable) <- c(paste("train", cl, sep = "_"), paste("test", cl, sep = "_"))
  #     }
  #     for(i in seq_len(k)){
  #       testSet <- which(speciesData$fold == i)
  #       trainSet <- which(speciesData$fold != i)
  #       foldNum[testSet] <- i
  #       foldList[[i]] <- assign(paste0("fold", i), list(trainSet, testSet))
  #       if(is.null(species)){
  #         trainTestTable$train[i] <- length(trainSet)
  #         trainTestTable$test[i] <- length(testSet)
  #       } else{
  #         countrain <- table(speciesData[trainSet ,species, drop = TRUE])
  #         countest <- table(speciesData[testSet ,species, drop = TRUE])
  #         trainTestTable[i, which(cl %in% names(countrain))] <- countrain
  #         trainTestTable[i, clen + which(cl %in% names(countest))] <- countest
  #       }
  #       if(biomod2Format==TRUE){ # creating a biomod2 DataSplitTable for validation
  #         colm <- paste0("RUN", i)
  #         biomodTable[,colm] <- FALSE
  #         biomodTable[trainSet, colm] <- TRUE
  #       }
  #     }
  #     if(any(trainTestTable <= numLimit)){
  #       zerofolds <- which(apply(trainTestTable, 1, function(x) any(x <= numLimit)))
  #       if(length(zerofolds) > 1){
  #         warning("The folds ", paste(zerofolds, collapse = ", "), " have class(es) with ", numLimit, " or less records")
  #       } else{
  #         warning("The fold ", zerofolds, " has class(es) with ", numLimit, " or less records")
  #       }
  #     }
  #     if(biomod2Format==TRUE){
  #       biomodTable <- as.matrix(biomodTable)
  #       theList <- list(folds=foldList, foldID=foldNum, biomodTable=biomodTable, k=k, species=species, records=trainTestTable)
  #     } else{
  #       theList <- list(folds=foldList, foldID=foldNum, biomodTable=NULL, k=k, species=species, records=trainTestTable)
  #     }
  #   } else stop("'The raster layer is empty!'")
  # } else stop('The input file is not a valid R raster file')
  # if(verbose) print(trainTestTable)
  # specify the output class
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
