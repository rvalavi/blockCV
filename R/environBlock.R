#' Use environmental clustering to separate train and test folds
#'
#' Environmental blocking for cross-validation. This function uses clustering methods to specify sets of similar environmental
#' conditions based on the input covariates. Species data corresponding to any of these groups or clusters are assigned to a
#' fold. This function does the clustering in raster space and species data. Clustering is done using \code{\link[stats]{kmeans}}
#' for both approaches (for raster using \pkg{RStoolbox} which use the same function internally). This function works on single or
#' multiple raster files; multiple rasters need to be in a brick or stack format.
#'
#' As k-means algorithms use Euclidean distance to estimate clusters, the input covariates should be quantitative variables. Since
#' variables with wider ranges of values might dominate the clusters and bias the environmental clustering (Hastie et al., 2009),
#' all the input rasters are first standardized within the function. This is done either by normalizing based on subtracting the
#' mean and dividing by the standard deviation of each raster (the default) or optionally by standardizing using linear scaling
#' to constrain all raster values between 0 and 1.
#'
#' By default, the clustering is done in the raster space. In this approach the clusters will be consistent throughout the region
#' and across species (in the same region). However, this may result in a cluster(s) that covers none of the species records,
#' espcially when species data is not dispersed throughout the region or the number of clusters (k or folds) is high. In this
#' case, the number of folds is less than specified \code{k}. If \code{rasterBlock = FALSE}, the clustering will be done in
#' species points and the number of the folds will be the same as \code{k}.
#'
#' Note that the input raster layer should cover all the species points, otherwise an error will rise. The records with no raster
#' value should be deleted prior to the analysis or another raster layer would be provided.
#'
#' @param rasterLayer RasterLayer, RasterBrick or RasterStack of covariates to identify environmental groups.
#' @inheritParams buffering
#' @inheritParams spatialBlock
#' @param rasterBlock Logical. If TRUE, the clustering is done in the raster layer rather than species data. See details for
#' more information.
#' @param standardization Standardize input raster layers. Three possible inputs are "normal" (the default), "standard" and "none".
#' See details for more information.
#' @param numLimit Integer value. The minimum number of points in each category of data (\emph{training-presence},
#' \emph{training-absence}, \emph{testing-presence} and \emph{testing-absence}). Shows a message if the number of points
#' in any of the folds happens to be less than this number.
#'
#' @seealso \code{\link{spatialBlock}} and \code{\link{buffering}} for alternative blocking strategies; \code{\link{foldExplorer}}
#' for visualisation of the generated folds.
#' @seealso For \emph{DataSplitTable} see \code{\link[biomod2]{BIOMOD_cv}} in \pkg{biomod2} package. \code{\link[RStoolbox]{unsuperClass}}
#' for clustering.
#'
#' @references Hastie, T., Tibshirani, R., & Friedman, J. (2009). The elements of statistical learning: Data Mining, Inference,
#' and Prediction (2nd ed., Vol. 1). Springer series in statistics New York.
#'
#' Roberts et al., 2017. Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. Ecography. 40: 913-929.
#'
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
#'
#' @examples
#' \dontrun{
#'
#' # load package data
#' awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))
#' # import presence-absence species data
#' PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
#' # make a sf object from data.frame
#' pa_data <- sf::st_as_sf(PA coords = c("x", "y"), crs = raster::crs(awt))
#'
#' # environmental clustering
#' eb <- envBlock(rasterLayer = awt,
#'                speciesData = pa_data,
#'                species = "Species", # name of the column with species data
#'                k = 5,
#'                standardization = "standard",
#'                rasterBlock = TRUE,
#'                numLimit = 50)
#' }
#'
envBlock <- function(rasterLayer,
                     speciesData,
                     species = NULL,
                     k = 5,
                     standardization = "normal",
                     rasterBlock = TRUE,
                     biomod2Format = TRUE,
                     numLimit = 0){
  ## change the sp objects to sf
  if(methods::is(speciesData, "SpatialPoints")){
    speciesData <- sf::st_as_sf(speciesData)
  } else if(!methods::is(speciesData, "sf")){
    stop("speciesData should be a sf or SpatialPoints object")
  }
  ## check if species is a col in speciesData
  if(!is.null(species)){
    if(species %in% colnames(speciesData) == FALSE){
      cat("There is no match between the columns name in 'speciesData' and 'species' argument (response variable).\n")
      species <- NULL
    }
  }
  if(methods::is(rasterLayer, 'Raster')){
    if(raster::nlayers(rasterLayer) >= 1){
      foldList <- list()
      foldNum <- rep(NA, nrow(speciesData))
      biomodTable <- data.frame(RUN1=rep(TRUE, nrow(speciesData)))
      if(standardization=="standard"){
        rasterLayer <- standardize(rasterLayer)
      } else if(standardization=="normal"){
        rasterLayer <- normalize(rasterLayer)
      } else if(standardization=="none" || is.null(standardization)){
        rasterLayer <- rasterLayer
      }
      if(rasterBlock==TRUE){
        unsClass <- RStoolbox::unsuperClass(rasterLayer, nSamples = 10000, nClasses = k, nIter = 500) # more samples or a new argument?!
        unsMap <- unsClass$map
        speciesData$fold <- raster::extract(unsMap, speciesData)
        if(anyNA(speciesData$fold)){
          stop("The input raster layer does not cover all the species points.")
        }
      } else{
        ext <- raster::extract(rasterLayer, speciesData)
        if(anyNA(ext)){
          stop("The input raster layer does not cover all the species points.")
        }
        kms <- stats::kmeans(ext, centers = k, iter.max = 500, nstart = 25)
        speciesData$fold <- kms$cluster
      }
      if(is.null(species)){
        trainTestTable <- base::data.frame(train = rep(0, k), test = 0)
      } else{
        cl <- sort(unique(speciesData[, species, drop = TRUE]))
        clen <- length(cl)
        trainTestTable <- as.data.frame(matrix(0, nrow = k, ncol = clen * 2))
        names(trainTestTable) <- c(paste("train", cl, sep = "_"), paste("test", cl, sep = "_"))
      }
      for(i in seq_len(k)){
        testSet <- which(speciesData$fold == i)
        trainSet <- which(speciesData$fold != i)
        foldNum[testSet] <- i
        foldList[[i]] <- assign(paste0("fold", i), list(trainSet, testSet))
        if(is.null(species)){
          trainTestTable$train[i] <- length(trainSet)
          trainTestTable$test[i] <- length(testSet)
        } else{
          countrain <- table(speciesData[trainSet ,species, drop = TRUE])
          countest <- table(speciesData[testSet ,species, drop = TRUE])
          trainTestTable[i, which(cl %in% names(countrain))] <- countrain
          trainTestTable[i, clen + which(cl %in% names(countest))] <- countest
        }
        if(biomod2Format==TRUE){ # creating a biomod2 DataSplitTable for validation
          colm <- paste0("RUN", i)
          biomodTable[,colm] <- FALSE
          biomodTable[trainSet, colm] <- TRUE
        }
      }
      if(any(trainTestTable <= numLimit)){
        zerofolds <- which(apply(trainTestTable, 1, function(x) any(x <= numLimit)))
        if(length(zerofolds) > 1){
          warning("The folds ", paste(zerofolds, collapse = ", "), " have class(es) with ", numLimit, " or less records")
        } else{
          warning("The fold ", zerofolds, " has class(es) with ", numLimit, " or less records")
        }
      }
      if(biomod2Format==TRUE){
        biomodTable <- as.matrix(biomodTable)
        theList <- list(folds=foldList, foldID=foldNum, biomodTable=biomodTable, k=k, species=species, records=trainTestTable)
      } else{
        theList <- list(folds=foldList, foldID=foldNum, biomodTable=NULL, k=k, species=species, records=trainTestTable)
      }
    } else stop("'The raster layer is empty!'")
  } else stop('The input file is not a valid R raster file')
  print(trainTestTable)
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
