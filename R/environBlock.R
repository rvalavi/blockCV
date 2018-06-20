standardize <- function(x){
  stzRaster <- raster::stack()
  for(i in 1:raster::nlayers(x)){
    stzRaster <- raster::stack(stzRaster, (x[[i]] - raster::minValue(x[[i]])) / (raster::maxValue(x[[i]]) - raster::minValue(x[[i]])))
  }
  return(stzRaster)
}

normalize <- function(x){
  stzRaster <- raster::stack()
  for(i in 1:raster::nlayers(x)){
    meanR <- mean(raster::values(x[[i]]), na.rm=TRUE)
    sdR <- sd(raster::values(x[[i]]), na.rm=TRUE)
    stzRaster <- raster::stack(stzRaster, ((x[[i]] - meanR) / sdR))
  }
  return(stzRaster)
}

#' Use environmental clustering to separate train and test folds
#'
#' Environmental blocking for cross-validation. This function uses clustering methods to specify sets of similar environmental
#' conditions based on the input covariates. Species data corresponding to any of these groups or clusters are assigned to a
#' fold. This function does the clustering in raster space and species data. Clustering is done using \code{\link[stats]{kmeans}}
#' for both approaches (for raster using \pkg{RStoolbox} which use the same function internally).
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
#' # make a SpatialPointsDataFrame object from data.frame
#' pa_data <- sp::SpatialPointsDataFrame(PA[,c("x", "y")], PA, proj4string=raster::crs(awt))
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
envBlock <- function(rasterLayer, speciesData, species=NULL, k=5, standardization="normal", rasterBlock=TRUE, biomod2Format=TRUE, numLimit=0){
  if(methods::is(rasterLayer, 'Raster')){
    if(raster::nlayers(rasterLayer) >= 1){
      foldList <- list()
      if(!is.null(species)){
        presences <- speciesData[speciesData@data[,species]==1,] # creating a layer of presence data
        trainTestTable <- base::data.frame(trainPr=rep(0, k), trainAb=0, testPr=0, testAb=0)
      } else{
        trainTestTable <- base::data.frame(train=rep(0, k), test=0)
      }
      biomodTable <- data.frame(RUN1=rep(TRUE, nrow(speciesData)))
      if(standardization=="standard"){
        rasterLayer <- standardize(rasterLayer)
      } else if(standardization=="normal"){
        rasterLayer <- normalize(rasterLayer)
      } else if(standardization=="none" || is.null(standardization)){
        rasterLayer <- rasterLayer
      }
      if(rasterBlock==TRUE){
        unsClass <- RStoolbox::unsuperClass(rasterLayer, nSamples = 10000, nClasses = k, nIter = 500)
        unsMap <- unsClass$map
        speciesData@data$fold <- raster::extract(unsMap, speciesData)
        if(any(is.na(speciesData@data$fold))){
          stop("The input raster layer does not cover all the species points.")
        }
      } else{
        ext <- raster::extract(rasterLayer, speciesData)
        if(any(is.na(ext))){
          stop("The input raster layer does not cover all the species points.")
        }
        kms <- stats::kmeans(ext, centers=k, iter.max = 500, nstart = 25)
        speciesData$fold <- kms$cluster
      }
      for(i in 1:k){
        testSet <- which(speciesData@data$fold == i)
        trainSet <- which(speciesData@data$fold != i)
        foldList[[i]] <- assign(paste0("fold", i), list(trainSet, testSet))
        if(!is.null(species)){
          lnPrsences <- length(presences)
          lnAbsences <- length(speciesData) - lnPrsences
          trainPoints <- speciesData[trainSet, ]
          trainTestTable$trainPr[i] <- length(trainPoints[trainPoints@data[,species]==1,])
          trainTestTable$trainAb[i] <- length(trainPoints[trainPoints@data[,species]!=1,])
          trainTestTable$testPr[i] <- lnPrsences - trainTestTable$trainPr[i]
          trainTestTable$testAb[i] <- lnAbsences - trainTestTable$trainAb[i]
        } else{
          trainTestTable$train[i] <- length(trainSet)
          trainTestTable$test[i] <- length(testSet)
        }
        if(biomod2Format==TRUE){ # creating a biomod2 DataSplitTable for validation
          colm <- paste0("RUN", i)
          biomodTable[,colm] <- FALSE
          biomodTable[trainSet, colm] <- TRUE
        }
      }
      if(!is.null(species)){
        for(j in 1:ncol(trainTestTable)){
          minNUm <- min(trainTestTable[,j])
          if(minNUm <= numLimit){
            colm <- names(trainTestTable)[j]
            if(colm == "trainPr"){
              colLimit <- "traning presences"
            } else if(colm == "trainAb"){
              colLimit <- "training absences"
            } else if(colm == "testPr"){
              colLimit <- "testing presences"
            } else if(colm == "testAb"){
              colLimit <- "testing absences"
            }
            message(paste("The number of records in", colLimit, "is less than or equal", numLimit))
          }
        }
      } else{
        for(j in 1:2){
          minNUm <- min(trainTestTable[,j])
          if(minNUm <= numLimit){
            colm <- names(trainTestTable)[j]
            if(colm == "train"){
              colLimit <- "traning folds"
            } else if(colm == "test"){
              colLimit <- "testing folds"
            }
            message(paste("The number of records in", colLimit, "is less than or equal", numLimit))
          }
        }
      }
      if(biomod2Format==TRUE){
        biomodTable <- as.matrix(biomodTable)
        theList <- list(folds=foldList, biomodTable=biomodTable, k=k, species=species, records=trainTestTable)
      } else{
        theList <- list(folds=foldList, biomodTable=NULL, k=k, species=species, records=trainTestTable)
      }
    } else stop("'The raster layer is empty!'")
  } else stop('The input file is not a valid R raster file')
  print(trainTestTable)
  # specify the output class
  class(theList) <- c("EnvironmentalBlock")
  return(theList)
}


#' @export
print.EnvironmentalBlock <- function(x, ...){
  print(class(x))
}

#' @export
summary.EnvironmentalBlock <- function(object, ...){
  print("Number of recoreds in each category")
  print(object$records)
}
