#' Use spatial blocks to separate train and test folds
#'
#' This function creates spatially separated folds based on a pre-specified distance. It assigns blocks to the training and
#' testing folds  \strong{randomly},  \strong{systematically} or in a  \strong{checkerboard pattern}. The distance (\code{theRange})
#' should be in \strong{metres}, regardless of the unit of the reference system of
#' the input data (for more information see the details section). By default,
#' the function creates blocks according to the extent and shape of the study area, assuming that the user has considered the
#' landscape for the given species and case study. Alternatively, blocks can solely be created based on species spatial data.
#' Blocks can also be offset so the origin is not at the outer
#' corner of the rasters. Instead of providing a distance, the blocks can also be created by specifying a number of rows and/or
#' columns and divide the study area into vertical or horizontal bins, as presented in Wenger & Olden (2012) and Bahn & McGill (2012).
#' Finally, the blocks can be specified by a user-defined spatial polygon layer.
#'
#'
#' To keep the consistency, all the functions use \strong{metres} as their unit. In this function, when the input map
#' has geographic coordinate system (decimal degrees), the block size is calculated based on deviding \code{theRange} by
#' 111325 (the standard distance of a degree in metres, on the Equator) to change the unit to degree. This value is optional
#' and can be changed by user via \code{degMetre} argument.
#'
#' The \code{xOffset} and \code{yOffset} can be used to change the spatial position of the blocks. It can also be used to
#' assess the sensitivity of analysis results to shifting in the blocking arrangements. These options are available when \code{theRange}
#' is defined. By default the region is located in the middle of the blocks and by setting the offsets, the blocks will shift.
#'
#' Roberts et. al. (2017) suggest that blocks should be substantially bigger than the range of spatial
#' autocorrelation (in model residual) to obtain realistic error estimates, while a buffer with the size of
#' the spatial autocorrelation range would result in a good estimation of error. This is because of the so-called
#' edge effect (O'Sullivan & Unwin, 2014), whereby points located on the edges of the blocks of opposite sets are
#' not separated spatially. Blocking with a buffering strategy overcomes this issue (see \code{\link{buffering}}).
#'
#'
#' @inheritParams buffering
#' @param species Character (optional). Indicating the name of the column in which species data (response variable e.g. 0s and 1s) is stored.
#' This argument is used \emph{to make folds with evenly distributed records}. \strong{This option only works by random fold selection and with binary or
#' multi-class responses} e.g. species presence-absence/background or land cover classes for remote sensing image classification.
#' If \code{speceis = NULL} the response classes will be treated the same and only training and testing records
#' will be counted and balanced.
#' @param blocks A sf or SpatialPolygons object to be used as the blocks (optional). This can be a user defined polygon and it must cover all
#' the species points.
#' @param theRange Numeric value of the specified range by which blocks are created and training/testing data are separated.
#' This distance should be in \strong{metres}. The range could be explored by \code{spatialAutoRange()} and \code{rangeExplorer()} functions.
#' @param k Integer value. The number of desired folds for cross-validation. The default is \code{k = 5}.
#' @param selection Type of assignment of blocks into  folds. Can be \strong{random} (default), \strong{systematic} or \strong{checkerboard}.
#' The checkerboard does not work with user-defined spatial blocks.
#' @param iteration Integer value. The number of attempts to create folds that fulfil the set requirement for minimum number
#' of points in each trainig and testing fold (for each response class e.g. \emph{train_0}, \emph{train_1}, \emph{test_0}
#' and \emph{test_1}), as specified by \code{species} and \code{numLimit} arguments.
#' @param numLimit Integer value. The minimum number of points in each training and testing folds.
#' If \code{numLimit = 0}, the most evenly dispersed number of records is chosen (given the number of iteration).
#' This option no longer accepts NULL as input. If it is set to NULL, 0 is used instead.
#' @param maskBySpecies Since version 1.1, this option is always set to \code{TRUE}.
#' @param degMetre Integer. The conversion rate of metres to degree. See the details section for more information.
#' @param rasterLayer A raster object for visualisation (optional). If provided, this will be used to specify the blocks covering the area.
#' @param border A sf or SpatialPolygons object to clip the block based on it (optional).
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
#' @param verbose Logical. To print the report of the recods per fold.
#'
#' @seealso \code{\link{spatialAutoRange}} and \code{\link{rangeExplorer}} for selecting block size; \code{\link{buffering}}
#' and \code{\link{envBlock}} for alternative blocking strategies; \code{\link{foldExplorer}} for visualisation of the generated folds.
#' @seealso For \emph{DataSplitTable} see \code{\link[biomod2]{BIOMOD_cv}} in \pkg{biomod2} package
#'
#' @references Bahn, V., & McGill, B. J. (2012). Testing the predictive performance of distribution models. Oikos, 122(3), 321-331.
#'
#' O'Sullivan, D., Unwin, D.J., (2010). Geographic Information Analysis, 2nd ed. John Wiley & Sons.
#'
#' Roberts et al., (2017). Cross-validation strategies for data with temporal, spatial, hierarchical,
#' or phylogenetic structure. Ecography. 40: 913-929.
#'
#' Wenger, S.J., Olden, J.D., (2012). Assessing transferability of ecological models: an underappreciated aspect of statistical
#' validation. Methods Ecol. Evol. 3, 260-267.
#'
#' @return An object of class S3. A list of objects including:
#'    \itemize{
#'     \item{folds - a list containing the folds. Each fold has two vectors with the training (first) and testing (second) indices}
#'     \item{foldID - a vector of values indicating the number of the fold for each observation (each number corresponds to the same point in species data)}
#'     \item{biomodTable - a matrix with the folds to be used in \pkg{biomod2} package}
#'     \item{k - number of the folds}
#'     \item{blocks - SpatialPolygon of the blocks}
#'     \item{range - the distance band of separating trainig and testing folds, if provided}
#'     \item{species - the name of the species (column), if provided}
#'     \item{plots - ggplot object}
#'     \item{records - a table with the number of points in each category of training and testing}
#'     }
#' @export
#'
#' @examples
#' \donttest{
#'
#' # load package data
#' library(sf)
#'
#' awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))
#' # import presence-absence species data
#' PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
#' # make a sf object from data.frame
#' pa_data <- sf::st_as_sf(PA, coords = c("x", "y"), crs = raster::crs(awt))
#'
#' # spatial blocking by specified range and random assignment
#' sb1 <- spatialBlock(speciesData = pa_data,
#'                     species = "Species",
#'                     theRange = 70000,
#'                     k = 5,
#'                     selection = "random",
#'                     iteration = 100,
#'                     numLimit = NULL,
#'                     biomod2Format = TRUE,
#'                     xOffset = 0.3, # shift the blocks horizontally
#'                     yOffset = 0)
#'
#' # spatial blocking by row/column and systematic fold assignment
#' sb2 <- spatialBlock(speciesData = pa_data,
#'                     species = "Species",
#'                     rasterLayer = awt,
#'                     rows = 5,
#'                     cols = 8,
#'                     k = 5,
#'                     selection = "systematic",
#'                     biomod2Format = TRUE)
#'
#' }
spatialBlock <- function(speciesData,
                         species = NULL,
                         blocks = NULL,
                         rasterLayer = NULL,
                         theRange = NULL,
                         rows = NULL,
                         cols = NULL,
                         k = 5L,
                         selection = "random",
                         iteration = 100L,
                         numLimit = 0L,
                         maskBySpecies = TRUE,
                         degMetre = 111325,
                         border = NULL,
                         showBlocks = TRUE,
                         biomod2Format = TRUE,
                         xOffset = 0,
                         yOffset = 0,
                         progress = TRUE,
                         verbose = TRUE){
  if(showBlocks){
    # check for availability of ggplot2
    pkg <- c("ggplot2")
    pkgna <- names(which(sapply(sapply(pkg, find.package, quiet = TRUE), length) == 0))
    if(length(pkgna) > 0){
      message("This function requires ", pkg, " package for plotting.", "\nWould you like to install it now?\n1: yes\n2: no")
      user <- readline(prompt = paste0("Selection: "))
      if(tolower(user) %in% c("1", "yes", "y")){
        utils::install.packages(pkgna)
      } else{
        stop("Please install ggplot2 package or set showBlocks = FALSE.")
      }
    }
  }
  if(!is.element(selection, c("systematic", "random", "checkerboard"))){
    stop("The selection argument must be 'systematic', 'random' or 'checkerboard'")
  }
  chpattern <- FALSE
  if(selection == "checkerboard"){
    chpattern <- TRUE
    k <- 2
  }
  if(is.null(numLimit)){
    numLimit <- 0
  }
  ## change the sp objects to sf
  if(methods::is(speciesData, "SpatialPoints")){
    speciesData <- sf::st_as_sf(speciesData)
  } else if(!methods::is(speciesData, "sf")){
    stop("speciesData should be a sf or SpatialPoints object")
  }
  if(!is.null(blocks)){
    if(methods::is(blocks, "SpatialPolygons")){
      blocks <- sf::st_as_sf(blocks)
    } else if(!methods::is(blocks, "sf")){
      stop("blocks, should be a spatial or sf object")
    }
  }
  if(!is.null(border)){
    if(methods::is(border, "SpatialPolygons")){
      border <- sf::st_as_sf(border)
    } else if(!methods::is(border, "sf")){
      stop("blocks, should be a spatial or sf object")
    }
  }
  if(maskBySpecies==FALSE){
    message("Since version 1.1, this option is always set to TRUE.\n")
  }
  ## check if species is a col in speciesData
  if(!is.null(species)){
    if(species %in% colnames(speciesData) == FALSE){
      warning("There is no match between the columns name in 'speciesData' and 'species' argument (response variable).\n")
      species <- NULL
    }
  }
  if(is.null(blocks)){
    if(is.null(rasterLayer)){
      net <- rasterNet(speciesData,
                       resolution = theRange,
                       xbin = cols,
                       ybin = rows,
                       degree = degMetre,
                       xOffset = xOffset,
                       yOffset = yOffset,
                       checkerboard = chpattern)
      subBlocks <- net[speciesData,] # subset the blocks
      if(!is.null(border)){
        subBlocks <- sf::st_crop(subBlocks, border)
      }
    } else{
      net <- rasterNet(rasterLayer[[1]],
                       resolution = theRange,
                       xbin = cols,
                       ybin = rows,
                       degree = degMetre,
                       xOffset = xOffset,
                       yOffset = yOffset,
                       checkerboard = chpattern)
      subBlocks <- net[speciesData,]
      if(!is.null(border)){
        subBlocks <- sf::st_crop(subBlocks, border)
      }
    }
  } else{
    subBlocks <- blocks[speciesData,]
    if(selection == "checkerboard"){
      selection <- "systematic"
      message("'checkerboard' fold selection does not work with user defined blocks. 'systematic' will be used instead.\n")
    }
  }
  iteration <- as.integer(iteration)
  if(iteration < 1 || !is.numeric(iteration)){
    iteration <- 1L
    message("The interation has been set to 1! \n", "Iteration must be a positive numeric value.\n")
  } else if(is.numeric(iteration) && iteration >= 10000){
    message("The process might take a while, due to the large number of iterations.\n")
  }
  if(progress==TRUE && numLimit == 0){
    pb <- progress::progress_bar$new(format = " Progress [:bar] :percent in :elapsed",
                                     total=iteration, clear=FALSE, width=75) # add progress bar
  }
  ## do the intersection once and outside of the loop
  # subBlocksDF <- sf:::as.data.frame.sgbp(sf::st_intersects(speciesData, subBlocks))
  subBlocksDF <- as.data.frame(sf::st_intersects(speciesData, subBlocks))
  names(subBlocksDF) <- c("records", "blocks")
  # randomly remove the repeated records occurred on the edges of blocks
  if(nrow(subBlocksDF) > nrow(speciesData)){
    subBlocksDF <- subBlocksDF[sample(nrow(subBlocksDF)), ]
    subBlocksDF <- subBlocksDF[!duplicated(subBlocksDF$records), ]
  } else if(nrow(subBlocksDF) < nrow(speciesData) || anyNA(subBlocksDF)){
    nonoverlap <- nrow(speciesData) - nrow(subBlocksDF)
    warning("At least ", nonoverlap, " of the points are not on the defined spatial blocks")
  }
  nrowBlocks <- length(unique(subBlocksDF$blocks))
  maxNumRecord <- 0
  maxSD <- Inf
  for(i in seq_len(iteration)){
    if(k > nrowBlocks){
      stop("'k' is bigger than the number of spatial blocks\n",
           "The number of spatial blocks is: ", nrowBlocks)
    } else if(k < 2){
      stop("'k' must be 2 or higher")
    }
    if(selection=='systematic'){
      foldDF <- data.frame(blocks = seq_len(nrowBlocks), folds = systematicNum(subBlocks, k))
      subBlocksDF <- merge(x = subBlocksDF, y = foldDF, by = "blocks", all.x = TRUE)
    } else if(selection=='checkerboard'){
      foldDF <- data.frame(blocks = seq_len(nrowBlocks), folds = subBlocks$layer)
      subBlocksDF <- merge(x = subBlocksDF, y = foldDF, by = "blocks", all.x = TRUE)
    } else if(selection=='random'){
      subBlocksDF <- subBlocksDF[, c("records", "blocks")] # to avoid repetition in iterations
      foldDF <- data.frame(blocks = seq_len(nrowBlocks), folds = 0)
      # create random folds
      num <- floor(nrowBlocks / k)
      foldDF$folds[seq_len(num * k)] <- sample(rep(seq_len(k), num), num * k)
      if(nrowBlocks %% k != 0){
        rest <- nrowBlocks %% k
        unfold <- which(foldDF$folds==0)
        foldDF$folds[unfold] <- sample(seq_len(k), rest, replace = FALSE)
      }
      subBlocksDF <- merge(x = subBlocksDF, y = foldDF, by = "blocks", all.x = TRUE)
    }
    if(is.null(species)){
      trainTestTable <- data.frame(train=rep(0, k), test=0)
    } else{
      cl <- sort(unique(speciesData[, species, drop = TRUE]))
      clen <- length(cl)
      trainTestTable <- as.data.frame(matrix(0, nrow = k, ncol = clen * 2))
      names(trainTestTable) <- c(paste("train", cl, sep = "_"), paste("test", cl, sep = "_"))
    }
    foldList <- list()
    foldNum <- rep(NA, nrow(subBlocksDF))
    biomodTable <- data.frame(RUN1=rep(TRUE, nrow(subBlocksDF)))
    for(p in seq_len(k)){
      trainSet <- subBlocksDF$records[which(subBlocksDF$folds != p)]
      testSet <- subBlocksDF$records[which(subBlocksDF$folds == p)]
      foldNum[testSet] <- p
      foldList[[p]] <- assign(paste0("fold", p), list(trainSet, testSet))
      if(is.null(species)){
        trainTestTable$train[p] <- length(trainSet)
        trainTestTable$test[p] <- length(testSet)
      } else{
        countrain <- table(speciesData[trainSet ,species, drop = TRUE])
        countest <- table(speciesData[testSet ,species, drop = TRUE])
        trainTestTable[p, which(cl %in% names(countrain))] <- countrain
        trainTestTable[p, clen + which(cl %in% names(countest))] <- countest
      }
      if(biomod2Format==TRUE){ # creating a biomod2 DataSplitTable for validation
        colm <- paste0("RUN", p)
        biomodTable[,colm] <- FALSE
        biomodTable[trainSet, colm] <- TRUE
      }
    }
    if(selection == "random"){
      if(is.numeric(numLimit) && numLimit > 0){
        if(any(trainTestTable < numLimit)==FALSE){ # exit the loop if meet the limit number
          break
        }
      } else if(numLimit == 0){ # find the highest minimum number in the table and store relevant objects
        if(min(trainTestTable) >= maxNumRecord && stats::sd(unlist(trainTestTable)) < maxSD){
          trainTestTable2 <- trainTestTable
          maxNumRecord <- min(trainTestTable2)
          maxSD <- stats::sd(unlist(trainTestTable))
          subBlocksDF2 <- subBlocksDF
          foldList2 <- foldList
          foldNum2 <- foldNum
          biomodTable2 <- biomodTable
          iter <- i
        }
      } else stop("numLimit argument should be a numeric value equal or hagher than 0 or be NULL")
      if(progress == TRUE && numLimit == 0){
        pb$tick() # update progress bar
      }
    } else{
      break
    }
  }
  if(numLimit == 0 && selection == "random"){ # return the best bloks, table etc.
    subBlocksDF <- subBlocksDF2
    trainTestTable <- trainTestTable2
    foldList <- foldList2
    foldNum <- foldNum2
    biomodTable <- biomodTable2
    if(verbose) cat(paste0("The best folds was in iteration ", iter, ":\n"))
  }
  if(verbose) print(trainTestTable)
  if(any(trainTestTable <= numLimit)){
    zerofolds <- which(apply(trainTestTable, 1, function(x) any(x == numLimit)))
    if(length(zerofolds) > 1){
      warning("The folds ", paste(zerofolds, collapse = ", "), " have class(es) with ", numLimit, " (or less) records")
    } else{
      warning("The fold ", zerofolds, " has class(es) with ", numLimit, " (or less) records")
    }
  }
  ## add the folds number to the blocks
  fold_of_block <- subBlocksDF[, c("blocks", "folds")]
  fold_of_block <- fold_of_block[!duplicated(fold_of_block), ]
  fold_of_block <- fold_of_block[order(fold_of_block$blocks), ]
  subBlocks$folds <- fold_of_block$folds
  if(showBlocks){
    if(is.null(rasterLayer)){
      p2 <- ggplot2::ggplot() +
        ggplot2::geom_sf(data = subBlocks,
                         color ="red",
                         fill ="orangered4",
                         alpha = 0.04,
                         size = 0.2) +
        ggplot2::geom_sf_text(ggplot2::aes_string(label = "folds"),
                              data = subBlocks) +
        ggplot2::labs(x = "", y = "") + # or set the axes labes to NULL
        ggplot2::ggtitle("Spatial blocks",
                         subtitle=paste("The", selection, "fold assignment")) +
        ggplot2::theme_bw()
    } else{
      if(methods::is(rasterLayer, "Raster")){
        samp <- raster::sampleRegular(rasterLayer[[1]], 5e+05, asRaster=TRUE)
        map_df <- raster::as.data.frame(samp, xy=TRUE, centroids=TRUE, na.rm=TRUE)
        colnames(map_df) <- c("Easting", "Northing", "MAP")
        mid <- stats::median(map_df$MAP)
        p2 <- ggplot2::ggplot() +
          ggplot2::geom_raster(data = map_df, ggplot2::aes_string(y="Northing", x="Easting", fill="MAP")) +
          ggplot2::scale_fill_gradient2(low="darkred", mid="yellow", high="darkgreen", midpoint=mid) +
          ggplot2::guides(fill = FALSE) +
          ggplot2::geom_sf(data = subBlocks,
                           color ="red",
                           fill ="orangered4",
                           alpha = 0.04,
                           size = 0.2) +
          ggplot2::geom_sf_text(ggplot2::aes_string(label = "folds"),
                                data = subBlocks) +
          ggplot2::labs(x = "", y = "") + # set the axes labes to NULL
          ggplot2::ggtitle("Spatial blocks",
                           subtitle=paste("The", selection, "fold assignment")) +
          ggplot2::theme_bw()
      }
    }
    plot(p2)
  } else{
    p2 <- NULL
  }
  # save the objects
  if(biomod2Format){
    biomodTable <- as.matrix(biomodTable)
    biomodTable2 <- biomodTable
  } else{
    biomodTable2 <- NULL
  }
  theList <- list(folds = foldList,
                  foldID = foldNum,
                  biomodTable = biomodTable2,
                  k = k,
                  blocks = sf::as_Spatial(subBlocks),
                  species = species,
                  range = theRange,
                  plots = p2,
                  records = trainTestTable)
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
