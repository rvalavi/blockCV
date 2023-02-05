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
#' of points in each trainig and testing fold (for each response class e.g. \emph{train_0}, \emph{train_1}, \emph{test_0}
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
                         iteration = 50L,
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
                         seed = 42,
                         progress = TRUE,
                         verbose = TRUE){

  message("This function is deprecated! Please use 'cv_spatial' instead.")

  # speciesData <- .check_x(speciesData, name = "speciesData")
  #
  # .check_pkgs("sf")
  #
  # if(!is.null(species)){
  #   if(!species %in% colnames(speciesData)){
  #     warning(sprintf("There is no column named '%s' in 'speciesData'.\n", species))
  #     species <- NULL
  #   }
  # }
  #
  # # checks for pre-defined folds
  # if(selection == "predefined"){
  #   if(is.null(foldsCol) || is.null(blocks)){
  #     stop("The 'blocks' and 'foldsCol' should be specified for 'predefined' selection")
  #   }
  #   if(!foldsCol %in% colnames(blocks)){
  #     stop(sprintf("There is no column named '%s' in 'blocks'.\n", foldsCol))
  #   }
  #   if(!is.numeric(blocks[,foldsCol, drop = TRUE])){
  #     stop("The fold numbers in 'foldsCol' must be integer numbers.")
  #   }
  # }
  #
  # # change the r to terra object
  # if(!is.null(rasterLayer)){
  #   rasterLayer <- .check_r(rasterLayer, name = "rasterLayer")
  # }
  #
  # out <- cv_spatial(
  #   x = speciesData,
  #   column = species,
  #   r = rasterLayer,
  #   k = k,
  #   hexagon = FALSE,
  #   flat_top = FALSE,
  #   size = theRange,
  #   rows_cols = c(rows, cols),
  #   selection = selection,
  #   iteration = iteration,
  #   user_blocks = blocks,
  #   folds_column = foldsCol,
  #   deg_to_metre = degMetre,
  #   biomod2 = biomod2Format,
  #   offset = c(xOffset, yOffset),
  #   seed = seed,
  #   progress = progress,
  #   report = verbose,
  #   plot = showBlocks
  # )
  #
  # theList <- list(folds = out$folds_list,
  #                 foldID = out$folds_ids,
  #                 biomodTable = out$biomod_table,
  #                 k = k,
  #                 blocks = sf::as_Spatial(out$blocks),
  #                 species = out$column,
  #                 range = out$size,
  #                 plots = if(showBlocks) cv_plot(out) else NULL,
  #                 records = out$records)




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
  if(!is.element(selection, c("systematic", "random", "checkerboard", "predefined"))){
    stop("The selection argument must be 'random', 'systematic', 'checkerboard', or 'predefined'.")
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
      stop("border, should be a spatial or sf object")
    }
  }
  # remove this for the next update!
  if(maskBySpecies==FALSE){
    message("Since version 1.1, this option is always set to TRUE.\n")
  }
  ## check if species is a col in speciesData
  if(!is.null(species)){
    if(species %in% colnames(speciesData) == FALSE){
      warning("There is no match between the column names in 'speciesData' and 'species' argument (response variable).\n")
      species <- NULL
    }
  }
  if(selection == "predefined"){
    if(is.null(foldsCol) || is.null(blocks)){
      stop("The blocks and foldsCol should be specified for 'predefined' fold selection")
    }
    if(foldsCol %in% colnames(blocks) == FALSE){
      stop("There is no match between the column names in the 'blocks' object and the specified 'foldsCol' argument.\n")
    }
    if(!is.numeric(blocks[,foldsCol, drop = TRUE])){
      stop("The fold values in the specified foldsCol should be integer numbers.")
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
    message("The interation has been set to 1! \n", "Iteration must be a positive integer value.\n")
  } else if(is.numeric(iteration) && iteration >= 1000){
    message("The process might take a while, due to the large number of iterations.\n")
  }
  # if(progress==TRUE && numLimit == 0){
  #   pb <- progress::progress_bar$new(format = " Progress [:bar] :percent in :elapsed",
  #                                    total=iteration, clear=FALSE, width=75) # add progress bar
  # }
  if(progress && numLimit == 0){
    pb <- utils::txtProgressBar(min = 0, max = iteration, style = 3)
  }
  ## do the intersection once and outside of the loop
  subBlocksDF <- as.data.frame(sf::st_intersects(sf::st_geometry(speciesData), sf::st_geometry(subBlocks)))
  names(subBlocksDF) <- c("records", "blocks")
  # randomly remove the repeated records occurred on the edges of blocks
  if(nrow(subBlocksDF) > nrow(speciesData)){
    if(!is.null(seed)){
      set.seed(seed)
    }
    # set.seed(42)
    subBlocksDF <- subBlocksDF[sample(nrow(subBlocksDF)), ]
    subBlocksDF <- subBlocksDF[!duplicated(subBlocksDF$records), ]
  } else if(nrow(subBlocksDF) < nrow(speciesData) || anyNA(subBlocksDF)){
    nonoverlap <- nrow(speciesData) - nrow(subBlocksDF)
    warning("At least ", nonoverlap, " of the points are not on the defined spatial blocks")
  }
  nrowBlocks <- length(unique(subBlocksDF$blocks))
  maxNumRecord <- 0
  maxSD <- Inf
  if(!is.null(seed)){
    set.seed(seed)
  }
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
    } else if(selection=='predefined'){
      foldDF <- data.frame(blocks = seq_len(nrowBlocks), folds = subBlocks[ , foldsCol, drop = TRUE])
      subBlocksDF <- merge(x = subBlocksDF, y = foldDF, by = "blocks", all.x = TRUE)
    }
    # creat records table
    if(is.null(species)){
      trainTestTable <- data.frame(train=rep(0, k), test=0)
    } else{
      cl <- sort(unique(speciesData[, species, drop = TRUE]))
      clen <- length(cl)
      trainTestTable <- as.data.frame(matrix(0, nrow = k, ncol = clen * 2))
      names(trainTestTable) <- c(paste("train", cl, sep = "_"), paste("test", cl, sep = "_"))
    }
    # count the number of points in each fold
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
      # if(progress == TRUE && numLimit == 0){
      #   pb$tick() # update progress bar
      # }
      if(progress && numLimit == 0){ # if iteration is higher than 5?
        utils::setTxtProgressBar(pb, i)
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
        ggplot2::geom_sf_text(ggplot2::aes(label = get("folds")),
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
          ggplot2::geom_tile(
            data = map_df,
            ggplot2::aes(y=get("Northing"), x=get("Easting"), fill=get("MAP"))) +
          ggplot2::scale_fill_gradient2(low="darkred", mid="yellow", high="darkgreen", midpoint=mid) +
          ggplot2::guides(fill = "none") +
          ggplot2::geom_sf(data = subBlocks,
                           color ="red",
                           fill ="orangered4",
                           alpha = 0.04,
                           size = 0.2) +
          ggplot2::geom_sf_text(ggplot2::aes(label = get("folds")),
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
                  # blocks = sf::as_Spatial(subBlocks),
                  blocks = subBlocks,
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



rasterNet <- function(x,
                      resolution = NULL,
                      xbin = NULL,
                      ybin = NULL,
                      mask = FALSE,
                      degree = 111325,
                      xOffset = NULL,
                      yOffset = NULL,
                      checkerboard = FALSE,
                      maxpixels = 1e5){
  if(methods::is(x, "sf")){
    x <- sf::as_Spatial(x)
  }
  ext <- raster::extent(x)
  extRef <- raster::extent(x)
  if(is.na(raster::projection(x))){
    mapext <- raster::extent(x)[1:4]
    if(all(mapext >= -180) && all(mapext <= 180)){
      resolution <- resolution / degree
      warning("The input layer has no CRS defined. Based on the extent of the input map it is assumed to have an un-projected reference system")
    } else {
      resolution <- resolution
      warning("The input layer has no CRS defined. Based on the extent of the input map it is assumed to have a projected reference system")
    }
  } else{
    if(raster::isLonLat(x)){
      resolution <- resolution / degree
    } else{
      resolution <- resolution
    }
  }
  if(!is.null(xbin) && is.null(ybin)){
    rasterNet <- raster::raster(ext, nrow=1, ncol=xbin, crs=raster::projection(x))
  } else if(is.null(xbin) && !is.null(ybin)){
    rasterNet <- raster::raster(ext, nrow=ybin, ncol=1, crs=raster::projection(x))
  } else if(!is.null(xbin) && !is.null(ybin)){
    rasterNet <- raster::raster(ext, nrow=ybin, ncol=xbin, crs=raster::projection(x))
  } else if(is.null(xbin) && is.null(ybin) && !is.null(resolution)){
    xrange <- raster::xmax(x) - raster::xmin(x) # number of columns
    yrange <- raster::ymax(x) - raster::ymin(x) # number of rows
    xPix <- ceiling(xrange / resolution)
    yPix <- ceiling(yrange / resolution)
    xdif <- ((xPix * resolution) - xrange) / 2 # the difference of extent divided by 2 to split on both sides
    ydif <- ((yPix * resolution) - yrange) / 2
    ext@xmin <- raster::xmin(x) - xdif
    ext@xmax <- raster::xmax(x) + xdif
    ext@ymin <- raster::ymin(x) - ydif
    ext@ymax <- raster::ymax(x) + ydif
    if(!is.null(xOffset)){
      if(xOffset > 1 || xOffset < 0){stop("xOffset should be between 0 and 1")}
      ext@xmin <- ext@xmin + (resolution * xOffset)
      ext@xmax <- ext@xmax + (resolution * xOffset)
    }
    if(!is.null(yOffset)){
      if(yOffset > 1 || yOffset < 0){stop("yOffset should be between 0 and 1")}
      ext@ymin <- ext@ymin + (resolution * yOffset)
      ext@ymax <- ext@ymax + (resolution * yOffset)
    }
    # adding cells if needed
    if(ext@xmin > extRef@xmin){ # add one column by increasing the extent and number of bins
      ext@xmin <- ext@xmin - resolution
      xPix <- xPix + 1
    }
    if(ext@ymin > extRef@ymin){
      ext@ymin <- ext@ymin - resolution
      yPix <- yPix + 1
    }
    rasterNet <- raster::raster(ext, nrow=yPix, ncol=xPix, crs=raster::projection(x))
  } else stop("A value should be specified for the block size")
  if(checkerboard){
    raster::values(rasterNet) <- seq_len(raster::ncell(rasterNet))
    m <- raster::as.matrix(rasterNet)
    for(i in seq_len(raster::ncol(rasterNet))){
      if(i %% 2 == 0){
        m[,i] <- rep(1:2, nrow(m))[seq_len(nrow(m))]
      } else{
        m[,i] <- rep(2:1, nrow(m))[seq_len(nrow(m))]
      }
    }
    raster::values(rasterNet) <- m # rasterNet[] <- m
  } else{
    raster::values(rasterNet) <- seq_len(raster::ncell(rasterNet))
  }
  rasterNet <- raster::rasterToPolygons(rasterNet)
  if(mask){
    if(methods::is(x, "Raster")){
      points <- raster::rasterToPoints(x[[1]], spatial=TRUE)
      if(nrow(points) > 75e4){
        maxpixels <- ifelse(maxpixels > 75e4, 75e4, maxpixels)
        points2 <- points[sample(nrow(points), maxpixels, replace=FALSE), ]
        rasterNet <- raster::intersect(rasterNet, points2)
      } else {
        rasterNet <- raster::intersect(rasterNet, points)
      }
    } else{
      rasterNet <- raster::intersect(rasterNet, x)
    }
  }
  return(sf::st_as_sf(rasterNet))
}

systematicNum <- function(layer, num=5){
  n <- nrow(layer)
  if(n %% num == 0){
    a <- n/num
    c <- rep(1:num, a)
  } else {
    a <- floor(n/num)
    b <- n %% num
    c <- c(rep(1:num, a), 1:b)
  }
  return(c)
}
