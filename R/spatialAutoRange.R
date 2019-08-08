rasterNet <- function(x, resolution=NULL, xbin=NULL, ybin=NULL, mask=FALSE, degree=111325, xOffset=NULL, yOffset=NULL,
                      checkerboard=FALSE, maxpixels=250000){
  ext <- raster::extent(x)
  extRef <- raster::extent(x)
  if(is.na(sp::proj4string(x))){
    mapext <- raster::extent(x)[1:4]
    if(mapext >= -180 && mapext <= 180){
      resolution <- resolution / degree
      warning("The input layer has no CRS defined. Based on the extent of the input map it is assumed to have an un-projected reference system")
    } else {
      resolution <- resolution
      warning("The input layer has no CRS defined. Based on the extent of the input map it is assumed to have a projected reference system")
    }
  } else{
    if(sp::is.projected(sp::SpatialPoints((matrix(1:10, 5, byrow=FALSE)), proj4string=crs(x)))){
      resolution <- resolution
    } else{
      resolution <- resolution / degree
    }
  }
  if(!is.null(xbin) && is.null(ybin)){
    rasterNet <- raster::raster(ext, nrow=1, ncol=xbin, crs=crs(x))
  } else if(is.null(xbin) && !is.null(ybin)){
    rasterNet <- raster::raster(ext, nrow=ybin, ncol=1, crs=crs(x))
  } else if(!is.null(xbin) && !is.null(ybin)){
    rasterNet <- raster::raster(ext, nrow=ybin, ncol=xbin, crs=crs(x))
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
    rasterNet <- raster::raster(ext, nrow=yPix, ncol=xPix, crs=crs(x))
  } else stop("A value should be specified for the block size")
  if(checkerboard == TRUE){
    values(rasterNet) <- 1:ncell(rasterNet)
    m <- as.matrix(rasterNet)
    for(i in 1:ncol(rasterNet)){
      if(i %% 2 == 0){
        m[,i] <- rep(1:2, nrow(m))[1:nrow(m)]
      } else{
        m[,i] <- rep(2:1, nrow(m))[1:nrow(m)]
      }
    }
    rasterNet[] <- m
  } else{
    values(rasterNet) <- 1:ncell(rasterNet)
  }
  rasterNet <- raster::rasterToPolygons(rasterNet)
  if(mask==TRUE){
    if(methods::is(x, 'Raster')){
      points <- raster::rasterToPoints(x[[1]], spatial=TRUE)
      if(nrow(points) > 750000){
        points2 <- points[sample(1:nrow(points), maxpixels, replace=FALSE), ]
        rasterNet <- raster::intersect(rasterNet, points2)
      } else  rasterNet <- raster::intersect(rasterNet, points)
    } else{
      rasterNet <- raster::intersect(rasterNet, x)
    }
  }
  return(rasterNet)
}



multiplot <- function(..., plotlist=NULL, file, cols=2, layout=NULL) {
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' Measure spatial autocorrelation in the predictor raster files
#'
#' This function provides a quantitative basis for choosing block size. The spatial autocorrelation in all continuous
#' predictor variables available as raster layers is assessed and reported. The function estimates spatial autocorrelation
#' ranges of all input raster layers. This is the range over which observations are independent and is determined by
#' constructing the empirical variogram, a fundamental geostatistical tool for measuring spatial autocorrelation.
#' The empirical variogram models the structure of spatial autocorrelation by measuring variability between all possible
#' pairs of points (O’Sullivan and Unwin, 2010). Results are plotted. See the details section for further information.
#'
#'
#' The input raster layers should be continuous for computing the variograms and estimating the range of spatial
#' autocorrelation. The input rasters should also have a specified coordinate reference system. However, if the reference
#' system is not specified, the function attempts to guess it based on the extent of the map. It assumes an unprojected
#' reference system for layers with extent lying between -180 and 180, and a projected reference system otherwise.
#'
#' Variograms are calculated based on the distances between pairs of points, so unprojected rasters (in degrees) will
#' not give an accurate result (especially over large latitudinal extents). For unprojected rasters, the great circle
#' distance (rather than Euclidian distance) is used to calculate the spatial distances between pairs of points. To
#' enable more accurate estimate, it is recommended to transform unprojected maps (geographic coordinate
#' system / latitude-longitude) to a projected metric reference system (e.g. UTM, Lambert) where it is possible.
#' See \code{\link[automap]{autofitVariogram}} from \pkg{automap} and \code{\link[gstat]{variogram}} from \pkg{gstat} packages
#' for further information.
#'
#'
#' @param rasterLayer RasterLayer, RasterBrick or RasterStack of covariates to find spatial autocorrelation range.
#' @param sampleNumber Integer. The number of sample points of each raster layer to fit variogram models. It is 5000 by default,
#' however it can be increased by user to represent their region well (relevant to the extent and resolution of rasters).
#' @param border A SpatialPolygons* or sf object for clipping output blocks. This increases the computation time slightly.
#' @param speciesData A spatial or sf object. If provided, the \code{sampleNumber} is ignored and
#' variograms are created based on species locations. This option is not recommended if the species data is not
#' evenly distributed across the whole study area and/or the number of records is low.
#' @param showPlots Logical. Show final plot of spatial blocks and autocorrelation ranges.
#' @param maxpixels Number of random pixels to select the blocks over the study area.
#' @param plotVariograms Logical. Plot fitted variograms. This can also be done after the analysis. Set to \code{FALSE} by default.
#' @param doParallel Logical. Run in parallel when more than one raster layer is available. Given multiple CPU cores, it is
#' recommended to set it to \code{TRUE} when there is a large number of rasters to process.
#' @param nCores Integer. Number of CPU cores to run in parallel. If \code{nCores = NULL} half of available cores in your
#' machine will be used.
#' @param progress Logical. Shows progress bar. It works only when \code{doParallel = FALSE}.
#' @param degMetre Integer. The conversion rate of metres to degree. This is for constructing spatial
#' blocks for visualisation. When the input map is in geographic coordinate system (decimal degrees), the block size is
#' calculated based on deviding the calculated \emph{range} by this value to convert to the input map's unit
#' (by default 111325; the standard distance of a degree in metres, on the Equator).
#'
#' @import automap
#' @import foreach
#' @import doParallel
#'
#' @references O’Sullivan, D., Unwin, D.J., 2010. Geographic Information Analysis, 2nd ed. John Wiley & Sons.
#'
#' Roberts et al., 2017. Cross-validation strategies for data with temporal, spatial, hierarchical,
#' or phylogenetic structure. Ecography. 40: 913-929.
#'
#' @return An object of class S3. A list object including:
#'     \itemize{
#'     \item{range - the suggested range, which is the median of all calculated ranges}
#'     \item{rangeTable - a table of input covariates names and their autocorrelation range}
#'     \item{plots - the output plot (the plot is shown by default)}
#'     \item{sampleNumber}
#'     \item{variograms - fitted variograms for all layers}
#'     }
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # load the example raster data
#' awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))
#' # import presence-absence species data
#' PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
#' # coordinate reference system
#' Zone55s <- "+proj=utm +zone=55 +south +ellps=GRS80 +units=m +no_defs"
#' # make a SpatialPointsDataFrame object from data.frame
#' pa_data <- sp::SpatialPointsDataFrame(PA[,c("x", "y")], PA, proj4string=CRS(Zone55s))
#'
#' # run the model in parallel
#' range1 <- spatialAutoRange(rasterLayer = awt,
#'                            sampleNumber = 5000, # number of cells to be used
#'                            doParallel = TRUE,
#'                            nCores = NULL, # use half of the CPU cores
#'                            plotVariograms = FALSE,
#'                            showPlots = TRUE)
#'
#' range2 <- spatialAutoRange(rasterLayer = awt,
#'                            speciesData = pa_data, # use species locations to create variogram(s)
#'                            doParallel = TRUE)
#'
#' # run the model with no parallel
#' range3 <- spatialAutoRange(rasterLayer = awt,
#'                            sampleNumber = 5000,
#'                            doParallel = FALSE,
#'                            showPlots = TRUE,
#'                            progress = TRUE)
#'
#' # show the result
#' summary(range1)
#' }
spatialAutoRange <- function(rasterLayer, sampleNumber=5000, border=NULL, speciesData=NULL,
                             doParallel=TRUE, nCores=NULL, showPlots=TRUE, degMetre=111325,
                             maxpixels=1e+05, plotVariograms=FALSE, progress=TRUE){
  if(!is.null(speciesData)){
    if((methods::is(speciesData, "SpatialPoints") || methods::is(speciesData, "sf"))==FALSE){
      stop("speciesData should be SpatialPoints* or sf object")
    }
  }
  if(methods::is(rasterLayer, 'Raster')){
    if(any(is.factor(rasterLayer))){warning("rasterLayer should not include any factor layer")}
    numLayer <- raster::nlayers(rasterLayer)
    if(is.na(sp::proj4string(rasterLayer))){
      mapext <- raster::extent(rasterLayer)[1:4]
      if(mapext >= -180 && mapext <= 180){
        raster::crs(rasterLayer) <- sp::CRS("+init=epsg:4326")
        warning("The input layer has no CRS defined. Based on the extent of the input map it is assumed to have geographic coordinate system")
      }
    }
    # reduce the sampleNumber if the raster does not have enough cells
    if(raster::ncell(rasterLayer) < 10 * sampleNumber){
      rp <- raster::rasterToPoints(rasterLayer[[1]])
      if(nrow(rp) < sampleNumber){
        sampleNumber <- nrow(rp)
        message("The sample number reduced to ", sampleNumber, ", the total number of available cells")
      }
    }
    if(numLayer==1){
      if(is.null(speciesData)){
        rasterPoints <- raster::rasterToPoints(rasterLayer, spatial=TRUE)
        set.seed(2017)
        points <- rasterPoints[sample(1:nrow(rasterPoints), sampleNumber, replace=FALSE), ]
        names(points) <- 'target'
      } else{
        points <- raster::extract(rasterLayer, speciesData, na.rm=TRUE, sp=TRUE)
        names(points)[ncol(points)] <- "target"
      }
      fittedVar = automap::autofitVariogram(target~1, points)
      theRange <- fittedVar$var_model[2,3]
      if(plotVariograms==TRUE){
        plot(fittedVar)
      }
    } else if(numLayer>1){
      df <- data.frame(layers=1:numLayer, range=1:numLayer, sill=1:numLayer)
      variogramList <- list()
      message(paste('There are', numLayer, 'raster layers'))
      if(doParallel==TRUE){
        if(is.null(nCores)){
          nCores <- ceiling((parallel::detectCores()) / 2)
        }
        cl <- parallel::makeCluster(nCores) # use snow clusters
        doParallel::registerDoParallel(cl) # set up a parallel backend for foreach package
        pp <- foreach::foreach(r = 1:numLayer, .inorder=TRUE, .packages=c('raster', 'automap')) %dopar% {
          # rasterPoints <- raster::rasterToPoints(rasterLayer[[r]], spatial=TRUE)
          # set.seed(2017)
          # points <- rasterPoints[sample(1:nrow(rasterPoints), sampleNumber, replace=FALSE), ]
          # names(points) <- 'target'
          if(is.null(speciesData)){
            rasterPoints <- raster::rasterToPoints(rasterLayer[[r]], spatial=TRUE)
            set.seed(2017)
            points <- rasterPoints[sample(1:nrow(rasterPoints), sampleNumber, replace=FALSE), ]
            names(points) <- 'target'
          } else{
            points <- raster::extract(rasterLayer[[r]], speciesData, na.rm=TRUE, sp=TRUE)
            names(points)[ncol(points)] <- "target"
          }
          fittedVar <- automap::autofitVariogram(target~1, points)
        }
        for(v in 1:length(pp)){
          df$range[v] <- pp[[v]]$var_model[2,3]
          df$sill[v] <- pp[[v]]$var_model[2,2]
          df$layers[v] <- names(rasterLayer)[v]
        }
        variogramList <- pp # save variogram of all layer
        parallel::stopCluster(cl)
        foreach::registerDoSEQ()
      } else{
        if(progress==TRUE){
          pb <- progress::progress_bar$new(format = " Progress [:bar] :percent in :elapsed",
                                 total=numLayer, clear=FALSE, width=75) # add progress bar
        }
        for(r in 1:numLayer){
          name <- names(rasterLayer[[r]])
          # rasterPoints <- raster::rasterToPoints(rasterLayer[[r]], spatial=TRUE)
          # set.seed(2017)
          # points <- rasterPoints[sample(1:nrow(rasterPoints), sampleNumber, replace=FALSE), ]
          # names(points) <- 'target'
          if(is.null(speciesData)){
            rasterPoints <- raster::rasterToPoints(rasterLayer[[r]], spatial=TRUE)
            set.seed(2017)
            points <- rasterPoints[sample(1:nrow(rasterPoints), sampleNumber, replace=FALSE), ]
            names(points) <- 'target'
          } else{
            points <- raster::extract(rasterLayer[[r]], speciesData, na.rm=TRUE, sp=TRUE)
            names(points)[ncol(points)] <- "target"
          }
          fittedVar <- automap::autofitVariogram(target~1, points)
          variogramList[[r]] <- fittedVar
          df$range[r] <- fittedVar$var_model[2,3]
          df$sill[r] <- fittedVar$var_model[2,2]
          df$layers[r] <- name;
          if(progress==TRUE){
            pb$tick() # update progress bar
          }
        }
        variogramList <- variogramList # save variogram of all layer
      }
      # calculate all ranges and mean them for block size
      theRange <- stats::median(df$range)
      modelInfo <- df[order(df$range),] # save range and sill of all layers
      if(plotVariograms==TRUE){
        for(v in 1:numLayer){
          plot(variogramList[[v]])
        }
      }
    } else stop('The raster layer is empty!')
  } else stop('The input file is not a valid R raster file')
  # creating the blocks based on calculated autocorrelation range
  # check the spatial reference sytem for specifiying block size
  if(is.na(sp::proj4string(rasterLayer))){
    mapext <- raster::extent(rasterLayer)[1:4]
    if(mapext >= -180 && mapext <= 180){
      theRange2 <- theRange * 1000
      if(numLayer>1){
        modelInfo$range <- modelInfo$range * 1000
      }
      xaxes <- "Longitude"
      yaxes <- "Latitude"
    } else{
      theRange2 <- theRange
      xaxes <- "Easting"
      yaxes <- "Northing"
    }
  } else{
    if(sp::is.projected(sp::SpatialPoints((matrix(1:10, 5, byrow=FALSE)), proj4string=raster::crs(rasterLayer)))){
      theRange2 <- theRange
      xaxes <- "Easting"
      yaxes <- "Northing"
    } else{
      theRange2 <- theRange * 1000
      if(numLayer>1){
        modelInfo$range <- modelInfo$range * 1000
      }
      xaxes <- "Longitude"
      yaxes <- "Latitude"
    }
  }
  if(is.null(border)){
    subBlocks <- rasterNet(rasterLayer[[1]], resolution=theRange2, degree=degMetre, mask=TRUE, maxpixels =maxpixels)
  } else{
    net <- rasterNet(rasterLayer[[1]], resolution=theRange2, degree=degMetre, mask=FALSE)
    if(methods::is(border, "sf")){
      border <- sf::as_Spatial(border)
    }
    subBlocks <- raster::crop(net, border)
  }
  if(numLayer>1){
    if(is.null(speciesData)){
      ptnum <- sampleNumber
    } else{
      ptnum <- nrow(speciesData)
    }
    p1 <- ggplot2::ggplot(data=modelInfo, aes(x=stats::reorder(factor(layers), range), y=range, fill=range))+
      geom_bar(stat="identity")+
      xlab("Layers") + ylab("Range (m)") +
      theme(axis.text.x = element_text(angle=75, hjust=1)) +
      ggtitle('Autocorrelation range',
              subtitle=paste('Based on', ptnum, 'sample points'))+
      guides(fill=FALSE)+
      geom_hline(yintercept=theRange2, color='red', size=0.9, linetype=2)+
      annotate("text", x=floor(nrow(modelInfo)/3), y= (theRange2 + (max(modelInfo$range)/20)),
               label="Block size", color='red')
  }
  # plot raster file in ggplot2
  samp <- raster::sampleRegular(rasterLayer[[1]], 5e+05, asRaster=TRUE)
  map.df <- raster::as.data.frame(samp, xy=TRUE, centroids=TRUE, na.rm=TRUE)
  colnames(map.df) <- c("Easting", "Northing", "MAP")
  mid <- mean(map.df$MAP)
  p2 <- ggplot2::ggplot(data=map.df, aes(y=Northing, x=Easting)) +
    geom_raster(aes(fill=MAP)) +
    coord_fixed() +
    scale_fill_gradient2(low="darkred", mid="yellow", high="darkgreen", midpoint=mid) +
    guides(fill=FALSE) +
    ggtitle('Spatial blocks', subtitle=paste("Based on", round(theRange2), "metres distance")) +
    geom_polygon(aes(x = long, y = lat, group=id),
                 data = subBlocks, color ="red",
                 fill ="orangered4",
                 alpha = 0.04,
                 size = 0.2) +
    xlab(xaxes) + ylab(yaxes)
  if(showPlots==TRUE){
    if(numLayer>1){
      multiplot(p1, p2)
    } else{
      plot(p2)
    }
  }
  if(numLayer>1){
    finalList <- list(range=theRange2, rangeTable=modelInfo, plots=list(barchart = p1, mapplot = p2),
                      sampleNumber=sampleNumber, variograms=variogramList)
  } else{
    finalList <- list(range=theRange2, plots=list(mapplot = p2), sampleNumber=sampleNumber, variograms=fittedVar)
  }
  # gc() # to release the occupied RAM in windows OS
  # specify the output class
  class(finalList) <- c("SpatialAutoRange", class(finalList))
  return(finalList)
}


#' @export
print.SpatialAutoRange <- function(x, ...){
  print(class(x))
}


#' @export
plot.SpatialAutoRange <- function(x, y, ...){
  if(length(x$plots) == 2){
    multiplot(x$plots$barchart, x$plots$mapplot)
  } else{
    plot(x$plots$mapplot)
  }
}

#' @export
summary.SpatialAutoRange <- function(object, ...){
  print("Summary statistics of spatial autocorrelation ranges of all input layers")
  print(summary(object$rangeTable$range))
  print(object$rangeTable[,1:2])
}
