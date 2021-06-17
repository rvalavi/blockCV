#' Measure spatial autocorrelation in the predictor raster files
#'
#' This function provides a quantitative basis for choosing block size. The spatial autocorrelation in all continuous
#' predictor variables available as raster layers is assessed and reported. The function estimates spatial autocorrelation
#' ranges of all input raster layers. This is the range over which observations are independent and is determined by
#' constructing the empirical variogram, a fundamental geostatistical tool for measuring spatial autocorrelation.
#' The empirical variogram models the structure of spatial autocorrelation by measuring variability between all possible
#' pairs of points (O'Sullivan and Unwin, 2010). Results are plotted. See the details section for further information.
#'
#'
#' The input raster layers should be continuous for computing the variograms and estimating the range of spatial
#' autocorrelation. The input rasters should also have a specified coordinate reference system. However, if the reference
#' system is not specified, the function attempts to guess it based on the extent of the map. It assumes an unprojected
#' reference system for layers with extent lying between -180 and 180, and a projected reference system otherwise.
#'
#' Variograms are calculated based on the distances between pairs of points, so unprojected rasters (in degrees) will
#' not give an accurate result (especially over large latitudinal extents). For unprojected rasters, \emph{the great circle distance} (rather
#' than Euclidian distance) is used to calculate the spatial distances between pairs of points. To
#' enable more accurate estimate, it is recommended to transform unprojected maps (geographic coordinate
#' system / latitude-longitude) to a projected metric reference system (e.g. UTM or Lambert) where it is possible.
#' See \code{\link[automap]{autofitVariogram}} from \pkg{automap} and \code{\link[gstat]{variogram}} from \pkg{gstat} packages
#' for further information.
#'
#'
#' @param rasterLayer A raster object of covariates to find spatial autocorrelation range.
#' @param sampleNumber Integer. The number of sample points of each raster layer to fit variogram models. It is 5000 by default,
#' however it can be increased by user to represent their region well (relevant to the extent and resolution of rasters).
#' @inheritParams spatialBlock
#' @param speciesData A spatial or sf object (optional). If provided, the \code{sampleNumber} is ignored and
#' variograms are created based on species locations. This option is not recommended if the species data is not
#' evenly distributed across the whole study area and/or the number of records is low.
#' @param showPlots Logical. Show final plot of spatial blocks and autocorrelation ranges.
#' @param maxpixels Number of random pixels to select the blocks over the study area.
#' @param plotVariograms Logical. Plot fitted variograms. This can also be done after the analysis.
#' It is \code{FALSE} by default.
#' @param doParallel Logical. Run in parallel when more than one raster layer is available. Given multiple CPU cores, it is
#' recommended to set it to \code{TRUE} when there is a large number of rasters to process.
#' @param nCores Integer. Number of CPU cores to run in parallel. If \code{nCores = NULL} half of available cores in your
#' machine will be used.
#' @param progress Logical. Shows progress bar. It works only when \code{doParallel = FALSE}.
#' @param degMetre Numeric. The conversion rate of metres to degree. This is for constructing spatial
#' blocks for visualisation. When the input map is in geographic coordinate system (decimal degrees), the block size is
#' calculated based on deviding the calculated \emph{range} by this value to convert to the input map's unit
#' (by default 111325; the standard distance of a degree in metres, on the Equator).
#'
#' @references O'Sullivan, D., Unwin, D.J., (2010). Geographic Information Analysis, 2nd ed. John Wiley & Sons.
#'
#' Roberts et al., (2017). Cross-validation strategies for data with temporal, spatial, hierarchical,
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
#' \donttest{
#'
#' # load the example raster data
#' awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))
#'
#' # run the model in parallel
#' range1 <- spatialAutoRange(rasterLayer = awt,
#'                            sampleNumber = 5000, # number of cells to be used
#'                            doParallel = TRUE,
#'                            nCores = 2, # if NULL, it uses half of the CPU cores
#'                            plotVariograms = FALSE,
#'                            showPlots = TRUE)
#'
#'
#' # run the model with no parallel
#' range3 <- spatialAutoRange(rasterLayer = awt,
#'                            sampleNumber = 5000,
#'                            doParallel = FALSE,
#'                            progress = TRUE)
#'
#' # show the result
#' summary(range1)
#' }
spatialAutoRange <- function(rasterLayer,
                             sampleNumber = 5000L,
                             border = NULL,
                             speciesData = NULL,
                             doParallel = FALSE,
                             nCores = NULL,
                             showPlots = TRUE,
                             degMetre = 111325,
                             maxpixels = 1e+05,
                             plotVariograms = FALSE,
                             progress = TRUE){
  # check for availability of required packages
  if(doParallel){
    pkg <- c("ggplot2", "cowplot", "automap", "future.apply")
  } else{
    pkg <- c("ggplot2", "cowplot", "automap")
  }
  pkgna <- names(which(sapply(sapply(pkg, find.package, quiet = TRUE), length) == 0))
  if(length(pkgna) > 0){
    nm <- paste(pkgna, collapse = ", ")
    message("This function requires these packages: ", nm, "\nWould you like to install them now?\n1: yes\n2: no")
    user <- readline(prompt = paste0("Selection: "))
    if(tolower(user) %in% c("1", "yes", "y")){
      utils::install.packages(pkgna)
    } else{
      stop("Please install these packages: ", nm)
    }
  }
  if(!is.null(speciesData)){
    if((methods::is(speciesData, "SpatialPoints") || methods::is(speciesData, "sf"))==FALSE){
      stop("speciesData should be a sf or SpatialPoints object")
    }
  }
  if(methods::is(rasterLayer, "Raster")){
    if(any(raster::is.factor(rasterLayer))){
      stop("rasterLayer should not include any factor layer")
    }
    numLayer <- raster::nlayers(rasterLayer)
    if(is.na(raster::projection(rasterLayer))){
      mapext <- raster::extent(rasterLayer)[1:4]
      if(all(mapext >= -180) && all(mapext <= 180)){
        raster::projection(rasterLayer) <- "+proj=longlat +datum=WGS84 +no_defs"
        warning("The input layer has no CRS defined. Based on the extent of the input map it is assumed to have geographic coordinate system")
      }
    }
    # reduce the sampleNumber if the raster does not have enough cells
    if(raster::ncell(rasterLayer) < 10 * sampleNumber && is.null(speciesData)){
      rp <- length(which(!is.na(raster::values(rasterLayer[[1]]))))
      if(rp < sampleNumber){
        sampleNumber <- rp
        message("The sampleNumber reduced to ", sampleNumber, "; the total number of available cells.\n")
      }
    }
    if(numLayer==1){
      if(is.null(speciesData)){
        rasterPoints <- raster::rasterToPoints(rasterLayer, spatial=TRUE)
        set.seed(2017)
        points <- rasterPoints[sample(nrow(rasterPoints), sampleNumber, replace=FALSE), ]
        names(points) <- "target"
      } else{
        points <- raster::extract(rasterLayer, speciesData, na.rm=TRUE, sp=TRUE)
        names(points)[ncol(points)] <- "target"
      }
      fittedVar <- automap::autofitVariogram(target~1, points)
      theRange <- fittedVar$var_model[2,3]
      if(plotVariograms){
        plot(fittedVar)
      }
    } else if(numLayer > 1){
      df <- data.frame(layers = seq_len(numLayer), range = 1, sill = 1)
      variogramList <- vector(mode = "list")
      message(paste("There are", numLayer, "raster layers\n"))
      if(doParallel){
        if(is.null(nCores)){
          nCores <- ceiling(future::availableCores() / 2)
        }
        suppressWarnings(future::plan("multiprocess", workers = nCores))
        pp <- future.apply::future_lapply(seq_len(numLayer),
                                          fitvario,
                                          spdata = speciesData,
                                          rdata = rasterLayer,
                                          sn = sampleNumber)
        future::plan("sequential")
        for(v in seq_len(length(pp))){
          df$range[v] <- pp[[v]]$var_model[2,3]
          df$sill[v] <- pp[[v]]$var_model[2,2]
          df$layers[v] <- names(rasterLayer)[v]
        }
        variogramList <- pp # save variogram of all layer
      } else{
        if(progress){
          pb <- progress::progress_bar$new(format = " Progress [:bar] :percent in :elapsed",
                                           total=numLayer, clear=FALSE, width=75) # add progress bar
        }
        for(r in seq_len(numLayer)){
          name <- names(rasterLayer[[r]])
          variogramList[[r]] <- fitvario(r, spdata = speciesData, rdata = rasterLayer, sn = sampleNumber)
          df$range[r] <- variogramList[[r]]$var_model[2,3]
          df$sill[r] <- variogramList[[r]]$var_model[2,2]
          df$layers[r] <- name;
          if(progress){
            pb$tick() # update progress bar
          }
        }
        variogramList <- variogramList # save variogram of all layer
      }
      # calculate all ranges and mean them for block size
      theRange <- stats::median(df$range)
      modelInfo <- df[order(df$range),] # save range and sill of all layers
      if(plotVariograms){
        for(v in seq_len(numLayer)){
          plot(variogramList[[v]])
        }
      }
    } else stop("The raster layer is empty!")
  } else stop("The input file is not a valid raster object")
  # creating the blocks based on calculated autocorrelation range
  if(is.na(raster::projection(rasterLayer))){
    mapext <- raster::extent(rasterLayer)[1:4]
    if(all(mapext >= -180) && all(mapext <= 180)){
      theRange2 <- theRange * 1000
      if(numLayer > 1){
        modelInfo$range <- modelInfo$range * 1000
      }
    } else{
      theRange2 <- theRange
    }
  } else{
    if(raster::isLonLat(rasterLayer)){
      theRange2 <- theRange * 1000
      if(numLayer > 1){
        modelInfo$range <- modelInfo$range * 1000
      }
    } else{
      theRange2 <- theRange
    }
  }
  if(is.null(border)){
    subBlocks <- rasterNet(rasterLayer[[1]], resolution=theRange2, degree=degMetre, mask=TRUE, maxpixels =maxpixels)
  } else{
    net <- rasterNet(rasterLayer[[1]], resolution=theRange2, degree=degMetre, mask=FALSE)
    if(!methods::is(border, "sf")){
      border <- sf::st_as_sf(border)
    }
    subBlocks <- sf::st_crop(net, border)
  }
  if(numLayer > 1){
    ptnum <- ifelse(is.null(speciesData), sampleNumber, nrow(speciesData))
    p1 <- ggplot2::ggplot() +
      ggplot2::geom_bar(ggplot2::aes(x=stats::reorder(factor(get("layers")), get("range")), y=get("range"), fill=get("range")),
                        stat="identity",data=modelInfo,) +
      ggplot2::labs(x = "Layers", y = "Range (m)") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle=75, hjust=1)) +
      ggplot2::ggtitle("Autocorrelation range", subtitle=paste("Based on", ptnum, "sample points"))+
      ggplot2::guides(fill = "none") +
      ggplot2::geom_hline(yintercept=theRange2, color='red', size=0.9, linetype=2) +
      ggplot2::annotate("text", x=floor(nrow(modelInfo)/3),
                        y =  (theRange2 + (max(modelInfo$range)/20)),
                        label="Block size", color='red')
  }
  # plot raster file in ggplot2
  samp <- raster::sampleRegular(rasterLayer[[1]], 5e+05, asRaster=TRUE)
  map_df <- raster::as.data.frame(samp, xy=TRUE, centroids=TRUE, na.rm=TRUE)
  colnames(map_df) <- c("Easting", "Northing", "MAP")
  mid <- stats::median(map_df$MAP)
  p2 <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = map_df, ggplot2::aes_string(y="Northing", x="Easting", fill="MAP")) +
    ggplot2::scale_fill_gradient2(low="darkred", mid="yellow", high="darkgreen", midpoint=mid) +
    ggplot2::guides(fill = "none") +
    ggplot2::geom_sf(data = subBlocks,
                     color ="red",
                     fill ="orangered4",
                     alpha = 0.04,
                     size = 0.2) +
    ggplot2::labs(x = "", y = "") + # set the axes labes to NULL
    ggplot2::ggtitle("Spatial blocks", subtitle = paste("Using", round(theRange2), "metres block size")) +
    ggplot2::theme_bw()
  if(showPlots){
    if(numLayer > 1){
      plot(cowplot::plot_grid(p1, p2))
    } else{
      plot(p2)
    }
  }
  if(numLayer > 1){
    finalList <- list(range = theRange2,
                      rangeTable = modelInfo,
                      plots = list(barchart = p1, mapplot = p2),
                      sampleNumber = sampleNumber,
                      variograms = variogramList)
  } else{
    finalList <- list(range = theRange2,
                      plots = list(mapplot = p2),
                      sampleNumber = sampleNumber,
                      variograms = fittedVar)
  }
  # specify the output class
  class(finalList) <- c("SpatialAutoRange")
  return(finalList)
}


#' @export
#' @method print SpatialAutoRange
print.SpatialAutoRange <- function(x, ...){
  print(class(x))
}


#' @export
#' @method plot SpatialAutoRange
plot.SpatialAutoRange <- function(x, y, ...){
  if(length(x$plots) == 2){
    plot(cowplot::plot_grid(x$plots$barchart, x$plots$mapplot))
  } else{
    plot(x$plots$mapplot)
  }
}

#' @export
#' @method summary SpatialAutoRange
summary.SpatialAutoRange <- function(object, ...){
  cat("Summary statistics of spatial autocorrelation ranges of all input layers:\n")
  print(summary(object$rangeTable$range))
  print(object$rangeTable[,1:2])
}
