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


fitvario <- function(r, spdata, rdata, sn){
  if(is.null(spdata)){
    rasterPoints <- raster::rasterToPoints(rdata[[r]], spatial = TRUE)
    set.seed(2017)
    points <- rasterPoints[sample(nrow(rasterPoints), sn, replace = FALSE),]
    names(points) <- "target"
  } else{
    points <- raster::extract(rdata[[r]], spdata, na.rm = TRUE, sp = TRUE)
    names(points)[ncol(points)] <- "target"
  }
  fittedVar <- automap::autofitVariogram(target~1, points)
  return(fittedVar)
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
    sdR <- stats::sd(raster::values(x[[i]]), na.rm=TRUE)
    stzRaster <- raster::stack(stzRaster, ((x[[i]] - meanR) / sdR))
  }
  return(stzRaster)
}
