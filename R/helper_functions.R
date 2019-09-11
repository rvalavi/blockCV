rasterNet <- function(x, resolution=NULL, xbin=NULL, ybin=NULL, mask=FALSE, degree=111325, xOffset=NULL, yOffset=NULL,
                      checkerboard=FALSE, maxpixels=250000){
  if(methods::is(x, "sf")){
    x <- sf::as_Spatial(x)
  }
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
  return(sf::st_as_sf(rasterNet))
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
    sdR <- sd(raster::values(x[[i]]), na.rm=TRUE)
    stzRaster <- raster::stack(stzRaster, ((x[[i]] - meanR) / sdR))
  }
  return(stzRaster)
}
