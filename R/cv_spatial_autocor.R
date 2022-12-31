
cv_spatial_autocor <- function(x,
                               r,
                               column = NULL,
                               num_sample = 5000L,
                               deg_to_metre = 111325,
                               # parallel = FALSE,
                               plot = TRUE,
                               progress = TRUE,
                               ... # extra arguments for cv_plot
){

  # check for availability of required packages
  pkg <- c("ggplot2", "cowplot", "automap", "terra")
  .pkg_checks(pkg)

  # check x is an sf object
  if(!missing(x)){
    if(!methods::is(x, "sf")){
      tryCatch(
        {
          x <- sf::st_as_sf(x)
        },
        error = function(cond) {
          message("'x' is not convertible to an sf object!")
          message("'x' must be an sf or spatial* object.")
        }
      )
    }
  }

  # x and column must be provided
  if(!missing(x) && is.null(column)){
    stop("When 'x' is provided, 'column' must also be provided. Otherwise, provide only 'r'.")
  }

  # is column in x?
  if(!missing(x) && !is.null(column)){
    if(!column %in% colnames(x)){
      stop(sprintf("There is no column named '%s' in 'x'.\n", column))
    }
  }

  # change the r to terra object
  if(!missing(r)){
    if(!methods::is(r, "SpatRaster")){
      tryCatch(
        {
          r <- terra::rast(r)
        },
        error = function(cond) {
          message("'r' is not convertible to a terra SpatRaster object!")
          message("'r' must be a SpatRaster, stars, Raster* object, or (multiple) path to raster files on disk.")
        }
      )
    }
  }

  if(!missing(r) && missing(x)){
    if(any(terra::is.factor(r))){
      stop("'r' should not include any factor layers!")
    }

    nlayer <- terra::nlyr(r)

    # assume wgs84 if crs is na
    if(is.na(sf::st_crs(r))){
      mapext <- terra::ext(r)[1:4]
      if(all(mapext >= -180) && all(mapext <= 180)){
        terra::crs(r) <- "epsg:4326"
        warning("The input layer CRS is undefined! Based on the extent of the input map it is assumed to have geographic coordinate system (epsg:4326)")
      } else{
        stop("Please provide a raster object with a defined CRS.")
      }
    }

    # reduce the num_sample if the raster does not have enough cells
    if(terra::ncell(r) < 5 * num_sample){
      rp <- length(terra::cells(r))
      if(rp < num_sample){
        num_sample <- rp
        message("The num_sample reduced to ", num_sample, "; the total number of available cells.\n")
      }
    }

    # to be used for filtering blocks
    samp_point <- terra::spatSample(r[[1]],
                                    size = 5e4,
                                    method = "regular",
                                    as.points = TRUE,
                                    na.rm = TRUE)
  }

  if(!missing(x)){
    if(is.na(sf::st_crs(x))){
      stop("The coordinate reference system of 'x' must be defined.")
    }
    # no raster is provided
    nlayer <- 1
  }

  # if(is.null(n_cores)){
  #   n_cores <- ceiling(future::availableCores() / 2)
  # }
  # future::plan("multisession", workers = n_cores)
  # vario_list <- future.apply::future_lapply(
  #   seq_along(nlayer),
  #   .fit_variogram,
  #   rr = rpoints,
  #   xx = NULL,
  #   column = NULL,
  #   rnames = names(r),
  #   rcrs=terra::crs(r),
  #   future.seed = NULL
  # )
  # # future::plan("sequential")

  if(missing(x) && progress) pb <- txtProgressBar(min = 0, max = nlayer, style = 3)

  # fitting wariogram models
  vario_list <- lapply(
    seq_len(nlayer),
    .fit_variogram,
    rr = switch(missing(x), r, NULL),
    xx = switch(!missing(x), x, NULL),
    column = switch(!missing(x), column, NULL),
    num_sample = num_sample,
    progress = ifelse(missing(x) && progress, TRUE, FALSE),
    pb = if(missing(x) && progress) pb else NULL
  )

  # make a dataframe from variograms data
  vario_data <- data.frame(layers = seq_len(nlayer), range = 1, sill = 1)
  for(v in seq_along(vario_list)){
    vario_data$layers[v] <- if(missing(x)) names(r)[v] else column
    vario_data$range[v] <- vario_list[[v]]$var_model[2,3]
    vario_data$sill[v] <- vario_list[[v]]$var_model[2,2]
  }

  # order them for plotting
  vario_data <- vario_data[order(vario_data$range), ] # save range and sill of all layers

  x_obj <- if(missing(x)) sf::st_as_sf(samp_point) else x

  size <- the_range <- stats::median(vario_data$range)

  if(sf::st_is_longlat(x_obj)){
    vario_data$range <- vario_data$range * 1000
    the_range <- the_range * 1000
    size <- the_range / deg_to_metre
  }

  # make an object for plotting
  vis_block <- sf::st_make_grid(x_obj, cellsize = round(size), what = "polygons")
  vis_block <- sf::st_sf(vis_block[x_obj])
  vis_block$folds <- 1:nrow(vis_block)
  plot_data <- list(blocks = vis_block)
  class(plot_data) <- "cv_spatial"

  if(nlayer > 1){
    ptnum <- ifelse(missing(x), num_sample, nrow(x))
    p1 <- .make_bar_plot(vario_data, the_range, ptnum)
  }

  # plot the spatial blocks
  p2 <- cv_plot(
    cv = plot_data,
    r = switch(!missing(r), r, NULL),
    label_size = 0,
    ...
  )
  p2 <- p2 + ggplot2::ggtitle(label = "Spatial blocks",
                     subtitle = paste("Using",
                                      round(the_range),
                                      "metres block size"))


  if(plot){
    if(nlayer > 1){
      plot(cowplot::plot_grid(p1, p2))
    } else{
      plot(p2)
    }
  }

  final_list <- list(range = the_range,
                     range_table = vario_data,
                     plots = if(nlayer > 1) list(barchart = p1, map_plot = p2) else p2,
                     num_sample = num_sample,
                     variograms = vario_list)

  # specify the output class
  class(final_list) <- c("cv_spatial_autocor")
  return(final_list)
}


#' @export
#' @method print cv_spatial_autocor
print.cv_spatial_autocor <- function(x, ...){
  print(class(x))
}


#' @export
#' @method plot cv_spatial_autocor
plot.cv_spatial_autocor <- function(x, y, ...){
  if(length(x$plots) == 2){
    plot(cowplot::plot_grid(x$plots$barchart, x$plots$mapplot))
  } else{
    plot(x$plots$mapplot)
  }
}

#' @export
#' @method summary cv_spatial_autocor
summary.cv_spatial_autocor <- function(object, ...){
  cat("Summary statistics of spatial autocorrelation ranges of all input layers:\n")
  print(summary(object$rangeTable$range))
  print(object$rangeTable[,1:2])
}
