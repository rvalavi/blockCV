# Author: Roozbeh Valavi
# contact: valavi.r@gmail.com
# Date : January 2023
# Version 0.1
# Licence GPL v3

# check for x
.x_check <- function(x, name = "x"){
  if(!methods::is(x, "sf")){
    tryCatch(
      {
        x <- sf::st_as_sf(x)
      },
      error = function(cond) {
        message(sprintf("'%s' is not convertible to an sf object!", name))
        message(sprintf("'%s' must be an sf or spatial* object.", name))
      }
    )
  }
  return(x)
}

# check for r
.r_check <- function(r, name = "r"){
  if(!methods::is(r, "SpatRaster")){
    tryCatch(
      {
        r <- terra::rast(r)
      },
      error = function(cond) {
        message(sprintf("'%s' is not convertible to a terra SpatRaster object!", name))
        message(sprintf("'%s' must be a SpatRaster, stars, Raster* object, or path to raster a file on disk.", name))
      }
    )
  }
  return(r)
}

# check for required packages
.pkg_check <- function(pkg){
  pkgna <- names(which(sapply(sapply(pkg, find.package, quiet = TRUE), length) == 0))
  if(length(pkgna) > 0){
    nm <- paste(pkgna, collapse = ", ")
    message("This function requires these packages: ", nm, "\nWould you like to install them now?\n1: yes\n2: no")
    user <- readline(prompt = paste0("Selection: "))
    if(tolower(user) %in% c("1", "yes", "y")){
      utils::install.packages(pkgna)
    } else{
      stop("Please install these packages for function to work: ", nm)
    }
  }
}


# generate fold number in checkerboard pattern or systematic
.fold_assign <- function(blocks, n, checkerboard=FALSE){
  # solve problems of digits in R
  old <- options("digits") # save away original options
  options(digits = 22) # change the option
  on.exit(options(old))
  # compute centroids
  cent <- sf::st_centroid(blocks)
  xy <- as.data.frame(sf::st_coordinates(cent))
  # to avoid problem of digits precision
  xy$X <- as.factor(xy$X)
  xy$Y <- as.factor(xy$Y)
  # get the dimension of blocks
  # xlev <- levels(xy$X)
  ylev <- levels(xy$Y)
  ny <- length(ylev)

  len <- nrow(xy)
  xy$ids <- seq_len(len)
  xy <- xy[order(xy$Y), ]
  xy$z <- 0

  if(checkerboard){
    for(i in rev(seq_len(ny))){
      wyi <- which(xy$Y == ylev[i])
      nx <- length(wyi)
      if(i %% 2){
        xy$z[wyi] <- rep(1:2, length.out = nx)
      } else{
        xy$z[wyi] <- rep(2:1, length.out = nx)
      }
    }
  } else{
    xy$z <- rep(1:n, length.out = len)
  }

  blocks <- sf::st_sf(blocks)
  blocks$block_id <- 1:nrow(blocks)
  xy <- xy[order(xy$ids), ]
  blocks$folds <- xy$z

  return(blocks)
}


# transform x and fold numbers for plotting
.x_to_long <- function(x, cv, num_plot=1:10){
  # get the folds list
  folds_list <- cv$folds_list
  # The iteration must be a natural number
  tryCatch(
    {
      num_plot <- abs(as.integer(num_plot))
      num_plot <- sort(num_plot)
    },
    error = function(cond) {
      message("'num_plot' must be natural numbers.")
    }
  )
  # length of the folds
  k <- length(folds_list)
  if(max(num_plot) > k){
    num_plot <- num_plot[num_plot <= k]
  }
  # get the length of unique ids
  if(methods::is(cv, "cv_buffer")){
    len <- length(unique(unlist(cv$folds_list)))
  } else{
    len <- length(unlist(folds_list[[1]]))
    if(len != nrow(x)){
      stop("Number of rows in 'x' does not match the folds in 'cv'!")
    }
  }
  # create a dataframe temp
  df <- data.frame(id = seq_len(len))
  # make the indices in x
  for(i in num_plot){
    df[, paste("Fold", i, sep = "")] <- NA
    test <- folds_list[[i]][[2]]
    train <- folds_list[[i]][[1]]
    df[test, paste("Fold", i, sep = "")] <- 0
    df[train, paste("Fold", i, sep = "")] <- 1
  }
  # get the geometry column name
  sf_colname <- attr(x, "sf_column")
  # cbind x and the df with fold ids
  xf <- cbind(x, df)
  # convert to dataframe for reshaping
  x_df <- as.data.frame(xf)
  # name of columns to rehspae long
  fold_names <- paste("Fold", num_plot, sep = "")
  # reshape x-df to long
  x_reshape <- stats::reshape(x_df,
                              direction = "long",
                              idvar = "id",
                              varying = fold_names,
                              times = fold_names,
                              v.names = "value",
                              timevar = "folds"
  )
  # convert back to sf
  x_long <- sf::st_as_sf(x_reshape, sf_column_name = sf_colname)
  # convert to factor for plotting
  x_long$value <- as.factor(x_long$value)
  levels(x_long$value) <- c("Test", "Train")

  return(x_long)
}


# make a bar plot for cv_spatial_autocor
.make_bar_plot <- function(vario_data, the_range, ptnum){
  # change the scale to km
  vario_data$range <- vario_data$range / 1000
  the_range <- the_range / 1000

  p <- ggplot2::ggplot(
    data = vario_data,
    ggplot2::aes(y = get("range"),
                 x = stats::reorder(factor(get("layers")),
                                    get("range"),
                                    decreasing = FALSE),
                 color = get("range"))
    ) +
    # ggplot2::geom_bar(
    #   ggplot2::aes(x = stats::reorder(factor(get("layers")), get("range")),
    #                y = get("range"),
    #                fill = get("range")),
    #   stat = "identity", data = vario_data,) +
    ggplot2::geom_point(size = 4) +
    ggplot2::geom_segment(
      ggplot2::aes_string(x = "layers",
                          xend = "layers",
                          y = 0,
                          yend = "range"),
      size = 1.5
    ) +
    ggplot2::labs(x = "Raster layers", y = "Range (km)") +
    ggplot2::theme_bw() +
    ggplot2::ggtitle("Autocorrelation range", subtitle = paste("Based on", ptnum, "sample points"))+
    ggplot2::guides(color = "none") +
    ggplot2::geom_hline(yintercept = the_range, color = 'red', size = 0.5, linetype = 2) +
    ggplot2::annotate("text", x = floor(nrow(vario_data) / 3),
                      y =  (the_range + (max(vario_data$range) / 20)),
                      angle = 270,
                      label = "Block size",
                      color = 'red') +
    ggplot2::coord_flip()

  return(p)
}


# auto-fit variogram models
.fit_variogram <- function(i,
                           xx = NULL,
                           rr = NULL,
                           column = NULL,
                           num_sample = 1e4,
                           progress = FALSE,
                           pb = NULL){
  if(is.null(xx)){
    points <- terra::spatSample(
      x = rr[[i]],
      size = num_sample,
      method = "random",
      as.points = TRUE,
      na.rm = TRUE
    )
    points <- sf::as_Spatial(sf::st_as_sf(points))
    names(points) <- "target"
  } else{
    points <- xx[column]
    points <- sf::as_Spatial(points)
    names(points) <- "target"
  }
  fit_vario <- automap::autofitVariogram(target~1, points)
  if(progress) utils::setTxtProgressBar(pb, i)

  return(fit_vario)
}

