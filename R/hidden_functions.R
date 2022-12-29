# generate fold number in checkerboard pattern or systematic
.fold_assign <- function(blocks, n, checkerboard=FALSE){
  cent <- sf::st_centroid(blocks)
  xy <- as.data.frame(sf::st_coordinates(cent))
  # get the dimension of blocks
  nx <- length(unique(xy$X))
  ny <- length(unique(xy$Y))
  len <- nrow(xy)
  # get the order of y dimension
  ids <- seq_len(len)
  ids_y <- ids[order(xy$Y)]
  # generate fold ids
  z <- c()
  if(checkerboard){
    if(nx %% 2){
      z <- rep(1:2, length.out = len)
    } else{
      for(i in 1:ny){
        if(i %% 2){
          z <- c(z, rep(1:2, length.out = nx))
        } else{
          z <- c(z, rep(2:1, length.out = nx))
        }
      }
    }
  } else{
    z <- rep(1:n, length.out = len)
  }
  # reorder folds by y dimension
  z <- z[order(ids_y)]
  # add them to the blocks
  blocks <- sf::st_sf(blocks)
  blocks$block_id <- 1:nrow(blocks)
  blocks$folds <- z

  return(blocks)
}

# transform x and fold numbers for plotting
.x_to_long <- function(x, cv, num_plot=1:10){

  folds_list <- cv$folds_list

  # The iteration must be a natural number
  tryCatch(
    {
      num_plot <- abs(as.integer(num_plot))
      num_plot <- sort(num_plot)
    },
    error = function(cond) {
      message("'num_plot' must be a natural numbers.")
    }
  )

  # length of the folds
  k <- length(folds_list)
  if(max(num_plot) > k){
    num_plot <- num_plot[num_plot <= k]
    # message("'max(num_plot)' is higher than k (the number of folds)! Only numbers lower than k are retianed.")
  }

  if(class(cv) != "cv_spatial_loo"){
    len <- length(unlist(folds_list[[1]]))
    if(len != nrow(x)){
      stop("Number of rows in `x` does not match the folds in `cv_object`!")
    }
  }


  df <- data.frame(id = seq_len(len))

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
                              v.names="value",
                              timevar="folds"
  )
  # convert back to sf
  x_long <- sf::st_as_sf(x_reshape, sf_column_name = sf_colname)
  # convert to factor for plotting
  x_long$value <- as.factor(x_long$value)
  levels(x_long$value) <- c("Test", "Train")

  return(x_long)
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
