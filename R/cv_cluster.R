cv_cluster <- function(x,
                       column = NULL,
                       r = NULL,
                       k = 5L,
                       scale = TRUE,
                       raster_cluster = FALSE,
                       num_sample = 10000L,
                       biomod2 = TRUE,
                       print = TRUE){

  # check x is an sf object
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

  # is column in x?
  if(!is.null(column)){
    if(!column %in% colnames(x)){
      warning(sprintf("There is no column named '%s' in 'x'.\n", column))
      column <- NULL
    }
  }

  # change the r to terra object
  if(!is.null(r)){
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

  if(!is.null(r)){
    if(terra::nlyr(r) < 1){
      stop("'r' is not a valid raster.")
    }
  }

  if(!is.null(r)){
    if(scale){
      tryCatch(
        {
          r <- terra::scale(r)
        },
        error = function(cond) {
          message("Scaling the raster failed!")
        }
      )
    }
  }

  # create train-test table
  if(is.null(column)){
    train_test_table <- data.frame(train = rep(0, k), test = 0)
  } else{
    cl <- sort(unique(x[, column, drop = TRUE]))
    clen <- length(cl)
    train_test_table <- as.data.frame(matrix(0, nrow = k, ncol = clen * 2))
    names(train_test_table) <- c(paste("train", cl, sep = "_"), paste("test", cl, sep = "_"))
  }

  fold_list <- list()
  fold_ids <- rep(NA, nrow(x))
  biomod_table <- data.frame(RUN1 = rep(TRUE, nrow(x)))

  if(is.null(r)){
    # spatial cluster on xy
    xy <- sf::st_coordinates(x)
    kms <- stats::kmeans(xy, centers = k, iter.max = 500, nstart = 25)
    x$fold <- as.integer(kms$cluster)

  } else{
    # a general check to make sure x is covered by r
    x_vals <- terra::extract(r, x, ID = FALSE)
    if(anyNA(x_vals)){
      stop("The input raster layer does not cover all the column points.")
    }

    if(raster_cluster){
      # check number of raster cells
      if(terra::ncell(r) < 5 * num_sample){
        rp <- length(terra::cells(r))
        if(rp < num_sample){
          num_sample <- rp
          message("The num_sample reduced to ", num_sample, "; the total number of available cells.\n")
        }
      }
      sampr <- terra::spatSample(r, size = num_sample, method = "random", na.rm = TRUE)
      sampr <- sampr[stats::complete.cases(sampr), ]
      sampr <- rbind(x_vals, sampr)
      kms <- stats::kmeans(sampr, centers = k, iter.max = 500, nstart = 25)
      x$fold <- kms$cluster[seq_len(nrow(x))]
    } else{
      kms <- stats::kmeans(x_vals, centers = k, iter.max = 500, nstart = 25)
      x$fold <- kms$cluster
    }

  }

  for(i in seq_len(k)){
    test_set <- which(x$fold == i)
    train_set <- which(x$fold != i)
    fold_ids[test_set] <- i
    fold_list[[i]] <- assign(paste0("fold", i), list(train_set, test_set))
    if(is.null(column)){
      train_test_table$train[i] <- length(train_set)
      train_test_table$test[i] <- length(test_set)
    } else{
      countrain <- table(x[train_set ,column, drop = TRUE])
      countest <- table(x[test_set ,column, drop = TRUE])
      train_test_table[i, which(cl %in% names(countrain))] <- countrain
      train_test_table[i, clen + which(cl %in% names(countest))] <- countest
    }
    if(biomod2){ # creating a biomod2 DataSplitTable for validation
      colm <- paste0("RUN", i)
      biomod_table[,colm] <- FALSE
      biomod_table[train_set, colm] <- TRUE
    }
  }

  # give a warning is any folds is empty
  zerofolds <- which(apply(train_test_table, 1, function(x) any(x < 1)))
  if(length(zerofolds) > 0){
    if(length(zerofolds) > 1){
      warning("Folds ", paste(zerofolds, collapse = ", "), " have class(es) with zero records")
    } else{
      warning("Fold ", zerofolds, " has class(es) with zero records")
    }
  }

  final_objs <- list(
    folds_list = fold_list,
    folds_ids = fold_ids,
    biomod_table = switch(biomod2, as.matrix(biomod_table), NULL),
    k = k,
    column = column,
    type = ifelse(is.null(r), "Spatial Cluster", "Covariate Cluster"),
    records = train_test_table
  )

  if(print) print(train_test_table)
  # specify the output class
  class(final_objs) <- c("cv_cluster")
  return(final_objs)
}


#' @export
#' @method print cv_cluster
print.cv_cluster <- function(x, ...){
  print(class(x))
}

#' @export
#' @method summary cv_cluster
summary.cv_cluster <- function(object, ...){
  print("Number of recoreds in each category")
  print(object$records)
}
