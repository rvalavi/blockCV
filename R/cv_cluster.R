#' Use environmental or spatial clustering to separate train and test folds
#'
#' This function uses clustering methods to specify sets of similar environmental
#' conditions based on the input covariates, or cluster of spatial coordinates of the sample data.
#' Sample data (i.e. species data) corresponding to any of
#' these groups or clusters are assigned to a fold. Clustering is done
#' using \code{\link[stats]{kmeans}} for both approaches. The only requirement is \code{x} that leads to
#' a clustering of the confidantes of sample data. Otherwise, by providing \code{r}, environmental
#' clustering is done.
#'
#' As k-means algorithms use Euclidean distance to estimate clusters, the input raster covariates should be quantitative variables.
#' Since variables with wider ranges of values might dominate the clusters and bias the environmental clustering (Hastie et al., 2009),
#' all the input rasters are first scaled and centred (\code{scale = TRUE}) within the function.
#'
#' If \code{raster_cluster = TRUE}, the clustering is done in the raster space. In this approach the clusters will be consistent throughout the region
#' and different sample datasets in the same region (for comparison). However, this may result in a cluster(s)
#' that covers none of the species records (the spatial location of response samples),
#' especially when species data is not dispersed throughout the region or the number of clusters (k or folds) is high. In this
#' case, the number of folds is less than specified \code{k}. If \code{raster_cluster = FALSE}, the clustering will be done in
#' species points and the number of the folds will be the same as \code{k}.
#'
#' Note that the input raster layer should cover all the species points, otherwise an error will rise. The records with no raster
#' value should be deleted prior to the analysis or another raster layer must be provided.
#'
#' @inheritParams cv_spatial
#' @param column character (optional). Indicating the name of the column in which response variable (e.g. species data as a binary
#'  response i.e. 0s and 1s) is stored. This is only used to see whether all the folds contain all the classes in the final report.
#' @param r a terra SpatRaster object of covariates to identify environmental groups. If provided, clustering will be done
#' in environmental space rather than spatial coordinates of sample points.
#' @param scale logical; whether to scale the input rasters (recommended) for clustering.
#' @param raster_cluster logical; if \code{TRUE}, the clustering is done over the entire raster layer,
#' otherwise it will be over the extracted raster values of the sample points. See details for more information.
#' @param num_sample integer; the number of samples from raster layers to build the clusters (when \code{raster_cluster = FALSE}).
#' @param ... additional arguments for \code{stats::kmeans} function, e.g. \code{algorithm = "MacQueen"}.
#'
#' @seealso \code{\link{cv_buffer}} and \code{\link{cv_spatial}}
#'
#' @references Hastie, T., Tibshirani, R., & Friedman, J. (2009). The elements of statistical learning: Data mining, inference, and prediction ( 2nd ed., Vol. 1).
#'
#' @return An object of class S3. A list of objects including:
#'    \itemize{
#'     \item{folds_list - a list containing the folds. Each fold has two vectors with the training (first) and testing (second) indices}
#'     \item{folds_ids - a vector of values indicating the number of the fold for each observation (each number corresponds to the same point in x)}
#'     \item{biomod_table - a matrix with the folds to be used in \pkg{biomod2} package}
#'     \item{k - number of the folds}
#'     \item{column - the name of the column if provided}
#'     \item{type - indicates whether spatial or environmental clustering was done.}
#'     \item{records - a table with the number of points in each category of training and testing}
#'     }
#' @export
#'
#' @examples
#' \donttest{
#' library(blockCV)
#'
#' # import presence-absence species data
#' points <- read.csv(system.file("extdata/", "species.csv", package = "blockCV"))
#' # make an sf object from data.frame
#' pa_data <- sf::st_as_sf(points, coords = c("x", "y"), crs = 7845)
#'
#' # load raster data
#' path <- system.file("extdata/au/", package = "blockCV")
#' files <- list.files(path, full.names = TRUE)
#' covars <- terra::rast(files)
#'
#' # spatial clustering
#' set.seed(6)
#' sc <- cv_cluster(x = pa_data,
#'                  column = "occ", # optional; name of the column with response
#'                  k = 5)
#'
#' # environmental clustering
#' set.seed(6)
#' ec <- cv_cluster(r = covars, # if provided will be used for environmental clustering
#'                  x = pa_data,
#'                  column = "occ", # optional; name of the column with response
#'                  k = 5,
#'                  scale = TRUE)
#'
#' }
cv_cluster <- function(
    x,
    column = NULL,
    r = NULL,
    k = 5L,
    scale = TRUE,
    raster_cluster = FALSE,
    num_sample = 10000L,
    biomod2 = TRUE,
    report = TRUE,
    ...
){

  # check x is an sf object
  x <- .check_x(x)
  # is column in x?
  column <- .check_column(column, x)

  if(!is.null(r)){
    # check for required packages
    pkg <- c("terra")
    .check_pkgs(pkg)
    # check r
    r <- .check_r(r)
    # check r layers
    if(terra::nlyr(r) < 1){
      stop("'r' is not a valid raster.")
    }
    # scale?
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


  if(is.null(r)){
    # spatial cluster on xy
    xy <- sf::st_coordinates(x)
    kms <- stats::kmeans(xy, centers = k, iter.max = 500, nstart = 25, ...)
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
      kms <- stats::kmeans(sampr, centers = k, iter.max = 500, nstart = 25, ...)
      x$fold <- kms$cluster[seq_len(nrow(x))]
    } else{
      kms <- stats::kmeans(x_vals, centers = k, iter.max = 500, nstart = 25, ...)
      x$fold <- kms$cluster
    }

  }

  # create train-test table
  if(is.null(column)){
    train_test_table <- data.frame(train = rep(0, k), test = 0)
  } else{
    cl <- sort(unique(x[, column, drop = TRUE]))
    clen <- length(cl)
    .check_classes(clen, column) # column should be binary or categorical
    train_test_table <- as.data.frame(matrix(0, nrow = k, ncol = clen * 2))
    names(train_test_table) <- c(paste("train", cl, sep = "_"), paste("test", cl, sep = "_"))
  }

  fold_list <- list()
  fold_ids <- rep(NA, nrow(x))
  biomod_table <- data.frame(RUN1 = rep(TRUE, nrow(x)))

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
    if(biomod2){ # creating a biomod2 CV.user.table for validation
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
    type = ifelse(is.null(r), "Spatial Cluster", "Environmental Cluster"),
    records = train_test_table
  )

  if(report) print(train_test_table)
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
  cat("Number of recoreds in each training and testing fold:\n")
  print(object$records)
}
