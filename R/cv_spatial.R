#' Use spatial blocks to separate train and test folds
#'
#' This function creates spatially separated folds based on a distance to number of row and/or column.
#' It assigns blocks to the training and testing folds \strong{randomly}, \strong{systematically} or
#' in a \strong{checkerboard pattern}. The distance (\code{size})
#' should be in \strong{metres}, regardless of the unit of the reference system of
#' the input data (for more information see the details section). By default,
#' the function creates blocks according to the extent and shape of the spatial sample data (\code{x} e.g. the species occurrence),
#' Alternatively, blocks can be created based on \code{r} assuming that the user has considered the
#' landscape for the given species and case study.
#' Blocks can also be offset so the origin is not at the outer corner of the rasters.
#' Instead of providing a distance, the blocks can also be created by specifying a number of rows and/or
#' columns and divide the study area into vertical or horizontal bins, as presented in Wenger & Olden (2012) and Bahn & McGill (2012).
#' Finally, the blocks can be specified by a user-defined spatial polygon layer.
#'
#' To keep the consistency, all the functions use \strong{metres} as their unit. In this function, when the input map
#' has geographic coordinate system (decimal degrees), the block size is calculated based on dividing \code{size} by
#' \code{deg_to_metre} (111325 as default, the standard distance of a degree in metres, on the Equator) to change
#' the unit to degree.
#'
#' The \code{offset} can be used to change the spatial position of the blocks. It can also be used to
#' assess the sensitivity of analysis results to shifting in the blocking arrangements.
#' These options are available when \code{size} is defined. By default the region is
#' located in the middle of the blocks and by setting the offsets, the blocks will shift.
#'
#' Roberts et. al. (2017) suggest that blocks should be substantially bigger than the range of spatial
#' autocorrelation (in model residual) to obtain realistic error estimates, while a buffer with the size of
#' the spatial autocorrelation range would result in a good estimation of error. This is because of the so-called
#' edge effect (O'Sullivan & Unwin, 2014), whereby points located on the edges of the blocks of opposite sets are
#' not separated spatially. Blocking with a buffering strategy overcomes this issue (see \code{\link{cv_buffer}}).
#'
#'
#' @param x a simple features (sf) or SpatialPoints object of spatial sample data (e.g., species data or ground truth sample for image classification).
#' @param column character (optional). Indicating the name of the column in which response variable (e.g. species data as a binary
#' response i.e. 0s and 1s) is stored to find balanced records in cross-validation folds. If \code{column = NULL}
#' the response variable classes will be treated the same and only training and testing records will be counted.
#' This is used for binary (e.g. presence-absence/background) or multi-class responses (e.g. land cover classes for
#' remote sensing image classification), and \emph{you can ignore it when the response variable is
#' continuous or count data}.
#' @param r a terra SpatRaster object (optional). If provided, its extent will be used to specify the blocks.
#' It also supports \emph{stars}, \emph{raster}, or path to a raster file on disk.
#' @param k integer value. The number of desired folds for cross-validation. The default is \code{k = 5}.
#' @param hexagon logical. Creates hexagonal (default) spatial blocks. If \code{FALSE}, square blocks is created.
#' @param flat_top logical. Creating hexagonal blocks with topped flat.
#' @param size numeric value of the specified range by which blocks are created and training/testing data are separated.
#' This distance should be in \strong{metres}. The range could be explored by \code{\link{cv_spatial_autocor}}
#' and \code{\link{cv_block_size}} functions.
#' @param rows_cols integer vector. Two integers to define the blocks based on row and
#' column e.g. \code{c(10, 10)} or \code{c(5, 1)}. Hexagonal blocks uses only the first one. This
#' option is ignored when \code{size} is provided.
#' @param selection type of assignment of blocks into folds. Can be \strong{random} (default), \strong{systematic}, \strong{checkerboard}, or \strong{predefined}.
#' The checkerboard does not work with hexagonal and user-defined spatial blocks. If the \code{selection = 'predefined'}, user-defined
#' blocks and \code{folds_column} must be supplied.
#' @param iteration integer value. The number of attempts to create folds with balanced records. Only works when \code{selection = "random"}.
#' @param user_blocks an sf or SpatialPolygons object to be used as the blocks (optional). This can be a user defined polygon and it must cover all
#' the species (response) points. If \code{selection = 'predefined'}, this argument and \strong{folds_column} must be supplied.
#' @param folds_column character. Indicating the name of the column (in \code{user_blocks}) in which the associated folds are stored.
#' This argument is necessary if you choose the 'predefined' selection.
#' @param deg_to_metre integer. The conversion rate of metres to degree. See the details section for more information.
#' @param biomod2 logical. Creates a matrix of folds that can be directly used in the \pkg{biomod2} package as
#' a \emph{DataSplitTable} for cross-validation.
#' @param offset two number between 0 and 1 to shift blocks by that proportion of block size.
#' This option only works when \code{size} is provided.
#' @param seed integer; a random seed for reproducibility.
#' @param progress logical; whether to shows a progress bar for random fold selection.
#' @param report logical; whether to print the report of the records per fold.
#' @param plot logical; whether to plot the final blocks with fold numbers in ggplot.
#' You can re-create this with \code{\link{cv_plot}}.
#' @param ... additional option for \code{\link{cv_plot}}.
#'
#'
#' @seealso \code{\link{cv_buffer}} and \code{\link{cv_cluster}}; \code{\link{cv_spatial_autocor}} and \code{\link{cv_block_size}} for selecting block size
#' @seealso For \emph{DataSplitTable} see \code{\link[biomod2]{BIOMOD_cv}} in \pkg{biomod2} package
#'
#' @references Bahn, V., & McGill, B. J. (2012). Testing the predictive performance of distribution models. Oikos, 122(3), 321-331.
#'
#' O'Sullivan, D., Unwin, D.J., (2010). Geographic Information Analysis, 2nd ed. John Wiley & Sons.
#'
#' Roberts et al., (2017). Cross-validation strategies for data with temporal, spatial, hierarchical,
#' or phylogenetic structure. Ecography. 40: 913-929.
#'
#' Wenger, S.J., Olden, J.D., (2012). Assessing transferability of ecological models: an underappreciated aspect of statistical
#' validation. Methods Ecol. Evol. 3, 260-267.
#'
#' @return An object of class S3. A list of objects including:
#'    \itemize{
#'     \item{folds_list - a list containing the folds. Each fold has two vectors with the training (first) and testing (second) indices}
#'     \item{folds_ids - a vector of values indicating the number of the fold for each observation (each number corresponds to the same point in species data)}
#'     \item{biomod_table - a matrix with the folds to be used in \pkg{biomod2} package}
#'     \item{k - number of the folds}
#'     \item{size - input size, if not null}
#'     \item{column - the name of the column if provided}
#'     \item{blocks - SpatialPolygon of the blocks}
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
#' # hexagonal spatial blocking by specified size and random assignment
#' sb1 <- cv_spatial(x = pa_data,
#'                   column = "occ",
#'                   size = 450000,
#'                   k = 5,
#'                   selection = "random",
#'                   iteration = 50)
#'
#' # spatial blocking by row/column and systematic fold assignment
#' sb2 <- cv_spatial(x = pa_data,
#'                   column = "occ",
#'                   rows_cols = c(8, 10),
#'                   k = 5,
#'                   hexagon = FALSE,
#'                   selection = "systematic")
#'
#' }
cv_spatial <- function(
    x,
    column = NULL,
    r = NULL,
    k = 5L,
    hexagon = TRUE,
    flat_top = FALSE,
    size = NULL,
    rows_cols = c(10, 10),
    selection = "random",
    iteration = 50L,
    user_blocks = NULL,
    folds_column = NULL,
    deg_to_metre = 111325,
    biomod2 = TRUE,
    offset = c(0, 0),
    seed = NULL,
    progress = TRUE,
    report = TRUE,
    plot = TRUE,
    ... # other arguments for cv_plot
){

  # pre-run checks ----------------------------------------------------------

  # check for availability of ggplot2
  if(plot){
    pkg <- c("ggplot2")
    .check_pkgs(pkg)
  }
  # check for selection arg
  selection <- match.arg(selection, choices = c("random", "systematic", "checkerboard", "predefined"))
  # check x is an sf object
  x <- .check_x(x)
  # is column in x?
  column <- .check_column(column, x)

  # check for user_blocks format
  if(!is.null(user_blocks)){
    user_blocks <- .check_x(user_blocks, name = user_blocks)
  }
  # checks for pre-defined folds
  if(selection == "predefined"){
    if(is.null(folds_column) || is.null(user_blocks)){
      stop("The 'user_blocks' and 'folds_column' should be specified for 'predefined' selection")
    }
    if(!folds_column %in% colnames(user_blocks)){
      stop(sprintf("There is no column named '%s' in 'user_blocks'.\n", folds_column))
    }
    if(!is.numeric(user_blocks[,folds_column, drop = TRUE])){
      stop("The fold numbers in 'folds_column' must be integer numbers.")
    }
  }

  # change the r to terra object
  if(!is.null(r)){
    r <- .check_r(r)
    r <- r[[1]]
  }

  # if hex; selection muse random or systematic
  if(hexagon && selection %in% c("checkerboard", "predefined")){
    selection <- "random"
    message("Hexagon blocks can only be used with random or systematic selection!")
    message("The random selection is used.")
  }

  if(selection == "checkerboard"){
    k <- 2
  }

  # iterations --------------------------------------------------------------

  # The iteration must be a natural number
  tryCatch(
    {
      iteration <- abs(as.integer(iteration))
      iteration <- max(1, iteration)
    },
    error = function(cond) {
      message("'iteration' must be a natural number.")
    }
  )

  # turn off progress if...
  if(selection != "random") progress <- FALSE
  if(iteration < 3) progress <- FALSE

  if(progress){
    pb <- utils::txtProgressBar(min = 0, max = iteration, style = 3)
  }

  # creating blocks ---------------------------------------------------------

  if(is.null(user_blocks)){
    # select the object to make grid
    x_obj <- if(is.null(r)) x else r

    if(is.null(size)){
      # select only row for hexagonal blocks
      if(hexagon) rows_cols <- rows_cols[1]
      blocks <- sf::st_make_grid(x_obj, n=rev(rows_cols), square=!hexagon, what="polygons", flat_topped=flat_top)
    } else{
      # convert metres to degrees
      if(sf::st_is_longlat(x_obj)) size <- size / deg_to_metre
      # prepare the offset values
      offset <- size * (abs(offset) %% 1)
      if(length(offset) < 2) offset[2] <- 0
      xm <- as.numeric(sf::st_bbox(x_obj)[1])
      ym <- as.numeric(sf::st_bbox(x_obj)[2])
      xoff <- xm - offset[1]
      yoff <- ym - offset[2]

      tryCatch(
        {
          blocks <- sf::st_make_grid(x_obj,
                                     cellsize = size,
                                     offset = c(xoff, yoff),
                                     square = !hexagon,
                                     what = "polygons",
                                     flat_topped = flat_top)
        },
        error = function(cond) {
          message("Could not create spatial blocks! possibly because of using a very small block size.")
          message("Remember, size is in metres not the unit of the CRS.")
        }
      )
    }
  } else{
    blocks <- sf::st_geometry(user_blocks)
  }

  ## subset the blocks by x
  sub_blocks <- blocks[x]
  blocks_len <- length(sub_blocks)

  # The iteration must be a natural number
  tryCatch(
    {
      k <- abs(as.integer(k))
    },
    error = function(cond) {
      message("'k' must be a natural number.")
    }
  )

  if(k > blocks_len){
    stop("'k' is bigger than the number of spatial blocks: ", blocks_len, ".\n")
  } else if(k < 2){
    stop("'k' must be a natural number equal or higher than 2.")
  }

  # x and block intersection ------------------------------------------------

  ## do the intersection once and outside of the loop
  blocks_df <- as.data.frame(sf::st_intersects(sf::st_geometry(x), sub_blocks))
  names(blocks_df) <- c("records", "block_id")
  # randomly remove the repeated records occurred on the edges of blocks
  if(nrow(blocks_df) > nrow(x)){
    if(!is.null(seed)){
      set.seed(seed)
    }
    blocks_df <- blocks_df[sample(nrow(blocks_df)), ]
    blocks_df <- blocks_df[!duplicated(blocks_df$records), ]
  } else if(nrow(blocks_df) < nrow(x) || anyNA(blocks_df)){
    nonoverlap <- nrow(x) - nrow(blocks_df)
    warning("At least ", nonoverlap, " of the points are not within the defined spatial blocks")
  }

  # create records table
  if(is.null(column)){
    train_test_table <- data.frame(train = rep(0, k), test = 0)
  } else{
    cl <- sort(unique(x[, column, drop = TRUE]))
    clen <- length(cl)
    train_test_table <- as.data.frame(matrix(0, nrow = k, ncol = clen * 2))
    names(train_test_table) <- c(paste("train", cl, sep = "_"), paste("test", cl, sep = "_"))
  }
  # create a table for biomod
  biomod_table <- data.frame(RUN1 = rep(TRUE, nrow(blocks_df)))

  # for selecting best iteration
  min_num <- 0
  max_sd <- Inf

  # iteration for creating folds --------------------------------------------

  if(!is.null(seed)){
    set.seed(seed)
  }
  # iteration if random selection, otherwise only 1 round
  for(i in seq_len(iteration)){

    if(selection=='checkerboard'){
      # use the full blocks for checkerboard pattern
      blocks <- .fold_assign(blocks, n = 2, checkerboard = TRUE)
      sub_blocks <- blocks[x, ]
      sub_blocks$block_id <- seq_len(nrow(sub_blocks))
      fold_df <- sf::st_drop_geometry(sub_blocks)
      blocks_df <- merge(x = blocks_df, y = fold_df, by = "block_id", all.x = TRUE)
    }

    if(selection=='systematic'){
      sub_blocks <- .fold_assign(sub_blocks, n = k)
      fold_df <- sf::st_drop_geometry(sub_blocks)
      blocks_df <- merge(x = blocks_df, y = fold_df, by = "block_id", all.x = TRUE)
    }

    if(selection=='predefined'){
      fold_df <- data.frame(block_id = seq_len(blocks_len), folds = user_blocks[, folds_column, drop = TRUE])
      blocks_df <- merge(x = blocks_df, y = fold_df, by = "block_id", all.x = TRUE)
    }

    if(selection=='random'){
      blocks_df <- blocks_df[, c("records", "block_id")] # to avoid repetition in iterations
      fold_df <- data.frame(block_id = seq_len(blocks_len), folds = 0)
      # create random folds with equal proportion
      num <- floor(blocks_len / k)
      fold_df$folds[seq_len(num * k)] <- sample(rep(seq_len(k), num), num * k)
      if(blocks_len %% k != 0){
        rest <- blocks_len %% k
        unfold <- which(fold_df$folds==0)
        fold_df$folds[unfold] <- sample(seq_len(k), rest, replace = FALSE)
      }
      blocks_df <- merge(x = blocks_df, y = fold_df, by = "block_id", all.x = TRUE)
    }

    # count the number of points in each fold
    fold_list <- list()
    fold_vect <- rep(NA, nrow(blocks_df))
    for(p in seq_len(k)){
      train_set <- blocks_df$records[which(blocks_df$folds != p)]
      test_set <- blocks_df$records[which(blocks_df$folds == p)]
      fold_vect[test_set] <- p
      fold_list[[p]] <- assign(paste0("fold", p), list(train_set, test_set))
      if(is.null(column)){
        train_test_table$train[p] <- length(train_set)
        train_test_table$test[p] <- length(test_set)
      } else{
        countrain <- table(x[train_set, column, drop = TRUE])
        countest <- table(x[test_set, column, drop = TRUE])
        train_test_table[p, which(cl %in% names(countrain))] <- countrain
        train_test_table[p, clen + which(cl %in% names(countest))] <- countest
      }
      if(biomod2){ # creating a biomod2 DataSplitTable for validation
        colm <- paste0("RUN", p)
        biomod_table[, colm] <- FALSE
        biomod_table[train_set, colm] <- TRUE
      }
    }

    # save the best folds in the iteration
    if(selection == "random"){
      if(min(train_test_table) >= min_num && stats::sd(unlist(train_test_table)) < max_sd){
        train_test_table2 <- train_test_table
        min_num <- min(train_test_table2)
        max_sd <- stats::sd(unlist(train_test_table))
        blocks_df2 <- blocks_df
        fold_list2 <- fold_list
        fold_vect2 <- fold_vect
        biomod_table2 <- biomod_table
        iter <- i
      }
      if(progress){ # if iteration is higher than 5?
        utils::setTxtProgressBar(pb, i)
      }
    } else{
      break
    }
  }

  if(selection == "random"){ # return the best bloks, table etc.
    sub_blocks <- sf::st_sf(sub_blocks)
    sub_blocks$block_id <- seq_len(nrow(sub_blocks))
    blocks_df_filter <- blocks_df2[,c("block_id","folds")]
    blocks_df_filter <- blocks_df_filter[!duplicated(blocks_df_filter), ]
    sub_blocks <- merge(x = sub_blocks, y = blocks_df_filter, by = "block_id", all.x = TRUE)
    blocks_df <- blocks_df2
    train_test_table <- train_test_table2
    fold_list <- fold_list2
    fold_vect <- fold_vect2
    biomod_table <- biomod_table2
  }
  cat("\n")
  if(report) print(train_test_table)

  if(any(train_test_table < 1)){
    zerofolds <- which(apply(train_test_table, 1, function(x) any(x < 1)))
    if(length(zerofolds) > 1){
      warning("The folds ", paste(zerofolds, collapse = ", "), " have class(es) with zero records")
    } else{
      warning("The fold ", zerofolds, " has class(es) with zero records")
    }
  }

  # remove the NA blocks
  sub_blocks <- sub_blocks[stats::complete.cases(sub_blocks$folds), ]

  # final objects for exporting
  final_objs <- list(
    folds_list = fold_list,
    folds_ids = fold_vect,
    biomod_table = switch(biomod2, as.matrix(biomod_table), NULL),
    k = k,
    size = size,
    column = column,
    blocks = sub_blocks,
    records = train_test_table
  )

  class(final_objs) <- c("cv_spatial")

  # plot with the cv_plot function
  if(plot){
    p1 <- cv_plot(
      cv = final_objs,
      r = switch(!is.null(r), r, NULL),
      ...
    )
    plot(p1)
  }

  return(final_objs)
}

#' @export
#' @method print cv_spatial
print.cv_spatial <- function(x, ...){
  print(class(x))
}


#' @export
#' @method plot cv_spatial
plot.cv_spatial <- function(x, y, ...){
  plot(x$blocks["folds"])
  message("Please use cv_plot function to plot each fold")
}


#' @export
#' @method summary cv_spatial
summary.cv_spatial <- function(object, ...){
  cat("Number of recoreds in each training and testing fold:\n")
  print(object$records)
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
