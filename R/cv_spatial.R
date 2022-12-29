cv_spatial <- function(
    x,
    column = NULL,
    r = NULL,
    k = 5L,

    hexagon = TRUE,
    flat_top = FALSE,

    cell_size = NULL,
    rows_cols = c(10, 10),

    selection = "random",
    iteration = 50L,

    user_blocks = NULL,
    folds_column = NULL,

    deg_to_metre = 111325,
    biomod2 = TRUE,

    # xOffset = 0,
    # yOffset = 0,
    offset = c(0, 0),

    seed = NULL,

    progress = TRUE,
    print = TRUE,
    plot = TRUE
){

  # pre-run checks ----------------------------------------------------------

  if(plot){
    # check for availability of ggplot2
    pkg <- c("ggplot2")
    pkgna <- names(which(sapply(sapply(pkg, find.package, quiet = TRUE), length) == 0))
    if(length(pkgna) > 0){
      message("This function requires ", pkg, " package for plotting.", "\nWould you like to install it now?\n1: yes\n2: no")
      user <- readline(prompt = paste0("Selection: "))
      if(tolower(user) %in% c("1", "yes", "y")){
        utils::install.packages(pkgna)
      } else{
        stop("Please install ggplot2 package or set plot = FALSE.")
      }
    }
  }

  if(!is.element(selection, c("systematic", "random", "checkerboard", "predefined"))){
    stop("The selection argument must be 'random', 'systematic', 'checkerboard', or 'predefined'.")
  }

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

  # check for user_blocks format
  if(!is.null(user_blocks)){
    if(!methods::is(user_blocks, "sf")){
      tryCatch(
        {
          user_blocks <- sf::st_as_sf(user_blocks)
        },
        error = function(cond) {
          message("'user_blocks' is not convertible to an sf object!")
          message("'user_blocks' must be an sf or spatial* object.")
        }
      )
    }
  }

  # is column in x?
  if(!is.null(column)){
    if(!column %in% colnames(x)){
      warning(sprintf("There is no column named '%s' in 'x'.\n", column))
      column <- NULL
    }
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
    if(!methods::is(r, "SpatRaster")){
      tryCatch(
        {
          r <- terra::rast(r)
          r <- r[[1]]
        },
        error = function(cond) {
          message("'r' is not convertible to a terra SpatRaster object!")
          message("'r' must be a SpatRaster, stars, Raster* object, or path to raster a file on disk.")
        }
      )
    }
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

  if(progress){
    pb <- progress::progress_bar$new(format = " Progress [:bar] :percent in :elapsed",
                                     total=iteration, clear=FALSE, width=75) # add progress bar
  }

  # creating blocks ---------------------------------------------------------

  if(is.null(user_blocks)){
    # select the object to make grid
    x_obj <- if(is.null(r)) x else r

    if(is.null(cell_size)){
      blocks <- sf::st_make_grid(x_obj, n=rev(rows_cols), square=!hexagon, what="polygons", flat_topped=flat_top)
    } else{
      # convert metres to degrees
      if(sf::st_is_longlat(x_obj)) cell_size <- cell_size / deg_to_metre
      blocks <- sf::st_make_grid(x_obj, cellsize=cell_size, square=!hexagon, what="polygons", flat_topped=flat_top)
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
        pb$tick() # update progress bar
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
  if(print) print(train_test_table)

  if(any(train_test_table < 1)){
    zerofolds <- which(apply(train_test_table, 1, function(x) any(x < 1)))
    if(length(zerofolds) > 1){
      warning("The folds ", paste(zerofolds, collapse = ", "), " have class(es) with zero records")
    } else{
      warning("The fold ", zerofolds, " has class(es) with zero records")
    }
  }

  # final objects for exporting
  final_objs <- list(folds_list = fold_list,
                     folds_ids = fold_vect,
                     biomod_table = switch(biomod2, as.matrix(biomod_table), NULL),
                     k = k,
                     column = column,
                     blocks = sub_blocks,
                     records = train_test_table)

  class(final_objs) <- c("cv_spatial")

  # plot with the cv_plot function
  if(plot){
    p1 <- cv_ggplot(
      cv = final_objs,
      r = switch(!is.null(r), r, NULL)
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
  message("Please use cv_plot function to plot each fold interactively.")
}


#' @export
#' @method summary cv_spatial
summary.cv_spatial <- function(object, ...){
  cat("Number of recoreds in each training and testing fold:\n")
  print(object$records)
}

