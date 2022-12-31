cv_buffer <- function(x,
                      column = NULL,
                      size,
                      presence_background = FALSE,
                      add_background = FALSE,
                      progress = TRUE){

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

  # x's CRS must be defined
  if(is.na(sf::st_crs(x))){
    stop("The coordinate reference system of x must be defined.")
  }
  # is column in x?
  if(!is.null(column)){
    if(!column %in% colnames(x)){
      warning(sprintf("There is no column named '%s' in 'x'.\n", column))
      column <- NULL
    }
  }

  dmatrix <- sf::st_distance(x)
  distuni <- dmatrix[1,1] # take the unit to avoid using units package
  distuni[1] <- size
  fold_list <- list()

  if(!is.null(column)){
    if(presence_background){
      unqsp <- unique(x[, column, drop = TRUE])
      if(!is.numeric(unqsp) || any(unqsp < 0) || any(unqsp > 1)){
        stop("Presence-background option is only for species data with 0s (backgrounds/pseudo-absences) and 1s (presences).\n", "The data should be numeric.\n")
      }
      prI <- which(x[, column, drop = TRUE] == 1) # presence indices to loop through
      n <- length(prI)
      cl <- sort(unique(x[, column, drop = TRUE]))
      clen <- length(cl)
      train_test_table <- as.data.frame(matrix(0, nrow = n, ncol = clen * 2))
      names(train_test_table) <- c(paste("train", cl, sep = "_"), paste("test", cl, sep = "_"))
      if(progress){
        pb <- txtProgressBar(min = 0, max = n, style = 3)
      }
      j <- 0
      for(i in prI){ # loop through presences
        j <- j + 1
        train_set <- which(dmatrix[i, ] > distuni)
        if(add_background){
          test <- which(dmatrix[i, ] <= distuni)
          inside <- x[test, ]
          test_set <- test[which(inside[, column, drop = TRUE] == 0)]
          test_set[length(test_set) + 1] <- i
        } else{
          test_set <- i
        }
        fold_list[[j]] <- assign(paste0("fold", j), list(train_set, test_set))
        countrain <- table(x[train_set ,column, drop = TRUE])
        countest <- table(x[test_set ,column, drop = TRUE])
        train_test_table[j, which(cl %in% names(countrain))] <- countrain
        train_test_table[j, clen + which(cl %in% names(countest))] <- countest
        if(progress){
          setTxtProgressBar(pb, i)
        }
      }
    } else{
      n <- nrow(x)
      cl <- sort(unique(x[, column, drop = TRUE]))
      clen <- length(cl)
      train_test_table <- as.data.frame(matrix(0, nrow = n, ncol = clen * 2))
      names(train_test_table) <- c(paste("train", cl, sep = "_"), paste("test", cl, sep = "_"))
      if(progress){
        pb <- txtProgressBar(min = 0, max = n, style = 3)
      }
      for(i in seq_len(n)){
        train_set <- which(dmatrix[i, ] > distuni)
        test_set <- i
        fold_list[[i]] <- assign(paste0("fold", i), list(train_set, test_set))
        countrain <- table(x[train_set ,column, drop = TRUE])
        countest <- table(x[test_set ,column, drop = TRUE])
        train_test_table[i, which(cl %in% names(countrain))] <- countrain
        train_test_table[i, clen + which(cl %in% names(countest))] <- countest
        if(progress){
          setTxtProgressBar(pb, i)
        }
      }
    }
  } else{ # data with no column column
    n <- nrow(x)
    train_test_table <- base::data.frame(train=rep(0, n), test=0)
    if(progress){
      pb <- txtProgressBar(min = 0, max = n, style = 3)
    }
    for(i in seq_len(n)){
      train_set <- which(dmatrix[i, ] > distuni)
      test_set <- i
      fold_list[[i]] <- assign(paste0("fold", i), list(train_set, test_set))
      train_test_table$train[i] <- length(train_set)
      train_test_table$test[i] <- length(test_set)
      if(progress){
        setTxtProgressBar(pb, i)
      }
    }
  }
  final_objs <- list(folds_list = fold_list,
                     k = n,
                     column = column,
                     size = size,
                     presence_background = presence_background,
                     records = train_test_table)
  class(final_objs) <- c("cv_buffer")
  return(final_objs)
}


#' @export
#' @method print cv_buffer
print.cv_buffer <- function(x, ...){
  print(class(x))
}

#' @export
#' @method summary cv_buffer
summary.cv_buffer <- function(object, ...){
  print("Number of recoreds in each category")
  print(object$records)
}
