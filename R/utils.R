# A collection of function to facilitate block generaion and reporting

# is it a LOO CV object?
.is_loo <- function(x){
    inherits(x, c("cv_buffer", "cv_nndm"))
}


# count the train and test records
.table_summary <- function(fold_list, x, column, n){
    if(is.null(column)){
        tt_count <- data.frame(train = rep(0, n), test = 0)
        for(i in seq_len(n)){
            train_set <- fold_list[[i]][[1]]
            test_set <- fold_list[[i]][[2]]
            tt_count$train[i] <- length(train_set)
            tt_count$test[i] <- length(test_set)
        }
    } else{
        cl <- sort(unique(x[, column, drop = TRUE]))
        clen <- length(cl)
        .check_classes(clen, column) # column should be binary or categorical
        tt_count <- as.data.frame(matrix(0, nrow = n, ncol = clen * 2))
        names(tt_count) <- c(paste("train", cl, sep = "_"), paste("test", cl, sep = "_"))
        for(i in seq_len(n)){
            train_set <- fold_list[[i]][[1]]
            test_set <- fold_list[[i]][[2]]
            countrain <- table(x[train_set, column, drop = TRUE])
            countest <- table(x[test_set, column, drop = TRUE])
            tt_count[i, which(cl %in% names(countrain))] <- countrain
            tt_count[i, clen + which(cl %in% names(countest))] <- countest
        }
    }

    return(tt_count)
}

