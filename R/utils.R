# A collection of function to facilitate block generation and reporting

# is it a LOO CV object?
.is_loo <- function(x){
    inherits(x, c("cv_buffer", "cv_nndm"))
}


# prepare the response classes/strata used for fold balance summaries
.column_response <- function(x, column, n_bins = NULL){
    if(is.null(column)){
        return(list(values = NULL, levels = NULL, bins = NULL))
    }

    values <- x[, column, drop = TRUE]
    if(!is.null(n_bins)){
        n_bins <- .check_n_bins(n_bins)
        if(.is_continuous_response(values, n_bins)){
            return(.quantile_response(values, n_bins))
        }
    }

    cl <- sort(unique(values))
    .check_classes(length(cl), column) # column should be binary or categorical
    list(values = values, levels = cl, bins = NULL)
}


.check_n_bins <- function(n_bins){
    if(length(n_bins) != 1){
        stop("'n_bins' must be a single integer value.")
    }
    n_bins <- suppressWarnings(as.integer(n_bins))
    if(is.na(n_bins) || n_bins < 2){
        stop("'n_bins' must be an integer value of 2 or higher.")
    }
    n_bins
}


.is_continuous_response <- function(values, n_bins){
    values <- values[!is.na(values)]
    if(!is.numeric(values) || length(values) == 0){
        return(FALSE)
    }
    n_unique <- length(unique(values))
    if(n_unique <= n_bins){
        return(FALSE)
    }
    any(values != floor(values)) || n_unique > 15L
}


.quantile_response <- function(values, n_bins){
    vals <- values[!is.na(values)]
    breaks <- unique(as.numeric(stats::quantile(
        vals,
        probs = seq(0, 1, length.out = n_bins + 1),
        na.rm = TRUE
    )))

    labels <- paste0("Q", seq_len(length(breaks) - 1))
    binned <- cut(
        values,
        breaks = breaks,
        include.lowest = TRUE,
        labels = labels,
        ordered_result = TRUE
    )
    bins <- data.frame(
        bin = labels,
        lower = breaks[-length(breaks)],
        upper = breaks[-1],
        stringsAsFactors = FALSE
    )
    attr(bins, "requested_bins") <- n_bins

    list(values = binned, levels = labels, bins = bins)
}


.records_table <- function(n, response){
    if(is.null(response$values)){
        return(data.frame(train = rep(0, n), test = 0))
    }

    cl <- response$levels
    tt_count <- as.data.frame(matrix(0, nrow = n, ncol = length(cl) * 2))
    names(tt_count) <- c(paste("train", cl, sep = "_"), paste("test", cl, sep = "_"))
    attr(tt_count, "column_bins") <- response$bins
    tt_count
}


.records_table_row <- function(tt_count, row, train_set, test_set, response){
    if(is.null(response$values)){
        tt_count$train[row] <- length(train_set)
        tt_count$test[row] <- length(test_set)
        return(tt_count)
    }

    cl <- as.character(response$levels)
    clen <- length(cl)
    countrain <- table(response$values[train_set])
    countest <- table(response$values[test_set])
    tt_count[row, match(names(countrain), cl)] <- as.integer(countrain)
    tt_count[row, clen + match(names(countest), cl)] <- as.integer(countest)
    attr(tt_count, "column_bins") <- response$bins
    tt_count
}


.print_column_bins <- function(response){
    if(!is.null(response$bins)){
        cat("Quantile bins used for the continuous 'column':\n")
        print(response$bins, row.names = FALSE)
    }
}


# count the train and test records
.table_summary <- function(fold_list, x, column, n, n_bins = NULL){
    response <- .column_response(x, column, n_bins = n_bins)
    tt_count <- .records_table(n, response)
    for(i in seq_len(n)){
        train_set <- fold_list[[i]][[1]]
        test_set <- fold_list[[i]][[2]]
        tt_count <- .records_table_row(tt_count, i, train_set, test_set, response)
    }

    return(tt_count)
}

