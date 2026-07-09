# A collection of function to facilitate block generation and reporting

# is it a LOO CV object?
.is_loo <- function(x){
    inherits(x, c("cv_buffer", "cv_nndm"))
}


# blend environmental covariates with the coordinates for spatially constrained
# environmental clustering. Each block is normalised by its number of columns and
# weighted so the k-means Euclidean distance on the result equals
# sqrt((1 - w) * d_env^2 + w * d_geo^2): w = 0 is pure environmental clustering and
# w = 1 clusters on the coordinates alone. The coordinates are always standardised
# (map units differ across axes and datasets); the covariates are standardised only
# when scale_env = TRUE, so scale = FALSE keeps the covariates on their native units.
# Constant (zero-variance) columns scale to NaN and are set to 0 so they contribute
# no distance.
.augment_geo <- function(env, xy, w, scale_env = TRUE){
    zscore <- function(m){
        m <- scale(as.matrix(m))
        m[is.nan(m)] <- 0
        m
    }
    env <- as.matrix(env)
    if(scale_env) env <- zscore(env)
    xy <- zscore(xy)
    env <- env * sqrt((1 - w) / NCOL(env))
    xy  <- xy  * sqrt(w / NCOL(xy))
    cbind(env, xy)
}


# prepare the response classes/strata used for fold balance summaries
.column_response <- function(x, column, num_bins = NULL){
    if(is.null(column)){
        return(list(values = NULL, levels = NULL, bins = NULL))
    }

    values <- x[, column, drop = TRUE]
    if(!is.null(num_bins)){
        num_bins <- .check_num_bins(num_bins)
        if(.is_continuous_response(values, num_bins)){
            return(.quantile_response(values, num_bins))
        }
    }

    cl <- sort(unique(values))
    .check_classes(length(cl), column) # column should be binary or categorical
    list(values = values, levels = cl, bins = NULL)
}


.check_num_bins <- function(num_bins){
    if(length(num_bins) != 1){
        stop("'num_bins' must be a single integer value.")
    }
    num_bins <- suppressWarnings(as.integer(num_bins))
    if(is.na(num_bins) || num_bins < 2){
        stop("'num_bins' must be an integer value of 2 or higher.")
    }
    num_bins
}


.is_continuous_response <- function(values, num_bins){
    values <- values[!is.na(values)]
    if(!is.numeric(values) || length(values) == 0){
        return(FALSE)
    }
    n_unique <- length(unique(values))
    if(n_unique <= num_bins){
        return(FALSE)
    }
    any(values != floor(values)) || n_unique > 15L
}


.quantile_response <- function(values, num_bins){
    vals <- values[!is.na(values)]
    breaks <- unique(as.numeric(stats::quantile(
        vals,
        probs = seq(0, 1, length.out = num_bins + 1),
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
    attr(bins, "requested_bins") <- num_bins

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


# shared pretty-printer for the fold-generating cv_* objects. Prints a one-line
# description, a set of "key: value" detail lines, and the per-fold train/test
# records table (capped for leave-one-out designs, which have one fold per
# point). Row names of the records table are the fold numbers, so they are kept.
# Returns invisible(x) so the object is not echoed a second time after the
# formatted output.
.print_cv_folds <- function(x, description, details = NULL, cap = 12L){
    cat(sprintf("blockCV %s: %s\n", class(x)[1], description))
    for(nm in names(details)){
        cat(sprintf("  %s: %s\n", nm, format(details[[nm]], scientific = FALSE)))
    }
    recs <- x$records
    if(!is.null(recs) && nrow(recs)){
        cat("\nNumber of records in each fold (train/test):\n")
        if(nrow(recs) > cap){
            print(utils::head(recs, cap))
            cat(sprintf("... %d more fold(s); use summary() for a full breakdown.\n",
                        nrow(recs) - cap))
        } else {
            print(recs)
        }
    }
    invisible(x)
}


# count the train and test records
.table_summary <- function(fold_list, x, column, n, num_bins = NULL){
    response <- .column_response(x, column, num_bins = num_bins)
    tt_count <- .records_table(n, response)
    for(i in seq_len(n)){
        train_set <- fold_list[[i]][[1]]
        test_set <- fold_list[[i]][[2]]
        tt_count <- .records_table_row(tt_count, i, train_set, test_set, response)
    }

    return(tt_count)
}


# test columns of the records table the balance search should optimise. For
# presence-background data only the presence (level 1) test cell is balanced so
# the many background records cannot dominate the objective; otherwise NULL is
# returned and every test column is used (see .balance_folds).
.balance_opt_cols <- function(response, presence_bg){
    if(!isTRUE(presence_bg) || is.null(response$values)){
        return(NULL)
    }
    if(!"1" %in% as.character(response$levels)){
        return(NULL)
    }
    "test_1"
}


# score a candidate fold assignment. The objective is evaluated on the test
# columns only: a train cell is the complement of its test cell
# (train = n_c - test), so per class the two have identical fold-to-fold
# variance -- balancing test balances train -- and dropping train stops the much
# larger train counts from dominating the objective. Returns a length-2 vector
# compared lexicographically (lower is better on each element):
#   [1] empty_cells -- empty test cells for classes with at least k records.
#       This is count-based; block layout may still prevent full coverage.
#   [2] imbalance  -- Pearson chi-square of the test counts against an equal
#       split (expected e_c = n_c / k), summed over classes. Dividing by e_c
#       normalises each class so an abundant class cannot swamp a rare one, and a
#       near-empty rare-class cell (large deviation relative to a small e_c) is
#       penalised hardest -- the failure balancing most needs to avoid.
.balance_score <- function(train_test_table, test_cols, k){
    Tt <- as.matrix(train_test_table[, test_cols, drop = FALSE])
    n_c <- colSums(Tt)
    fillable <- n_c >= k
    empty_cells <- if(any(fillable)) sum(Tt[, fillable, drop = FALSE] == 0) else 0
    keep <- n_c > 0
    if(any(keep)){
        e_c <- n_c[keep] / k
        # column-recycle rather than %*% diag(): diag() of a length-1 vector is
        # an identity matrix, which would break the single-column (PB) case.
        dev <- sweep(Tt[, keep, drop = FALSE], 2, e_c)
        imbalance <- sum(sweep(dev^2, 2, e_c, "/"))
    } else {
        imbalance <- 0
    }
    c(empty_cells, imbalance)
}


# lexicographic 'is score a strictly better than best b' on (empty cells, imbalance)
.balance_better <- function(a, b){
    a[1] < b[1] || (a[1] == b[1] && a[2] < b[2])
}


# assign spatial pieces (blocks or clusters) to k folds using a random search
# that keeps the most balanced split, scored by .balance_score (fewest eligible
# empty test cells, then lowest per-class test imbalance). Shared by cv_spatial
# (random selection) and cv_cluster (balance = TRUE). blocks_df must have a
# 'records' column (the index of each point) and a 'block_id' column (the piece
# each point belongs to). 'opt_cols' names the test columns to optimise (used
# for presence-background balancing); NULL uses every test column.
.balance_folds <- function(blocks_df, blocks_len, k, iteration, response,
                           seed = NULL, biomod2 = TRUE, progress = FALSE,
                           opt_cols = NULL){

    # keep only the columns needed for the assignment
    blocks_df <- blocks_df[, c("records", "block_id")]
    train_test_table <- .records_table(k, response)
    biomod_table <- data.frame(RUN1 = rep(TRUE, nrow(blocks_df)))

    # test columns the objective is scored on, and the best score so far
    test_cols <- if(is.null(opt_cols)){
        grep("^test", names(train_test_table), value = TRUE)
    } else {
        opt_cols
    }
    best <- c(Inf, Inf)

    if(progress){
        pb <- utils::txtProgressBar(min = 0, max = iteration, style = 3)
        on.exit(close(pb), add = TRUE)
    }

    if(!is.null(seed)){
        set.seed(seed)
    }

    for(i in seq_len(iteration)){
        # create random folds with equal proportion of pieces
        fold_df <- data.frame(block_id = seq_len(blocks_len), folds = 0)
        num <- floor(blocks_len / k)
        fold_df$folds[seq_len(num * k)] <- sample(rep(seq_len(k), num), num * k)
        if(blocks_len %% k != 0){
            rest <- blocks_len %% k
            unfold <- which(fold_df$folds == 0)
            fold_df$folds[unfold] <- sample(seq_len(k), rest, replace = FALSE)
        }

        # map the piece folds back to the points
        blocks_df_i <- merge(x = blocks_df, y = fold_df, by = "block_id", all.x = TRUE)

        # count the number of points in each fold
        train_test_table[] <- 0
        fold_list <- list()
        fold_vect <- rep(NA, nrow(blocks_df_i))
        for(p in seq_len(k)){
            train_set <- blocks_df_i$records[which(blocks_df_i$folds != p)]
            test_set <- blocks_df_i$records[which(blocks_df_i$folds == p)]
            fold_vect[test_set] <- p
            fold_list[[p]] <- assign(paste0("fold", p), list(train_set, test_set))
            train_test_table <- .records_table_row(train_test_table, p, train_set, test_set, response)
            if(biomod2){ # creating a biomod2 CV.user.table for validation
                colm <- paste0("RUN", p)
                biomod_table[, colm] <- FALSE
                biomod_table[train_set, colm] <- TRUE
            }
        }

        # score this assignment and keep it if it beats the best so far
        score <- .balance_score(train_test_table, test_cols, k)
        if(.balance_better(score, best)){
            best <- score
            train_test_table2 <- train_test_table
            blocks_df2 <- blocks_df_i
            fold_list2 <- fold_list
            fold_vect2 <- fold_vect
            biomod_table2 <- biomod_table
        }
        if(progress){
            utils::setTxtProgressBar(pb, i)
        }
    }

    list(
        folds_list = fold_list2,
        folds_ids = fold_vect2,
        biomod_table = biomod_table2,
        records = train_test_table2,
        blocks_df = blocks_df2
    )
}

