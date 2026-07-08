#' Leave-group-out cross-validation using an existing grouping factor
#'
#' This function creates cross-validation folds from a grouping column that is already
#' present in the data -- for example a \emph{site}, \emph{plot}, \emph{campaign}, or
#' \emph{individual} identifier. All records that share a group are always kept together,
#' so a group is never split across the training and testing sets. This is the
#' hierarchical / grouped ("leave-group-out") blocking that is not covered by the
#' clustering functions: unlike \code{\link{cv_cluster}} (which \emph{builds} groups with
#' k-means) or \code{\link{cv_knndm}} (whose \code{"hierarchical"} option is agglomerative
#' clustering), \code{cv_group} takes the groups as given.
#'
#' The number of folds is controlled by \code{k}:
#' \itemize{
#'  \item{\code{k = NULL} (default) or \code{k} greater than or equal to the number of groups
#'   -- \strong{leave-group-out}: every group forms its own fold and is left out once. In this
#'   mode \code{balance} has no effect because the folds are fixed by the groups.}
#'  \item{\code{k} less than the number of groups -- the groups are \strong{merged into \code{k}
#'   folds}. With \code{balance = FALSE} (default) the groups are distributed across the folds
#'   deterministically (keeping a similar number of groups per fold). With \code{balance = TRUE}
#'   the groups are assigned to the \code{k} folds over \code{iteration} random attempts, keeping
#'   the split that best balances the training/testing records (or the classes/bins of
#'   \code{column} when it is provided), mirroring the balancing used by \code{\link{cv_cluster}}.}
#' }
#'
#' In every case whole groups move together, so \code{k < n_groups} still yields folds that
#' respect the grouping structure while giving fewer, larger folds than leave-group-out.
#'
#' @inheritParams cv_spatial
#' @param group_col character. The name of the column in \code{x} that holds the grouping factor
#' (e.g. site, plot, campaign, or individual ID). Records sharing a value are kept together in the
#' same fold. Missing values are not allowed.
#' @param k integer (optional). The number of desired folds. When \code{NULL} (default) or greater
#' than or equal to the number of groups, leave-group-out is used (one fold per group). When smaller
#' than the number of groups, the groups are merged into \code{k} folds. See \sQuote{Details}.
#' @param column character (optional). Indicating the name of the column in which response variable
#' (e.g. species data as a binary response i.e. 0s and 1s) is stored. It is used to report whether all
#' the folds contain all the classes and, when \code{balance = TRUE}, to balance those classes across
#' the folds. Continuous numeric responses are binned into quantiles using \code{num_bins} before
#' records are counted.
#' @param balance logical. Only used when \code{k} is smaller than the number of groups. If \code{TRUE},
#' the groups are assigned to the \code{k} folds over \code{iteration} random attempts to balance the
#' training/testing records (or the classes/bins of \code{column} when it is provided). If \code{FALSE}
#' (default), the groups are distributed across the folds deterministically. Ignored for leave-group-out.
#' @param iteration integer value. The number of random attempts to assign the groups to folds when
#' \code{balance = TRUE}.
#' @param seed integer; a random seed for reproducibility of the balancing search.
#'
#' @seealso \code{\link{cv_cluster}} and \code{\link{cv_spatial}};
#' \code{\link{cv_plot}} to visualise, and \code{\link{cv_distance}} and \code{\link{cv_similarity}} to evaluate, the folds
#'
#' @return An object of class S3. A list of objects including:
#'    \itemize{
#'     \item{folds_list - a list containing the folds. Each fold has two vectors with the training (first) and testing (second) indices}
#'     \item{folds_ids - a vector of values indicating the number of the fold for each observation (each number corresponds to the same point in x)}
#'     \item{biomod_table - a matrix with the folds to be used in \pkg{biomod2} package}
#'     \item{k - number of the folds}
#'     \item{column - the name of the column if provided}
#'     \item{group_col - the name of the grouping column}
#'     \item{type - indicates whether leave-group-out or merged grouping was used}
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
#' # add an example grouping column (e.g. survey site)
#' pa_data$site <- sample(paste0("site_", 1:8), nrow(pa_data), replace = TRUE)
#'
#' # leave-group-out: one fold per site
#' lgo <- cv_group(x = pa_data,
#'                 group_col = "site",
#'                 column = "occ") # optional; name of the column with response
#'
#' # merge the sites into 4 balanced folds
#' set.seed(6)
#' gm <- cv_group(x = pa_data,
#'                group_col = "site",
#'                column = "occ",
#'                k = 4,
#'                balance = TRUE)
#'
#' }
cv_group <- function(
        x,
        group_col,
        k = NULL,
        column = NULL,
        balance = FALSE,
        iteration = 100L,
        seed = NULL,
        biomod2 = TRUE,
        num_bins = 4L,
        report = interactive(),
        progress = interactive()
){

    # check x is an sf object
    x <- .check_x(x)
    # group_col is required and must exist in x
    group_col <- .check_group_col(group_col, x)
    # is column in x?
    column <- .check_column(column, x)
    # iteration must be a natural number
    iteration <- suppressWarnings(as.integer(iteration))
    if(length(iteration) != 1L || is.na(iteration)){
        stop("'iteration' must be a single integer value.")
    }
    iteration <- max(1L, abs(iteration))

    # map the grouping factor to integer codes 1..n_groups
    groups <- x[, group_col, drop = TRUE]
    if(anyNA(groups)){
        stop("'group_col' contains missing values; every record must belong to a group.")
    }
    grp <- factor(groups)
    group_id <- as.integer(grp)
    n_groups <- nlevels(grp)
    if(n_groups < 2L){
        stop("'group_col' must contain at least two distinct groups.")
    }

    # resolve k and the cross-validation mode
    logo <- is.null(k)
    if(!logo){
        k <- suppressWarnings(as.integer(k))
        if(length(k) != 1L || is.na(k)){
            stop("'k' must be a single integer value or NULL.")
        }
        k <- max(1L, abs(k))
        if(k >= n_groups){
            if(k > n_groups){
                message("'k' (", k, ") is >= the number of groups (", n_groups,
                        "); using leave-group-out with k = ", n_groups, ".\n")
            }
            logo <- TRUE
        } else if(k < 2L){
            stop("'k' must be 2 or higher.")
        }
    }
    if(logo){
        k <- n_groups
        if(isTRUE(balance)){
            message("'balance' is ignored for leave-group-out (folds are fixed by the groups).\n")
            balance <- FALSE
        }
    }

    # create train-test table
    response <- .column_response(x, column, num_bins = num_bins)

    # progress bar only makes sense for the balancing search
    if(!balance || iteration < 3) progress <- FALSE

    if(balance){
        # assign the groups to k folds, balancing the records
        blocks_df <- data.frame(records = seq_len(nrow(x)), block_id = group_id)
        res <- .balance_folds(
            blocks_df = blocks_df,
            blocks_len = n_groups,
            k = k,
            iteration = iteration,
            response = response,
            seed = seed,
            biomod2 = biomod2,
            progress = progress,
            opt_cols = NULL
        )
        fold_list <- res$folds_list
        fold_ids <- res$folds_ids
        biomod_table <- res$biomod_table
        train_test_table <- res$records

    } else{
        # deterministic fold id per group (one per fold for leave-group-out,
        # otherwise a round-robin that keeps a similar number of groups per fold)
        group_fold <- if(logo) seq_len(n_groups) else rep(seq_len(k), length.out = n_groups)
        point_fold <- group_fold[group_id]

        train_test_table <- .records_table(k, response)
        fold_list <- list()
        fold_ids <- rep(NA, nrow(x))
        biomod_table <- data.frame(RUN1 = rep(TRUE, nrow(x)))

        for(i in seq_len(k)){
            test_set <- which(point_fold == i)
            train_set <- which(point_fold != i)
            fold_ids[test_set] <- i
            fold_list[[i]] <- assign(paste0("fold", i), list(train_set, test_set))
            train_test_table <- .records_table_row(train_test_table, i, train_set, test_set, response)
            if(biomod2){ # creating a biomod2 CV.user.table for validation
                colm <- paste0("RUN", i)
                biomod_table[,colm] <- FALSE
                biomod_table[train_set, colm] <- TRUE
            }
        }
    }

    # give a warning is any folds is empty
    zerofolds <- which(apply(train_test_table, 1, function(x) any(x < 1)))
    if(length(zerofolds) > 0){
        zero_text <- "class(es)"
        if(!is.null(response$bins)) zero_text <- "class/bin(es)"
        if(length(zerofolds) > 1){
            warning("Folds ", paste(zerofolds, collapse = ", "), " have ", zero_text, " with zero records")
        } else{
            warning("Fold ", zerofolds, " has ", zero_text, " with zero records")
        }
    }

    final_objs <- list(
        folds_list = fold_list,
        folds_ids = fold_ids,
        biomod_table = switch(biomod2, as.matrix(biomod_table)),
        k = k,
        column = column,
        presence_bg = FALSE,
        group_col = group_col,
        n_groups = n_groups,
        type = ifelse(logo, "Leave-group-out", "Grouped (merged to k folds)"),
        records = train_test_table
    )

    if(report){
        .print_column_bins(response)
        print(train_test_table)
    }
    # specify the output class
    class(final_objs) <- c("cv_group")

    return(final_objs)
}


#' @export
#' @method print cv_group
print.cv_group <- function(x, ...){
    print(class(x))
}

#' @export
#' @method plot cv_group
plot.cv_group <- function(x, y, data = NULL, ...){
    .plot_cv_fold_map(
        cv = x,
        y = if(!missing(y)) y else NULL,
        data = data,
        has_y = !missing(y),
        ...
    )
}

#' @export
#' @method summary cv_group
summary.cv_group <- function(object, ...){
    cat("Number of recoreds in each training and testing fold:\n")
    print(object$records)
}
