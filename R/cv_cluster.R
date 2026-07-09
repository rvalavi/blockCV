#' Use environmental or spatial clustering to separate train and test folds
#'
#' This function uses clustering methods to specify sets of similar environmental
#' conditions based on the input covariates, or cluster of spatial coordinates of the sample data.
#' Sample data (i.e. species data) corresponding to any of
#' these groups or clusters are assigned to a fold. Clustering is done
#' using \code{\link[stats]{kmeans}} for both approaches. The only requirement is \code{x} that leads to
#' a clustering of the coordinates of sample data. Otherwise, by providing \code{r}, environmental
#' clustering is done.
#' Environmental clustering can also be spatially restricted with \code{spatial_weight},
#' and fold sizes or classes can be balanced with \code{balance = TRUE}; see Details.
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
#' By default (\code{balance = FALSE}) the points are clustered directly into \code{k} groups, so the folds are compact but their
#' sizes -- and the number of records of each class -- can be very uneven. When \code{balance = TRUE}, the points are instead
#' grouped into a larger number of clusters (\code{k * k_multiplier}) that are then assigned to \code{k} folds over \code{iteration}
#' random attempts, keeping the split that best balances the training/testing records (or the classes/bins of \code{column} when
#' it is provided). This mirrors the balancing used by \code{\link{cv_spatial}} with \code{selection = "random"}, and \code{k_multiplier}
#' controls the trade-off between fold balance (higher values) and fold compactness (lower values).
#'
#' For environmental clustering (when \code{r} is provided), \code{spatial_weight} adds a soft spatial compactness
#' pressure by blending the coordinates into the covariates before clustering, so the k-means Euclidean distance
#' becomes \eqn{(1 - w) d^2_{env} + w \, d^2_{geo}}: \code{spatial_weight = 0} (default) is pure environmental
#' clustering, larger values pull the environmental clusters to be more geographically compact, and
#' \code{spatial_weight = 1} clusters on the coordinates alone. This is a soft pressure -- it favours, but does not
#' guarantee, spatially contiguous folds. The coordinates are always standardised (so a
#' \code{spatial_weight = 1} clustering is \emph{not} identical to the coordinate-only clustering obtained with
#' \code{r = NULL}, which uses the raw coordinates), while the covariates are standardised only when
#' \code{scale = TRUE}; \code{spatial_weight} is best calibrated with the default \code{scale = TRUE}. Coordinates
#' should be projected (a warning is issued for lon/lat data).
#'
#' @details
#' For presence-background data (\code{presence_bg = TRUE}), \code{column} holds \code{1} for presences and
#' \code{0} for \emph{background} points -- locations sampled across the study area to represent the available
#' conditions rather than confirmed absences. When \code{balance = TRUE} the balancing then equalises only the
#' presence records across folds, so the abundant background cannot dominate the split.
#'
#' @inheritParams cv_spatial
#' @param column character (optional). Indicating the name of the column in which response variable (e.g. species data as a binary
#'  response i.e. 0s and 1s) is stored. It is used to report whether all the folds contain all the classes and, when
#'  \code{balance = TRUE}, to balance those classes across the folds.
#'  Continuous numeric responses are binned into quantiles using \code{num_bins} before records are counted.
#' @param r a terra SpatRaster object of covariates to identify environmental groups. If provided, clustering will be done
#' in environmental space rather than spatial coordinates of sample points. Only numeric (quantitative) covariates are
#' supported; categorical (factor) layers are rejected because k-means relies on Euclidean distance.
#' @param scale logical; whether to scale the input rasters (recommended) for clustering.
#' @param raster_cluster logical; if \code{TRUE}, the clustering is done over the entire raster layer,
#' otherwise it will be over the extracted raster values of the sample points. See details for more information.
#' @param num_sample integer; the number of samples from raster layers to build the clusters (when \code{raster_cluster = FALSE}).
#' @param spatial_weight numeric in \code{[0, 1]}; only used (and only validated) for environmental clustering (when \code{r} is
#' provided). It adds a soft spatial compactness pressure by blending the geographic coordinates into the environmental clustering.
#' \code{0} (the default) does pure environmental clustering; larger values (e.g. \code{0.3}--\code{0.5}) give environmentally
#' coherent folds that are also geographically separated; \code{1} clusters on the coordinates alone. See \sQuote{Details}.
#' @param balance logical. If \code{TRUE}, the sample points are first grouped into \code{k * k_multiplier} clusters (see \code{k_multiplier})
#' and these small clusters are then assigned to \code{k} folds over \code{iteration} random attempts to balance the training/testing
#' records (or the classes/bins of \code{column} when it is provided). If \code{FALSE} (default), the points are clustered directly
#' into \code{k} folds, which keeps the folds compact but does not control their size or class balance.
#' @param presence_bg logical; whether to treat \code{column} as species presence-background data (0s for
#' background points and 1s for presences; see \sQuote{Details}). When \code{TRUE} (and \code{balance = TRUE}), the balancing
#' search equalises only the presence (1s) records across folds so the many background points cannot dominate
#' the objective; the background points are still clustered spatially but ignored when scoring the balance.
#' Requires a binary numeric \code{column}. The default is \code{FALSE}.
#' @param k_multiplier integer. The multiplier controlling how many clusters are created before they are merged into folds
#' (i.e. \code{k * k_multiplier}). Only used when \code{balance = TRUE}. Larger values give more balanced folds at the cost of less
#' compact folds; smaller values keep the folds compact but less balanced. The default is \code{5}.
#' @param iteration integer value. The number of random attempts to assign the clusters to folds when \code{balance = TRUE}.
#' @param seed integer; a random seed for reproducibility of the balancing search.
#' @param ... additional arguments for \code{stats::kmeans} function, e.g. \code{algorithm = "MacQueen"}.
#'
#' @seealso \code{\link{cv_buffer}} and \code{\link{cv_spatial}};
#' \code{\link{cv_plot}} to visualise, and \code{\link{cv_distance}} and \code{\link{cv_similarity}} to evaluate, the folds
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
#'     \item{spatial_weight - the spatial weight used for environmental clustering (\code{NA} for spatial clustering).}
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
#' # spatially constrained environmental clustering
#' set.seed(6)
#' sec <- cv_cluster(r = covars,
#'                   x = pa_data,
#'                   column = "occ",
#'                   k = 5,
#'                   scale = TRUE,
#'                   spatial_weight = 0.4) # blend geography into the environmental clusters
#'
#' # spatial clustering with balanced folds
#' set.seed(6)
#' bc <- cv_cluster(x = pa_data,
#'                  column = "occ",
#'                  k = 5,
#'                  balance = TRUE, # balance the records across the folds
#'                  k_multiplier = 3) # cluster into k * k_multiplier groups, then merge into k folds
#'
#' # presence-background data: balance folds based on presences
#' points_pb <- read.csv(system.file("extdata/", "species_pb.csv", package = "blockCV"))
#' pb_data <- sf::st_as_sf(points_pb, coords = c("x", "y"), crs = 7845)
#' # pb_data <- pb_data[c(which(pb_data$occ == 1)[1:50], which(pb_data$occ == 0)[1:200]), ]
#'
#' bc_pb <- cv_cluster(x = pb_data,
#'                     column = "occ",
#'                     k = 5,
#'                     presence_bg = TRUE,
#'                     balance = TRUE,
#'                     k_multiplier = 5)
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
        spatial_weight = 0,
        balance = FALSE,
        presence_bg = FALSE,
        k_multiplier = 5L,
        iteration = 100L,
        seed = NULL,
        biomod2 = TRUE,
        num_bins = 4L,
        report = interactive(),
        progress = interactive(),
        ...
){

    # check x is an sf object
    x <- .check_x(x)
    # is column in x?
    column <- .check_column(column, x)
    # validate presence-background data (0/1 column) when requested
    invisible(.presence_index(x, column, presence_bg))
    # k_multiplier and iteration must be natural numbers
    k_multiplier <- max(1L, abs(as.integer(k_multiplier)))
    iteration <- max(1L, abs(as.integer(iteration)))
    # spatial_weight only affects environmental clustering; warn if it is set without a raster
    if(is.null(r) && !isTRUE(spatial_weight == 0)){
        warning("'spatial_weight' is only used for environmental clustering (when 'r' is provided); ignored.")
    }

    if(!is.null(r)){
        # spatial_weight blends geography into environmental clustering; must be in [0, 1]
        if(length(spatial_weight) != 1 || !is.numeric(spatial_weight) ||
           is.na(spatial_weight) || spatial_weight < 0 || spatial_weight > 1){
            stop("'spatial_weight' must be a single number between 0 and 1.")
        }
        # check for required packages
        pkg <- c("terra")
        .check_pkgs(pkg)
        # check r
        r <- .check_r(r)
        # check points fall within the raster extent
        .check_within(x, r)
        # check r layers
        if(terra::nlyr(r) < 1){
            stop("'r' is not a valid raster.")
        }
        # k-means needs numeric covariates; reject categorical (factor) layers
        .check_r_numeric(r)
        # scale?
        if (scale){
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

    # over-cluster into many small groups when balancing, otherwise one per fold
    n_centers <- k
    if(balance){
        n_centers <- min(k * k_multiplier, nrow(x))
        if(n_centers < k * k_multiplier){
            message("The number of clusters reduced to ", n_centers, "; the total number of sample points.\n")
        }
    }

    if (is.null(r)){
        # spatial cluster on xy; k-means uses Euclidean distance on the coordinates,
        # which is only approximate for geographic (lon/lat) data
        if(isTRUE(sf::st_is_longlat(x))){
            warning("k-means clustering uses Euclidean distance on coordinates; geographic (lon/lat) coordinates are approximate -- consider projecting 'x'.")
        }
        xy <- sf::st_coordinates(x)
        .check_kmeans_k(n_centers, nrow(xy), what = "sample points")
        kms <- stats::kmeans(xy, centers = n_centers, iter.max = 500, nstart = 25, ...)
        cluster_ids <- as.integer(kms$cluster)

    } else{
        # a general check to make sure x is covered by r
        x_vals <- terra::extract(r, x, ID = FALSE)
        if (anyNA(x_vals)){
            stop("The input raster layer does not cover all the column points.")
        }
        # blending geography (spatial_weight > 0) uses Euclidean distance on coordinates
        if(spatial_weight > 0 && isTRUE(sf::st_is_longlat(x))){
            warning("'spatial_weight' uses Euclidean distance on coordinates; geographic (lon/lat) coordinates are approximate -- consider projecting 'x'.")
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
            if(spatial_weight > 0){
                # sample cells with their coordinates (xy = TRUE prepends x, y) so the
                # sampled cells and the sample points share both covariates and geography
                sampr <- terra::spatSample(r, size = num_sample, method = "random", na.rm = TRUE, xy = TRUE)
                sampr <- sampr[stats::complete.cases(sampr[, -(1:2), drop = FALSE]), , drop = FALSE]
                env <- rbind(x_vals, sampr[, -(1:2), drop = FALSE])
                xy  <- rbind(sf::st_coordinates(x), as.matrix(sampr[, 1:2]))
                feat <- .augment_geo(env, xy, spatial_weight, scale_env = scale)
            } else{
                # original path (unchanged) for exact backward compatibility
                sampr <- terra::spatSample(r, size = num_sample, method = "random", na.rm = TRUE)
                sampr <- sampr[stats::complete.cases(sampr), ]
                feat <- rbind(x_vals, sampr)
            }
            .check_kmeans_k(n_centers, nrow(feat), what = "sample points and raster cells")
            kms <- stats::kmeans(feat, centers = n_centers, iter.max = 500, nstart = 25, ...)
            cluster_ids <- kms$cluster[seq_len(nrow(x))]
        } else{
            if(spatial_weight > 0){
                feat <- .augment_geo(x_vals, sf::st_coordinates(x), spatial_weight, scale_env = scale)
            } else{
                feat <- x_vals # original path (unchanged) for exact backward compatibility
            }
            .check_kmeans_k(n_centers, nrow(feat), what = "sample points")
            kms <- stats::kmeans(feat, centers = n_centers, iter.max = 500, nstart = 25, ...)
            cluster_ids <- kms$cluster
        }

    }

    # create train-test table
    response <- .column_response(x, column, num_bins = num_bins)

    # progress bar only makes sense for the balancing search
    if(!balance || iteration < 3) progress <- FALSE

    if(balance){
        # assign the over-clustered groups to k folds, balancing the records
        blocks_df <- data.frame(records = seq_len(nrow(x)), block_id = cluster_ids)
        res <- .balance_folds(
            blocks_df = blocks_df,
            blocks_len = n_centers,
            k = k,
            iteration = iteration,
            response = response,
            seed = seed,
            biomod2 = biomod2,
            progress = progress,
            opt_cols = .balance_opt_cols(response, presence_bg)
        )
        fold_list <- res$folds_list
        fold_ids <- res$folds_ids
        biomod_table <- res$biomod_table
        train_test_table <- res$records

    } else{
        train_test_table <- .records_table(k, response)
        fold_list <- list()
        fold_ids <- rep(NA, nrow(x))
        biomod_table <- data.frame(RUN1 = rep(TRUE, nrow(x)))

        for(i in seq_len(k)){
            test_set <- which(cluster_ids == i)
            train_set <- which(cluster_ids != i)
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
        biomod_table = if(biomod2) as.matrix(biomod_table) else NULL,
        k = k,
        column = column,
        presence_bg = presence_bg,
        type = ifelse(is.null(r), "Spatial Cluster", "Environmental Cluster"),
        spatial_weight = if(is.null(r)) NA_real_ else spatial_weight,
        records = train_test_table
    )

    if(report){
        .print_column_bins(response)
        print(train_test_table)
    }
    # specify the output class
    class(final_objs) <- c("cv_cluster")

    return(final_objs)
}


#' @export
#' @method print cv_cluster
print.cv_cluster <- function(x, ...){
    w <- x$spatial_weight
    desc <- if(is.null(x$type) || identical(x$type, "Spatial Cluster")){
        "spatial clustering (geographical coordinates)"
    } else if(is.null(w) || is.na(w) || w <= 0){
        "environmental clustering (feature space)"
    } else {
        "spatially-constrained environmental clustering (feature space + geography)"
    }
    details <- list("Folds" = x$k)
    if(!is.null(w) && !is.na(w) && w > 0) details[["Spatial weight"]] <- w
    if(!is.null(x$column)) details[["Balancing column"]] <- x$column
    details[["Presence-background"]] <- if(isTRUE(x$presence_bg)) "yes" else "no"
    .print_cv_folds(x, desc, details)
}

#' @export
#' @method plot cv_cluster
plot.cv_cluster <- function(x, y, data = NULL, ...){
    .plot_cv_fold_map(
        cv = x,
        y = if(!missing(y)) y else NULL,
        data = data,
        has_y = !missing(y),
        ...
    )
}

#' @export
#' @method summary cv_cluster
summary.cv_cluster <- function(object, ...){
    cat("Number of records in each training and testing fold:\n")
    print(object$records)
}
