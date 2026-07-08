#' Compute similarity measures to evaluate possible extrapolation in testing folds
#'
#' This function evaluates environmental similarity between training and testing folds,
#' helping to detect potential extrapolation in the testing data. It supports three
#' similarity outputs: MESS, L1, and L2. The L1 and L2 options are distance-based
#' similarity scores.
#'
#' The MESS is calculated as described in Elith et al. (2010). MESS represents
#' how similar a point in a testing fold is to a training fold (as a reference
#' set of points), with respect to a set of predictor variables in \code{r}.
#' The negative values are the sites where at least one variable has a value that is outside
#' the range of environments over the reference set, so these are novel environments.
#'
#' When using the L1 (Manhattan) or L2 (Euclidean) score options (experimental), the
#' function performs the following steps for each test sample:
#'
#' \itemize{
#' \item{1. Calculates the minimum distance between each test sample and all training samples
#'    in the same fold using the selected metric (L1 or L2).}
#' \item{2. Calculates a baseline distance: the average of the minimum distances between a set
#'    of random background samples (defined by \code{num_sample}) from the raster and all training/test
#'    samples combined.}
#' \item{3. Computes a similarity score by subtracting the test sample’s minimum distance from
#'    the baseline average. A higher score indicates the test sample is more similar to
#'    the training data, while lower or negative scores indicate novelty.}
#' }
#'
#' This provides a simple, distance-based similarity score, not a raw distance: values
#' below zero indicate test samples that are farther from their training data than the
#' random-background baseline. Note that this approach is experimental.
#'
#' When the supplied \code{cv} object was built with \code{presence_bg = TRUE}, the similarity
#' is computed on the \emph{presences} only: the background points (locations sampled to represent
#' the available conditions rather than confirmed absences) are excluded from both the training and
#' testing sets of every fold, matching the framing used in \code{\link{cv_distance}}. This is read
#' automatically from \code{cv}.
#'
#' @inheritParams cv_plot
#' @param x a simple features (sf) object of the spatial sample points used for creating
#' the \code{cv} object.
#' @param r a terra SpatRaster object of environmental predictor that are going to be used for modelling. This
#' is used to calculate similarity between the training and testing points.
#' @param method the similarity method: MESS, L1, or L2. Read the details section.
#' @param type character; \code{"distribution"} (default) draws the per-fold similarity distributions,
#' while \code{"map"} plots the sample points in geographical space and colours each test point by its
#' similarity value, showing \emph{where} extrapolation occurs.
#' @param num_sample number of random raster samples used for the L1/L2 baseline.
#' @param seed integer; an optional random seed. The L1/L2 baseline is built from a random raster
#' sample (\code{num_sample}); set \code{seed} to make the result reproducible.
#' @param num_plots a vector of indices of folds for plotting (default uses all).
#' @param jitter_width numeric; the width of jitter points.
#' @param points_size numeric; the size of points.
#' @param points_alpha numeric; the opacity of points
#' @param points_colors character; a character vector of colours for the diverging value scale.
#' Defaults to a red-grey-blue ramp (red = novel/negative, blue = similar/positive) centred at zero.
#' @inheritParams cv_spatial
#'
#' @seealso \code{\link{cv_spatial}}, \code{\link{cv_cluster}}, \code{\link{cv_buffer}}, and \code{\link{cv_nndm}};
#' \code{\link{cv_plot}} to visualise, and \code{\link{cv_distance}} to evaluate, the folds
#'
#' @references Elith, J., Kearney, M., & Phillips, S. (2010). The art of modelling range-shifting species: The art of modelling range-shifting species. Methods in Ecology and Evolution, 1(4), 330–342.
#'
#' @return a ggplot object. A per-fold extrapolation summary is attached as
#' \code{attr(p, "extrapolation")}: a data.frame with the number of test points, the percentage
#' flagged as novel (similarity value below zero), the minimum and median similarity, and -- for
#' \code{method = "MESS"} -- the most limiting variable (the covariate most often driving the
#' novelty). The overall novelty rate is also shown in the plot subtitle.
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
#' # hexagonal spatial blocking by specified size and random assignment
#' sb <- cv_spatial(x = pa_data,
#'                  column = "occ",
#'                  size = 450000,
#'                  k = 5,
#'                  iteration = 1)
#'
#' # compute extrapolation
#' cv_similarity(cv = sb, r = covars, x = pa_data, method = "MESS")
#'
#' }
cv_similarity <- function(
        cv,
        x,
        r,
        num_plots = seq_along(cv$folds_list),
        method = "MESS",
        type = "distribution",
        num_sample = 10000L,
        seed = NULL,
        jitter_width = 0.1,
        points_size = 2,
        points_alpha = 0.7,
        points_colors = NULL,
        progress = TRUE
){

    # check required packages
    pkg <- c("ggplot2")
    .check_pkgs(pkg)
    # check x is an sf object
    x <- .check_x(x)
    # change the r to terra object
    r <- .check_r(r)
    # check points fall within the raster extent
    .check_within(x, r)

    # check cv obj
    .check_cv(cv)

    method <- match.arg(tolower(method), choices = c("mess", "l1", "l2"))
    MESS <- method == "mess"
    L1 <- method == "l1"
    type <- match.arg(tolower(type), choices = c("distribution", "map"))

    # The iteration must be a natural number
    tryCatch(
        {
            num_plots <- abs(as.integer(num_plots))
            num_plots <- sort(num_plots)
        },
        error = function(cond) {
            message("'num_plots' must be a natural numbers.")
        }
    )

    # get the folds list
    folds_list <- cv$folds_list
    # length of the folds
    k <- length(folds_list)
    if(max(num_plots) > k){
        num_plots <- num_plots[num_plots <= k]
    }
    # presence-background: assess similarity on the presences only (mirror cv_distance),
    # i.e. drop the background points from every fold's train/test sets
    pbg <- isTRUE(cv$presence_bg) && !is.null(cv$column) && cv$column %in% colnames(x)
    presence_idx <- if(pbg) which(x[, cv$column, drop = TRUE] == 1) else NULL
    # extract the raster values
    points <- terra::extract(r, x, ID = FALSE)
    # to set as nrow for df; cv_buffer has only one target points unless P-BG
    n <- nrow(points)
    if(.is_loo(cv)){
        n <- ifelse(cv$presence_bg, nrow(points), 1)
    }
    # number of predictors
    m <- ncol(points)
    df <- data.frame(id = seq_len(n))
    # add progress bar
    if(progress) pb <- utils::txtProgressBar(min = 0, max = length(num_plots), style = 3)

    if (!MESS) {
        # scale a dataset based on params of another scaled dataset
        f <- function(x, s) {
            mns <- attr(s,"scaled:center")
            sds <- attr(s,"scaled:scale")
            sds[sds == 0 | is.na(sds)] <- 1 # guard constant predictors (avoid divide-by-zero)
            for (i in names(mns)) {
                x[[i]] <- (x[[i]] - mns[i]) / sds[i]
            }
            attr(x, "scaled:center") <- NULL
            attr(x, "scaled:scale") <- NULL
            return(x)
        }
        # the baseline is a random raster sample; seed it so results are reproducible
        if(!is.null(seed)) set.seed(seed)
        rand_points <- terra::spatSample(r, size = num_sample, na.rm = TRUE)
        rand_points <- rand_points[stats::complete.cases(rand_points), ] # just make sure..
        # scale and centre the random points
        rand_points <- scale(rand_points)
        # constant predictors become NaN after scaling; set them to the centred value (0)
        zero_var <- which(attr(rand_points, "scaled:scale") == 0 |
                              is.na(attr(rand_points, "scaled:scale")))
        if(length(zero_var)) rand_points[, zero_var] <- 0
        # scale and centre the samples points based on the same stats
        points <- as.matrix(f(x = points, s = rand_points))
        rand_points <- as.matrix(rand_points)
        # remove the attributes
        attr(rand_points, "scaled:center") <- NULL
        attr(rand_points, "scaled:scale") <- NULL
        # samples with missing covariates cannot be used by the distance routine
        complete_idx <- which(stats::complete.cases(points))
    }

    # collect a per-fold extrapolation summary (a test point is "novel" when its value < 0)
    vn <- colnames(points)
    extr <- data.frame()
    total_n <- 0L; total_novel <- 0L
    # per-point values (with their original row index in 'x') for the spatial map
    map_vals <- data.frame()

    # calculate similarity for testing data
    for(i in num_plots){
        df[, paste("Fold", i, sep = "")] <- NA
        train <- folds_list[[i]][[1]]
        test <- folds_list[[i]][[2]]
        # presence-background: keep the presences only (drop background points)
        if(pbg){
            train <- intersect(train, presence_idx)
            test  <- intersect(test,  presence_idx)
        }
        # L1/L2 cannot use samples with missing covariate values
        if(!MESS){
            train <- intersect(train, complete_idx)
            test  <- intersect(test,  complete_idx)
        }
        # nothing left to compare in this fold
        if(!length(train) || !length(test)){
            if(progress) utils::setTxtProgressBar(pb, i)
            next
        }
        if (MESS) {
            mes <- sapply(seq_len(m), function(j) .messi3(points[test, j], points[train, j]))
            if(.is_loo(cv)){
                mmes <- min(mes)
            } else{
                # a single test point collapses 'mes' to a vector; keep it a 1-row matrix
                if(is.null(dim(mes))) mes <- matrix(mes, nrow = length(test))
                mmes <- apply(mes, 1, min, na.rm = TRUE)
            }
        } else {
            mmes <- similarity_cpp(
                train_mat = points[train, , drop = FALSE],
                test_mat  = points[test,  , drop = FALSE],
                rand_mat  = rand_points,
                L1 = L1
            )
        }
        df[seq_along(mmes), paste("Fold", i, sep = "")] <- mmes
        # keep each test point's value with its original index in 'x' (for type = "map")
        map_vals <- rbind(map_vals, data.frame(idx = test, value = mmes))

        # per-fold extrapolation summary
        finite <- is.finite(mmes)
        novel  <- finite & (mmes < 0)
        n_i    <- sum(finite)
        # most limiting (most dissimilar) variable among novel points; MESS only
        mod_i  <- NA_character_
        if(MESS && any(novel)){
            var_idx <- if(is.null(dim(mes))) which.min(mes) else apply(mes[novel, , drop = FALSE], 1, which.min)
            mod_i <- names(sort(table(vn[var_idx]), decreasing = TRUE))[1]
        }
        extr <- rbind(extr, data.frame(
            fold = i,
            n_test = n_i,
            pct_novel = if(n_i) round(100 * sum(novel) / n_i, 1) else NA_real_,
            min = if(n_i) round(min(mmes[finite]), 2) else NA_real_,
            median = if(n_i) round(stats::median(mmes[finite]), 2) else NA_real_,
            limiting_var = mod_i,
            stringsAsFactors = FALSE
        ))
        total_n <- total_n + n_i
        total_novel <- total_novel + sum(novel)

        if(progress) utils::setTxtProgressBar(pb, i)
    }
    # overall extrapolation rate across the plotted folds
    overall_pct <- if(total_n > 0) 100 * total_novel / total_n else NA_real_

    # ---- labels/colours shared by both plot types ----
    col_name <- switch(method, mess = "MESS", l1 = "L1 similarity score", l2 = "L2 similarity score")
    novel_label <- if(MESS) "novel (MESS < 0)" else "beyond the baseline (score < 0)"
    sub_txt <- if(is.na(overall_pct)) NULL else
        sprintf("Extrapolating test points %s: %.1f%% overall", novel_label, overall_pct)
    pbg_caption <- if(pbg) "Presence-background object: similarity computed on presences only" else NULL
    cols <- c("#B2182B", "#D6604D", "#F4A582", "#F7F7F7", "#92C5DE", "#4393C3", "#2166AC")
    points_colors <- if(is.null(points_colors)) cols else points_colors

    # ---- spatial map: colour each test point by its similarity to show *where* extrapolation occurs ----
    if(type == "map"){
        map_sf <- x
        map_sf$.sim <- NA_real_
        map_sf$.sim[map_vals$idx] <- map_vals$value
        map_sf <- map_sf[is.finite(map_sf$.sim), ]
        if(!nrow(map_sf)) stop("No finite similarity values to map.")
        maxabs <- max(abs(map_sf$.sim))
        p1 <- ggplot2::ggplot(map_sf) +
            ggplot2::geom_sf(ggplot2::aes(colour = get(".sim")),
                             size = points_size, alpha = points_alpha) +
            ggplot2::scale_color_gradientn(colours = points_colors,
                                           limits = c(-maxabs, maxabs),
                                           na.value = "#BEBEBE") +
            ggplot2::labs(colour = col_name, subtitle = sub_txt, caption = pbg_caption) +
            ggplot2::theme_bw()
        attr(p1, "extrapolation") <- extr
        return(p1)
    }

    fold_names <- paste("Fold", num_plots, sep = "")
    # reshape for plotting
    mes_reshp <- stats::reshape(
        df,
        direction = "long",
        idvar = "id",
        varying = fold_names,
        times = fold_names,
        v.names =  "value",
        timevar = "folds"
    )

    # remove NAs
    mes_reshp <- mes_reshp[stats::complete.cases(mes_reshp), ]
    if(.is_loo(cv)) mes_reshp$folds <- as.numeric(substr(mes_reshp$folds, 5, 25))
    # get the max value for color legend
    maxabs <- max(abs(mes_reshp$value))
    # define point colors
    cols <- c("#B2182B", "#D6604D", "#F4A582", "#F7F7F7", "#92C5DE", "#4393C3", "#2166AC")
    points_colors <- if(is.null(points_colors)) cols else points_colors

    # provide alternatives for class(cv)
    goem_buffer <- ggplot2::geom_point(size = points_size, alpha = points_alpha)
    geom_other <- ggplot2::geom_jitter(width = jitter_width, size = points_size, alpha = points_alpha)
    geom_vio <- ggplot2::geom_violin(
        ggplot2::aes(x = get("folds"), y = get("value"), group = get("folds")),
        inherit.aes = FALSE,
        fill = NA
    )
    # which geom to choose
    geom_exta <- if(.is_loo(cv)) goem_buffer else geom_other
    mes_reshp$folds <- factor(mes_reshp$folds)

    col_name <- switch(
        method,
        mess = "MESS",
        l1 = "L1 similarity score",
        l2 = "L2 similarity score"
    )
    y_name <- switch(
        method,
        mess = "MESS Values",
        l1 = "Similarity score (baseline - distance)",
        l2 = "Similarity score (baseline - distance)"
    )

    # headline extrapolation rate for the subtitle
    novel_label <- if(MESS) "novel (MESS < 0)" else "beyond the baseline (score < 0)"
    sub_txt <- if(is.na(overall_pct)) NULL else
        sprintf("Extrapolating test points %s: %.1f%% overall", novel_label, overall_pct)

    # shade the novel region (value < 0)
    novel_shade <- ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
                                     fill = "#D53E4F", alpha = 0.06)
    # per-fold "% novel" labels (only meaningful for k-fold designs, not leave-one-out)
    geom_pct <- NULL
    if(!.is_loo(cv) && nrow(extr)){
        pct_lab <- data.frame(
            folds = factor(paste0("Fold", extr$fold), levels = levels(mes_reshp$folds)),
            pct = extr$pct_novel
        )
        pct_lab <- pct_lab[!is.na(pct_lab$folds) & !is.na(pct_lab$pct), ]
        # show a decimal only below 1% so small-but-nonzero folds never read as "0%"
        pct_lab$label <- ifelse(pct_lab$pct == 0, "0%",
                         ifelse(pct_lab$pct < 1, sprintf("%.1f%%", pct_lab$pct),
                                                 sprintf("%.0f%%", pct_lab$pct)))
        if(nrow(pct_lab)){
            geom_pct <- ggplot2::geom_text(
                data = pct_lab,
                ggplot2::aes(x = get("folds"), y = Inf, label = get("label")),
                inherit.aes = FALSE, vjust = 1.4, size = 4.2, fontface = "bold",
                # warn: red for folds with any extrapolation, neutral for the rest
                colour = ifelse(pct_lab$pct > 0, "#B2182B", "grey35")
            )
        }
    }

    p1 <- ggplot2::ggplot(
        data = mes_reshp,
        ggplot2::aes(x = get("folds"), y = get("value"), colour = get("value"))) +
        novel_shade +
        ggplot2::geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
        geom_exta +
        switch(!.is_loo(cv), geom_vio) +
        geom_pct +
        ggplot2::scale_color_gradientn(colours = points_colors,
                                       limits = c(-maxabs, maxabs),
                                       na.value = "#BEBEBE03") +
        ggplot2::labs(x = "Folds", y = y_name, col = col_name, subtitle = sub_txt,
                      caption = if(pbg) "Presence-background object: similarity computed on presences only" else NULL) +
        ggplot2::theme_bw()

    # per-fold extrapolation summary table
    attr(p1, "extrapolation") <- extr

    return(p1)
}


# calculating mess for cv_similarity
# function borrowed from dismo:::.messi3
.messi3 <- function(p, v) {
    # seems 2-3 times faster than messi2
    v <- stats::na.omit(v)
    f <- 100*findInterval(p, sort(v)) / length(v)
    minv <- min(v)
    maxv <- max(v)
    res <- 2*f
    f[is.na(f)] <- -99
    i <- f>50 & f<100
    res[i] <- 200-res[i]

    i <- f==0
    res[i] <- 100*(p[i]-minv)/(maxv-minv)
    i <- f==100
    res[i] <- 100*(maxv-p[i])/(maxv-minv)
    res
}
