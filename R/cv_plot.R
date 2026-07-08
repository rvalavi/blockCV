#' Visualising folds created by blockCV in ggplot
#'
#' This function visualises the folds create by blockCV. It also accepts a raster
#' layer to be used as background in the output plot.
#' The output of \code{cv_plot()} is a \code{ggplot} object, so it can be
#' customised with standard \code{ggplot2} layers. For example, you can add a
#' title, change the theme, or modify colours and scales.
#'
#' @details
#' The \code{bg_alpha} argument only applies to objects created with \code{presence_bg = TRUE}. There the
#' \code{cv} response holds \code{1} for presences and \code{0} for \emph{background} points -- locations
#' sampled across the study area to represent the available conditions rather than confirmed absences -- and
#' the background points are drawn faded so the presences stand out. This point-level \dQuote{background} is
#' unrelated to the raster \code{r} used as a map backdrop.
#'
#' @param cv a blockCV cv_* object; a \code{cv_spatial}, \code{cv_cluster}, \code{cv_group},
#' \code{cv_buffer}, \code{cv_nndm}, or \code{cv_knndm}
#' @param x a simple features (sf) or SpatialPoints object of the spatial sample data used for creating
#' the \code{cv} object. This is required for point-based objects such as \code{cv_cluster},
#' \code{cv_group}, \code{cv_buffer}, \code{cv_nndm}, and \code{cv_knndm}; it can be omitted for
#' \code{cv_spatial} objects.
#' @param r a terra SpatRaster object (optional). If provided, it will be used as background of the plots.
#' It also supports \emph{stars}, \emph{raster}, or path to a raster file on disk.
#' @param nrow integer; number of rows for facet plot
#' @param ncol integer; number of columns for facet plot
#' @param num_plots a vector of indices of folds; by default the first 10 are shown (if available).
#' You can choose any of the folds to be shown e.g. \code{1:3} or \code{c(2, 7, 16, 22)}
#' @param max_pixels integer; maximum number of pixels used for plotting \code{r}
#' @param raster_colors character; a character vector of colours for raster background e.g. \code{terrain.colors(20)}
#' @param points_colors character; two colours to be used for train and test points
#' @param points_alpha numeric; the opacity of points
#' @param bg_alpha numeric; opacity of the \emph{background} points (response \code{0}) when the \code{cv}
#' object was built with \code{presence_bg = TRUE} (see \sQuote{Details}). Lower values fade them so the
#' presences (response \code{1}) stand out; set \code{bg_alpha = points_alpha} to disable the fading. Has no
#' effect for objects that are not presence-background. The default is \code{0.1}.
#' @param label_size integer; size of fold labels when a \code{cv_spatial} object is used.
#' @param remove_na logical; whether to remove excluded points in \code{cv_buffer} from the plot
#' @param combine_folds logical; if \code{TRUE}, all folds are shown in a single map with points
#' coloured by their fold ID instead of separate train/test facets. Only available for
#' \code{cv_spatial}, \code{cv_cluster}, \code{cv_group} and \code{cv_knndm} objects.
#' @param fold_colors character; a vector of colours for the folds when \code{combine_folds = TRUE};
#' by default a qualitative palette is generated for the number of folds.
#'
#' @seealso \code{\link{cv_distance}} and \code{\link{cv_similarity}} to evaluate the folds;
#' \code{\link{cv_spatial}}, \code{\link{cv_cluster}}, \code{\link{cv_buffer}}, \code{\link{cv_nndm}} and
#' \code{\link{cv_knndm}} to create them
#' @importFrom grDevices gray.colors
#' @return a ggplot object
#' @export
#'
#' @examples
#' \donttest{
#' library(blockCV)
#'
#' # import presence-absence species data
#' points <- read.csv(system.file("extdata/", "species.csv", package = "blockCV"))
#' pa_data <- sf::st_as_sf(points, coords = c("x", "y"), crs = 7845)
#'
#' # spatial clustering
#' sc <- cv_cluster(x = pa_data, k = 5)
#'
#' # now plot the create folds
#' cv_plot(cv = sc,
#'         x = pa_data, # sample points
#'         nrow = 2,
#'         points_alpha = 0.5)
#'
#' }
cv_plot <- function(
        cv,
        x,
        r = NULL,
        nrow = NULL,
        ncol = NULL,
        num_plots = 1:10,
        max_pixels = 3e5,
        remove_na = TRUE,
        raster_colors = gray.colors(10, alpha = 1),
        points_colors = c("#E69F00", "#56B4E9"),
        points_alpha = 0.7,
        bg_alpha = 0.1,
        label_size = 4,
        combine_folds = FALSE,
        fold_colors = NULL
){
    # check for availability of ggplot2
    .check_pkgs(c("ggplot2"))
    # check it's a cv obj
    .check_cv(cv)

    # combined fold map is only for objects with a unique fold per point, and needs points
    if(combine_folds){
        if(.is_loo(cv)) stop("'combine_folds = TRUE' is not supported for cv_buffer or cv_nndm (leave-one-out).")
        if(missing(x)) stop(.cv_plot_missing_x_message(cv), call. = FALSE)
    }

    # check x is an sf object
    if(!missing(x)){
        x <- .check_x(x)
    }

    # change the r to terra object
    if(!is.null(r)){
        r <- .check_r(r)
        r <- r[[1]]
    }
    # is it a cv_spatial object?
    is_spatial <- inherits(cv, "cv_spatial")
    # does the object carry spatial blocks to draw? (cv_spatial always; cv_knndm when kept)
    has_blocks <- !is.null(cv$blocks)
    # presence-background objects: fade the background (0s) so the presences stand out
    pbg_column <- if(isTRUE(cv$presence_bg)) cv$column else NULL
    pbg_plot <- !is.null(pbg_column) && !missing(x) && pbg_column %in% names(x)

    # make geom_tile for raster plots
    if(!is.null(r)){
        map_df <- terra::as.data.frame(
            x = terra::spatSample(
                x = r,
                size = max_pixels,
                method = "regular",
                as.raster = TRUE
            ),
            xy = TRUE,
            na.rm = TRUE
        )
        colnames(map_df) <- c("x", "y", "value")

        geom_rast <- ggplot2::geom_tile(
            data = map_df,
            ggplot2::aes(x = get("x"), y = get("y"), fill = get("value")))

        geom_rast_col <- ggplot2::scale_fill_gradientn(colours = raster_colors)
    }
    # make geom_sf for spatial blocks
    if(has_blocks){
        blocks <- cv$blocks
        geom_poly <- ggplot2::geom_sf(
            data = sf::st_geometry(blocks),
            inherit.aes = FALSE,
            colour = "red",
            fill = "orangered4",
            alpha = 0.04,
            linewidth = 0.2
        )
    }

    if(!missing(x)){
        if(combine_folds){
            # single map: colour each point by its fold ID
            if(length(cv$folds_ids) != nrow(x)){
                stop("Number of rows in 'x' does not match the folds in 'cv'!")
            }
            x$folds <- as.factor(cv$folds_ids)
            if(is.null(fold_colors)){
                fold_colors <- grDevices::hcl.colors(length(cv$folds_list), "Dark 3")
            } else if(!is.character(fold_colors)){
                stop("'fold_colors' must be a character vector.")
            } else if(length(fold_colors) < length(levels(x$folds))){
                stop(sprintf(
                    "'fold_colors' must provide at least %s colours for the folds in 'cv'.",
                    length(levels(x$folds))
                ))
            }
        } else{
            x_long <- .x_to_long(x, cv, num_plots = num_plots)
            # exclude NAs from cv_buffer
            if(.is_loo(cv) && remove_na){
                x_long <- x_long[which(stats::complete.cases(x_long$value)), ]
            }
        }
    } else{
        # point-based CV objects store fold indices, so the original sample data are needed
        if(!is_spatial) stop(.cv_plot_missing_x_message(cv), call. = FALSE)
    }

    geom_sftext <- if (label_size > 0) {
        ggplot2::geom_sf_text(
            ggplot2::aes(label = get("folds")),
            size = label_size, fun.geometry = sf::st_centroid
        )
    } else NULL

    if(missing(x)){
        if(is_spatial){
            p1 <- ggplot2::ggplot(data = blocks) +
                switch(!is.null(r), geom_rast) + # only switch works with ggplot
                switch(!is.null(r), geom_rast_col) +
                ggplot2::geom_sf(colour = "red",
                                 fill = "orangered4",
                                 alpha = 0.04,
                                 linewidth = 0.2) +
                switch(!is.null(geom_sftext), geom_sftext) +
                ggplot2::labs(x = "", y = "") + # or set the axes labes to NULL
                ggplot2::scale_x_continuous(guide = ggplot2::guide_axis(check.overlap = TRUE)) +
                ggplot2::theme_minimal() +
                ggplot2::guides(fill = "none")
        }

    } else if(combine_folds){
        bg_pts <- if(pbg_plot){ v <- x[[pbg_column]]; !is.na(v) & v == 0 } else NULL
        p1 <- ggplot2::ggplot(data = x) +
            switch(!is.null(r), geom_rast) + # only switch works with ggplot
            switch(!is.null(r), geom_rast_col) +
            switch(has_blocks, geom_poly) +
            .cv_point_layers(x, ggplot2::aes(col = get("folds")), points_alpha, bg_alpha, bg_pts) +
            ggplot2::scale_color_manual(values = fold_colors) +
            ggplot2::labs(x = "", y = "", col = "Folds",
                          caption = if(pbg_plot) "Presence-background: background points (0) shown faded" else NULL) + # set the axes labes to NULL
            ggplot2::theme_bw() +
            ggplot2::guides(fill = "none")

    } else{
        bg_long <- if(pbg_plot){ v <- x_long[[pbg_column]]; !is.na(v) & v == 0 } else NULL
        p1 <- ggplot2::ggplot(data = x_long) +
            switch(!is.null(r), geom_rast) + # only switch works with ggplot
            switch(!is.null(r), geom_rast_col) +
            switch(has_blocks, geom_poly) +
            .cv_point_layers(x_long, ggplot2::aes(col = get("value")), points_alpha, bg_alpha, bg_long) +
            ggplot2::scale_color_manual(values = points_colors, na.value = "#BEBEBE03") +
            ggplot2::facet_wrap(~get("folds"), nrow = nrow, ncol = ncol) +
            ggplot2::labs(x = "", y = "", col = "",
                          caption = if(pbg_plot) "Presence-background: background points (0) shown faded" else NULL) + # set the axes labes to NULL
            ggplot2::theme_bw() +
            ggplot2::guides(fill = "none")
    }

    return(p1)
}


# point layer(s) for cv_plot. When 'bg' is supplied, the background is drawn first
# at 'bg_alpha' and the presences on top at 'points_alpha', so presence-background 
# maps highlight the presences. Both layers keep the same 'mapping' (fold or train/test colour).
.cv_point_layers <- function(data, mapping, points_alpha, bg_alpha, bg = NULL){
    if(is.null(bg) || !any(bg)){
        return(ggplot2::geom_sf(mapping, alpha = points_alpha))
    }
    list(
        ggplot2::geom_sf(data = data[bg, ], mapping = mapping, alpha = bg_alpha),
        ggplot2::geom_sf(data = data[!bg, ], mapping = mapping, alpha = points_alpha)
    )
}


.cv_plot_missing_x_message <- function(cv, caller = "cv_plot"){
    cls <- class(cv)[1]
    how <- if(identical(caller, "plot")){
        paste(
            "Supply it as the second argument to plot(), or use the 'data' argument;",
            "e.g. plot(cv, samples) or plot(cv, data = samples).",
            "You can also call cv_plot(cv = cv, x = samples)."
        )
    } else{
        "Supply it with the 'x' argument; e.g. cv_plot(cv = cv, x = samples)."
    }

    paste(
        "The original sample data are required to plot a", cls, "object.",
        "Use the same sf or Spatial object that was passed to the CV function.",
        how
    )
}


.plot_cv_fold_map <- function(cv, y = NULL, data = NULL, has_y = FALSE, ...){
    if(has_y && !is.null(data)){
        stop("Use only one of 'y' or 'data' when providing sample data.", call. = FALSE)
    }

    sample_data <- if(has_y) y else data
    if(is.null(sample_data)){
        stop(.cv_plot_missing_x_message(cv, caller = "plot"), call. = FALSE)
    }

    p1 <- cv_plot(cv = cv, x = sample_data, ...)
    plot(p1)
    invisible(p1)
}


# transform x and fold numbers for plotting
.x_to_long <- function(x, cv, num_plots=1:10){
    # get the folds list
    folds_list <- cv$folds_list
    # The iteration must be a natural number
    tryCatch(
        {
            num_plots <- abs(as.integer(num_plots))
            num_plots <- sort(num_plots)
        },
        error = function(cond) {
            message("'num_plots' must be natural numbers.")
        }
    )
    # length of the folds
    k <- length(folds_list)
    if(max(num_plots) > k){
        num_plots <- num_plots[num_plots <= k]
    }
    # number of original sample points (folds index into these)
    len <- .cv_n_points(cv)
    if(len != nrow(x)){
        stop("Number of rows in 'x' does not match the folds in 'cv'!")
    }
    # create a dataframe temp
    df <- data.frame(id = seq_len(len))
    # make the indices in x
    for (i in num_plots){
        df[, paste("Fold", i, sep = "")] <- NA
        test <- folds_list[[i]][[2]]
        train <- folds_list[[i]][[1]]
        df[test, paste("Fold", i, sep = "")] <- 0
        df[train, paste("Fold", i, sep = "")] <- 1
    }
    # get the geometry column name
    sf_colname <- attr(x, "sf_column")
    # cbind x and the df with fold ids
    xf <- cbind(x, df)
    # convert to dataframe for reshaping
    x_df <- as.data.frame(xf)
    # name of columns to rehspae long
    fold_names <- paste("Fold", num_plots, sep = "")
    # reshape x-df to long
    x_reshape <- stats::reshape(
        x_df,
        direction = "long",
        idvar = "id",
        varying = fold_names,
        times = fold_names,
        v.names = "value",
        timevar = "folds"
    )
    # convert back to sf
    x_long <- sf::st_as_sf(x_reshape, sf_column_name = sf_colname)
    # convert to factor for plotting
    x_long$value <- as.factor(x_long$value)
    levels(x_long$value) <- c("Test", "Train")

    return(x_long)
}

