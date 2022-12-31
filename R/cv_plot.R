
cv_plot <- function(
    cv,
    x,
    r = NULL,
    nrow = NULL,
    ncol = NULL,
    num_plots = sample(length(cv$folds_list), min(length(cv$folds_list), 10)),
    max_pixels = 3e5,
    raster_colors = gray.colors(10, alpha = 1),
    points_colors = c("#E69F00", "#56B4E9"),
    points_alpha = 0.7,
    label_size = 4
){
  # check for availability of ggplot2
  pkg <- c("ggplot2")
  .pkg_checks(pkg)

  if(!class(cv) %in% c("cv_spatial", "cv_cluster", "cv_buffer")){
    stop("'cv' must be a blockCV cv_* object.")
  }

  # check x is an sf object
  if(!missing(x)){
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
          message("'r' must be a SpatRaster, stars, Raster* object, or path to a raster file on disk.")
        }
      )
    }
  }

  # make geom_tile for raster plots
  if(!is.null(r)){
    map_df <- terra::spatSample(r,
                                size = max_pixels,
                                method="regular",
                                xy=TRUE,
                                na.rm=TRUE)
    colnames(map_df) <- c("x", "y", "value")

    geom_rast <- ggplot2::geom_tile(data = map_df,
                                    ggplot2::aes_string(x="x", y="y", fill="value"))
    geom_rast_col <- ggplot2::scale_fill_gradientn(colours = raster_colors)
  }
  # make geom_sf for spatial blocks
  if(class(cv) == "cv_spatial"){
    blocks <- cv$blocks
    geom_poly <- ggplot2::geom_sf(data = sf::st_geometry(blocks),
                                  inherit.aes = FALSE,
                                  colour = "red",
                                  fill = "orangered4",
                                  alpha = 0.04,
                                  size = 0.2)
  }

  if(!missing(x)){
    x_long <- .x_to_long(x, cv, num_plot = num_plots)
  } else{
    # stop if x is missing for buffer and cluster
    if(class(cv) != "cv_spatial") stop("'x' is required for plotting cv_cluster and cv_buffer.")
  }

  if(missing(x)){
    if(class(cv) == "cv_spatial"){

      p1 <- ggplot2::ggplot(data = blocks) +
        switch(!is.null(r), geom_rast, NULL) +
        switch(!is.null(r), geom_rast_col, NULL) +
        ggplot2::geom_sf(colour = "red",
                         fill = "orangered4",
                         alpha = 0.04,
                         size = 0.2) +
        ggplot2::geom_sf_text(ggplot2::aes_string(label = "folds"),
                              size = label_size, fun.geometry = sf::st_centroid) +
        ggplot2::labs(x = "", y = "") + # or set the axes labes to NULL
        ggplot2::scale_x_continuous(guide = ggplot2::guide_axis(check.overlap = TRUE)) +
        ggplot2::theme_minimal() +
        ggplot2::guides(fill = "none")

    }

  } else{

    p1 <- ggplot2::ggplot(data = x_long) +
      switch(!is.null(r), geom_rast, NULL) +
      switch(!is.null(r), geom_rast_col, NULL) +
      switch(class(cv) == "cv_spatial", geom_poly, NULL) +
      ggplot2::geom_sf(ggplot2::aes_string(col = "value"),
                       alpha = points_alpha) +
      ggplot2::scale_color_manual(values = points_colors, na.value = "#BEBEBE03") +
      ggplot2::facet_wrap(~folds, nrow = nrow, ncol = ncol) +
      ggplot2::labs(x = "", y = "", col = "") + # set the axes labes to NULL
      ggplot2::theme_bw() +
      ggplot2::guides(fill = "none")

  }

  return(p1)
}
