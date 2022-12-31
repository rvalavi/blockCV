#' Explore the generated folds
#'
#' This function is deprecated! Please use \code{\link{cv_plot}} function for plotting the folds.
#'
foldExplorer <- function(blocks, rasterLayer, speciesData){
  stop(
    "This function is deprecated! Please use `cv_plot` function for plotting the folds."
  )
}


cv_block_size <- function(r,
                          x = NULL,
                          column = NULL,
                          min_size = NULL,
                          max_size = NULL){
  # check for required packages
  pkg <- c("ggplot2", "shiny")
  .pkg_checks(pkg)

  # check x is an sf object
  if(!is.null(x)){
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
  # is column in x?
  if(!is.null(x) && !is.null(column)){
    if(!column %in% colnames(x)){
      stop(sprintf("There is no column named '%s' in 'x'.\n", column))
    }
  }

  # change the r to terra object
  if(!missing(r)){
    if(!methods::is(r, "SpatRaster")){
      tryCatch(
        {
          r <- terra::rast(r)
        },
        error = function(cond) {
          message("'r' is not convertible to a terra SpatRaster object!")
          message("'r' must be a SpatRaster, stars, Raster* object, or (multiple) path to raster files on disk.")
        }
      )
    }
  }

  # define the object for analysis
  x_obj <- if(missing(r)) x else r
  bbox <- sf::st_bbox(x_obj)
  xrange <- (bbox[3] - bbox[1])
  yrange <- (bbox[4] - bbox[2])

  if(is.na(sf::st_crs(x_obj))){
    stop("'r' and/or 'x' must have defined coordinate refernce systems.")
  }

  size_short <- min(xrange, yrange)
  size_long <- min(xrange, yrange)
  min_size <- ifelse(is.null(min_size), size_short / 20, min_size)
  max_size <- ifelse(is.null(max_size), size_short, max_size)

  if(sf::st_is_longlat(x_obj)){
    min_size <- min_size * 111325
    max_size <- max_size * 111325
    size_long <- size_long * 111325
    size_short <- size_short * 111325
  }

  if(max_size > size_long){ # limit the maximum range to the maximum extent of the map
    message("Selected range is bigger than the extent of the base map. The range is set to the biggest extent (X or Y).\n")
    max_size <- size_long
  }
  min_size <- round(min_size)
  max_size <- round(max_size)
  val <- round(size_short / 10)

  if(!missing(r)){
    map_df <- terra::spatSample(r[[1]],
                                size = 2e5,
                                method = "regular",
                                xy = TRUE,
                                na.rm = TRUE)
    colnames(map_df) <- c("x", "y", "value")

    base_plot <- ggplot2::ggplot() +
      ggplot2::geom_tile(data = map_df,
                         ggplot2::aes_string(x = "x",
                                             y = "y",
                                             fill = "value")) +
      ggplot2::scale_fill_gradientn(colours = gray.colors(20, alpha = 1)) +
      ggplot2::guides(fill = "none")

  } else{
    base_plot <- ggplot2::ggplot()
  }

  if(!is.null(x)){
    geom_x <- ggplot2::geom_sf(
      data = x,
      ggplot2::aes_string(colour = switch(!is.null(column), column, NULL)),
      inherit.aes = FALSE,
      alpha = 0.5
    )
  }

  # create UI for shniy app
  ui <- shiny::fluidPage(
    shiny::wellPanel(
      shiny::fluidRow(
        shiny::column(6, offset = 3,
                      shiny::sliderInput(inputId = "num",
                                         label = "Choose a range (m) to create spatial blocks",
                                         value = val,
                                         min = min_size,
                                         max = max_size)
        )
      ),
      shiny::HTML("<br/>"),
      shiny::fluidRow(
        shiny::submitButton()

      ),
      shiny::HTML("<br/>"),

      shiny::fluidRow(
        shiny::column(12, shiny::plotOutput(outputId = "ggplot"))
      )
    )
  )

  # create shiny server and main code
  server <- function(input, output){
    output$ggplot <- shiny::renderPlot({

      plot_size <- if(sf::st_is_longlat(x_obj)) round(input$num) / 111325 else round(input$num)
      vis_block <- sf::st_make_grid(x_obj, cellsize = plot_size, what = "polygons")

      p1 <- base_plot +
        ggplot2::geom_sf(data = vis_block,
                         color = "red",
                         fill = "orangered4",
                         alpha = 0.04,
                         size = 0.2) +
        switch(!is.null(x), geom_x, NULL) +
        ggplot2::ggtitle("Spatial blocks",
                         subtitle=paste("Using",
                                        input$num,
                                        "(m) block size")) +
        ggplot2::labs(x = "", y = "", color = column)

      # plot ggplot
      plot(p1)
    })
  }
  # starting the shiny app
  shiny::shinyApp(ui = ui, server = server)
}
