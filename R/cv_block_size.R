#' Explore spatial block size
#'
#' This function assists selection of block size. It allows the user to visualise the blocks
#' interactively, viewing the impact of block size on number and arrangement of blocks in
#' the landscape (and optionally on the distribution of species data in those blocks).
#' Slide to the selected block size, and click \emph{Apply Changes} to change the block size.
#'
#' @inheritParams cv_spatial
#' @param x a simple features (sf) or SpatialPoints object of spatial sample data. If \code{r} is supplied, this
#' is only added to the plot. Otherwise, the extent of \code{x} is used for creating the blocks.
#' @param column character (optional). Indicating the name of the column in which response variable (e.g.
#' species data as a binary response i.e. 0s and 1s) is stored to be shown on the plot.
#' @param min_size numeric; the minimum size of the blocks (in metres) to explore.
#' @param max_size numeric; the maximum size of the blocks (in metres) to explore.
#'
#' @return an interactive shiny session
#' @export
#'
#' @examples
#' \donttest{
#' if(interactive()){
#' library(blockCV)
#'
#' # import presence-absence species data
#' points <- read.csv(system.file("extdata/", "species.csv", package = "blockCV"))
#' pa_data <- sf::st_as_sf(points, coords = c("x", "y"), crs = 7845)
#'
#' # manually choose the size of spatial blocks
#' cv_block_size(x = pa_data,
#'               column = "occ",
#'               min_size = 2e5,
#'               max_size = 9e5)
#'
#' }
#' }
#'
cv_block_size <- function(r, # priority
                          x = NULL,
                          column = NULL,
                          min_size = NULL,
                          max_size = NULL){
  # check for required packages
  pkg <- c("ggplot2", "shiny")
  .check_pkgs(pkg)

  # check x is an sf object
  if(!is.null(x)){
    x <- .check_x(x)
    column <- .check_column(column, x)
  }
  # change the r to terra object
  if(!missing(r)){
    r <- .check_r(r)
    r <- r[[1]]
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
      ggplot2::geom_tile(
        data = map_df,
        ggplot2::aes(x = get("x"), y = get("y"), fill = get("value"))) +
      ggplot2::scale_fill_gradientn(colours = gray.colors(20, alpha = 1)) +
      ggplot2::guides(fill = "none")

  } else{
    base_plot <- ggplot2::ggplot()
  }

  if(!is.null(x)){
    geom_x <- ggplot2::geom_sf(
      data = x,
      switch(!is.null(column), ggplot2::aes(colour = {{ column }}), NULL),
      # ggplot2::aes(colour = {{ switch(!is.null(column), column, NULL) }}),
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
  server <- function(input, output, session){
    output$ggplot <- shiny::renderPlot({

      # stop app after session ends
      session$onSessionEnded(function() {
        shiny::stopApp()
      })

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
