#' Explore the generated folds
#'
#' A function for visualising the generated folds on a map, and allowing interactive exploration of the data in the folds,
#' using the \pkg{RStudio Shiny} app.
#'
#' @param blocks An SpatialBlock, EnvironmentalBlock or BufferedBlock object.
#' @param rasterLayer A raster object as background map for visualisation.
#' @inheritParams buffering
#'
#' @seealso \code{\link{spatialBlock}}, \code{\link{buffering}} and \code{\link{envBlock}}
#'
#' @return An interactive map showing folds and the species data, that can be used  to explore folds. Note that this can also
#' be opened in a web browser window. When you return to the R console, press "Esc" to return to the prompt.
#' @export
#'
#' @examples
#' \donttest{
#' if(interactive()){
#'
#' # load package data
#' awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))
#' # import presence-absence species data
#' PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
#' # make a sf object from data.frame
#' pa_data <- sf::st_as_sf(PA, coords = c("x", "y"), crs = raster::crs(awt))
#'
#' # spatial blocking by specified range and random assignment
#' sb <- spatialBlock(speciesData = pa_data,
#'                    species = "Species",
#'                    rasterLayer = awt,
#'                    theRange = 70000,
#'                    k = 5,
#'                    selection = "random",
#'                    iteration = 100)
#'
#' foldExplorer(sb, awt, pa_data)
#'
#' # buffering with presence-absence data
#' bf <- buffering(speciesData= pa_data,
#'                 species= "Species", # to count the number of presences and absences
#'                 theRange= 70000,
#'                 spDataType = "PA",
#'                 progress = TRUE)
#'
#' foldExplorer(bf, awt, pa_data)
#'
#' # environmental clustering
#' eb <- envBlock(rasterLayer = awt,
#'                speciesData = pa_data,
#'                species = "Species",
#'                k = 5)
#'
#' foldExplorer(eb, awt, pa_data)
#'
#' }
#' }
#'
foldExplorer <- function(blocks, rasterLayer, speciesData){
  # check for required packages
  pkg <- c("ggplot2", "cowplot", "shiny", "shinydashboard")
  pkgna <- names(which(sapply(sapply(pkg, find.package, quiet = TRUE), length) == 0))
  if(length(pkgna) > 0){
    nm <- paste(pkgna, collapse = ", ")
    message("This function requires these packages: ", nm,
            "\nWould you like to install them now?\n1: yes\n2: no")
    user <- readline(prompt = paste0("Selection: "))
    if(tolower(user) %in% c("1", "yes", "y")){
      utils::install.packages(pkgna)
    } else{
      stop("Please install these packages: ", nm)
    }
  }
  # testing the input arguments
  if(is.null(rasterLayer)){
    stop("A raster layer should be provided")
  } else if(is.null(speciesData)){
    stop("Species data should be provided")
  } else if(is.null(blocks)){
    stop("An object of SpatialBlock, EnvironmentalBlock or BufferedBlock is needed")
  }
  # select the default value
  if(class(blocks) == "SpatialBlock"){
    polyObj <- blocks$blocks
  } else{
    polyObj <- NULL
  }
  folds <- blocks$folds
  kmax <- length(folds)
  species <- blocks$species
  # set x and y coordinates
  if(methods::is(speciesData, "SpatialPoints")){
    speciesData <- sf::st_as_sf(speciesData)
  } else if(!methods::is(speciesData, "sf")){
    stop("speciesData should be a sf or SpatialPoints object")
  }
  # plot raster file in ggplot2
  samp <- raster::sampleRegular(rasterLayer[[1]], 5e+05, asRaster=TRUE)
  map_df <- raster::as.data.frame(samp, xy=TRUE, centroids=TRUE, na.rm=TRUE)
  colnames(map_df) <- c("Easting", "Northing", "MAP")
  mid <- stats::median(map_df$MAP)
  basePlot <- ggplot2::ggplot() +
    ggplot2::geom_raster(data=map_df, ggplot2::aes_string(y="Northing", x="Easting", fill="MAP")) +
    ggplot2::scale_fill_gradient2(low="darkred", mid="yellow", high="darkgreen", midpoint=mid) +
    ggplot2::guides(fill = FALSE) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "", y = "")
  # create UI for shniy app
  ui <- shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(title = "Fold Explorer"),
    shinydashboard::dashboardSidebar(disable = TRUE),
    shinydashboard::dashboardBody(
      shiny::fluidRow(
        shiny::sidebarLayout(
          shiny::sidebarPanel(shiny::fluidRow(
            shiny::wellPanel(shiny::sliderInput(inputId = "num",
                                                label = "Choose a fold to display",
                                                value = 2, min = 1, max = kmax, step = 1L)
            )
          ),
          shiny::fluidRow("Record's Statistics"),
          if(is.null(species)){
            shiny::fluidRow(shiny::textOutput(outputId = "Tr"),
                            shiny::textOutput(outputId = "Ts"))
          } else{
            shiny::fluidRow(shiny::tableOutput(outputId = "Tb"))
          }
          ),
          shiny::mainPanel(shiny::plotOutput(outputId = "ggplot"))
        )
      )
    )
  )
  # create shiny server and main code
  server <- function(input, output){
    output$ggplot <- shiny::renderPlot({
      trainSet <- unlist(folds[[input$num]][1])
      testSet <- unlist(folds[[input$num]][2])
      training <- speciesData[trainSet, ]
      testing <- speciesData[testSet, ]
      if(class(blocks) == "SpatialBlock"){
        plotPoly <- polyObj[polyObj$folds==input$num, ]
        plotPoly <- sf::st_as_sf(plotPoly)
      }
      if(is.null(species)){
        if(class(blocks) == "SpatialBlock"){
          ptr <- basePlot + ggplot2::geom_sf(data = plotPoly, color ="red", fill ="orangered4",
                                             alpha = 0.04, size = 0.2) +
            ggplot2::geom_sf(data = training, alpha = 0.7, color = "blue", size = 2) +
            ggplot2::ggtitle('Training set')
          # ploting test data
          pts <- basePlot + ggplot2::geom_sf(data = plotPoly, color ="red", fill ="orangered4",
                                             alpha = 0.04, size = 0.2) +
            ggplot2::geom_sf(data = testing, alpha = 0.7, color = "blue", size = 2) +
            ggplot2::ggtitle('Testing set')
        } else{
          ptr <- basePlot + ggplot2::geom_sf(data = training,  alpha = 0.7, color = "blue", size = 2) +
            ggplot2::ggtitle('Training set')
          # ploting test data
          pts <- basePlot + ggplot2::geom_sf(data = testing,  alpha = 0.7, color = "blue", size = 2) +
            ggplot2::ggtitle("Testing set")
        }
      } else{
        if(class(blocks) == "SpatialBlock"){
          ptr <- basePlot + ggplot2::geom_sf(data = plotPoly, color ="red", fill ="orangered4",
                                             alpha = 0.04, size = 0.2) +
            ggplot2::geom_sf(data = training, ggplot2::aes(color = get(species)), show.legend = "point",
                             alpha = 0.7, size = 2) +
            ggplot2::labs(color = species) +
            ggplot2::ggtitle('Training set')
          # ploting test data
          pts <- basePlot + ggplot2::geom_sf(data = plotPoly, color ="red", fill ="orangered4",
                                             alpha = 0.04, size = 0.2) +
            ggplot2::geom_sf(data = testing, ggplot2::aes(color = get(species)), show.legend = "point",
                             alpha = 0.7, size = 2) +
            ggplot2::labs(color = species) +
            ggplot2::ggtitle('Testing set')

        } else{
          ptr <- basePlot + ggplot2::geom_sf(data = training, ggplot2::aes(color = get(species)), show.legend = "point",
                                             alpha = 0.7, size = 2) +
            ggplot2::labs(color = species) +
            ggplot2::ggtitle('Training set')
          # ploting test data
          pts <- basePlot + ggplot2::geom_sf(data = testing, ggplot2::aes(color = get(species)), show.legend = "point",
                                             alpha = 0.7, size = 2) +
            ggplot2::labs(color = species) +
            ggplot2::ggtitle("Testing set")
        }
      }
      plot(cowplot::plot_grid(ptr, pts))
    })
    if(is.null(species)){
      output$Tr <- shiny::renderText({paste("Number of training records: ", blocks$records[input$num,1])})
      output$Ts <- shiny::renderText({paste("Number of testing records:   ", blocks$records[input$num,2])})
    } else{
      output$Tb <- shiny::renderTable({blocks$records[input$num,]})
    }
  }
  # starting the shiny app
  shiny::shinyApp(ui = ui, server = server)
}



#' Explore spatial block size
#'
#' This function assists selection of block size. It allows the user to visualise the blocks interactively,
#' viewing the impact of block size on number and arrangement of blocks in the landscape (and optionally on
#' the distribution of species data in those blocks). \strong{Slide} to the slected block size, and click \strong{Apply Changes}
#' to change the block size.
#'
#' The only required argument for this function is the \code{rasterLayer}. The rest are optional. If the \code{rangeTable}
#' is provided, the minimum, maximum and initial ranges for searching the size of spatial blocks will be selected
#' based on the spatial autocorrelation range of covariates. It is also possible to restrict the allowable range
#' of block sizes by using the \code{minRange} and \code{maxRanege} arguments.
#'
#' @inheritParams foldExplorer
#' @param speciesData A simple features (sf) or SpatialPoints object containing species data (response variable). If provided, the species data will be shown on the map.
#' @param species Character value indicating the name of the field in which the species data (response variable e.g. 0s and 1s) are stored.
#' If provided, species presence and absence data will be shown in different colours.
#' @param rangeTable A data.frame created by \code{spatialAutoRange} function containing spatial autocorrelation parameters of all covariates.
#' @param minRange A numeric value to set the minimum possible range for creating spatial blocks. It is used to limit the searching domain of
#' spatial block size.
#' @param maxRange A numeric value to set the maximum possible range for creating spatial blocks. It is used to limit the searching
#' domain of spatial block size.
#'
#' @seealso \code{\link{spatialBlock}}; \code{\link{spatialAutoRange}} for the \code{rangeTable}
#'
#' @return An interactive map with blocks (and optionally species data) superimposed. Note that this can also be opened in a
#' web browser window. When you return to the R console, press "Esc" to return to the prompt.
#' @export
#'
#' @examples
#' \donttest{
#' if(interactive()){
#'
#' # load package data
#' awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))
#' # import presence-absence species data
#' PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
#' # make a sf object from data.frame
#' pa_data <- sf::st_as_sf(PA, coords = c("x", "y"), crs = raster::crs(awt))
#'
#' rangeExplorer(rasterLayer = awt) # the only mandatory input
#'
#' # add species data to add them on the map
#' rangeExplorer(rasterLayer = awt,
#'               speciesData = pa_data,
#'               species = "Species",
#'               rangeTable = NULL,
#'               minRange = 30000, # limit the search domain
#'               maxRange = 100000)
#' }
#' }
rangeExplorer <- function(rasterLayer,
                          speciesData=NULL,
                          species=NULL,
                          rangeTable=NULL,
                          minRange=NULL,
                          maxRange=NULL){
  # check for required packages
  pkg <- c("ggplot2", "shiny", "shinydashboard", "geosphere")
  pkgna <- names(which(sapply(sapply(pkg, find.package, quiet = TRUE), length) == 0))
  if(length(pkgna) > 0){
    nm <- paste(pkgna, collapse = ", ")
    message("This function requires these packages: ", nm,
            "\nWould you like to install them now?\n1: yes\n2: no")
    user <- readline(prompt = paste0("Selection: "))
    if(tolower(user) %in% c("1", "yes", "y")){
      utils::install.packages(pkgna)
    } else{
      stop("Please install these packages: ", nm)
    }
  }
  if(!is.null(speciesData)){
    if(methods::is(speciesData, "SpatialPoints")){
      speciesData <- sf::st_as_sf(speciesData)
    } else if(!methods::is(speciesData, "sf")){
      stop("speciesData should be a sf or SpatialPoints object")
    }
  }
  # plot raster file in ggplot2
  Xmx <- raster::xmax(rasterLayer)
  Xmn <- raster::xmin(rasterLayer)
  Ymx <- raster::ymax(rasterLayer)
  Ymn <- raster::ymin(rasterLayer)
  Ymean <- (Ymx - Ymn)/2
  if(is.na(raster::projection(rasterLayer))){
    mapext <- raster::extent(rasterLayer)[1:4]
    if(mapext >= -180 && mapext <= 180){
      xrange <- geosphere::distGeo(c(Xmx, Ymean), c(Xmn, Ymean))
      yrange <- geosphere::distGeo(c(Xmn, Ymx), c(Xmn, Ymn))
      maxy <- max(c(xrange, yrange))
      resol <- raster::res(rasterLayer) * 111000
      xaxes <- "Longitude"
      yaxes <- "Latitude"
      message("The input layer has no CRS defined. Based on the extent of the input map it is assumed to have an un-projected reference system.\n")
    } else {
      xrange <- Xmx - Xmn
      yrange <- Ymx - Ymn
      maxy <- max(c(xrange, yrange))
      resol <- raster::res(rasterLayer)
      xaxes <- "Easting"
      yaxes <- "Northing"
    }
  } else{
    if(raster::isLonLat(rasterLayer)){
      xrange <- geosphere::distGeo(c(Xmx, Ymean), c(Xmn, Ymean))
      yrange <- geosphere::distGeo(c(Xmn, Ymx), c(Xmn, Ymn))
      maxy <- max(c(xrange, yrange))
      resol <- raster::res(rasterLayer) * 111000
      xaxes <- "Longitude"
      yaxes <- "Latitude"
    } else{
      xrange <- Xmx - Xmn
      yrange <- Ymx - Ymn
      maxy <- max(c(xrange, yrange))
      resol <- raster::res(rasterLayer)
      xaxes <- "Easting"
      yaxes <- "Northing"
    }
  }
  samp <- raster::sampleRegular(rasterLayer[[1]], 5e+05, asRaster=TRUE)
  map_df <- raster::as.data.frame(samp, xy=TRUE, centroids=TRUE, na.rm=TRUE)
  colnames(map_df) <- c("Easting", "Northing", "MAP")
  mid <- stats::median(map_df$MAP)
  basepl <- ggplot2::ggplot() +
    ggplot2::geom_raster(data=map_df, ggplot2::aes_string(y="Northing", x="Easting", fill="MAP")) +
    ggplot2::scale_fill_gradient2(low="darkred", mid="yellow", high="darkgreen", midpoint=mid) +
    ggplot2::guides(fill=FALSE) +
    ggplot2::theme_bw()
  # set the base plot
  if(!is.null(speciesData) && is.null(species)){
    basePlot <- basepl + ggplot2::geom_sf(data = speciesData, alpha = 0.6, color="blue", size = 2)
  } else if(!is.null(speciesData) && !is.null(species)){
    basePlot <- basepl + ggplot2::geom_sf(data = speciesData, ggplot2::aes(color = get(species)), show.legend = "point", alpha = 0.6, size = 2)
  } else if(is.null(speciesData)){
    basePlot <- basepl
  }
  # define min and max range based on the parameters
  if(!is.null(rangeTable)){
    minR <- round(min(rangeTable$range))
    maxR <- round(max(rangeTable$range))
    val <- floor(stats::median(rangeTable$range))
  } else if(is.null(rangeTable) && is.null(maxRange) && !is.null(minRange)){
    minR <- minRange
    maxR <- maxy
    val <- mean(c(xrange, yrange)) / 10
  } else if(is.null(rangeTable) && !is.null(maxRange) && is.null(minRange)){
    minR <- resol[1] * 100
    maxR <- maxRange
    val <- mean(c(xrange, yrange)) / 10
  } else if(is.null(rangeTable) && is.null(maxRange) && is.null(minRange)){
    minR <- resol[1] * 100
    maxR <- maxy
    val <- mean(c(xrange, yrange)) / 10
  } else if(is.null(rangeTable) && !is.null(maxRange) && !is.null(minRange)){
    minR <- minRange
    maxR <- maxRange
    val <- mean(c(xrange, yrange)) / 10
  }
  if(maxR > maxy){ # limit the maximum range to the maximum extent of the map
    message("Selected range is bigger than the extent of the base map. The range is set to the biggest extent (X or Y).\n")
    maxR <- maxy
  }
  minR <- round(minR)
  maxR <- round(maxR)
  val <- round(val)
  # create UI for shniy app
  ui <- shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(title = "Range Explorer"),
    shinydashboard::dashboardSidebar(disable = TRUE),
    shinydashboard::dashboardBody(
      shiny::wellPanel(
        shiny::fluidRow(
          shiny::column(6, offset = 3,
                        shiny::sliderInput(inputId = "num",
                                           label = "Choose a range (m) to create spatial blocks",
                                           value = val, min = minR, max = maxR)
          )
        ),
        shiny::submitButton()
      ),
      shiny::fluidRow(
        shiny::column(12, shiny::plotOutput(outputId = "ggplot")
        )
      )
    )
  )
  # create shiny server and main code
  server <- function(input, output){
    output$ggplot <- shiny::renderPlot({
      subBlocks <- rasterNet(rasterLayer[[1]], resolution=input$num, mask=TRUE, maxpixels=150000)
      p2 <- basePlot + ggplot2::geom_sf(data = subBlocks, color ="red",
                                        fill ="orangered4", alpha = 0.04, size = 0.2) +
        ggplot2::ggtitle("Spatial blocks", subtitle=paste("Using", input$num, "(m) block size")) +
        ggplot2::labs(x = "", y = "", color = species)
      # plot ggplot
      plot(p2)
    })
  }
  # starting the shiny app
  shiny::shinyApp(ui = ui, server = server)
}
