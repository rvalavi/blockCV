#' Explore the generated folds
#'
#' A function for visualising the generated folds on a map, and allowing interactive exploration of the data in the folds, using the RStudio Shiny app.
#'
#' @param blocks An SpatialBlock, EnvironmentalBlock or BufferedBlock object.
#' @param rasterLayer A RasterLayer, RasterBrick or RasterStack object as background map for visualisation.
#' @inheritParams buffering
#'
#' @import ggplot2
#' @import shiny
#' @import shinydashboard
#' @seealso \code{\link{spatialBlock}}, \code{\link{buffering}} and \code{\link{envBlock}}
#'
#' @return An interactive map showing folds and the species data, that can be used  to explore folds. Note that this can also
#' be opened in a web browser window. When you return to the R console, press “Esc” to return to the prompt.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # load package data
#' awt <- raster::brick(system.file("extdata", "awt.tif", package = "blockCV"))
#' # import presence-absence species data
#' PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
#' # make a SpatialPointsDataFrame object from data.frame
#' pa_data <- sp::SpatialPointsDataFrame(PA[,c("x", "y")], PA, proj4string=raster::crs(awt))
#'
#' # spatial blocking by specified range and random assignment
#' sb <- spatialBlock(speciesData = pa_data,
#'                    species = "Species",
#'                    rasterLayer = awt,
#'                    theRange = 66000,
#'                    k = 5,
#'                    selection = 'random',
#'                    iteration = 250,
#'                    numLimit = NULL,
#'                    maskBySpecies = FALSE,
#'                    biomod2Format = TRUE)
#'
#' foldExplorer(sb, awt, pa_data)
#'
#' # buffering with presence-absence data
#' bf <- buffering(speciesData= pa_data,
#'                 species= "Species", # to count the number of presences and absences
#'                 theRange= 66500,
#'                 spDataType = "PA",
#'                 progress = T)
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
#'
foldExplorer <- function(blocks, rasterLayer, speciesData){
  # testing the input arguments
  if(is.null(rasterLayer)){
    stop("A raster layer should be provided")
  } else if(is.null(speciesData)){
    stop("Species data should be provided")
  } else if(is.null(blocks)){
    stop("An object of SpatialBlock, EnvironmentalBlock or BufferedBlock is needed")
  }
  # select the default value
  palpha <- 0.6
  if(class(blocks) == "BufferedBlock"){
    polyObj <- NULL
    if(blocks$dataType == "PB"){
      palpha <- 0.5
    }
  } else if(class(blocks) == "SpatialBlock"){
    polyObj <- blocks$blocks
  } else if(class(blocks) == "EnvironmentalBlock"){
    polyObj <- NULL
  } else{
    stop("blocks must be an object of SpatialBlock, EnvironmentalBlock or BufferedBlock")
  }
  folds <- blocks$folds
  kmax <- length(folds)
  species <- blocks$species
  # set x and y coordinates
  coor <- sp::coordinates(speciesData)
  coor <- as.data.frame(coor)
  speciesData@data <- cbind(speciesData@data, coor)
  # plot raster file in ggplot2
  samp <- raster::sampleRegular(rasterLayer[[1]], 5e+05, asRaster=TRUE)
  map.df <- raster::as.data.frame(samp, xy=TRUE, centroids=TRUE, na.rm=TRUE)
  colnames(map.df) <- c("Easting", "Northing", "MAP")
  mid <- stats::median(map.df$MAP)
  basePlot <- ggplot2::ggplot(data=map.df, aes(y=Northing, x=Easting)) + geom_raster(aes(fill=MAP)) + coord_fixed() +
    scale_fill_gradient2(low="darkred", mid="yellow", high="darkgreen", midpoint=mid) + guides(fill=FALSE)
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
                                                value = 2, min = 1, max = kmax, step = 1)
            )
          ),
          shiny::fluidRow("Record's Statistics"),
          if(is.null(species)){
            shiny::fluidRow(shiny::textOutput(outputId = "Tr"),
                            shiny::textOutput(outputId = "Ts"))
          } else{
            shiny::fluidRow(shiny::textOutput(outputId = "TrP"),
                            shiny::textOutput(outputId = "TrA"),
                            shiny::textOutput(outputId = "TsP"),
                            shiny::textOutput(outputId = "TsA"))
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
      training2 <- sp::SpatialPoints(training, proj4string = crs(training))
      testing <- speciesData[testSet, ]
      testing2 <- sp::SpatialPoints(testing, proj4string = crs(testing))
      if(is.null(species)){
        trdf <- data.frame(ID=1:length(training2))
        training <- sp::SpatialPointsDataFrame(training2, trdf)
        coor <- sp::coordinates(training)
        coor <- as.data.frame(coor)
        training@data <- cbind(training@data, coor)
        names(training)[2:3] <- c("Easting", "Northing")
        tsdf <- data.frame(ID=1:length(testing2))
        testing <- sp::SpatialPointsDataFrame(testing2, tsdf)
        coor <- sp::coordinates(testing)
        coor <- as.data.frame(coor)
        testing@data <- cbind(testing@data, coor)
        names(testing)[2:3] <- c("Easting", "Northing")
        rm(training2, testing2, trdf, tsdf, coor)
      } else{
        trdf <- data.frame(ID=1:length(training2))
        trdf$Species <- training@data[,species]
        training <- sp::SpatialPointsDataFrame(training2, trdf)
        coor <- sp::coordinates(training)
        coor <- as.data.frame(coor)
        training@data <- cbind(training@data, coor)
        names(training)[3:4] <- c("Easting", "Northing")
        training$Species <- as.factor(training$Species)
        tsdf <- data.frame(ID=1:length(testing2))
        tsdf$Species <- testing@data[,species]
        testing <- sp::SpatialPointsDataFrame(testing2, tsdf)
        coor <- sp::coordinates(testing)
        coor <- as.data.frame(coor)
        testing@data <- cbind(testing@data, coor)
        names(testing)[3:4] <- c("Easting", "Northing")
        testing$Species <- as.factor(testing$Species)
        rm(training2, testing2, trdf, tsdf, coor)
      }
      thePoly <- polyObj[polyObj$folds==input$num, ]
      plotPoly <- ggplot2::fortify(thePoly)
      if(is.null(species)){
        if(class(blocks) == "SpatialBlock"){
          ptr <- basePlot + geom_polygon(aes(x = long, y = lat, group=group), data = plotPoly,
                                         color ="red", fill ="orangered4", alpha = 0.04, size = 0.2) +
            geom_point(aes(x = Easting, y = Northing), data = training@data, alpha = 0.6, color="blue", size = 2) +
            ggtitle('Training set')
          # ploting test data
          pts <- basePlot + geom_polygon(aes(x = long, y = lat, group=group), data = plotPoly,
                                         color ="red", fill ="orangered4", alpha = 0.04, size = 0.2) +
            geom_point(aes(x = Easting, y = Northing), data = testing@data, alpha = 0.6, color="blue", size = 2) +
            ggtitle('Testing set')
        } else if(class(blocks) == "BufferedBlock" || class(blocks) == "EnvironmentalBlock"){
          ptr <- basePlot + geom_point(aes(x = Easting, y = Northing), data = training@data, alpha = 0.6, color="blue", size = 2) +
            ggtitle('Training set')
          # ploting test data
          pts <- basePlot + geom_point(aes(x = Easting, y = Northing), data = testing@data, alpha = 0.6, color="blue", size = 2) +
            ggtitle('Testing set')
        }
      } else{
        if(class(blocks) == "SpatialBlock"){
          ptr <- basePlot + geom_polygon(aes(x = long, y = lat, group=group), data = plotPoly,
                                         color ="red", fill ="orangered4", alpha = 0.04, size = 0.2) +
            geom_point(aes(x = Easting, y = Northing, colour=Species), data = training@data, alpha = 0.6, size = 2) +
            scale_colour_manual(values = c("blue", "red"),  labels = c("Absence/Background", "Presence")) +
            ggtitle('Training set')
          # ploting test data
          if(blocks$records[input$num,3] == 0){
            pts <- basePlot + geom_polygon(aes(x = long, y = lat, group=group), data = plotPoly,
                                           color ="red", fill ="orangered4", alpha = 0.04, size = 0.2) +
              geom_point(aes(x = Easting, y = Northing, colour=Species), data = testing@data, alpha = 0.6, size = 2) +
              scale_colour_manual(values = c("blue"),  labels = c("Absence/Background")) +
              ggtitle('Testing set')
          } else if(blocks$records[input$num,4] == 0){
            pts <- basePlot + geom_polygon(aes(x = long, y = lat, group=group), data = plotPoly,
                                           color ="red", fill ="orangered4", alpha = 0.04, size = 0.2) +
              geom_point(aes(x = Easting, y = Northing, colour=Species), data = testing@data, alpha = 0.6, size = 2) +
              scale_colour_manual(values = c("red"),  labels = c("Presence")) +
              ggtitle('Testing set')
          } else{
            pts <- basePlot + geom_polygon(aes(x = long, y = lat, group=group), data = plotPoly,
                                           color ="red", fill ="orangered4", alpha = 0.04, size = 0.2) +
              geom_point(aes(x = Easting, y = Northing, colour=Species), data = testing@data, alpha = 0.6, size = 2) +
              scale_colour_manual(values = c("blue", "red"),  labels = c("Absence/Background", "Presence")) +
              ggtitle('Testing set')
          }
        } else if(class(blocks) == "EnvironmentalBlock"){
          ptr <- basePlot + geom_point(aes(x = Easting, y = Northing, colour=Species), data = training@data, alpha = palpha, size = 2) +
            scale_colour_manual(values = c("blue", "red"),  labels = c("Absence/Background", "Presence")) +
            ggtitle('Training set')
          # ploting test data
          if(blocks$records[input$num,3] == 0){
            pts <- basePlot + geom_point(aes(x = Easting, y = Northing, colour=Species), data = testing@data, alpha = 0.6, size = 2) +
              scale_colour_manual(values = c("blue"),  labels = c("Absence/Background")) +
              ggtitle('Testing set')
          } else if(blocks$records[input$num,4] == 0){
            pts <- basePlot + geom_point(aes(x = Easting, y = Northing, colour=Species), data = testing@data, alpha = 0.6, size = 2) +
              scale_colour_manual(values = c("red"),  labels = c("Presence")) +
              ggtitle('Testing set')
          } else{
            pts <- basePlot + geom_point(aes(x = Easting, y = Northing, colour=Species), data = testing@data, alpha = 0.6, size = 2) +
              scale_colour_manual(values = c("blue", "red"),  labels = c("Absence/Background", "Presence")) +
              ggtitle('Testing set')
          }
        } else if(class(blocks) == "BufferedBlock"){
          # ploting test data
          if(blocks$dataType == "PA"){
            ptr <- basePlot + geom_point(aes(x = Easting, y = Northing, colour=Species), data = training@data, alpha = palpha, size = 2) +
              scale_colour_manual(values = c("blue", "red"),  labels = c("Absence", "Presence")) +
              ggtitle('Training set')
            if(blocks$records[input$num,3] == 0){
              pts <- basePlot + geom_point(aes(x = Easting, y = Northing, colour=Species), data = testing@data, alpha = 0.6, size = 2) +
                scale_colour_manual(values = c("blue"),  labels = c("Absence")) +
                ggtitle('Testing set')
            } else{
              pts <- basePlot + geom_point(aes(x = Easting, y = Northing, colour=Species), data = testing@data, alpha = 0.6, size = 2) +
                scale_colour_manual(values = c("red"),  labels = c("Presence")) +
                ggtitle('Testing set')
            }
          } else if(blocks$dataType == "PB"){
            ptr <- basePlot + geom_point(aes(x = Easting, y = Northing, colour=Species), data = training@data, alpha = palpha, size = 2) +
              scale_colour_manual(values = c("blue", "red"),  labels = c("Background", "Presence")) +
              ggtitle('Training set')
            if(blocks$records[input$num,4] == 0){
              pts <- basePlot + geom_point(aes(x = Easting, y = Northing, colour=Species), data = testing@data, alpha = 0.6, size = 2) +
                scale_colour_manual(values = c("red"),  labels = c("Presence")) +
                ggtitle('Testing set')
            } else{
              pts <- basePlot + geom_point(aes(x = Easting, y = Northing, colour=Species), data = testing@data, alpha = 0.6, size = 2) +
                scale_colour_manual(values = c("blue", "red"),  labels = c("Background", "Presence")) +
                ggtitle('Testing set')
            }
          }
        }
      }
      multiplot(ptr, pts)
    })
    if(is.null(species)){
      output$Tr <- shiny::renderText({paste("No. training records: ", blocks$records[input$num,1])})
      output$Ts <- shiny::renderText({paste("No. testing records:   ", blocks$records[input$num,2])})
    } else{
      if(class(blocks) == "BufferedBlock"){
        if(blocks$dataType == "PB"){
          output$TrP <- shiny::renderText({paste("No. training presence records:   ", blocks$records[input$num,1])})
          output$TrA <- shiny::renderText({paste("No. training background records:", blocks$records[input$num,2])})
          output$TsP <- shiny::renderText({paste("No. testing presence records:   ", blocks$records[input$num,3])})
          output$TsA <- shiny::renderText({paste("No. testing background records:", blocks$records[input$num,4])})
        } else{
          output$TrP <- shiny::renderText({paste("No. training presence records:   ", blocks$records[input$num,1])})
          output$TrA <- shiny::renderText({paste("No. training absence records:", blocks$records[input$num,2])})
          output$TsP <- shiny::renderText({paste("No. testing presence records:   ", blocks$records[input$num,3])})
          output$TsA <- shiny::renderText({paste("No. testing absence records:", blocks$records[input$num,4])})
        }
      } else{
        output$TrP <- shiny::renderText({paste("No. training presence records:   ", blocks$records[input$num,1])})
        output$TrA <- shiny::renderText({paste("No. training absence/background records:", blocks$records[input$num,2])})
        output$TsP <- shiny::renderText({paste("No. testing presence records:   ", blocks$records[input$num,3])})
        output$TsA <- shiny::renderText({paste("No. testing absence/background records:", blocks$records[input$num,4])})
      }
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
#' @param speciesData A SpatialPoints* object containing species data. If provided, the species data will be shown on the map.
#' @param species Character value indicating the name of the field in which the species presence-absence/backgrouns data (0s and 1s) are stored.
#' If provided, species presence and absence data will be shown in different colours.
#' @param rangeTable A data.frame created by \code{spatialAutoRange} function containing spatial autocorrelation parameters of all covariates.
#' @param minRange A numeric value to set the minimum possible range for creating spatial blocks. It is used to limit the searching domain of spatial block size.
#' @param maxRange A numeric value to set the maximum possible range for creating spatial blocks. It is used to limit the searching domain of spatial block size.
#'
#' @import ggplot2
#' @import shiny
#' @import shinydashboard
#'
#' @seealso \code{\link{spatialBlock}}; \code{\link{spatialAutoRange}} for the \code{rangeTable}
#'
#' @return An interactive map with blocks (and optionally species data) superimposed. Note that this can also be opened in a web browser
#' window. When you return to the R console, press “Esc” to return to the prompt.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # load package data
#' awt <- raster::brick(system.file("extdata", "awt.tif", package = "blockCV"))
#' # import presence-absence species data
#' PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
#' # make a SpatialPointsDataFrame object from data.frame
#' pa_data <- sp::SpatialPointsDataFrame(PA[,c("x", "y")], PA, proj4string=raster::crs(awt))
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
#'
#' }
rangeExplorer <- function(rasterLayer, speciesData=NULL, species=NULL, rangeTable=NULL, minRange=NULL, maxRange=NULL){
  # plot raster file in ggplot2
  samp <- raster::sampleRegular(rasterLayer[[1]], 5e+05, asRaster=TRUE)
  map.df <- raster::as.data.frame(samp, xy=TRUE, centroids=TRUE, na.rm=TRUE)
  colnames(map.df) <- c("Easting", "Northing", "MAP")
  mid <- mean(map.df$MAP)
  basepl <- ggplot(data=map.df, aes(y=Northing, x=Easting)) + geom_raster(aes(fill=MAP)) + coord_fixed() +
    scale_fill_gradient2(low="darkred", mid="yellow", high="darkgreen", midpoint=mid) + guides(fill=FALSE)
  if(!is.null(speciesData)){
    coor <- sp::coordinates(speciesData)
    coor <- as.data.frame(coor)
    speciesData@data <- cbind(speciesData@data, coor)
  }
  # set the base plot
  if(!is.null(speciesData) && is.null(species)){
    speciesXY <- sp::SpatialPoints(speciesData)
    speciesXY <- as.data.frame(speciesXY)
    names(speciesXY) <- c("Easting", "Northing")
    basePlot <- basepl + geom_point(aes(x = Easting, y = Northing),
                                    data = speciesXY, alpha = 0.6, color="blue", size = 2)
  } else if(!is.null(speciesData) && !is.null(species)){
    plotData2 <- sp::SpatialPoints(speciesData, proj4string = crs(speciesData))
    spdf <- data.frame(ID=1:length(plotData2))
    spdf$Species <- speciesData@data[,species]
    plotData <- sp::SpatialPointsDataFrame(plotData2, spdf)
    coor <- sp::coordinates(plotData)
    coor <- as.data.frame(coor)
    plotData@data <- cbind(plotData@data, coor)
    names(plotData)[3:4] <- c("Easting", "Northing")
    plotData$Species <- as.factor(plotData$Species)
    basePlot <- basepl + geom_point(aes(x = Easting, y = Northing, colour=Species), data = plotData@data, alpha = 0.6, size = 2) +
      scale_colour_manual(values = c("blue", "red"),  labels = c("Absence/Background", "Presence"))
  } else if(is.null(speciesData) && is.null(species)){
    basePlot <- basepl
  } else if(is.null(speciesData) && !is.null(species)){
    basePlot <- basepl
  }
  Xmx <- raster::xmax(rasterLayer)
  Xmn <- raster::xmin(rasterLayer)
  Ymx <- raster::ymax(rasterLayer)
  Ymn <- raster::ymin(rasterLayer)
  Ymean <- (Ymx - Ymn)/2
  if(is.na(sp::proj4string(rasterLayer))){
    mapext <- raster::extent(rasterLayer)[1:4]
    if(mapext >= -180 && mapext <= 180){
      xrange <- geosphere::distGeo(c(Xmx, Ymean), c(Xmn, Ymean))
      yrange <- geosphere::distGeo(c(Xmn, Ymx), c(Xmn, Ymn))
      maxy <- max(c(xrange, yrange))
      resol <- raster::res(rasterLayer) * 111000
      warning("The input layer has no CRS defined. Based on the extent of the input map it is assumed to have an un-projected reference system")
    } else {
      xrange <- Xmx - Xmn
      yrange <- Ymx - Ymn
      maxy <- max(c(xrange, yrange))
      resol <- raster::res(rasterLayer)
    }
  } else{
    if(sp::is.projected(sp::SpatialPoints((matrix(1:10, 5, byrow=FALSE)), proj4string=crs(rasterLayer)))){
      xrange <- Xmx - Xmn
      yrange <- Ymx - Ymn
      maxy <- max(c(xrange, yrange))
      resol <- raster::res(rasterLayer)
    } else{
      xrange <- geosphere::distGeo(c(Xmx, Ymean), c(Xmn, Ymean))
      yrange <- geosphere::distGeo(c(Xmn, Ymx), c(Xmn, Ymn))
      maxy <- max(c(xrange, yrange))
      resol <- raster::res(rasterLayer) * 111000
    }
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
    message("Selected range was bigger than the extent of base map. It has been set to the biggest extent (X or Y")
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
      rasterNet1 <- rasterNet(rasterLayer[[1]], resolution=input$num)
      net <- raster::rasterToPolygons(rasterNet1)
      points <- raster::rasterToPoints(rasterLayer[[1]], spatial=TRUE)
      if(nrow(points) > 1000000){
        points2 <- points[sample(1:nrow(points), 150000, replace=FALSE), ]
        subBlocks <- raster::intersect(net, points2)
      } else  subBlocks <- raster::intersect(net, points)
      p2 <- basePlot + ggtitle('Spatial blocks', subtitle=paste('Based on', input$num, '(m) distance')) +
        geom_polygon(aes(x = long, y = lat, group=id),
                     data = subBlocks, color ="red",
                     fill ="orangered4",
                     alpha = 0.04,
                     size = 0.2)
      # plot ggplot
      plot(p2)
    })
  }
  # starting the shiny app
  shiny::shinyApp(ui = ui, server = server)
}
