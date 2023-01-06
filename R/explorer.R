#' Explore spatial block size
#'
#' This function is deprecated and will be removed in future updates! Please use \code{\link{cv_block_size}} instead!
#'
#' @param raster layer for make plot
#' @param speciesData a simple features (sf) or SpatialPoints object containing species data (response variable). If provided, the species data will be shown on the map.
#' @param species character value indicating the name of the field in which the species data (response variable e.g. 0s and 1s) are stored.
#' If provided, species presence and absence data will be shown in different colours.
#' @param rangeTable deprecated option!
#' @param minRange a numeric value to set the minimum possible range for creating spatial blocks. It is used to limit the searching domain of
#' spatial block size.
#' @param maxRange a numeric value to set the maximum possible range for creating spatial blocks. It is used to limit the searching
#' domain of spatial block size.
#'
#' @seealso \code{\link{cv_block_size}}
#'
#' @export
rangeExplorer <- function(rasterLayer,
                          speciesData=NULL,
                          species=NULL,
                          rangeTable=NULL,
                          minRange=NULL,
                          maxRange=NULL){

  message("This function is deprecated! Please use 'cv_block_size' instead.")

  cv_block_size(r = rasterLayer, # priority
                x = speciesData,
                column = species,
                min_size = minRange,
                max_size = maxRange)

  # # check for required packages
  # pkg <- c("ggplot2", "shiny", "shinydashboard", "geosphere")
  # pkgna <- names(which(sapply(sapply(pkg, find.package, quiet = TRUE), length) == 0))
  # if(length(pkgna) > 0){
  #   nm <- paste(pkgna, collapse = ", ")
  #   message("This function requires these packages: ", nm,
  #           "\nWould you like to install them now?\n1: yes\n2: no")
  #   user <- readline(prompt = paste0("Selection: "))
  #   if(tolower(user) %in% c("1", "yes", "y")){
  #     utils::install.packages(pkgna)
  #   } else{
  #     stop("Please install these packages: ", nm)
  #   }
  # }
  # if(!is.null(speciesData)){
  #   if(methods::is(speciesData, "SpatialPoints")){
  #     speciesData <- sf::st_as_sf(speciesData)
  #   } else if(!methods::is(speciesData, "sf")){
  #     stop("speciesData should be a sf or SpatialPoints object")
  #   }
  # }
  # # plot raster file in ggplot2
  # Xmx <- raster::xmax(rasterLayer)
  # Xmn <- raster::xmin(rasterLayer)
  # Ymx <- raster::ymax(rasterLayer)
  # Ymn <- raster::ymin(rasterLayer)
  # Ymean <- (Ymx - Ymn)/2
  # if(is.na(raster::projection(rasterLayer))){
  #   mapext <- raster::extent(rasterLayer)[1:4]
  #   if(mapext >= -180 && mapext <= 180){
  #     xrange <- geosphere::distGeo(c(Xmx, Ymean), c(Xmn, Ymean))
  #     yrange <- geosphere::distGeo(c(Xmn, Ymx), c(Xmn, Ymn))
  #     maxy <- max(c(xrange, yrange))
  #     resol <- raster::res(rasterLayer) * 111000
  #     xaxes <- "Longitude"
  #     yaxes <- "Latitude"
  #     message("The input layer has no CRS defined. Based on the extent of the input map it is assumed to have an un-projected reference system.\n")
  #   } else {
  #     xrange <- Xmx - Xmn
  #     yrange <- Ymx - Ymn
  #     maxy <- max(c(xrange, yrange))
  #     resol <- raster::res(rasterLayer)
  #     xaxes <- "Easting"
  #     yaxes <- "Northing"
  #   }
  # } else{
  #   if(raster::isLonLat(rasterLayer)){
  #     xrange <- geosphere::distGeo(c(Xmx, Ymean), c(Xmn, Ymean))
  #     yrange <- geosphere::distGeo(c(Xmn, Ymx), c(Xmn, Ymn))
  #     maxy <- max(c(xrange, yrange))
  #     resol <- raster::res(rasterLayer) * 111000
  #     xaxes <- "Longitude"
  #     yaxes <- "Latitude"
  #   } else{
  #     xrange <- Xmx - Xmn
  #     yrange <- Ymx - Ymn
  #     maxy <- max(c(xrange, yrange))
  #     resol <- raster::res(rasterLayer)
  #     xaxes <- "Easting"
  #     yaxes <- "Northing"
  #   }
  # }
  # samp <- raster::sampleRegular(rasterLayer[[1]], 5e+05, asRaster=TRUE)
  # map_df <- raster::as.data.frame(samp, xy=TRUE, centroids=TRUE, na.rm=TRUE)
  # colnames(map_df) <- c("Easting", "Northing", "MAP")
  # mid <- stats::median(map_df$MAP)
  # basepl <- ggplot2::ggplot() +
  #   ggplot2::geom_raster(data=map_df, ggplot2::aes_string(y="Northing", x="Easting", fill="MAP")) +
  #   ggplot2::scale_fill_gradient2(low="darkred", mid="yellow", high="darkgreen", midpoint=mid) +
  #   ggplot2::guides(fill=FALSE) +
  #   ggplot2::theme_bw()
  # # set the base plot
  # if(!is.null(speciesData) && is.null(species)){
  #   basePlot <- basepl + ggplot2::geom_sf(data = speciesData, alpha = 0.6, color="blue", size = 2)
  # } else if(!is.null(speciesData) && !is.null(species)){
  #   basePlot <- basepl + ggplot2::geom_sf(data = speciesData, ggplot2::aes(color = get(species)), show.legend = "point", alpha = 0.6, size = 2)
  # } else if(is.null(speciesData)){
  #   basePlot <- basepl
  # }
  # # define min and max range based on the parameters
  # if(!is.null(rangeTable)){
  #   minR <- round(min(rangeTable$range))
  #   maxR <- round(max(rangeTable$range))
  #   val <- floor(stats::median(rangeTable$range))
  # } else if(is.null(rangeTable) && is.null(maxRange) && !is.null(minRange)){
  #   minR <- minRange
  #   maxR <- maxy
  #   val <- mean(c(xrange, yrange)) / 10
  # } else if(is.null(rangeTable) && !is.null(maxRange) && is.null(minRange)){
  #   minR <- resol[1] * 100
  #   maxR <- maxRange
  #   val <- mean(c(xrange, yrange)) / 10
  # } else if(is.null(rangeTable) && is.null(maxRange) && is.null(minRange)){
  #   minR <- resol[1] * 100
  #   maxR <- maxy
  #   val <- mean(c(xrange, yrange)) / 10
  # } else if(is.null(rangeTable) && !is.null(maxRange) && !is.null(minRange)){
  #   minR <- minRange
  #   maxR <- maxRange
  #   val <- mean(c(xrange, yrange)) / 10
  # }
  # if(maxR > maxy){ # limit the maximum range to the maximum extent of the map
  #   message("Selected range is bigger than the extent of the base map. The range is set to the biggest extent (X or Y).\n")
  #   maxR <- maxy
  # }
  # minR <- round(minR)
  # maxR <- round(maxR)
  # val <- round(val)
  # # create UI for shniy app
  # ui <- shinydashboard::dashboardPage(
  #   shinydashboard::dashboardHeader(title = "Range Explorer"),
  #   shinydashboard::dashboardSidebar(disable = TRUE),
  #   shinydashboard::dashboardBody(
  #     shiny::wellPanel(
  #       shiny::fluidRow(
  #         shiny::column(6, offset = 3,
  #                       shiny::sliderInput(inputId = "num",
  #                                          label = "Choose a range (m) to create spatial blocks",
  #                                          value = val, min = minR, max = maxR)
  #         )
  #       ),
  #       shiny::submitButton()
  #     ),
  #     shiny::fluidRow(
  #       shiny::column(12, shiny::plotOutput(outputId = "ggplot")
  #       )
  #     )
  #   )
  # )
  # # create shiny server and main code
  # server <- function(input, output){
  #   output$ggplot <- shiny::renderPlot({
  #     subBlocks <- rasterNet(rasterLayer[[1]], resolution=input$num, mask=TRUE, maxpixels=150000)
  #     p2 <- basePlot + ggplot2::geom_sf(data = subBlocks, color ="red",
  #                                       fill ="orangered4", alpha = 0.04, size = 0.2) +
  #       ggplot2::ggtitle("Spatial blocks", subtitle=paste("Using", input$num, "(m) block size")) +
  #       ggplot2::labs(x = "", y = "", color = species)
  #     # plot ggplot
  #     plot(p2)
  #   })
  # }
  # # starting the shiny app
  # shiny::shinyApp(ui = ui, server = server)
}
