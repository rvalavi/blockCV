## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github("rvalavi/blockCV")

## ---- results="hide", warning=FALSE, message=FALSE-----------------------
# loading the package
library(blockCV)

## ---- fig.height=5, fig.width=7.2, warning=FALSE, message=FALSE----------
# loading raster library
library(raster)
# import raster data
awt <- raster::brick(system.file("extdata", "awt.tif", package = "blockCV"))
# assign the names of the covariates
names(awt) <- c("bc01",  "bc04",  "bc05",  "bc06",  "bc12",  "bc15",  "bc17",  "bc20",  "bc33", "slope", "topo")


## ---- fig.height=4.5, fig.width=7.1--------------------------------------
# import presence-absence species data
PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
# make a SpatialPointsDataFrame object from data.frame
pa_data <- SpatialPointsDataFrame(PA[,c("x", "y")], PA, proj4string=crs(awt))
# see the first few rows
head(pa_data)
# plot species data on the map
plot(awt[[1]]) # plot raster data
points(pa_data[which(pa_data$Species==1), ], col="red") # add presence points
points(pa_data[which(pa_data$Species==0), ], col="blue") # add absence points
legend(x=500000, y=8250000, legend=c("Presence","Absence"), col=c(2, 4), pch=c(1,1), bty="n")


## ------------------------------------------------------------------------
# import presence-background species data
PB <- read.csv(system.file("extdata", "PB.csv", package = "blockCV"))
# make a SpatialPointsDataFrame object from data.frame
pb_data <- SpatialPointsDataFrame(PB[,c("x", "y")], PB, proj4string=crs(awt))
# number of presence and background records
table(pb_data$Species)


