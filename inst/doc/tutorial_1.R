## ---- eval=FALSE--------------------------------------------------------------
#  # install stable version from CRAN
#  install.packages("blockCV", dependencies = TRUE)
#  
#  # install latest update from GitHub
#  remotes::install_github("rvalavi/blockCV", dependencies = TRUE)
#  

## ----message=TRUE, warning=TRUE-----------------------------------------------
# loading the package
library(blockCV)


## ---- fig.height=5, fig.width=7.2, warning=FALSE, message=FALSE---------------
library(sf) # working with spatial vector data
library(terra) # working with spatial raster data
library(tmap) # plotting maps

# load raster data
# the pipe operator |> is available for R version 4.1 or higher
rasters <- system.file("extdata/au/", package = "blockCV") |>
  list.files(full.names = TRUE) |>
  terra::rast()


## ---- fig.height=4.5, fig.width=7.1-------------------------------------------
# load species presence-absence data and convert to sf
points <- read.csv(system.file("extdata/", "species.csv", package = "blockCV"))
head(points)


## ---- fig.height=4.5, fig.width=7.1-------------------------------------------
pa_data <- sf::st_as_sf(points, coords = c("x", "y"), crs = 7845)

## ---- fig.height=4.5, fig.width=7.1-------------------------------------------
tm_shape(rasters[[1]]) +
  tm_raster(legend.show = FALSE, n = 30, palette = gray.colors(10)) +
  tm_shape(pa_data) +
  tm_dots(col = "occ", style = "cat", size = 0.1)


## ---- results='hide', fig.keep='all', warning=FALSE, message=FALSE, fig.height=5, fig.width=7----
sb1 <- cv_spatial(x = pa_data,
                  column = "occ", # the response column (binary or multi-class)
                  k = 5, # number of folds
                  size = 350000, # size of the blocks in metres
                  selection = "random", # random blocks-to-fold
                  iteration = 50, # find evenly dispersed folds
                  biomod2 = TRUE) # also create folds for biomod2


## ---- warning=FALSE, message=FALSE, fig.height=5, fig.width=7-----------------
sb2 <- cv_spatial(x = pa_data,
                  column = "occ",
                  r = rasters, # optionally add a raster layer
                  k = 5, 
                  size = 350000, 
                  hexagon = FALSE, # use square blocks
                  selection = "random",
                  progress = FALSE, # turn off progress bar for vignette
                  iteration = 50, 
                  biomod2 = TRUE)


## ----warning=FALSE, message=FALSE, fig.height=5, fig.width=7------------------
# systematic fold assignment 
# and also use row/column for creating blocks instead of size
sb3 <- cv_spatial(x = pa_data,
                  column = "occ",
                  rows_cols = c(12, 10),
                  hexagon = FALSE,
                  selection = "systematic")


## ----warning=FALSE, message=FALSE, fig.height=5, fig.width=7------------------
# checkerboard block to CV fold assignment
sb4 <- cv_spatial(x = pa_data,
                  column = "occ",
                  size = 350000,
                  hexagon = FALSE,
                  selection = "checkerboard")


## ----warning=FALSE, message=FALSE, fig.height=5, fig.width=7------------------
tm_shape(sb4$blocks) +
  tm_fill(col = "folds", style = "cat")


## -----------------------------------------------------------------------------
# spatial clustering
set.seed(6)
scv <- cv_cluster(x = pa_data,
                  column = "occ", # optional: counting number of train/test records
                  k = 5)

## ----warning=FALSE, message=FALSE---------------------------------------------
# environmental clustering
set.seed(6)
ecv <- cv_cluster(x = pa_data,
                  column = "occ",
                  r = rasters,
                  k = 5, 
                  scale = TRUE)


## ----results='hide', fig.keep='all'-------------------------------------------
bloo <- cv_buffer(x = pa_data,
                  column = "occ",
                  size = 350000)


## ----fig.height=5, fig.width=7------------------------------------------------
nncv <- cv_nndm(x = pa_data,
                column = "occ",
                r = rasters,
                size = 350000,
                num_sample = 5000, 
                sampling = "regular",
                min_train = 0.1,
                plot = TRUE)


## ----warning=FALSE, message=FALSE, fig.height=6, fig.width=8------------------
cv_plot(cv = scv, 
        x = pa_data)


## ----warning=FALSE, message=FALSE, fig.height=5, fig.width=8------------------
cv_plot(cv = bloo,
        x = pa_data,
        num_plots = c(1, 50, 100)) # only show folds 1, 50 and 100


## ----warning=FALSE, message=FALSE, fig.height=5, fig.width=7------------------
cv_plot(cv = sb1,
        r = rasters,
        raster_colors = terrain.colors(10, alpha = 0.5),
        label_size = 4) 


## ----fig.height=4, fig.width=6------------------------------------------------
cv_similarity(cv = ecv, # the environmental clustering
              x = pa_data, 
              r = rasters, 
              progress = FALSE)


## ---- results='hide', fig.keep='all', warning=FALSE, message=FALSE, fig.height=5, fig.width=7.2----
sac1 <- cv_spatial_autocor(r = rasters, 
                           num_sample = 5000)


## -----------------------------------------------------------------------------
# class of the output result
class(sac1)

## -----------------------------------------------------------------------------
# summary statistics of the output
summary(sac1)

## ---- warning=FALSE, message=FALSE, fig.height=5, fig.width=7.2---------------
sac2 <- cv_spatial_autocor(x = pa_data, 
                           column =  "occ")


## ---- eval=TRUE, fig.height=4, fig.width=7------------------------------------
library(automap)

plot(sac2$variograms[[1]])


## ---- eval=FALSE--------------------------------------------------------------
#  cv_block_size(r = rasters)
#  

## ---- eval=FALSE--------------------------------------------------------------
#  cv_block_size(x = pa_data,
#                column = "occ") # optionally add the response
#  

## ---- eval=FALSE--------------------------------------------------------------
#  cv_block_size(x = pa_data,
#                column = "occ",
#                r = rasters,
#                min_size = 2e5,
#                max_size = 9e5)
#  

