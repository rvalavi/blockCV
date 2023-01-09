## ---- eval=FALSE--------------------------------------------------------------
#  # install stable version from CRAN
#  install.packages("blockCV", dependencies = TRUE)
#  
#  # install latest update from GitHub
#  remotes::install_github("rvalavi/blockCV", dependencies = TRUE)
#  

## ---- results="hide", warning=FALSE, message=FALSE----------------------------
# loading the package
library(blockCV)


## ---- fig.height=5, fig.width=7.2, warning=FALSE, message=FALSE---------------
library(ggplot2) # plotting
library(rasterVis) # plotting raster data
library(sf) # working with spatial vector data
library(terra) # working with saptial raster data
library(tmap) # plotting spatial data

# load raster data
rasters <- system.file("extdata/au/", package = "blockCV") |>
  list.files(full.names = TRUE) |>
  terra::rast()


## ---- fig.height=4.5, fig.width=7.1-------------------------------------------
# load species presence-asence data and convert to sf
points <- read.csv(system.file("extdata/", "species.csv", package = "blockCV"))
head(points)


## ---- fig.height=4.5, fig.width=7.1-------------------------------------------
pa_data <- sf::st_as_sf(points, coords = c("x", "y"), crs = 7845)

## ---- fig.height=4.5, fig.width=7.1-------------------------------------------
tm_shape(rasters[[1]]) +
  tm_raster(legend.show = FALSE, n = 30, palette = gray.colors(10)) +
  tm_shape(pa_data) +
  tm_dots(col = "occ", style = "cat", size = 0.1)


## ---- fig.height=4.5, fig.width=7.1-------------------------------------------
# plot species data on the map
rasterVis::gplot(rasters[[1]]) +
  geom_tile(aes(fill = value)) +
  scale_fill_distiller(palette = 4, na.value = NA) +
  guides(fill = "none") +
  geom_sf(data = pa_data, 
          aes(color = as.factor(occ)),
          inherit.aes = FALSE) +
  scale_color_manual(values = c("red", "yellow")) +
  theme_minimal()


## ---- results='hide', fig.keep='all', warning=FALSE, message=FALSE, fig.height=5, fig.width=7----
# spatial blocking by specified range with random assignment
sb <- cv_spatial(x = pa_data,
                 column = "occ",
                 k = 5,
                 size = 350000,
                 selection = "random",
                 iteration = 50)


## ----warning=FALSE, message=FALSE, fig.height=5, fig.width=7------------------
sb2 <- cv_spatial(x = pa_data,
                  column = "occ",
                  k = 5,
                  size = 350000,
                  hexagon = FALSE,
                  selection = "random",
                  iteration = 50,
                  progress = FALSE)

## ----warning=FALSE, message=FALSE, fig.height=5, fig.width=7------------------
# spatial blocking by rows with systematic assignment
sb3 <- cv_spatial(x = pa_data,
                  column = "occ",
                  k = 5,
                  size = 350000,
                  hexagon = FALSE,
                  selection = "checkerboard",
                  iteration = 50)


## ----warning=FALSE, message=FALSE, fig.height=5, fig.width=7------------------
tm_shape(sb3$blocks) +
  tm_fill(col = "folds", style = "cat")


## -----------------------------------------------------------------------------
# environmental clustering
sc <- cv_cluster(x = pa_data,
                 column = "occ",
                 k = 5) # here number is low for test

## ----eval=FALSE, warning=FALSE, message=FALSE---------------------------------
#  # environmental clustering
#  ec <- cv_cluster(x = pa_data,
#                   column = "occ",
#                   r = rasters,
#                   k = 5,
#                   scale = TRUE,
#                   num_sample = 500) # here number is low for test

## ----eval=FALSE---------------------------------------------------------------
#  # buffering with presence-background data
#  bf2 <- cv_buffer(x = pa_data,
#                   column = "occ",
#                   size = 350000)
#  

## ----warning=FALSE, message=FALSE, fig.height=5, fig.width=7------------------
cv_plot(cv = sc,
        x = pa_data, 
        r = rasters) # optionally add a raster background


## -----------------------------------------------------------------------------
cv_similarity(cv = sb, x = pa_data, r = rasters, progress = FALSE)


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

