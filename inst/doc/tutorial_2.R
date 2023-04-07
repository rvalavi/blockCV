## ----echo=FALSE---------------------------------------------------------------
options(scipen = 10)

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


## ----message=TRUE, warning=TRUE-----------------------------------------------
library(blockCV)


## ---- fig.keep='all', warning=FALSE, message=FALSE, fig.height=5, fig.width=7----
scv1 <- cv_spatial(
  x = pa_data,
  column = "occ", # the response column (binary or multi-class)
  r = rasters,
  k = 5, # number of folds
  size = 360000, # size of the blocks in metres
  selection = "random", # random blocks-to-fold
  iteration = 50, # find evenly dispersed folds
  progress = FALSE, # turn off progress bar
  biomod2 = TRUE, # also create folds for biomod2
  raster_colors = terrain.colors(10, rev = TRUE) # options from cv_plot for a better colour contrast
) 


## -----------------------------------------------------------------------------
range <- cv_spatial_autocor(
  x = pa_data, # species data
  column = "occ", # column storing presence-absence records (0s and 1s)
  plot = FALSE
)

range$range

## ----fig.height=5, fig.width=7------------------------------------------------
scv2 <- cv_nndm(
  x = pa_data,
  column = "occ",
  r = rasters,
  size = 360000, # range of spatial autocorrelation
  num_sample = 10000, # number of samples of prediction points
  sampling = "regular", # sampling methods; it can be random as well
  min_train = 0.1, # minimum portion to keep in each train fold
  plot = TRUE
)

## -----------------------------------------------------------------------------
# see the number of folds in scv2 object
scv2$k


## ----fig.height=5, fig.width=8------------------------------------------------
cv_plot(
  cv = scv2, # cv object
  x = pa_data, # species spatial data
  num_plots = c(1, 10, 100) # three of folds to plot
)

## ---- eval=FALSE--------------------------------------------------------------
#  # loading the libraries
#  library(randomForest)
#  library(precrec)
#  
#  # extract the raster values for the species points as a dataframe
#  model_data <- terra::extract(rasters, pa_data, df = TRUE, ID = FALSE)
#  # adding species column to the dataframe
#  model_data$occ <- as.factor(pa_data$occ)
#  head(model_data)
#  
#  # extract the fold indices from buffering object
#  # created in the previous section
#  # the folds_list works for all three blocking strategies
#  folds <- scv2$folds_list
#  
#  # create a data.frame to store the prediction of each fold (record)
#  test_table <- pa_data
#  test_table$preds <- NA
#  
#  for(k in seq_len(length(folds))){
#    # extracting the training and testing indices
#    # this way works with folds_list list (but not folds_ids)
#    trainSet <- unlist(folds[[k]][1]) # training set indices; first element
#    testSet <- unlist(folds[[k]][2]) # testing set indices; second element
#    rf <- randomForest(occ ~ ., model_data[trainSet, ], ntree = 500) # model fitting on training set
#    test_table$preds[testSet] <- predict(rf, model_data[testSet, ], type = "prob")[,2] # predict the test set
#  }
#  
#  # calculate Area Under the ROC and PR curves and plot the result
#  precrec_obj <- evalmod(scores = test_table$preds, labels = test_table$occ)
#  auc(precrec_obj)
#  

## ----echo=FALSE---------------------------------------------------------------
# to not run the model and reduce run time; result are calculated and loaded
read.csv("../man/figures/roc_rf.csv") 


## ---- eval=FALSE, fig.height=3.7, fig.width=7---------------------------------
#  library(ggplot2)
#  
#  autoplot(precrec_obj)
#  

## ----eval=FALSE---------------------------------------------------------------
#  # loading the library
#  library(biomod2)
#  
#  # extract the raster values for the species points as a dataframe
#  raster_values <- terra::extract(rasters, pa_data, df = TRUE, ID = FALSE)
#  
#  # 1. Formatting Data
#  biomod_data <- BIOMOD_FormatingData(resp.var = pa_data$occ,
#                                      expl.var = raster_values,
#                                      resp.xy = sf::st_coordinates(pa_data),
#                                      resp.name = "occ",
#                                      na.rm = TRUE)
#  
#  # 2. Defining the folds for data.split.table
#  # note that biomod_table should be used here not folds
#  # use generated folds from cv_spatial in previous section
#  spatial_cv_folds <- scv1$biomod_table
#  
#  # 3. Defining Models Options; using default options here.
#  biomod_options <- BIOMOD_ModelingOptions()
#  
#  # 4. Model fitting
#  biomod_model_out <- BIOMOD_Modeling(biomod_data,
#                                      models = c('GLM','MARS','GBM'),
#                                      bm.options = biomod_options,
#                                      data.split.table = spatial_cv_folds,
#                                      var.import = 0,
#                                      metric.eval = c('ROC'),
#                                      do.full.models = TRUE)
#  

## ---- eval=FALSE--------------------------------------------------------------
#  # 5. Model evaluation
#  biomod_model_eval <- get_evaluations(biomod_model_out)
#  biomod_model_eval[c("run", "algo", "metric.eval", "calibration", "validation")]
#  

## ----echo=FALSE---------------------------------------------------------------
# to not run the model and reduce run time; result are calculated and loaded
read.csv("../man/figures/evl_biomod.csv") 


