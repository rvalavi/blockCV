# version 2.1.0
* snowfall package for parallel processing is replaces by future.apply package; #7
* future.apply, shiny, shinydashboard, geosphere and ggplot2 packages moved to SUGGESTION packages. These are not install by default and the user is asked if needed. #7
* RStoolbox is no longer used for clustering
* an argument is added to envBlock function for sampling from raster layers
* no dependency on sp package any more
* doParallel = FALSE by default in spatialAutoRange function

# version 2.0.1
* print() and cat() are removed from the functions and verbose argument is added instead

# version 2.0.0
* most of the underlying functions are migrated to `sf` functions;
* the parallel processing changed from `foreach` to `snowfall`;
* the `species` argument in `spatialBlock` function accepts multi-class responses to find evenly distributed records in train and test folds;

# version 1.1.0
* change `spatialAutoRange` function to accepts rasters with low number of pixels; #2
* the `maskBySpecies = FALSE` in `spatialBlock` function is no longer supported;
* the `numLimit` argument in `spatialBlock` function is only accepts numeric values, and 0 means searching for evenly distributed folds;

# version 1.0.1
* add `speciesData` to `spatialAutoRange` function;
