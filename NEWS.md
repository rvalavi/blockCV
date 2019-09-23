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
