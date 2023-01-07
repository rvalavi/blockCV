# version 3.0
* function names have changed, all starting with `cv_*` now
* all old functions (v2.x) are deprecate now and replaces with the new function; although they should work fine
* better error handling is utilised
* dependencies on `raster`, `rgdal`, `rgeos`, `progress`, `future`, and `future.apply` are removed
* the `terra` package is now used for all raster processing with support for `stars` and `raster` formats
* several dependencies are removed and `blockCV` only imports `sf` package now
* spatial blocks now support hexagonal (default), rectangular, and user-defined blocks
* clustering function now works both on *environmental* rasters, and *spatial coordinates* of the sample points
* the `cv_spatial_autocor` function now calculates spatial autocorrelation range for either the response (i.e. the binary or continuous data) or a set of continuous raster covariates (as before)
* the new `cv_plot` function can be used to plot the folds of all blocking strategy with ggplot facets
* the newly developed function ******** is implemented
* ...


# version 2.1.4
* fixed CRAN error for ggplot guide 
* added rgdal as a suggest
* changed the crs of raster data in the package to avoid datum warnings

# version 2.1.3
* fix the warning for spatialBlock function on geographic coordinate system

# version 2.1.2
* predefined folds from user-defined blocks are noe accepted
* add seed argument to spatialBlock to have consistent results where needed

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
