# version 3.1.4
* Just the `biomod2` example is updated in vignettes; and the link in help file

# version 3.1.3
* the `biomod2` and `gstat` packages are added to the Suggests section
* Some minor edits in messages

# version 3.1.2
* The *iteration* in `cv_spatial` and `spatialBlock` is increased to 100 to make the result matches with v2.1.4
* Removed the requirement of C++11
* Some warnings are added for the miss use of the `column` argument

# version 3.1.1
* some internal fix.
* the `extend` parameter is now added to `spatialBlock` and the function now uses `cv_spatial` internally.
* the `user_blocks` in `cv_spatial` is restricted to random and predefined and systematic selection.
* no `raster` package dependency

# version 3.1.0
* the result of the `cv_spatial` function for square blocks now matches the one of version 2 function `spatialBlock` (i.e. fold assignment starts from top-right corner; this is not the case for hexagon blocks)
* square spatial blocks can be expanded to ensure no points fall outside the border of the blocks. This can be controlled by `extend` parameter now.

# version 3.0.3
* fixing a bug in counting records in the reporting of `cv_spatial`

# version 3.0.2
* fixing fold numbering of `cv_spatial` to reproducibility of earlier versions

# version 3.0.1
* Massive performance improvement in the C++ code of `cv_nndm` function for large datasets

# version 3.0
* Dependency to rgdal and rgeos are removed, and overall less dependency
* Function names have been changed, with all functions now starting with `cv_`
* The old functions (v2.x) still work to allow appropriate time for adapting the new code
* The CV blocking functions are now: `cv_spatial`, `cv_cluster`, `cv_buffer`, and `cv_nndm`
* Spatial blocks now support **hexagonal** (default), rectangular, and user-defined blocks
* A fast C++ implementation of **Nearest Neighbour Distance Matching (NNDM)** algorithm (Milà et al. 2022) is now added
* The NNDM algorithm can handle species presence-background data and other types of data
* The `cv_cluster` function generates blocks based on kmeans clustering. It now works on both environmental rasters and the **spatial coordinates of sample points**
* The `cv_spatial_autocor` function now calculates the spatial autocorrelation range for both the **response (i.e. the binary or continuous data)** and a set of continuous raster covariates
* The new `cv_plot` function allows for visualization of folds from all blocking strategies using ggplot facets
* The `terra` package is now used for all raster processing and supports both `stars` and `raster` objects, as well as files on disk.
* The new `cv_similarity` provides measures on possible extrapolation to testing folds

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
