# version 3.3.0

- Added a `balance` argument to make fold balancing explicit and consistent across functions. When `balance = TRUE`, the folds are chosen to even out the training/testing records (or the classes/bins of `column` when supplied); when `balance = FALSE`, `column` is used only for the records report and does not affect the folds. `cv_cluster` gains this option (with `balance = FALSE` keeping the previous single k-means behaviour): setting `balance = TRUE` over-clusters the points into `k * k_multiplier` groups (new `k_multiplier` argument, default `3`) and assigns them to folds over `iteration` random attempts, reusing the same balancing search as `cv_spatial`. This provides a clustering-based strategy for keeping a comparable number of presences/classes in each fold, avoiding empty partitions; issue [#34](https://github.com/rvalavi/blockCV/issues/34). `cv_spatial` (`balance = TRUE` by default, reproducing the previous random-selection behaviour) and `cv_knndm` (`balance = TRUE` by default) also expose the argument.
- Added `cv_knndm`, a k-fold implementation of Nearest Neighbour Distance Matching (Linnenbrink et al., 2024). The `"hierarchical"` and `"kmeans"` methods follow the paper (Kolmogorov-Smirnov random-CV gate, PC1-based merging of clusters into folds, and a `maxp` constraint), while the default `clustering = "blocks"` is a `blockCV`-specific variant that builds the groups from hexagonal (or square) spatial blocks. It supports geographical and feature-space matching (`space`), and accepts prediction points via `r`, `pred_points`, or `model_domain`; issue [#63](https://github.com/rvalavi/blockCV/issues/63).
- Added `cv_distance`, a diagnostic that works on any `blockCV` cross-validation object (`cv_spatial`, `cv_cluster`, `cv_buffer`, `cv_nndm`, or `cv_knndm`). It draws the nearest-neighbour distance distribution plot produced internally by `cv_nndm`/`cv_knndm` for a fold configuration that has already been generated, so the test-to-train distances (`CV` curve) can be compared against the prediction-to-sample distances (`Prediction` curve) and a nearest-neighbour/leave-one-out reference (`LOO` curve). With `add_random = TRUE` it overlays a random `k`-fold baseline (mean over `n_random` splits with a 10-90% band), and the Wasserstein-1 distance of each curve to the prediction distribution is reported and attached to the returned plot as `attr(p, "W")`. It supports geographical and feature-space matching (`space`), and accepts prediction points via `r`, `pred_points`, or `model_domain`.
- `cv_nndm` now accepts prediction points via `pred_points` or `model_domain` (in addition to `r`, which is no longer mandatory), and falls back to Euclidean distances with a warning when the sample data (`x`) has an undefined CRS, so it can be used with simulated or non-georeferenced grids; issue [#58](https://github.com/rvalavi/blockCV/issues/58).
- `column` now supports continuous (or count) responses: numeric columns are split into quantile bins via the new `n_bins` argument (default `4`) so that folds can be balanced and reported across the range of the response, rather than only across discrete classes. In `cv_spatial`, random selection uses these bins in its iterative balance search. Binary and low-cardinality categorical columns are still treated as classes, tied quantile breaks fall back to fewer bins, and `n_bins = NULL` restores the previous per-value behaviour. Available in `cv_spatial`, `cv_cluster`, `cv_buffer`, `cv_nndm`, and `cv_knndm` (where it also feeds the class-balance validity gate).
- Default reports, plots, and progress bars now run only in interactive sessions for `cv_spatial`, `cv_buffer`, `cv_cluster`, and `cv_nndm` where those output controls are available; issue [#57](https://github.com/rvalavi/blockCV/issues/57).
- Added `combine_folds` argument to `cv_plot` to show all folds in a single map with points coloured by their fold ID (with `fold_colors` to set the palette), as an alternative to the faceted train/test plots. Available for `cv_spatial`, `cv_cluster`, and `cv_knndm`.
- Improved plotting for point-based CV objects: `plot()` methods now accept the original sample data and call `cv_plot()` internally, while missing sample data produce a clearer error message; issues [#59](https://github.com/rvalavi/blockCV/issues/59) and [#60](https://github.com/rvalavi/blockCV/issues/60).
- Added a third tutorial showing how to use `blockCV` folds with `caret::trainControl()` via `index` and `indexOut`; issue [#48](https://github.com/rvalavi/blockCV/issues/48).

# version 3.2.0

- Added two new methods, L1 and L2 distances, to the `cv_similarity` function
- Fixed a warning in `cv_similarity` for colour aesthetics with ggplot
- Fixed the summary method and plotting for `cv_spatial_autocor`

# version 3.1.7

- Temporarily added `sp` package dependency to avoid CRAN error as required by `automap` package [#55].

# version 3.1.6

- Resolved unclear error messages; issue [#52](https://github.com/rvalavi/blockCV/issues/52) by A. Márcia Barbosa
- Resolved ggplot testing failure; issue [#54](https://github.com/rvalavi/blockCV/issues/54) by Teun van den Brand

# version 3.1.5

- Resolved background pattern artefacts in raster plotting; issue [#50](https://github.com/rvalavi/blockCV/issues/50) by Camila Neder.

# version 3.1.4

- Just the `biomod2` example is updated in vignettes; and the link in help file

# version 3.1.3

- the `biomod2` and `gstat` packages are added to the Suggests section
- Some minor edits in messages

# version 3.1.2

- The *iteration* in `cv_spatial` and `spatialBlock` is increased to 100 to make the result matches with v2.1.4
- Removed the requirement of C++11
- Some warnings are added for the miss use of the `column` argument

# version 3.1.1

- some internal fix.
- the `extend` parameter is now added to `spatialBlock` and the function now uses `cv_spatial` internally.
- the `user_blocks` in `cv_spatial` is restricted to random and predefined and systematic selection.
- no `raster` package dependency

# version 3.1.0

- the result of the `cv_spatial` function for square blocks now matches the one of version 2 function `spatialBlock` (i.e. fold assignment starts from top-right corner; this is not the case for hexagon blocks)
- square spatial blocks can be expanded to ensure no points fall outside the border of the blocks. This can be controlled by `extend` parameter now.

# version 3.0.3

- fixing a bug in counting records in the reporting of `cv_spatial`

# version 3.0.2

- fixing fold numbering of `cv_spatial` to reproducibility of earlier versions

# version 3.0.1

- Massive performance improvement in the C++ code of `cv_nndm` function for large datasets

# version 3.0

- Dependency to `rgdal` and `rgeos` are removed, and overall less dependency
- Function names have been changed, with all functions now starting with `cv_`
- The old functions (v2.x) still work to allow appropriate time for adapting the new code
- The CV blocking functions are now: `cv_spatial`, `cv_cluster`, `cv_buffer`, and `cv_nndm`
- Spatial blocks now support **hexagonal** (default), rectangular, and user-defined blocks
- A fast C++ implementation of **Nearest Neighbour Distance Matching (NNDM)** algorithm (Milà et al. 2022) is now added
- The NNDM algorithm can handle species presence-background data and other types of data
- The `cv_cluster` function generates blocks based on kmeans clustering. It now works on both environmental rasters and the **spatial coordinates of sample points**
- The `cv_spatial_autocor` function now calculates the spatial autocorrelation range for both the **response (i.e. the binary or continuous data)** and a set of continuous raster covariates
- The new `cv_plot` function allows for visualization of folds from all blocking strategies using ggplot facets
- The `terra` package is now used for all raster processing and supports both `stars` and `raster` objects, as well as files on disk.
- The new `cv_similarity` provides measures on possible extrapolation to testing folds

# version 2.1.4

- fixed CRAN error for ggplot guide
- added rgdal as a suggest
- changed the crs of raster data in the package to avoid datum warnings

# version 2.1.3

- fix the warning for spatialBlock function on geographic coordinate system

# version 2.1.2

- predefined folds from user-defined blocks are noe accepted
- add seed argument to spatialBlock to have consistent results where needed

# version 2.1.0

- snowfall package for parallel processing is replaces by future.apply package; #7
- future.apply, shiny, shinydashboard, geosphere and ggplot2 packages moved to SUGGESTION packages. These are not install by default and the user is asked if needed. #7
- RStoolbox is no longer used for clustering
- an argument is added to envBlock function for sampling from raster layers
- no dependency on sp package any more
- doParallel = FALSE by default in spatialAutoRange function

# version 2.0.1

- print() and cat() are removed from the functions and verbose argument is added instead

# version 2.0.0

- most of the underlying functions are migrated to `sf` functions;
- the parallel processing changed from `foreach` to `snowfall`;
- the `species` argument in `spatialBlock` function accepts multi-class responses to find evenly distributed records in train and test folds;

# version 1.1.0

- change `spatialAutoRange` function to accepts rasters with low number of pixels; #2
- the `maskBySpecies = FALSE` in `spatialBlock` function is no longer supported;
- the `numLimit` argument in `spatialBlock` function is only accepts numeric values, and 0 means searching for evenly distributed folds;

# version 1.0.1

- add `speciesData` to `spatialAutoRange` function;
