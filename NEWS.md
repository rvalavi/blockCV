# version 4.0.0

- Removed the legacy v2.x functions: `buffering`, `envBlock`, `foldExplorer`, `rangeExplorer`, `spatialAutoRange`, and `spatialBlock`.
- Renamed the `num_plot` argument of `cv_similarity` to `num_plots`, for consistency with `cv_plot` (a breaking change for existing code that set `num_plot`).
- Added `cv_knndm`, a k-fold Nearest Neighbour Distance Matching method with `"blocks"`, `"hierarchical"`, and `"kmeans"` clustering options. It supports geographical and feature-space matching via `space`, and accepts prediction points via `r`, `pred_points`, or `model_domain`; issue [#63](https://github.com/rvalavi/blockCV/issues/63).
- Added `cv_distance`, a diagnostic for comparing an existing `blockCV` fold design with the nearest-neighbour distance distribution expected in the prediction domain.
- Added `cv_group`, a leave-group-out cross-validation using an existing grouping column (e.g. site, plot, campaign, or individual) so records that share a group are never split across train and test. `k = NULL` gives one fold per group; a smaller `k` merges the groups into `k` folds, optionally balanced (`balance = TRUE`).
- Added `balance` to make fold balancing explicit in `cv_spatial`, `cv_cluster`, and `cv_knndm`. `cv_cluster(balance = FALSE)` keeps the previous single k-means behaviour; `balance = TRUE` uses the new `k_multiplier` argument to create candidate clusters for more even folds; issue [#34](https://github.com/rvalavi/blockCV/issues/34).
- Added presence-background balancing support so `cv_spatial`, `cv_cluster`, and `cv_knndm` can base fold balancing or matching on presences, preventing background points from dominating the split.
- Fixed `plot.cv_spatial_autocor()` so single-layer autocorrelation objects plot correctly and multi-layer plots include the stored map panel.
- Removed the `sp` dependency (a retiring package) from `cv_spatial_autocor`. Longlat inputs now pass an explicit great-circle cutoff to `automap::autofitVariogram` instead of being converted to a `sp` object, so `automap (>= 1.1-20)` is now required.
- Reworked the fold-balancing search behind `cv_spatial` (random selection) and now `cv_cluster(balance = TRUE)`. The objective is now scored on the test folds only and normalised per class (a Pearson chi-square against an equal split), and it minimises empty test cells for classes with at least `k` records before minimising imbalance. This replaces the previous rule (raise the smallest cell, then lower the standard deviation of all train and test cells at once), which could keep a split that leaves a class missing from a test fold and let an abundant class dominate the score. Because a different candidate assignment now wins, folds for a given `seed` may differ from earlier versions.
- Added support for continuous or count `column` values in `cv_spatial`, `cv_cluster`, `cv_buffer`, `cv_nndm`, and `cv_knndm`. Numeric columns are binned with `num_bins` (default `4`) for fold balancing and records reports; `num_bins = NULL` restores per-value handling.
- `cv_nndm` no longer requires `r` when `pred_points` or `model_domain` are supplied, and now falls back to Euclidean distances with a warning when `x` has no CRS; issue [#58](https://github.com/rvalavi/blockCV/issues/58).
- Fixed a train/test leakage bug in `cv_nndm(presence_bg = TRUE)` where the test fold used the presence's position instead of its row index, holding out the wrong point and leaving it in both train and test.
- `cv_similarity` now reports extrapolation with a per-fold summary attached as `attr(p, "extrapolation")`, the overall novelty rate in the subtitle, and shading for the novel region. A new `type = "map"` plots the sample points in space, colouring each test point by its similarity to show where extrapolation occurs.
- `cv_similarity` now computes presence-background similarity from presences only, matching `cv_distance`, and adds `seed` for reproducible random baselines. L1/L2 distance-based scores now handle single-layer rasters, single-point folds, missing covariates, and constant predictors more robustly.
- Reports, plots, and progress bars now run by default only in interactive sessions for `cv_spatial`, `cv_buffer`, `cv_cluster`, and `cv_nndm`; issue [#57](https://github.com/rvalavi/blockCV/issues/57).
- `cv_spatial` avoids duplicate block-to-point intersection work when subsetting and assigning blocks; issue [#40](https://github.com/rvalavi/blockCV/issues/40).
- Added `combine_folds` and `fold_colors` to `cv_plot` for single-map fold plots in `cv_spatial`, `cv_cluster`, and `cv_knndm`.
- Added `bg_alpha` to `cv_plot` to control the opacity of background points in presence-background fold plots.
- Improved point-based `plot()` methods so they accept the original sample data and report missing data more clearly; issues [#59](https://github.com/rvalavi/blockCV/issues/59) and [#60](https://github.com/rvalavi/blockCV/issues/60).
- Added `species_pb.csv`, a presence-background example dataset used in the `cv_buffer()` and `cv_cluster()` examples.
- Added a third tutorial showing how to use `blockCV` folds with `caret::trainControl()` via `index` and `indexOut`; issue [#48](https://github.com/rvalavi/blockCV/issues/48).

# version 3.2.0

- Added two new distance-based similarity scores, L1 and L2, to the `cv_similarity` function
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
