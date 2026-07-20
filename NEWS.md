# version 4.0-0

## Added

- Added `cv_knndm`, a k-fold Nearest Neighbour Distance Matching method with `"blocks"`, `"hierarchical"`, and `"kmeans"` clustering options. It supports geographical and feature-space matching via `space`, and accepts prediction points via `r`, `pred_points`, or `model_domain`; issue [#63](https://github.com/rvalavi/blockCV/issues/63).
- Added `cv_distance`, a diagnostic for comparing an existing `blockCV` fold design with the nearest-neighbour distance distribution expected in the prediction domain. It returns a `cv_distance` object carrying the per-fold test-to-nearest-train distance table in `$distances`, the Wasserstein-1 distances in `$W`, and the plot in `$plot`; a `plot` argument (default `TRUE`) controls only whether the plot is drawn.
- Added `cv_summary`, a one-call fold-quality diagnostic that gathers fold sizes and class prevalence, the `cv_distance` and `cv_similarity` diagnostics (when a raster or prediction domain is supplied), and a data.frame of warnings about degenerate folds (empty or single-class test folds, classes missing from training, tiny folds, severe class imbalance, and high leakage).
- Added `cv_group`, a leave-group-out cross-validation using an existing grouping column (e.g. site, plot, campaign, or individual) so records that share a group are never split across train and test. `k = NULL` gives one fold per group; a smaller `k` merges the groups into `k` folds, optionally balanced (`balance = TRUE`).
- Added `balance` to make fold balancing explicit in `cv_spatial`, `cv_cluster`, and `cv_knndm`. `cv_cluster(balance = FALSE)` keeps the previous single k-means behaviour; `balance = TRUE` uses the new `k_multiplier` argument to create candidate clusters for more even folds; issue [#34](https://github.com/rvalavi/blockCV/issues/34).
- Added presence-background balancing support so `cv_spatial`, `cv_cluster`, and `cv_knndm` can base fold balancing or matching on presences, preventing background points from dominating the split.
- Added `spatial_weight` to `cv_cluster` for spatially-constrained environmental clustering. It adds a soft spatial compactness pressure by blending the (standardised) coordinates into the covariates before k-means, so the Euclidean distance becomes `(1 - w) * d_env^2 + w * d_geo^2`. `spatial_weight = 0` (default) keeps the previous pure-environmental behaviour exactly, larger values give environmentally coherent folds that are also geographically separated, and `1` clusters on the coordinates alone.
- Added support for continuous or count `column` values in `cv_spatial`, `cv_cluster`, `cv_buffer`, `cv_nndm`, and `cv_knndm`. Numeric columns are binned with `num_bins` (default `4`) for fold balancing and records reports; `num_bins = NULL` restores per-value handling.
- Added `combine_folds` and `fold_colors` to `cv_plot` for single-map fold plots in `cv_spatial`, `cv_cluster`, `cv_group`, and `cv_knndm`.
- Added `bg_alpha` to `cv_plot` to control the opacity of background points in presence-background fold plots; it is capped at `points_alpha` (with a warning) so the background is never drawn more prominent than the presences, and the "background shown faded" caption appears only when the background is actually more transparent.
- Added `species_pb.csv`, a presence-background example dataset used in the `cv_buffer()` and `cv_cluster()` examples.
- Added a second tutorial covering `cv_summary`, `cv_similarity`, `cv_distance`, and block-size selection with `cv_spatial_autocor` and `cv_block_size`.
- Added a fourth tutorial showing how to use `blockCV` folds with `caret::trainControl()` via `index` and `indexOut`; issue [#48](https://github.com/rvalavi/blockCV/issues/48).

## Changed

- The minimum R version is now 3.6.0.
- Renamed the `num_plot` argument of `cv_similarity` to `num_plots`, for consistency with `cv_plot` (a breaking change for existing code that set `num_plot`).
- Removed the `sp` dependency (a retiring package) from `cv_spatial_autocor`. Longlat inputs now pass an explicit great-circle cutoff to `automap::autofitVariogram` instead of being converted to a `sp` object, so `automap (>= 1.1-20)` is now required.
- Reworked the fold-balancing search behind `cv_spatial` (random selection) and now `cv_cluster(balance = TRUE)`. The objective is now scored on the test folds only and normalised per class (a Pearson chi-square against an equal split), and it minimises empty test cells for classes with at least `k` records before minimising imbalance. This replaces the previous rule (raise the smallest cell, then lower the standard deviation of all train and test cells at once), which could keep a split that leaves a class missing from a test fold and let an abundant class dominate the score. Because a different candidate assignment now wins, folds for a given `seed` may differ from earlier versions.
- `cv_similarity` now computes presence-background similarity from presences only, matching `cv_distance`, and adds `seed` for reproducible random baselines. L1/L2 distance-based scores now handle single-layer rasters, single-point folds, missing covariates, and constant predictors more robustly, and the MESS score now guards constant (zero-range) predictors. `cv_similarity` also validates that the rows of `x` line up with the fold indices in `cv`, matching `cv_distance`.
- Reports, plots, and progress bars now run by default only in interactive sessions for `cv_spatial`, `cv_buffer`, `cv_cluster`, and `cv_nndm`; `cv_similarity` progress bars now do the same. Report flags now suppress printing only; fold record tables are still returned. Issue [#57](https://github.com/rvalavi/blockCV/issues/57).
- The `print` methods for the fold objects (`cv_spatial`, `cv_buffer`, `cv_cluster`, `cv_group`, `cv_knndm`, `cv_nndm`) and for `cv_spatial_autocor` now show an informative summary (method description, key settings, and the per-fold train/test record table) instead of only echoing the class name, matching the diagnostic `print` methods of `cv_distance`, `cv_similarity`, and `cv_summary`. Each returns the object invisibly. The `cv_cluster` print distinguishes spatial, environmental (feature-space), and spatially-constrained environmental clustering, and the returned object now carries the `spatial_weight` used. The `cv_spatial` object now carries the `block_shape` (`"hexagon"`, `"square"`, or `"user-defined"`) and the `selection` used, both shown in its print. The `cv_distance` and `cv_similarity` prints now note when the object is presence-background, making it clear the distance/similarity/novelty statistics are computed on presence points only (background samples excluded).
- Functions that need a suggested package now stop with a clear installation message instead of an interactive install prompt, so they behave consistently in scripts and non-interactive sessions.
- `cv_cluster` now rejects categorical (factor) layers in `r` with a clear error instead of silently clustering their integer codes, since environmental clustering uses k-means (Euclidean distance) and supports only numeric covariates. It also warns when spatial (coordinate-based) clustering runs on geographic (lon/lat) coordinates, since k-means uses Euclidean distance; consider projecting `x`. It now also validates `k` against the number of available sample points (or sampled raster cells) before calling `stats::kmeans`, stopping with a clearer message instead of the cryptic "more cluster centers than distinct data points" error.
- The `cv_*` functions now raise a clear error when `x` cannot be coerced to an `sf` object, instead of silently continuing with the invalid input and failing later with an unrelated message.
- `cv_nndm` no longer requires `r` when `pred_points` or `model_domain` are supplied, now falls back to Euclidean distances with a warning when `x` has no CRS, and errors when `x` and `r` have mismatched coordinate reference systems; issue [#58](https://github.com/rvalavi/blockCV/issues/58).
- `cv_nndm` now returns `exclusion`, a table with one row per fold giving the fold number, the row index of its test point in `x` (`test_id`), and the `exclusion_distance` matched to that point, so it is clear which record is held out at what distance. Unlike `cv_buffer`, whose exclusion distance is the constant `size`, NNDM matches a distinct radius to each test point; the radius is not capped by `size`, as points whose nearest neighbour already lies beyond `size` are left unthinned and keep their larger natural nearest-neighbour distance.
- `cv_similarity` now returns a `cv_similarity` object (a list with the per-fold extrapolation summary in `$extrapolation`, the overall novelty rate in `$overall`, and the plot in `$plot`) instead of a ggplot with the summary attached as an attribute. It gains a `plot` argument (default `TRUE`) that controls only whether the plot is drawn; the plot object is always built and returned. The overall novelty rate is shown in the title and the novel region is shaded, and a new `type = "map"` plots the sample points in space, colouring each test point by its similarity to show where extrapolation occurs.
- Improved the fold-object `plot()` methods, now including `cv_spatial`, so they accept the original sample data and report missing data more clearly; issues [#59](https://github.com/rvalavi/blockCV/issues/59) and [#60](https://github.com/rvalavi/blockCV/issues/60).
- `cv_spatial` avoids duplicate block-to-point intersection work when subsetting and assigning blocks; issue [#40](https://github.com/rvalavi/blockCV/issues/40).

## Removed

- Removed the legacy v2.x functions: `buffering`, `envBlock`, `foldExplorer`, `rangeExplorer`, `spatialAutoRange`, and `spatialBlock`; use the `cv_*` equivalents in new code.

## Fixed

- Fixed a train/test leakage bug in `cv_nndm(presence_bg = TRUE)` where the test fold used the presence's position instead of its row index, holding out the wrong point and leaving it in both train and test.
- Fixed `plot.cv_spatial_autocor()` so single-layer autocorrelation objects plot correctly and multi-layer plots include the stored map panel.
- `cv_plot` and `cv_similarity` now stop with a clear message when none of the requested `num_plots` folds exist (e.g. `num_plots = 99` for a 5-fold object), instead of silently reducing the selection to nothing and failing later with a confusing "undefined columns selected" error; out-of-range values are still dropped when at least one requested fold is valid.

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
