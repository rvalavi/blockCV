# blockCV <img src="man/figures/logo.png" align="right" width="120"/>

[![R build
status](https://github.com/rvalavi/blockCV/workflows/R-CMD-check/badge.svg)](https://github.com/rvalavi/blockCV/actions)
[![codecov](https://codecov.io/gh/rvalavi/blockCV/branch/master/graph/badge.svg)](https://codecov.io/gh/rvalavi/blockCV)
[![CRAN
version](https://www.r-pkg.org/badges/version/blockCV)](https://CRAN.R-project.org/package=blockCV)
[![total](https://cranlogs.r-pkg.org/badges/grand-total/blockCV)](https://CRAN.R-project.org/package=blockCV)
[![License](https://img.shields.io/badge/license-GPL%20(%3E=%203)-lightgrey.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![Methods in Ecology & Evolution](https://img.shields.io/badge/Methods%20in%20Ecology%20%26%20Evolution-10,%20225--232-blue.svg)](https://doi.org/10.1111/2041-210X.13107)

### Spatially and environmentally separated folds for cross-validation

The `blockCV` package creates spatially or environmentally separated
training and testing folds for **k-fold**, **leave-group-out**, and
**leave-one-out (LOO)** cross-validation. These folds support more
realistic evaluation of models fitted to spatially structured data,
including remote-sensing classification, soil mapping, and species
distribution modelling.

Alongside several fold-construction strategies, `blockCV` provides tools
for checking fold balance, comparing fold separation with the prediction
domain, and identifying environmental extrapolation. It can also estimate
spatial autocorrelation ranges in point data or continuous raster
covariates, providing an initial distance scale to investigate when
designing spatial folds.
 
## Main features

-   Six fold-construction strategies: spatial blocks (`cv_spatial`),
    spatial or environmental clustering (`cv_cluster`), existing grouping
    factors (`cv_group`), buffering (`cv_buffer`), leave-one-out nearest
    neighbour distance matching (`cv_nndm`), and k-fold nearest neighbour
    distance matching (`cv_knndm`)
-   Hexagonal (default), rectangular, or user-defined spatial blocks,
    assigned to folds using random, systematic, checkerboard, or predefined
    selection. Random assignment can search for balanced folds
-   Environmental clustering with optional spatial constraints through
    `spatial_weight`, and optional over-clustering to improve fold balance
-   Geographical or feature-space kNNDM using prediction locations supplied
    by a raster (`r`), prediction points (`pred_points`), or a polygon
    (`model_domain`), with block, hierarchical, or k-means grouping
-   Response-aware fold summaries and, where supported, balancing for
    binary, multi-class, continuous, count, and presence-background data.
    Continuous and count responses can be grouped into quantile bins using
    `num_bins`
-   Fold diagnostics through `cv_summary`, `cv_distance`, and
    `cv_similarity` to assess balance, train-test separation, agreement
    with prediction-domain distances, and environmental extrapolation
-   Fold visualisation with `cv_plot`, including faceted train-test maps
    for every strategy and combined-fold maps for k-fold methods
-   Spatial autocorrelation and interactive block-size tools
    (`cv_spatial_autocor` and `cv_block_size`) for exploring an initial
    separation distance
-   Raster processing with `terra`, including support for `stars`,
    `raster`, and raster files on disk

## What's new in v4.0

-   Added `cv_knndm` for k-fold nearest neighbour distance matching, with
    geographical and feature-space matching and block, hierarchical, or
    k-means grouping
-   Added `cv_group` for leave-group-out cross-validation based on an
    existing site, plot, campaign, individual, or other grouping factor
-   Added `cv_summary` for one-call fold-quality summaries and warnings,
    and `cv_distance` for comparing fold separation with nearest-neighbour
    distances in the prediction domain
-   Expanded `cv_similarity` with per-fold extrapolation summaries,
    overall novelty rates, and spatial map visualisation
-   Made fold balancing explicit in `cv_spatial`, `cv_cluster`, and
    `cv_knndm`; added presence-background balancing and quantile binning
    (`num_bins`) for continuous or count responses
-   Added spatially constrained environmental clustering through
    `cv_cluster(spatial_weight = ...)`
-   Added combined-fold maps to `cv_plot` and informative print methods
    for fold and diagnostic objects
-   Expanded `cv_nndm` and `cv_knndm` to accept prediction rasters,
    prediction points, or model-domain polygons
-   Removed the legacy v2.x function names. Other breaking changes include
    the new structured return value from `cv_similarity`, the rename of
    `num_plot` to `num_plots`, and interactive-only defaults for several
    automatic plots, reports, and progress bars

See [NEWS.md](NEWS.md) for the full changelog.

## Installation

To install the latest update of the package from GitHub use:

``` r
remotes::install_github("rvalavi/blockCV", build_vignettes = TRUE, dependencies = TRUE)
```

Or installing from CRAN:

``` r
install.packages("blockCV", dependencies = TRUE)
```

## Vignettes

To see the practical examples of the package see:

1.  [blockCV introduction: how to create block cross-validation
    folds](https://htmlpreview.github.io/?https://github.com/rvalavi/blockCV/blob/master/inst/doc/tutorial_1.html)
2.  [Choosing and diagnosing spatial
    folds](https://htmlpreview.github.io/?https://github.com/rvalavi/blockCV/blob/master/inst/doc/tutorial_2.html)
3.  [Block cross-validation for species distribution
    modelling](https://htmlpreview.github.io/?https://github.com/rvalavi/blockCV/blob/master/inst/doc/tutorial_3.html)
4.  [Using blockCV with
    `caret`](https://htmlpreview.github.io/?https://github.com/rvalavi/blockCV/blob/master/inst/doc/tutorial_4.html)

## Basic usage

This code snippet showcases some of the package's functionalities, but for more comprehensive tutorials, please refer to the vignette included with the package (and above).

``` r
# loading the package
library(blockCV)
library(sf) # working with spatial vector data
library(terra) # working with spatial raster data
```

``` r
# load raster data; the pipe operator |> is available in R v4.1 or higher
myrasters <- system.file("extdata/au/", package = "blockCV") |>
  list.files(full.names = TRUE) |>
  terra::rast()

# load species presence-absence data and convert to sf
pa_data <- read.csv(system.file("extdata/", "species.csv", package = "blockCV")) |>
  sf::st_as_sf(coords = c("x", "y"), crs = 7845)

```


``` r
# spatial blocking by specified range and random assignment
sb <- cv_spatial(
    x = pa_data, # sf or SpatialPoints of sample data (e.g. species data)
    column = "occ", # optional response column for fold records/balancing
    r = myrasters, # a raster for background (optional)
    size = 450000, # size of the blocks in metres
    k = 5, # number of folds
    hexagon = TRUE, # use hexagonal blocks - default
    selection = "random", # random blocks-to-fold
    iteration = 100, # search for balanced folds
    biomod2 = TRUE # also create folds for biomod2
)
```

![](https://i.ibb.co/WGfrF7B/Rplot1.png)

Or create spatial clusters for k-fold cross-validation:

``` r
# create spatial clusters
set.seed(6)
sc <- cv_cluster(
    x = pa_data, 
    column = "occ", # optionally count classes/bins in folds
    k = 5
)
```

``` r
# now plot the created folds
cv_plot(
    cv = sc, # a blockCV object
    x = pa_data, # sample points
    r = myrasters[[1]], # optionally add a raster background
    points_alpha = 0.5,
    nrow = 2
)
```

![](https://i.ibb.co/dGrF9xp/Rplot02.png)

Create k-fold NNDM folds:

``` r
# k-fold nearest neighbour distance matching
knn <- cv_knndm(
    x = pa_data,
    column = "occ", # optionally prefer class-complete folds
    r = myrasters[[1]], # prediction area, or use pred_points/model_domain
    k = 5,
    num_sample = 5000
)
```

``` r
# compare an existing fold design with prediction-domain distances
cv_distance(
    cv = sb,
    x = pa_data,
    r = myrasters[[1]],
    num_sample = 5000
)
```

Investigate spatial autocorrelation in the landscape to choose a
suitable size for spatial blocks:

``` r
# exploring the effective range of spatial autocorrelation in raster covariates or sample data
cv_spatial_autocor(
    r = myrasters, # a SpatRaster object or path to files
    num_sample = 5000, # number of cells to be used
    plot = TRUE
)
```

For the residual-based block-size guidance in Roberts et al. (2017), fit
the model first, add its residuals to the sample points, and pass that
residual column to `cv_spatial_autocor(x = ..., column = ...)`. Ranges
estimated from the raw response or raster covariates are exploratory
proxies and may mis-size blocks for residual autocorrelation.

Alternatively, you can manually choose the size of spatial blocks in an
interactive session using a Shiny app.

``` r
# a shiny interactive app to aid selecting a size for spatial blocks
cv_block_size(
    r = myrasters[[1]],
    x = pa_data, # optionally add sample points
    column = "occ",
    min_size = 2e5,
    max_size = 9e5
)
```

## Reporting issues

Please report issues at: <https://github.com/rvalavi/blockCV/issues>

## Citation

To cite package **blockCV** in publications, please use:

Valavi R, Elith J, Lahoz-Monfort JJ, Guillera-Arroita G. **blockCV: An R
package for generating spatially or environmentally separated folds for
k-fold cross-validation of species distribution models**. *Methods Ecol
Evol*. 2019; 10:225--232. <https://doi.org/10.1111/2041-210X.13107>
