# blockCV <img src="man/figures/logo.png" align="right" width="120"/>

[![R build
status](https://github.com/rvalavi/blockCV/workflows/R-CMD-check/badge.svg)](https://github.com/rvalavi/blockCV/actions)
[![codecov](https://codecov.io/gh/rvalavi/blockCV/branch/master/graph/badge.svg)](https://codecov.io/gh/rvalavi/blockCV)
[![CRAN
version](https://www.r-pkg.org/badges/version/blockCV)](https://CRAN.R-project.org/package=blockCV)
[![total](https://cranlogs.r-pkg.org/badges/grand-total/blockCV)](https://CRAN.R-project.org/package=blockCV)
[![License](https://img.shields.io/badge/license-GPL%20(%3E=%203)-lightgrey.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![DOI](https://zenodo.org/badge/116337503.svg)](https://zenodo.org/badge/latestdoi/116337503)
[![Methods in Ecology & Evolution](https://img.shields.io/badge/Methods%20in%20Ecology%20&%20Evolution-10,%20225-232-blue.svg)](https://doi.org/10.1111/2041-210X.13107)

### Spatial and environmental blocking for k-fold and LOO cross-validation

The package `blockCV` offers a range of functions for generating train
and test folds for **k-fold** and **leave-one-out (LOO)**
cross-validation (CV). It allows for separation of data spatially and
environmentally, with various options for block construction.
Additionally, it includes a function for assessing the level of spatial
autocorrelation in response or raster covariates, to aid in selecting an
appropriate distance band for data separation. The `blockCV` package is
suitable for the evaluation of a variety of spatial modelling
applications, including classification of remote sensing imagery, soil
mapping, and species distribution modelling (SDM). It also provides
support for different SDM scenarios, including presence-absence and
presence-background species data, rare and common species, and raster
data for predictor variables.

## Main features

-   There are four blocking methods: **spatial**, **clustering**,
    **buffers**, and **NNDM** (Nearest Neighbour Distance Matching)
    blocks
-   Several ways to construct spatial blocks
-   The assignment of the spatial blocks to cross-validation folds can
    be done in three different ways: **random**, **systematic** and
    **checkerboard pattern**
-   The spatial blocks can be assigned to cross-validation folds to have
    *evenly distributed records* for *binary* (e.g. species
    presence-absence/background) or *multi-class* responses (e.g. land
    cover classes for remote sensing image classification)
-   The buffering and NNDM functions can account for *presence-absence*
    and *presence-background* data types
-   Using geostatistical techniques to inform the choice of a suitable
    distance band by which to separate the data sets

## New updates of the version 3.0

The latest version `blockCV` (v3.0) features significant updates and changes. All function names have been revised to more general names, beginning with `cv_*`. Although the previous functions (version 2.x) will continue to work, they will be removed in future updates after being available for an extended period. It is highly recommended to update your code with the new functions provided below.

Some new updates:

-   Function names have been changed, with all functions now starting
    with `cv_`
-   The CV blocking functions are now: `cv_spatial`, `cv_cluster`,
    `cv_buffer`, and `cv_nndm`
-   Spatial blocks now support **hexagonal** (now, default),
    rectangular, and user-defined blocks
-   A fast C++ implementation of **Nearest Neighbour Distance Matching
    (NNDM)** algorithm (Milà et al. 2022) is now added
-   The NNDM algorithm can handle species presence-background data and
    other types of data
-   The `cv_cluster` function generates blocks based on kmeans
    clustering. It now works on both environmental rasters and the
    **spatial coordinates of sample points**
-   The `cv_spatial_autocor` function now calculates the spatial
    autocorrelation range for both the **response (i.e. binary or
    continuous data)** and a set of continuous raster covariates
-   The new `cv_plot` function allows for visualization of folds from
    all blocking strategies using ggplot facets
-   The `terra` package is now used for all raster processing and
    supports both `stars` and `raster` objects, as well as files on
    disk.
-   The new `cv_similarity` provides measures on possible extrapolation
    to testing folds

## Installation

To install the latest update of the package from GitHub use:

``` r
remotes::install_github("rvalavi/blockCV", dependencies = TRUE)
```

Or installing from CRAN:

``` r
install.packages("blockCV", dependencies = TRUE)
```

## Vignettes

To see the practical examples of the package see:

1.  [blockCV introduction: how to create block cross-validation
    folds](https://htmlpreview.github.io/?https://github.com/rvalavi/blockCV/blob/master/inst/doc/tutorial_1.html)
2.  [Block cross-validation for species distribution
    modelling](https://htmlpreview.github.io/?https://github.com/rvalavi/blockCV/blob/master/inst/doc/tutorial_2.html)
3.  Using blockCV with the `caret` and `tidymodels` (coming soon!)

## Basic usage

This code snippet showcases some of the package's functionalities, but for more comprehensive tutorials, please refer to the vignette included with the package (and above).

``` r
# loading the package
library(blockCV)
library(sf) # working with spatial vector data
library(terra) # working with spatial raster data
```

``` r
# load raster data; the pipe operator |> is available for R v4.1 or higher
myrasters <- system.file("extdata/au/", package = "blockCV") |>
  list.files(full.names = TRUE) |>
  terra::rast()

# load species presence-absence data and convert to sf
pa_data <- read.csv(system.file("extdata/", "species.csv", package = "blockCV")) |>
  sf::st_as_sf(coords = c("x", "y"), crs = 7845)

```


``` r
# spatial blocking by specified range and random assignment
sb <- cv_spatial(x = pa_data, # sf or SpatialPoints of sample data (e.g. species data)
                 column = "occ", # the response column (binary or multi-class)
                 r = myrasters, # a raster for background (optional)
                 size = 450000, # size of the blocks in metres
                 k = 5, # number of folds
                 hexagon = TRUE, # use hexagonal blocks - defualt
                 selection = "random", # random blocks-to-fold
                 iteration = 100, # to find evenly dispersed folds
                 biomod2 = TRUE) # also create folds for biomod2
```

![](https://i.ibb.co/WGfrF7B/Rplot1.png)

Or create spatial clusters for k-fold cross-validation:

``` r
# create spatial clusters
set.seed(6)
sc <- cv_cluster(x = pa_data, 
                 column = "occ", # optionally count data in folds (binary or multi-class)
                 k = 5)
```

``` r
# now plot the created folds
cv_plot(cv = sc, # a blockCV object
        x = pa_data, # sample points
        r = myrasters[[1]], # optionally add a raster background
        points_alpha = 0.5,
        nrow = 2)
```

![](https://i.ibb.co/dGrF9xp/Rplot02.png)

Investigate spatial autocorrelation in the landscape to choose a
suitable size for spatial blocks:

``` r
# exploring the effective range of spatial autocorrelation in raster covariates or sample data
cv_spatial_autocor(r = myrasters, # a SpatRaster object or path to files
                   num_sample = 5000, # number of cells to be used
                   plot = TRUE)
```

Alternatively, you can manually choose the size of spatial blocks in an
interactive session using a Shiny app.

``` r
# shiny app to aid selecting a size for spatial blocks
cv_block_size(r = myrasters[[1]],
              x = pa_data, # optionally add sample points
              column = "occ",
              min_size = 2e5,
              max_size = 9e5)
```

## Reporting issues

Please report issues at: <https://github.com/rvalavi/blockCV/issues>

## Citation

To cite package **blockCV** in publications, please use:

Valavi R, Elith J, Lahoz-Monfort JJ, Guillera-Arroita G. **blockCV: An R
package for generating spatially or environmentally separated folds for
k-fold cross-validation of species distribution models**. *Methods Ecol
Evol*. 2019; 10:225--232. <https://doi.org/10.1111/2041-210X.13107>
