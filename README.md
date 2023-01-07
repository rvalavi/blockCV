# blockCV <img src="man/figures/logo.png" align="right" width="120" />

[![Build Status](https://travis-ci.org/rvalavi/blockCV.svg?branch=master)](https://travis-ci.org/rvalavi/blockCV)
[![codecov](https://codecov.io/gh/rvalavi/blockCV/branch/master/graph/badge.svg)](https://codecov.io/gh/rvalavi/blockCV)
[![CRAN version](https://www.r-pkg.org/badges/version/blockCV)](https://CRAN.R-project.org/package=blockCV)
[![total](http://cranlogs.r-pkg.org/badges/grand-total/blockCV)](https://www.rpackages.io/package/blockCV)
[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%203%29-lightgrey.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![DOI](https://zenodo.org/badge/116337503.svg)](https://zenodo.org/badge/latestdoi/116337503)


### Spatial and environmental blocking for k-fold and LOO cross-validation   
   
In a nutshell, package **blockCV** provides functions to build train and test data sets using three general strategies: *buffers*, *spatial* and *environmental* blocks. It offers several options for how those blocks are constructed and how blocks are allocated to cross-validation folds. It includes a function that applies geostatistical techniques to investigate the existing level of spatial autocorrelation in the chosen predictor variables or the response variable (e.g. binary or continuous responses) to inform the choice of the block and buffer size. In addition, visualization tools further aid the selection of block size and provide an understanding of the spread of sample points across generated folds. 


## Main features
* There are three blocking methods: **buffers**, **spatial** and **clustering** blocks
* Several ways to construct spatial blocks
* The assignment of the spatial blocks to cross-validation folds can be done in three different ways: **random**, **systematic** and **checkerboard pattern**
* The spatial blocks can be assigned to cross-validation folds to have *evenly distributed records* for *binary* (e.g. species presence-absence/background) or *multi-class* responses (e.g. land cover classes for remote sensing image classification) 
* The buffering function can account for *presence-absence* and *presence-background* data types 
* Using geostatistical techniques to inform the choice of a suitable distance band by which to separate the data sets 

## New updates of version 3.x
* Function names have changed, all starting with `cv_` now
* The three main blocking functions are now: `cv_spatial`, `cv_cluster`, and `cv_buffer`
* Spatial blocks now support hexagonal (default), rectangular, and user-defined blocks
* Clustering function now works both on *environmental* rasters, and *spatial coordinates* of the sample points
* The `cv_spatial_autocor` function now calculates spatial autocorrelation range for either the response (i.e. the binary or continuous data) or a set of continuous raster covariates (as before)
* The new `cv_plot` function can be used to plot the folds of all blocking strategy with ggplot facets
* The newly developed function ******** is implemented
* The `terra` package is now used for all raster processing with support for `stars` and `raster` objects (or files on disk)
* The new `cv_extrapolate` provides measures on possible extrapolation to testing folds

**Note**: All function names have changed to more general names starting with `cv_*`. The old functions (v2.x) still work, but they will be removed in future versions. Please update your code with the new functions.


## Installation
To install the latest update of the package from GitHub use:

```r
remotes::install_github("rvalavi/blockCV", dependencies = TRUE)
```
Or installing from CRAN:

```r
install.packages("blockCV", dependencies = TRUE)
```

## Vignettes
To see the vignettes of the package use:

1- blockCV introduction: how to create block corss-validation folds

2- Block cross-validation for species distribution modelling

3- Using blockCV with `caret` and `tidymodels`

The vignette is also available via this [link](http://htmlpreview.github.io/?https://github.com/rvalavi/blockCV/blob/master/vignettes/BlockCV_for_SDM.html).


## Basic usage
The following is an example of using spatial block cross-validation for evaluation of species distribution modelling. You can find a comprehensive tutorial in the vignette of the package.

```r
# loading the package
library(blockCV)

# spatial blocking by specified range and random assignment
sb <- cv_spatial(x = pa_data, # sf or SpatialPoints of sample data
                 column = "occ", # the response column (binary or multi-class)
                 r = myrasters, # a raster for background (optional)
                 size = 450000, # size of the blocks in metres
                 k = 5, # number of folds
                 hexagon = TRUE, # use hexagonal blocks - defualt
                 selection = "random", # random blocks-to-fold
                 iteration = 100, # find evenly dispersed folds
                 biomod2 = TRUE) # also create folds for biomod2

```
![](https://i.ibb.co/WGfrF7B/Rplot1.png)

Or create spatial clusters for k-fold cross-validation:

```r
# create spatial clusters
set.seed(6)
sc <- cv_cluster(x = pa_data, 
                 column = "occ", # optionally count data in folds (binary or multi-class)
                 k = 5)

```

```r
# now plot the create folds
cv_plot(cv = sc, # a blockCV object
        x = pa_data, # sample points
        r = myrasters[[1]], # optionally add a raster background
        points_alpha = 0.5,
        nrow = 2)

```
![](https://i.ibb.co/dGrF9xp/Rplot02.png)



```r
# investigate spatial autocorrelation in the landscape
# this helps to choose a suitable size for spatial blocks
cv_spatial_autocor(r = myrasters, # a SpatRaster object or path to files
                   num_sample = 5000, # number of cells to be used
                   plot = TRUE)
```



```r
# alternatively, you can manually choose the size of spatial blocks 
cv_block_size(r = myrasters[[1]],
              x = pa_data, # optionally add sample points
              column = "occ",
              min_size = 2e5,
              max_size = 9e5)

```

## Issues
Please report issues at: 
[https://github.com/rvalavi/blockCV/issues](https://github.com/rvalavi/blockCV/issues)

## Citation
To cite package **blockCV** in publications, please use:

Valavi R, Elith J, Lahoz-Monfort JJ, Guillera-Arroita G. **blockCV: An R package for generating spatially or environmentally separated folds for k-fold cross-validation of species distribution models**. *Methods Ecol Evol*. 2019; 10:225–232. [https://doi.org/10.1111/2041-210X.13107](https://doi.org/10.1111/2041-210X.13107)

