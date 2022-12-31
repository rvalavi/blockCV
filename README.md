# blockCV <img src="man/figures/logo.png" align="right" width="120" />

[![Build Status](https://travis-ci.org/rvalavi/blockCV.svg?branch=master)](https://travis-ci.org/rvalavi/blockCV)
[![codecov](https://codecov.io/gh/rvalavi/blockCV/branch/master/graph/badge.svg)](https://codecov.io/gh/rvalavi/blockCV)
[![CRAN version](https://www.r-pkg.org/badges/version/blockCV)](https://CRAN.R-project.org/package=blockCV)
[![total](http://cranlogs.r-pkg.org/badges/grand-total/blockCV)](https://www.rpackages.io/package/blockCV)
[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%203%29-lightgrey.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![DOI](https://zenodo.org/badge/116337503.svg)](https://zenodo.org/badge/latestdoi/116337503)


### Spatial and environmental blocking for k-fold and LOO cross-validation   
   
In a nutshell, package **blockCV** provides functions to build train and test data sets using three general strategies: *buffers*, *spatial* and *environmental* blocks. It offers several options for how those blocks are constructed and how blocks are allocated to cross-validation folds. It includes a function that applies geostatistical techniques to investigate the existing level of spatial autocorrelation in the chosen predictor variables or the response variable (e.g. binary or continuous responses) to inform the choice of the block and buffer size. In addition, visualization tools further aid the selection of block size and provide an understanding of the spread of sample points across generated folds. 


## Features
Compared to other available packages, **blockCV** provides more strategies and control over fold generation including:

* There are three blocking methods: **buffers**, **spatial** and **clustering** blocks
* The assignment of the spatial blocks to cross-validation folds can be done in three different ways: **random**, **systematic** and **checkerboard** pattern
* The spatial blocks can be assigned to cross-validation folds to have **evenly distributed records** for **binary** (e.g. species presence-absence/background) or **multi-class** responses (e.g. land cover classes for remote sensing image classification) 
* The position of the spatial blocks can be modified 
* The buffering function can account for *presence-absence* and *presence-background* data types 
* Using geostatistical techniques to inform the choice of a suitable distance band by which to separate the data sets 

## New updates of version 3
* Function names have changed, all starting with `cv_*` now
* Spatial blocks now support hexagonal (default now), square, and user-defined blocks
* Clustering function now works both on *environmental* rasters and *spatial coordinates* of the sample points
* The `cv_spatial_autocor` function now calculates spatial autocorrelation range for either the response (i.e. the binary or continuous data) or a set of continuous raster covariates (as before)
* The new `cv_plot` function can be used to plot the folds of all blocking strategy
* The newly developed function ********.
* 


**Note**: All function names have changed to more general names starting with `cv_*`. The old functions (v2.x) are still working but they will be removed in future versions. Please update your code with the new naming.

## Installation
To install the package from GitHub use:

```r
remotes::install_github("rvalavi/blockCV", dependencies = TRUE)
```
Or installing from CRAN:

```r
install.packages("blockCV", dependencies = TRUE)
```

## Vignette
To see the vignette of the package use:

```r
browseVignettes("blockCV")
```
The vignette is also available via this [link](http://htmlpreview.github.io/?https://github.com/rvalavi/blockCV/blob/master/vignettes/BlockCV_for_SDM.html).


## Basic usage
The following is an example of using spatial block cross-validation for evaluation of species distribution modelling. You can find a comprehensive tutorial in the vignette of the package.

```r
# loading the package
library(blockCV)

# spatial blocking by specified range and random assignment
sb <- cv_spatial(x = pa_data, # sf or SpatialPoints of sample data
                 column = "occ", # the response column (binomial or multi-class)
                 r = myrasters, # a raster for background (optional)
                 size = 450000, # size of the blocks in metres
                 k = 5, # number of folds
                 selection = "random", # random blocks-to-fold
                 iteration = 100, # find evenly dispersed folds
                 biomod2 = TRUE) # also create folds for biomod2

```
![](https://i.ibb.co/WGfrF7B/Rplot1.png)

```r
# create spatial clusters
set.seed(1)
sc <- cv_cluster(x = pa_data, 
                 column = "occ", # optionally count data in folds
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
# ![](https://i.ibb.co/1MYWj8n/Rplot01.png)



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
              min_size = 100000,
              max_size = 500000)

```



## Citation
To cite package **blockCV** in publications, please use:

Valavi R, Elith J, Lahoz-Monfort JJ, Guillera-Arroita G. **blockCV: An R package for generating spatially or environmentally separated folds for k-fold cross-validation of species distribution models**. *Methods Ecol Evol*. 2019; 10:225–232. [https://doi.org/10.1111/2041-210X.13107](https://doi.org/10.1111/2041-210X.13107)

