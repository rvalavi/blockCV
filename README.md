# blockCV <img src="man/figures/logo.png" align="right" width="120" />

[![Build Status](https://travis-ci.org/rvalavi/blockCV.svg?branch=master)](https://travis-ci.org/rvalavi/blockCV)
[![codecov](https://codecov.io/gh/rvalavi/blockCV/branch/master/graph/badge.svg)](https://codecov.io/gh/rvalavi/blockCV)
[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%203%29-lightgrey.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
[![DOI](https://zenodo.org/badge/116337503.svg)](https://zenodo.org/badge/latestdoi/116337503)
[![CRAN version](https://www.r-pkg.org/badges/version/blockCV)](https://CRAN.R-project.org/package=blockCV)

In a nutshell, package **blockCV** provides functions to build train and test data sets using three general strategies: *buffers*, *spatial* and *environmental* blocks. It offers several options for how those blocks are constructed and how blocks are allocated to cross-validation folds. It includes a function that applies geostatistical techniques to investigate the existing level of spatial autocorrelation in the chosen predictor variables to inform the choice of block and buffer size. In addition, visualization tools further aid selection of block size and provide understanding of the spread of species data across generated folds. 

### Installation
To install the package from GitHub use:

```r
# install.packages("devtools")
devtools::install_github("rvalavi/blockCV")
```

### Bsic usage
The following is an example of using spatial block cross-validation for evaluation of species distribution modelling. You can find a comprehensive tutorial in the vignette of the package.

```r
# loading the package
library(blockCV)

# investigate spatial autocorrelation in raster covariates
spatialAutoRange(rasterLayer = awt, # raster file
                 sampleNumber = 5000, # number of cells to be used
                 doParallel = TRUE,
                 showPlots = TRUE)
```
![](https://image.ibb.co/eQGJZ8/spatial_Auto.jpg)

```r
# spatial blocking by specified range and random assignment
sb <- spatialBlock(speciesData = pa_data,
                   species = "Species", # the response column (binomial or multi-class)
                   rasterLayer = awt,
                   theRange = 70000, # size of the blocks
                   k = 5,
                   selection = "random",
                   iteration = 100, # find evenly dispersed folds
                   biomod2Format = TRUE)

```
![](https://image.ibb.co/dpDMnT/Rplot01.jpg)

```r
# explor the generated folds
foldExplorer(sb, awt, pa_data)

```

### Vignette
To see the vignette of the package use:

```r
browseVignettes("blockCV")
```
The vignette is also available via this [link](http://htmlpreview.github.io/?https://github.com/rvalavi/blockCV/blob/master/vignettes/BlockCV_for_SDM.html).

### Features
Compared to other available packages, **blockCV** provides more strategies and control over fold generation including:

* The assignment of the spatial blocks to cross-validation folds can be done in three different ways.
* The position of the spatial blocks can be modified. 
* The buffering function can account for presence-absence and presence-background data types.
* The variables are standardized to avoid wide range variables to dominate the environmental blocks
* Using geostatistical techniques to inform the choice of a suitable distance band by which to separate the data sets 


### Citation
To cite package ‘blockCV’ in publications use:

Valavi R, Elith J, Lahoz-Monfort JJ, Guillera-Arroita G. **blockCV: An r package for generating spatially or environmentally separated folds for k-fold cross-validation of species distribution models**. *Methods Ecol Evol*. 2019; 10:225–232. [https://doi.org/10.1111/2041-210X.13107](https://doi.org/10.1111/2041-210X.13107)

