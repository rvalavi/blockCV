[![Build Status](https://travis-ci.org/rvalavi/blockCV.svg)](https://travis-ci.org/rvalavi/blockCV)
[![codecov.io](https://codecov.io/github/rvalavi/blockCV/coverage.svg?branch=master)](https://codecov.io/github/rvalavi/blockCV?branch=master)

# blockCV

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

# spatial blocking by specified range and random assignment
sb <- spatialBlock(speciesData = pa_data,
                   species = "Species",
                   rasterLayer = awt,
                   theRange = 68000, # size of the blocks
                   k = 5,
                   selection = "random",
                   iteration = 250, # find evenly dispersed folds
                   biomod2Format = TRUE)

# explor the generated folds
foldExplorer(sb, awt, pa_data)

```

### Features
Compared to other available packages, **blockCV** provides more strategies and control over fold generation including:

* The assignment of the spatial blocks to cross-validation folds can be done in two different ways.
* The position of the spatial blocks can be modified. 
* The buffering function can account for presence-absence and presence-background data types.
* The variables are standardized to avoid wide range variables to dominate the environmental blocks
* Using geostatistical techniques to inform the choice of a suitable distance band by which to separate the data sets 


### Citation
To cite package ‘blockCV’ in publications use:

  Roozbeh Valavi, Gurutzeta Guillera-Arroita, José Lahoz-Monfort and Jane Elith
  (2018). blockCV: Spatial and environmental blocking for k-fold cross-validation.
  R package version 0.1.0.

Citation will be updated.
