% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blocking.R
\name{spatialBlock}
\alias{spatialBlock}
\title{Use spatial blocks to separate train and test folds}
\usage{
spatialBlock(
  speciesData,
  species = NULL,
  rasterLayer = NULL,
  theRange = NULL,
  rows = NULL,
  cols = NULL,
  k = 5L,
  selection = "random",
  iteration = 100L,
  blocks = NULL,
  foldsCol = NULL,
  numLimit = 0L,
  maskBySpecies = TRUE,
  degMetre = 111325,
  border = NULL,
  showBlocks = TRUE,
  biomod2Format = TRUE,
  xOffset = 0,
  yOffset = 0,
  seed = NULL,
  progress = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{speciesData}{A simple features (sf) or SpatialPoints object containing species data (response variable).}

\item{species}{Character (optional). Indicating the name of the column in which species data (response variable e.g. 0s and 1s) is stored.
This argument is used \emph{to make folds with evenly distributed records}. \strong{This option only works by random fold selection and with binary or
multi-class responses} e.g. species presence-absence/background or land cover classes for remote sensing image classification.
If \code{speceis = NULL} the response classes will be treated the same and only training and testing records
will be counted and balanced.}

\item{rasterLayer}{A raster object for visualisation (optional). If provided, this will be used to specify the blocks covering the area.}

\item{theRange}{Numeric value of the specified range by which blocks are created and training/testing data are separated.
This distance should be in \strong{metres}. The range could be explored by \code{spatialAutoRange()} and \code{rangeExplorer()} functions.}

\item{rows}{Integer value by which the area is divided into latitudinal bins.}

\item{cols}{Integer value by which the area is divided into longitudinal bins.}

\item{k}{Integer value. The number of desired folds for cross-validation. The default is \code{k = 5}.}

\item{selection}{Type of assignment of blocks into folds. Can be \strong{random} (default), \strong{systematic}, \strong{checkerboard}, or \strong{predefined}.
The checkerboard does not work with user-defined spatial blocks. If the selection = 'predefined', user-defined blocks and foldsCol must be supplied.}

\item{iteration}{Integer value. The number of attempts to create folds that fulfil the set requirement for minimum number
of points in each trainig and testing fold (for each response class e.g. \emph{train_0}, \emph{train_1}, \emph{test_0}
and \emph{test_1}), as specified by \code{species} and \code{numLimit} arguments.}

\item{blocks}{A sf or SpatialPolygons object to be used as the blocks (optional). This can be a user defined polygon and it must cover all
the species (response) points. If the selection = 'predefined', this argument (and foldsCol) must be supplied.}

\item{foldsCol}{Character. Indicating the name of the column (in user-defined blocks) in which the associated folds are stored.
This argument is necessary if you choose the 'predefined' selection.}

\item{numLimit}{Integer value. The minimum number of points in each training and testing folds.
If \code{numLimit = 0}, the most evenly dispersed number of records is chosen (given the number of iteration).
This option no longer accepts NULL as input. If it is set to NULL, 0 is used instead.}

\item{maskBySpecies}{Since version 1.1, this option is always set to \code{TRUE}.}

\item{degMetre}{Integer. The conversion rate of metres to degree. See the details section for more information.}

\item{border}{A sf or SpatialPolygons object to clip the block based on it (optional).}

\item{showBlocks}{Logical. If TRUE the final blocks with fold numbers will be created with ggplot and plotted. A raster layer could be specified
in \code{rasterlayer} argument to be as background.}

\item{biomod2Format}{Logical. Creates a matrix of folds that can be directly used in the \pkg{biomod2} package as
a \emph{DataSplitTable} for cross-validation.}

\item{xOffset}{Numeric value between \strong{0} and \strong{1} for shifting the blocks horizontally.
The value is the proportion of block size.}

\item{yOffset}{Numeric value between \strong{0} and \strong{1} for shifting the blocks vertically. The value is the proportion of block size.}

\item{seed}{Integer. A random seed generator for reproducibility.}

\item{progress}{Logical. If TRUE shows a progress bar when \code{numLimit = NULL} in random fold selection.}

\item{verbose}{Logical. To print the report of the recods per fold.}
}
\value{
An object of class S3. A list of objects including:
   \itemize{
    \item{folds - a list containing the folds. Each fold has two vectors with the training (first) and testing (second) indices}
    \item{foldID - a vector of values indicating the number of the fold for each observation (each number corresponds to the same point in species data)}
    \item{biomodTable - a matrix with the folds to be used in \pkg{biomod2} package}
    \item{k - number of the folds}
    \item{blocks - SpatialPolygon of the blocks}
    \item{range - the distance band of separating trainig and testing folds, if provided}
    \item{species - the name of the species (column), if provided}
    \item{plots - ggplot object}
    \item{records - a table with the number of points in each category of training and testing}
    }
}
\description{
This function creates spatially separated folds based on a pre-specified distance. It assigns blocks to the training and
testing folds  \strong{randomly},  \strong{systematically} or in a  \strong{checkerboard pattern}. The distance (\code{theRange})
should be in \strong{metres}, regardless of the unit of the reference system of
the input data (for more information see the details section). By default,
the function creates blocks according to the extent and shape of the study area, assuming that the user has considered the
landscape for the given species and case study. Alternatively, blocks can solely be created based on species spatial data.
Blocks can also be offset so the origin is not at the outer
corner of the rasters. Instead of providing a distance, the blocks can also be created by specifying a number of rows and/or
columns and divide the study area into vertical or horizontal bins, as presented in Wenger & Olden (2012) and Bahn & McGill (2012).
Finally, the blocks can be specified by a user-defined spatial polygon layer.
}
\details{
To keep the consistency, all the functions use \strong{metres} as their unit. In this function, when the input map
has geographic coordinate system (decimal degrees), the block size is calculated based on deviding \code{theRange} by
111325 (the standard distance of a degree in metres, on the Equator) to change the unit to degree. This value is optional
and can be changed by user via \code{degMetre} argument.

The \code{xOffset} and \code{yOffset} can be used to change the spatial position of the blocks. It can also be used to
assess the sensitivity of analysis results to shifting in the blocking arrangements. These options are available when \code{theRange}
is defined. By default the region is located in the middle of the blocks and by setting the offsets, the blocks will shift.

Roberts et. al. (2017) suggest that blocks should be substantially bigger than the range of spatial
autocorrelation (in model residual) to obtain realistic error estimates, while a buffer with the size of
the spatial autocorrelation range would result in a good estimation of error. This is because of the so-called
edge effect (O'Sullivan & Unwin, 2014), whereby points located on the edges of the blocks of opposite sets are
not separated spatially. Blocking with a buffering strategy overcomes this issue (see \code{\link{buffering}}).
}
\examples{
\donttest{

# load package data
library(sf)

awt <- raster::brick(system.file("extdata", "awt.grd", package = "blockCV"))
# import presence-absence species data
PA <- read.csv(system.file("extdata", "PA.csv", package = "blockCV"))
# make a sf object from data.frame
pa_data <- sf::st_as_sf(PA, coords = c("x", "y"), crs = raster::crs(awt))

# spatial blocking by specified range and random assignment
sb1 <- spatialBlock(speciesData = pa_data,
                    species = "Species",
                    theRange = 70000,
                    k = 5,
                    selection = "random",
                    iteration = 100,
                    numLimit = NULL,
                    biomod2Format = TRUE,
                    xOffset = 0.3, # shift the blocks horizontally
                    yOffset = 0)

# spatial blocking by row/column and systematic fold assignment
sb2 <- spatialBlock(speciesData = pa_data,
                    species = "Species",
                    rasterLayer = awt,
                    rows = 5,
                    cols = 8,
                    k = 5,
                    selection = "systematic",
                    biomod2Format = TRUE)

}
}
\references{
Bahn, V., & McGill, B. J. (2012). Testing the predictive performance of distribution models. Oikos, 122(3), 321-331.

O'Sullivan, D., Unwin, D.J., (2010). Geographic Information Analysis, 2nd ed. John Wiley & Sons.

Roberts et al., (2017). Cross-validation strategies for data with temporal, spatial, hierarchical,
or phylogenetic structure. Ecography. 40: 913-929.

Wenger, S.J., Olden, J.D., (2012). Assessing transferability of ecological models: an underappreciated aspect of statistical
validation. Methods Ecol. Evol. 3, 260-267.
}
\seealso{
\code{\link{spatialAutoRange}} and \code{\link{rangeExplorer}} for selecting block size; \code{\link{buffering}}
and \code{\link{envBlock}} for alternative blocking strategies; \code{\link{foldExplorer}} for visualisation of the generated folds.

For \emph{DataSplitTable} see \code{\link[biomod2]{BIOMOD_cv}} in \pkg{biomod2} package
}
