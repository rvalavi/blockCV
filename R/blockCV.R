#' blockCV: Spatial and Environmental Blocking for K-Fold and LOO Cross-Validation
#'
#' Simple random selection of training and testing folds in the structured environment leads to
#' an underestimation of error in the evaluation of spatial
#' predictions and may result in inappropriate model selection (Telford and Birks, 2009; Roberts et al., 2017). The use of spatial and
#' environmental blocks to separate training and testing sets has been suggested as a good strategy for realistic error estimation in datasets
#' with dependence structures, and more generally as a robust method for estimating the predictive performance of models used to predict mapped
#' distributions (Roberts et al., 2017). The package \code{blockCV} offers
#' a range of functions for generating train and test folds
#' for \strong{k-fold} and \strong{leave-one-out (LOO)} cross-validation (CV). It allows for separation
#' of data spatially and environmentally, with various options for block construction.
#' Additionally, it includes a function for assessing the level of spatial autocorrelation
#' in response or raster covariates, to aid in selecting an appropriate distance band for
#' data separation. The \code{blockCV} package is suitable for the evaluation of a variety of
#' spatial modelling applications, including classification of remote sensing imagery,
#' soil mapping, and species distribution modelling (SDM). It also provides support for
#' different SDM scenarios, including presence-absence and presence-background species
#' data, rare and common species, and raster data for predictor variables.
#'
#' @seealso \code{\link{cv_spatial}}, \code{\link{cv_cluster}}, \code{\link{cv_buffer}}, and \code{\link{cv_nndm}} for blocking strategies.
#'
#' @references Valavi, R., Elith, J., Lahoz-Monfort, J. J., & Guillera-Arroita, G. (2019). blockCV: An R package for generating spatially or environmentally separated folds for k-fold cross-validation of species distribution models. Methods in Ecology and Evolution, 10(2), 225-232. doi:10.1111/2041-210X.13107.
#'
#' @name blockCV
#' @author Roozbeh Valavi, Jane Elith, José Lahoz-Monfort, Ian Flint, and Gurutzeta Guillera-Arroita
#' @import sf
"_PACKAGE"
