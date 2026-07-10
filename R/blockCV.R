#' blockCV: Spatial, Environmental, and Grouped Cross-Validation
#'
#' Random assignment of spatially structured or grouped observations to training and
#' testing folds can underestimate prediction error and lead to inappropriate model
#' selection (see Valavi et al., 2019, and references therein). The \code{blockCV}
#' package creates training and testing folds for \strong{k-fold},
#' \strong{leave-group-out}, and \strong{leave-one-out (LOO)} cross-validation.
#' Strategies include spatial blocks, spatial or environmental clustering, existing
#' grouping factors, spatial buffers, and nearest neighbour distance matching. The
#' package also provides tools to visualise and diagnose fold designs, compare
#' train-test separation with the prediction domain, assess environmental novelty,
#' and investigate spatial autocorrelation and candidate block sizes. It supports
#' spatial modelling applications such as remote-sensing classification, soil
#' mapping, and species distribution modelling, including presence-absence and
#' presence-background data.
#'
#' @seealso \code{vignette("tutorial_1", package = "blockCV")} for examples of all fold-construction strategies, and \code{vignette("tutorial_2", package = "blockCV")} for fold assessment and design.
#'
#' @references Valavi, R., Elith, J., Lahoz-Monfort, J. J., & Guillera-Arroita, G. (2019). blockCV: An R package for generating spatially or environmentally separated folds for k-fold cross-validation of species distribution models. Methods in Ecology and Evolution, 10(2), 225-232. doi:10.1111/2041-210X.13107.
#'
#' @name blockCV
#' @author Roozbeh Valavi, Jane Elith, José Lahoz-Monfort, Ian Flint, and Gurutzeta Guillera-Arroita
#' @import sf
"_PACKAGE"
