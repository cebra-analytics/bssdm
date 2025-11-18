#' @description
#' The bssdm package provides tools for species distribution modelling
#' with a focus on biosecurity applications. It includes implementations
#' of range bagging and Climatch algorithms, along with similarity
#' analysis methods (MESS, ExDet) and model validation tools.
#'
#' @section Main functions:
#' * [rangebag()]: Range bagging species distribution modelling
#' * [climatch()]: Climatch species distribution modelling
#' * [mess()]: Multivariate Environmental Similarity Surfaces
#' * [exdet()]: Extrapolation detection in environmental space
#' * [boyce()]: Boyce index for model validation
#'
#' @section Model prediction:
#' * [predict.Rangebag()]: Predict method for Rangebag models
#' * [predict.Climatch()]: Predict method for Climatch models
#'
#' @references
#' Drake, J. M. (2015). Range bagging: a new method for ecological
#' niche modelling from presence-only data.
#' _Journal of the Royal Society Interface_, 12(107), 20150086.
#' \doi{10.1098/rsif.2015.0086}
#'
#' @references
#' ABARES (2020). Climatch v2.0 User Manual. Canberra.
#' \url{https://climatch.cp1.agriculture.gov.au/}
#'
#' @references
#' Elith, J., Kearney, M., & Phillips, S. (2010). The art of modelling
#' range-shifting species. _Methods in Ecology and Evolution_, 1(4), 330-342.
#' \doi{10.1111/j.2041-210X.2010.00036.x}
#'
#' @references
#' Mesgaran, M. B., Cousens, R. D., & Webber, B. L. (2014). Here be dragons:
#' a tool for quantifying novelty due to covariate range and correlation change
#' when projecting species distribution models. _Diversity and Distributions_,
#' 20(10), 1147-1159. \doi{10.1111/ddi.12209}
#'
#' @references
#' Hirzel, A. H., Le Lay, G., Helfer, V., Randin, C., & Guisan, A. (2006).
#' Evaluating the ability of habitat suitability models to predict species
#' presences. _Ecological Modelling_, 199(2), 142-152.
#' \doi{10.1016/j.ecolmodel.2006.05.017}
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
