#' S4 class to represent a range bagging SDM.
#'
#' The model object (S4) class component for an implementation of the range
#' bagging species distribution modelling (SDM) method (Drake, 2015).
#'
#' @slot method SDM method: "rangebag".
#' @slot variables List of climate (or environmental) variable names.
#' @slot presence The selected climate data corresponding to occurrences
#'   points.
#' @slot coordinates The coordinates for the selected climate data.
#' @slot ch_models A list of convex hull models (vertices).
#' @references Drake, J. M. (2015). Range bagging: a new method for ecological
#'   niche modelling from presence-only data.
#'   \emph{Journal of the Royal Society Interface}, 12(107), 20150086.
#'   \doi{10.1098/rsif.2015.0086}
methods::setClass("Rangebag",
                  slots = c(method = "character",
                            variables = "character",
                            presence = "data.frame",
                            coordinates = "data.frame",
                            ch_models = "list"))
