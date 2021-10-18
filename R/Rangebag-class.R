#' S4 class to represent a Range Bagging SDM.
#'
#' Description of Rangebag class ...
#'
#' @slot method SDM method: "rangebag".
#' @slot variables List of climate (or environmental) variable names.
#' @slot presence The selected climate data corresponding to occurrences points.
#' @slot coordinates The coordinates for the selected climate data.
#' @slot ch_models A list of convex hull models (vertices).
setClass("Rangebag",
         slots = c(method = "character",
                   variables = "character",
                   presence = "data.frame",
                   coordinates = "data.frame",
                   ch_models = "list"))
