#' S4 class to represent a Climatch SDM.
#'
#' Description of Climatch class ...
#'
#' @slot method SDM method: "climatch".
#' @slot algorithm Algorithm: "euclidean" or "closest_standard_score".
#' @slot variables List of climate (or environmental) variable names.
#' @slot sd The standard deviation of each variable calculated via the climate data (\emph{x}) or the \emph{sd_data} when provided.
#' @slot presence The selected (nearest within range) climate data for each occurrence point.
#' @slot coordinates The coordinates for the selected climate data.
#' @slot as_score Indication of whether to generate a score 0-10 or values 0-1.
setClass("Climatch",
         slots = c(method = "character",
                   algorithm = "character",
                   variables = "character",
                   sd = "numeric",
                   presence = "data.frame",
                   coordinates = "data.frame",
                   as_score = "logical"))
