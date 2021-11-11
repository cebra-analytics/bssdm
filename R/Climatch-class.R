#' S4 class to represent a climatch SDM.
#'
#' The model object (S4) class component for an implementation of the ABARES
#' Climatch species distribution modelling (SDM) method (ABARES, 2020).
#'
#' @slot method SDM method: "climatch".
#' @slot algorithm Algorithm: "euclidean" or "closest_standard_score".
#' @slot variables List of climate (or environmental) variable names.
#' @slot sd The standard deviation of each variable calculated via the climate
#'   data (\emph{x}) or the \emph{sd_data} when provided.
#' @slot presence The selected (nearest within range) climate data for each
#'   occurrence point.
#' @slot coordinates The coordinates for the selected climate data.
#' @slot as_score Indication of whether to generate a score 0-10 or values 0-1.
#' @references ABARES (2020). Climatch v2.0 User Manual. Canberra.
#'   \url{https://climatch.cp1.agriculture.gov.au/} Accessed: November 2021.
methods::setClass("Climatch",
                  slots = c(method = "character",
                            algorithm = "character",
                            variables = "character",
                            sd = "numeric",
                            presence = "data.frame",
                            coordinates = "data.frame",
                            as_score = "logical"))
