#' Climatch SDM model building method
#'
#' The model building component of an implementation of the ABARES Climatch
#' species distribution modelling (SDM) method (ABARES, 2020).
#'
#' @param x Climate (or environmental) data as a \code{raster::Raster*} or
#'   \code{terra::SpatRaster} (with any CRS), or a \code{data.frame} with
#'   WGS84 \emph{lon} and \emph{lat} columns.
#' @param p Species occurrence data as a \code{data.frame} (or \code{matrix})
#'   with WGS84 \emph{lon} and \emph{lat} columns.
#' @param algorithm Climatch method algorithm selected from "euclidean"
#'   (default) or "closest_standard_score".
#' @param d_max Maximum range distance (default = 50 km) used when matching
#'   occurrence points to nearest climate data points/cells.
#' @param sd_data Optional \code{data.frame} for calculating the standard
#'   deviation for climate variable, or a \code{vector} of pre-calculated
#'   values.
#' @param as_score Logical to indicate whether to generate a score 0-10
#'   (default = TRUE) or values 0-1 (FALSE).
#' @param ... Additional parameters.
#' @return A "Climatch" model S4 object containing slots:
#'   \describe{
#'     \item{\code{method}}{SDM method: "climatch".}
#'     \item{\code{algorithm}}{Algorithm: "euclidean" or
#'       "closest_standard_score").}
#'     \item{\code{variables}}{List of climate (or environmental) variable
#'       names.}
#'     \item{\code{sd}}{The standard deviation of each variable calculated via
#'       the climate data (\emph{x}) or the \emph{sd_data} when provided.}
#'     \item{\code{presence}}{The selected (nearest within range) climate data
#'       for each occurrence point.}
#'     \item{\code{coordinates}}{The coordinates for the selected climate
#'       data.}
#'     \item{\code{as_score}}{Indication of whether to generate a score 0-10
#'       or values 0-1.}
#'   }
#' @references ABARES (2020). Climatch v2.0 User Manual. Canberra.
#'   \url{https://climatch.cp1.agriculture.gov.au/} Accessed: November 2021.
#' @include Climatch-class.R
#' @export
climatch <- function(x, p,
                     algorithm = "euclidean",
                     d_max = 50, # km
                     sd_data = NULL,
                     as_score = TRUE, ...) {
  UseMethod("climatch")
}

#' @name climatch
#' @export
climatch.Raster <- function(x, p,
                            algorithm = "euclidean",
                            d_max = 50, # km
                            sd_data = NULL,
                            as_score = TRUE, ...) {

  # Project to lon/lat if necessary (needed for distance calculations)
  if (!raster::isLonLat(x)) {
    x <- raster::projectRaster(x, crs = "EPSG:4326")
  }

  # Convert x to a data frame with lon, lat, and variables
  x <- raster::as.data.frame(x, xy = TRUE, na.rm = TRUE)
  names(x)[1:2] <- c("lon", "lat")

  # Call the data frame version of the function
  climatch(x, p,
           algorithm = algorithm,
           d_max = d_max,
           sd_data = sd_data,
           as_score = as_score, ...)
}

#' @name climatch
#' @export
climatch.SpatRaster <- function(x, p,
                                algorithm = "euclidean",
                                d_max = 50, # km
                                sd_data = NULL,
                                as_score = TRUE, ...) {

  # Project to lon/lat if necessary (needed for distance calculations)
  if (!terra::is.lonlat(x)) {
    x <- terra::project(x, crs = "EPSG:4326")
  }

  # Convert x to a data frame with lon, lat, and variables
  x <- terra::as.data.frame(x, xy = TRUE, na.rm = TRUE)
  names(x)[1:2] <- c("lon", "lat")

  # Call the data frame version of the function
  climatch(x, p,
           algorithm = algorithm,
           d_max = d_max,
           sd_data = sd_data,
           as_score = as_score, ...)
}

#' @name climatch
#' @export
climatch.data.frame <- function(x, p,
                                algorithm = "euclidean",
                                d_max = 50, # km
                                sd_data = NULL,
                                as_score = TRUE, ...) {

  # Check that x and p have sufficient columns
  if (ncol(x) < 3) {
    stop("Climatch x data should have at least 3 columns.", call. = FALSE)
  }
  p <- as.data.frame(p)
  if (ncol(p) < 2) {
    stop("Climatch p data should have at least 2 columns.", call. = FALSE)
  }

  # Make sure x and p coordinate columns are 'lon' and 'lat'
  if (!all(c("lon", "lat") %in% names(x))) {
    stop("Climatch x coordinates should be 'lon' and 'lat'.", call. = FALSE)
  }
  if (!all(c("lon", "lat") %in% names(p))) {
    stop("Climatch p coordinates should be 'lon' and 'lat'.", call. = FALSE)
  }

  # Check algorithm
  if (!algorithm %in% c("euclidean", "closest_standard_score")) {
    stop(paste("Climatch algorithm should be 'euclidean' or",
               "'closest_standard_score'."),
         call. = FALSE)
  }

  # Calculate standard deviations when required
  first_variable_col <- (max(which(names(x) %in% c("lon", "lat"))) + 1)
  variables <- names(x)[first_variable_col:ncol(x)]
  if (!is.null(sd_data)) {
    if (is.data.frame(sd_data)) {
      if (!all(variables %in% names(sd_data))) {
        stop("Climatch sd data should have the same variables as x.",
             call. = FALSE)
      }
      sd_data <- unlist(lapply(sd_data[, variables], sd))
    } else if (is.numeric(sd_data) &&
               length(sd_data) != length(variables)) {
      stop("Climatch predict sd data should match x variables.",
           call. = FALSE)
    }
  } else {
    sd_data <- unlist(lapply(x[, variables], sd))
  }

  # Get data points (x) within a fixed range (d_max) of the occurrence
  # points (p). Select up to one data point per occurrence point.
  source_data <- cbind(id = 1:nrow(x), x[, c("lon", "lat")])
  selected_idx <- c()
  for (i in 1:nrow(p)) {

    # Use cheap distance calculation to find close points
    suppressMessages(
      distances <- geodist::geodist(p[i, c("lon", "lat")],
                                    source_data[, c("lon", "lat")],
                                    measure = "cheap")/1000) # km
    close_idx <- which(distances <= d_max*1.5) # allow for inaccuracy

    # Recalculate close points with accurate method
    if (length(close_idx)) {
      distances[close_idx] <-
        geosphere::distGeo(p[i, c("lon", "lat")],
                           source_data[close_idx, c("lon", "lat")])/1000 # km
    }

    # Add closest point when in range
    closest <- which.min(distances)
    if (distances[closest] <= d_max) {
      selected_idx <- c(selected_idx, source_data$id[closest])
      source_data <- source_data[-closest,]
    }
  }

  # Return a "Climatch" object with configured, calculated and selected data
  return(methods::new("Climatch",
                      method = "climatch",
                      algorithm = algorithm,
                      variables = variables,
                      sd = sd_data,
                      presence = x[selected_idx, variables],
                      coordinates = x[selected_idx, c("lon", "lat")],
                      as_score = as_score))
}
