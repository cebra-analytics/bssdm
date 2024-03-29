#' Range bagging SDM predict method
#'
#' The model prediction component of an implementation of the range bagging
#' species distribution modelling (SDM) method (Drake, 2015).
#'
#' @param object A "Rangebag" model S4 object containing slots:
#'   \describe{
#'     \item{\code{method}}{SDM method: "rangebag".}
#'     \item{\code{variables}}{List of climate (or environmental) variable
#'       names.}
#'     \item{\code{presence}}{The selected climate data corresponding to
#'       occurrences points.}
#'     \item{\code{coordinates}}{The coordinates for the selected climate
#'       data.}
#'     \item{\code{ch_models}}{A list of convex hull models (vertices).}
#'   }
#' @param x Climate (or environmental) data with corresponding model variables
#'   as a \code{terra::SpatRaster}, \code{raster::Raster*}, \code{data.frame},
#'   or \code{matrix}.
#' @param raw_output Logical to indicate whether to return raw predicted
#'   values (TRUE) or as an object (as per \emph{x}: FALSE). Default is NULL,
#'   returning either raw values or a spatial raster (as per \emph{x}).
#' @param filename Optional filename for writing spatial raster output (only).
#'   Default is "".
#' @param ... Additional parameters.
#' @return Predicted values as a raw vector or a \code{terra::SpatRaster},
#'   \code{raster::Raster*}, \code{data.frame}, or \code{matrix} (as per
#'   \emph{x}).
#' @references Drake, J. M. (2015). Range bagging: a new method for ecological
#'   niche modelling from presence-only data.
#'   \emph{Journal of the Royal Society Interface}, 12(107), 20150086.
#'   \doi{10.1098/rsif.2015.0086}
#' @export
predict.Rangebag <- function(object, x,
                             raw_output = NULL,
                             filename = "", ...) {

  # Ensure x data has matching object variables
  if (!(all(object@variables %in% names(x)) ||
        all(object@variables %in% colnames(x)))) {
    stop("Rangebag predict x should have same variables as object.",
         call. = FALSE)
  }

  # Extract data values for object variables from x
  if (class(x)[1] %in% c("Raster", "RasterStack", "RasterBrick")) {
    x_data <- raster::as.data.frame(x, xy = TRUE, na.rm = TRUE)
    x_coords <- x_data[, c("x", "y")]
    x_data <- as.matrix(x_data[, object@variables, drop = FALSE])
  } else if (class(x)[1] == "SpatRaster") {
    x_data <- terra::as.data.frame(x, xy = TRUE, na.rm = TRUE)
    x_coords <- x_data[, c("x", "y")]
    x_data <- as.matrix(x_data[, object@variables, drop = FALSE])
  } else if (is.data.frame(x) || is.matrix(x)) {
    x_data <- as.matrix(x[, object@variables, drop = FALSE])
  }

  # Count the number of convex hull model fits for each x data row/cell
  n_models <- length(object@ch_models)
  n_dim <- ncol(object@ch_models[[1]])
  ch_counts <- numeric(nrow(x_data))
  for (i in 1:n_models) {
    vars <- colnames(object@ch_models[[i]])
    if (n_dim == 1) {
      data_in_ch <- (x_data[, vars] >= object@ch_models[[i]][1, 1] &
                       x_data[, vars] <= object@ch_models[[i]][2, 1])
    } else {
      data_in_ch <- geometry::inhulln(
        geometry::convhulln(object@ch_models[[i]], options='Pp'),
        x_data[, vars])
    }
    ch_counts <- ch_counts + data_in_ch
  }

  # Return the count fraction as a raw vector, raster or data frame
  if (is.null(raw_output)) {
    if (is.data.frame(x) || is.matrix(x)) {
      raw_output <- TRUE
    } else {
      raw_output <- FALSE
    }
  }
  if (raw_output) {
    return(ch_counts/n_models)
  } else if (class(x)[1] %in% c("Raster", "RasterStack", "RasterBrick")) {
    output_rast <- raster::extend(
      raster::rasterFromXYZ(cbind(x_coords, predicted = ch_counts/n_models),
                            res = raster::res(x), crs = raster::crs(x)),
      raster::extent(x))
    if (filename != "") {
      output_rast <- raster::writeRaster(output_rast, filename)
    }
    return(output_rast)
  } else if (class(x)[1] == "SpatRaster") {
    return(terra::extend(
      terra::rast(cbind(x_coords, predicted = ch_counts/n_models),
                  type = "xyz", crs = terra::crs(x)), terra::ext(x),
      filename = filename))
  } else if (is.data.frame(x) || is.matrix(x)) {
    return(cbind(x[, which(!names(x) %in% object@variables)],
                 predicted = ch_counts/n_models))
  }
}
