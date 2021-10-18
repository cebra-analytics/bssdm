#' Climatch SDM predict method
#'
#' Description of Climatch SDM predict method...
#'
#' @param object A "Climatch" model S4 object containing slots:
#'   \describe{
#'     \item{\code{method}}{SDM method: "rangebag".}
#'     \item{\code{variables}}{List of climate (or environmental) variable names.}
#'     \item{\code{presence}}{The selected climate data corresponding to occurrences points.}
#'     \item{\code{coordinates}}{The coordinates for the selected climate data.}
#'     \item{\code{ch_models}}{A list of convex hull models (vertices).}
#'   }
#' @param x Climate (or environmental) data with corresponding model variables as a \code{Raster*}, \code{data.frame}, or \code{matrix}.
#' @param raw_output Logical to indicate whether to return raw predicted values (default = TRUE) or as an object (as per \emph{x}: FALSE).
#' @param ... Additional parameters.
#' @return Predicted values as a raw vector or a \code{Raster*}, \code{data.frame}, or \code{matrix} (as per \emph{x}).
#' @export
predict.Rangebag <- function(object, x,
                             raw_output = TRUE, ...) {

  # Ensure x data has matching object variables
  if (!(all(object@variables %in% names(x)) ||
        all(object@variables %in% colnames(x)))) {
    stop("Rangebag predict x should have same variables as object.",
         call. = FALSE)
  }

  # Extract data values for object variables from x
  if (class(x)[1] %in% c("Raster", "RasterStack", "RasterBrick")) {
    x_data <- raster::as.data.frame(x, xy=TRUE, na.rm = TRUE)
    x_coords <- x_data[, c("x", "y")]
    x_data <- as.matrix(x_data[, object@variables])
  } else if (is.data.frame(x) || is.matrix(x)) {
    x_data <- as.matrix(x[, object@variables])
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
  if (raw_output) {
    return(ch_counts/n_models)
  } else if (class(x)[1] %in% c("Raster", "RasterStack", "RasterBrick")) {
    return(raster::rasterFromXYZ(cbind(x_coords,
                                       predicted = ch_counts/n_models),
                                 res = raster::res(x),
                                 crs = raster::crs(x)))
  } else if (is.data.frame(x) || is.matrix(x)) {
    return(cbind(x[, which(!names(x) %in% object@variables)],
                 predicted = ch_counts/n_models))
  }
}
