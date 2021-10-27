#' Range bagging SDM model building method
#'
#' Description of Rangebag SDM method...
#'
#' @param x Climate (or environmental) data as a \code{raster::Raster*} or
#'   \code{terra::SpatRaster} (with any CRS), or a \code{data.frame} with
#'   WGS84 \emph{lon} and \emph{lat} columns.
#' @param p Species occurrence data as a \code{data.frame} (or \code{matrix})
#'   with WGS84 \emph{lon} and \emph{lat} columns.
#' @param n_models Number of convex hull models to build in sampled environment
#'   space (default = 100).
#' @param n_dim Number of dimensions (variables) of sampled convex hull models
#'   (default = 2).
#' @param sample_prop Proportion of environment data rows sampled for fitting
#'   each convex hull model (default = 0.5).
#' @param limit_occur Logical to indicate whether to limit occurrence data to
#'   one per environment data cell (default = TRUE).
#' @param ... Additional parameters.
#' @return A "Rangebag" model S4 object containing slots:
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
#' @include Rangebag-class.R
#' @export
rangebag <- function(x, p,
                     n_models = 100,
                     n_dim = 2,
                     sample_prop = 0.5,
                     limit_occur = TRUE, ...) {
  UseMethod("rangebag")
}

#' @name rangebag
#' @export
rangebag.Raster <- function(x, p,
                     n_models = 100,
                     n_dim = 2,
                     sample_prop = 0.5,
                     limit_occur = TRUE, ...) {

  # Check that p has sufficient columns
  p <- as.data.frame(p)
  if (ncol(p) < 2) {
    stop("Rangebag p data should have at least 2 columns.", call. = FALSE)
  }

  # Make sure p coordinate columns are 'lon' and 'lat'
  if (!all(c("lon", "lat") %in% names(p))) {
    stop("Rangebag p coordinates should be 'lon' and 'lat'.", call. = FALSE)
  }

  # Project to lon/lat if necessary
  if (!raster::isLonLat(x)) {
    x <- raster::projectRaster(x, crs = "EPSG:4326")
  }

  # Select data from cells in x corresponding to each occurrence point in p
  fit_idx <- raster::cellFromXY(x, p[, c("lon", "lat")])
  if (limit_occur) { # limit to one occurrence per cell
    fit_idx <- unique(fit_idx)
  }
  fit_data <- stats::na.omit(as.data.frame(x[fit_idx]))
  fit_coords <- as.data.frame(raster::xyFromCell(x, fit_idx))
  names(fit_coords) <- c("lon", "lat")

  # Fit {n_models} convex hull models to data samples
  ch_models <- list()
  for (i in 1:n_models) {

    # Sample {sample_prop} data rows and {n_dim} variable columns
    rows <- sample(nrow(fit_data), ceiling(sample_prop*nrow(fit_data)),
                   replace = FALSE)
    vars <- colnames(fit_data)[sample.int(ncol(fit_data), size = n_dim,
                                          replace = FALSE)]
    sample_data <- fit_data[rows, vars, drop = FALSE]

    # Fit convex hulls
    if(n_dim == 1) {
      ch_models[[i]] <- sample_data[c(which.min(sample_data[,1]),
                                      which.max(sample_data[,1])),,
                                    drop = FALSE]
    } else {
      ch_models[[i]] <- sample_data[unique(as.vector(
        geometry::convhulln(sample_data, options = 'Pp'))),]
    }
  }

  # Return a "Rangebag" object with selected data and convex hull models
  return(methods::new("Rangebag",
                      method = "rangebag",
                      variables = names(x),
                      presence = fit_data,
                      coordinates = fit_coords,
                      ch_models = ch_models))
}

#' @name rangebag
#' @export
rangebag.SpatRaster <- function(x, p,
                                n_models = 100,
                                n_dim = 2,
                                sample_prop = 0.5,
                                limit_occur = TRUE, ...) {

  # Call the Raster version of the function (swap later)
  rangebag(raster::stack(x), p,
           n_models = n_models,
           n_dim = n_dim,
           sample_prop = sample_prop,
           limit_occur = limit_occur, ...)
}

#' @name rangebag
#' @export
rangebag.data.frame <- function(x, p,
                                n_models = 100,
                                n_dim = 2,
                                sample_prop = 0.5,
                                limit_occur = TRUE, ...) {

  # Check that x has sufficient columns
  if (ncol(x) < 3) {
    stop("Rangebag x data should have at least 3 columns.", call. = FALSE)
  }

  # Make sure x coordinate columns are 'lon' and 'lat'
  if (!all(c("lon", "lat") %in% names(x))) {
    stop("Rangebag x coordinates should be 'lon' and 'lat'.", call. = FALSE)
  }

  # Convert to a Raster*
  ordered_idx <- c(which(c("lon", "lat") %in% names(x)),
                   which(!names(x) %in% c("lon", "lat")))
  x <- raster::rasterFromXYZ(x[, ordered_idx], crs = "EPSG:4326")

  # Call the Raster version of the function
  rangebag(x, p,
           n_models = n_models,
           n_dim = n_dim,
           sample_prop = sample_prop,
           limit_occur = limit_occur, ...)
}
