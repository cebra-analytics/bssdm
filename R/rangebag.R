#' Range bagging SDM model building method
#'
#' The model building component of an implementation of the range bagging
#' species distribution modelling (SDM) method (Drake, 2015).
#'
#' @param x Climate (or environmental) data as a \code{terra::SpatRaster} or
#'   \code{raster::Raster*} (with any CRS), or a \code{data.frame} with
#'   WGS84 \emph{lon} and \emph{lat} columns.
#' @param p Species occurrence data as a \code{data.frame} (or \code{matrix})
#'   with WGS84 \emph{lon} and \emph{lat} columns. Any points outside the extent
#'   of \code{x} will be ignored.
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
#' @references Drake, J. M. (2015). Range bagging: a new method for ecological
#'   niche modelling from presence-only data.
#'   \emph{Journal of the Royal Society Interface}, 12(107), 20150086.
#'   \doi{10.1098/rsif.2015.0086}
#' @include Rangebag-class.R
#' @export
rangebag <- function(x, p,
                     n_models = 100,
                     n_dim = 2,
                     sample_prop = 0.5,
                     limit_occur = TRUE, ...) {
  if (n_dim <= 0) {
    stop("n_dim must be greater than 0.", call. = FALSE)
  }
  if (n_models <= 0) {
    stop("n_models must be greater than 0.", call. = FALSE)
  }
  if (sample_prop <= 0 || sample_prop > 1) {
    stop(
      "sample_prop must be greater than 0 and less than or equal to 1.", 
      call. = FALSE
    )
  }
  UseMethod("rangebag")
}

#' @name rangebag
#' @export
rangebag.Raster <- function(x, p,
                     n_models = 100,
                     n_dim = 2,
                     sample_prop = 0.5,
                     limit_occur = TRUE, ...) {

  # Call the terra version of the function
  rangebag(terra::rast(x), p,
           n_models = n_models,
           n_dim = n_dim,
           sample_prop = sample_prop,
           limit_occur = limit_occur, ...)
}

#' @name rangebag
#' @export
rangebag.SpatRaster <- function(x, p,
                                n_models = 100,
                                n_dim = 2,
                                sample_prop = 0.5,
                                limit_occur = TRUE, ...) {

  # Ensure the number of layers are consistent with the number of dimensions
  if (terra::nlyr(x) < n_dim) {
    warning("Rangebag x data has fewer variables than n_dim.", call. = FALSE)
    n_dim <- terra::nlyr(x) # adjust
  }

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
  if (!terra::is.lonlat(x)) {
    x <- terra::project(x, "EPSG:4326")
  }

  # Select data from cells in x corresponding to each occurrence point in p
  fit_idx <- terra::cellFromXY(x, p[, c("lon", "lat")])
  if(any(is.na(fit_idx))) { # Remove any points outside the extent of the data
    n_outside <- length(which(is.na(fit_idx)))
    warning(sprintf(
      "%i point%s outside the extent of the environmental data and will be ignored",
      n_outside, if(n_outside > 1) "s are" else " is"
    ), call. = FALSE)
    fit_idx <- fit_idx[!is.na(fit_idx)]
  }
  if (limit_occur) { # limit to one occurrence per cell
    fit_idx <- unique(fit_idx)
  }
  fit_data <- stats::na.omit(
    cbind(as.data.frame(terra::xyFromCell(x, fit_idx)),
          as.data.frame(x[fit_idx])))
  fit_coords <- fit_data[,1:2]
  names(fit_coords) <- c("lon", "lat")
  fit_data <- fit_data[, -(1:2), drop = FALSE]

  # Fit {n_models} convex hull models to data samples
  ch_models <- list()
  for (i in 1:n_models) {

    # Sample {sample_prop} data rows and {n_dim} variable columns
    vars <- colnames(fit_data)[sample.int(ncol(fit_data), size = n_dim,
                                          replace = FALSE)]
    rows <- sample(nrow(fit_data), ceiling(sample_prop*nrow(fit_data)),
                   replace = FALSE)
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

  # Convert to a terra::SpatRaster
  ordered_idx <- c(match(c("lon", "lat"), names(x)),
                   which(!names(x) %in% c("lon", "lat")))
  x <- terra::rast(x[, ordered_idx], crs = "EPSG:4326")

  # Call the terra version of the function
  rangebag(x, p,
           n_models = n_models,
           n_dim = n_dim,
           sample_prop = sample_prop,
           limit_occur = limit_occur, ...)
}
