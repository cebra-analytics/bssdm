#' Climatch SDM predict method
#'
#' The model prediction component of an implementation of the ABARES Climatch
#' species distribution modelling (SDM) method (ABARES, 2020).
#'
#' @param object A "Climatch" model S4 object containing slots:
#'   \describe{
#'     \item{\code{method}}{SDM method: "climatch".}
#'     \item{\code{algorithm}}{Algorithm: "euclidean" or
#'       "closest_standard_score".}
#'     \item{\code{variables}}{List of climate (or environmental) variable
#'       names.}
#'     \item{\code{sd}}{The standard deviation of each variable calculated via
#'       the climate data (\emph{x}) or the \emph{sd_data} when provided.}
#'     \item{\code{presence}}{The selected (nearest within range) climate data
#'       for each occurrence point.}
#'     \item{\code{coordinates}}{The coordinates for the selected climate
#'       data.}
#'     \item{\code{as_score}}{Indication of whether to generate a score 0-10
#'       (TRUE) or values 0-1 (FALSE).}
#'   }
#' @param x Climate (or environmental) data with corresponding model variables
#'   as a \code{terra::SpatRaster}, \code{raster::Raster*}, \code{data.frame},
#'   or \code{matrix}.
#' @param algorithm Optional (overriding) Climatch method algorithm selected
#'   from "euclidean" or "closest_standard_score".
#' @param sd_data Optional (overriding) \code{data.frame} for calculating the
#'   standard deviation for climate variable, or a \code{vector} of
#'   pre-calculated values.
#' @param as_score Optional (overriding) logical to indicate whether to
#'   generate a score 0-10 (TRUE) or values 0-1 (FALSE).
#' @param raw_output Logical to indicate whether to return raw predicted
#'   values (TRUE) or as an object (as per \emph{x}: FALSE). Default is NULL,
#'   returning either raw values or a spatial raster (as per \emph{x}).
#' @param filename Optional filename for writing spatial raster output (only).
#'   Default is "".
#' @param parallel_cores Optional number of cores available for parallel
#'   processing, thus enable parallel processing. Default is NULL (serial).
#' @param ... Additional parameters.
#' @return Predicted values as a raw vector or a \code{terra::SpatRaster},
#'   \code{raster::Raster*}, \code{data.frame}, or \code{matrix} (as per
#'   \emph{x}).
#' @references ABARES (2020). Climatch v2.0 User Manual. Canberra.
#'   \url{https://climatch.cp1.agriculture.gov.au/} Accessed: November 2021.
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @export
predict.Climatch <- function(object, x,
                             algorithm = NULL,
                             sd_data = NULL,
                             as_score = NULL,
                             raw_output = NULL,
                             filename = "",
                             parallel_cores = NULL, ...) {

  # Transpose presence data (for performance)
  t_y_source <- t(as.matrix(object@presence))

  # Ensure x data has matching object variables
  if (!(all(object@variables %in% names(x)) ||
        all(object@variables %in% colnames(x)))) {
    stop("Climatch predict x should have same variables as object.",
         call. = FALSE)
  }

  # Extract data values for object variables from x
  if (class(x)[1] %in% c("Raster", "RasterStack", "RasterBrick")) {
    y_target <- as.matrix(raster::as.data.frame(x, xy = TRUE, na.rm = TRUE))
    y_coords <- y_target[, c("x", "y")]
    y_target <- y_target[, object@variables]
  } else if (class(x)[1] == "SpatRaster") {
    y_target <- as.matrix(terra::as.data.frame(x, xy = TRUE, na.rm = TRUE))
    y_coords <- y_target[, c("x", "y")]
    y_target <- y_target[, object@variables]
  } else if (is.data.frame(x) || is.matrix(x)) {
    y_target <- as.matrix(x[, object@variables])
  }

  # Calculate standard deviations when required
  if (!is.null(sd_data)) {
    if (is.data.frame(sd_data)) {
      if (!all(object@variables %in% names(sd_data))) {
        stop("Climatch predict sd data should have same variables as object.",
             call. = FALSE)
      }
      sd_k = unlist(lapply(sd_data[, sd_data], sd))
    } else if (is.numeric(sd_data) &&
               length(sd_data) == length(object@variables)) {
      sd_k <- sd_data
    } else {
      stop("Climatch predict sd data should match object variables.",
           call. = FALSE)
    }
  } else {
    sd_k <- object@sd
  }

  # Select the matching algorithm
  if (is.null(algorithm)) {
    algorithm <- object@algorithm
  }
  if (algorithm == "closest_standard_score") {
    cut_values = sort(c(Inf, 3.9, 1.645, 1.285, 1.004, 0.845, 0.675, 0.525,
                        0.385, 0.255, 0.125, -Inf))
  } else if (algorithm != "euclidean") {
    stop(paste("Climatch predict algorithm should be 'euclidean' or",
               "'closest_standard_score'."),
         call. = FALSE)
  }

  # Run matching algorithm (from Climatch v2.0 User Manual)
  # for source site i, target site j, climate values y for k variables
  if (!is.null(parallel_cores) && parallel_cores > 1) {

    doParallel::registerDoParallel(cores = parallel_cores)
    d_j <- foreach(
      j = 1:nrow(y_target),
      .combine = "c",
      .errorhandling = c("stop"),
      .noexport = c("object", "x", "y_coords")) %dopar% {
        if (algorithm == "euclidean") {

          # d_j = floor{[1 - min_i(sqrt(1/k*sum_k((y_ik - y_jk)^2/sd_k^2))]*10}
          max(floor((1 - min(sqrt(
            colMeans((t_y_source - y_target[j,])^2/sd_k^2)
          )))*10 + 1e-6), 0) # compensate for float inaccuracies when flooring

        } else if (algorithm == "closest_standard_score") {

          # d_j = 11 - min_i(max_k(cut(sqrt((y_ik - y_jk)^2/sd_k^2),
          #                            {cut_values})))
          11 - min(matrixStats::colMaxs(
            array(.bincode(abs(t_y_source - y_target[j,])/sd_k, cut_values),
                  dim(t_y_source))))
        }
      }
    doParallel::stopImplicitCluster()

  } else { # serial

    d_j <- rep(0, nrow(y_target))
    for (j in 1:nrow(y_target)) {
      if (algorithm == "euclidean") {

        # d_j = floor{[1 - min_i(sqrt(1/k*sum_k((y_ik - y_jk)^2/sd_k^2))]*10}
        d_j[j] <- max(floor((1 - min(sqrt(
          colMeans((t_y_source - y_target[j,])^2/sd_k^2)
        )))*10 + 1e-6), 0) # compensate for float inaccuracies when flooring

      } else if (algorithm == "closest_standard_score") {

        # d_j = 11 - min_i(max_k(cut(sqrt((y_ik - y_jk)^2/sd_k^2),
        #                            {cut_values})))
        d_j[j] <- 11 - min(matrixStats::colMaxs(
          array(.bincode(abs(t_y_source - y_target[j,])/sd_k, cut_values),
                dim(t_y_source))))
      }
    }
  }

  # Return as score 0-10 (as per Climatch site) or 0-1
  if (is.null(as_score)) {
    as_score <- object@as_score
  }
  if (!as_score) {
    d_j <- d_j/10
  }

  # Return a raw vector, raster or data frame
  if (is.null(raw_output)) {
    if (is.data.frame(x) || is.matrix(x)) {
      raw_output <- TRUE
    } else {
      raw_output <- FALSE
    }
  }
  if (raw_output) {
    return(d_j)
  } else if (class(x)[1] %in% c("Raster", "RasterStack", "RasterBrick")) {
    output_rast <- raster::extend(
      raster::rasterFromXYZ(cbind(y_coords, predicted = d_j),
                            res = raster::res(x), crs = raster::crs(x)),
      raster::extent(x))
    if (filename != "") raster::writeRaster(output_rast, filename)
    return(output_rast)
  } else if (class(x)[1] == "SpatRaster") {
    return(terra::extend(
      terra::rast(cbind(y_coords, predicted = d_j),
                  type = "xyz", crs = raster::crs(x)), terra::ext(x),
      filename = filename))
  } else if (is.data.frame(x) || is.matrix(x)) {
    return(cbind(x[, which(!names(x) %in% object@variables)],
                 predicted = d_j))
  }
}
