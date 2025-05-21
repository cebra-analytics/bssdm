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
#' @param parallel_cores Optional number of cores available for parallel
#'   processing, thus enable parallel processing. Default is NULL (serial).
#' @param debug Output additional debug information (memory block info, etc.).
#'   Default is FALSE.
#' @param ... Additional parameters.
#' @return Predicted values as a raw vector or a \code{terra::SpatRaster},
#'   \code{raster::Raster*}, \code{data.frame}, or \code{matrix} (as per
#'   \emph{x}).
#' @references Drake, J. M. (2015). Range bagging: a new method for ecological
#'   niche modelling from presence-only data.
#'   \emph{Journal of the Royal Society Interface}, 12(107), 20150086.
#'   \doi{10.1098/rsif.2015.0086}
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @export
predict.Rangebag <- function(object, x,
                             raw_output = NULL,
                             filename = "",
                             parallel_cores = NULL,
                             debug = FALSE, ...) {

  # Ensure x data has matching object variables
  if (!(all(object@variables %in% names(x)) ||
        all(object@variables %in% colnames(x)))) {
    stop("Rangebag predict x should have same variables as object.",
         call. = FALSE)
  }

  # Convert x raster to terra
  x_is_raster <- FALSE
  if (class(x)[1] %in% c("Raster", "RasterStack", "RasterBrick")) {
    x_is_raster <- TRUE
    x <- terra::rast(x)
  }

  # Handle terra raster data in blocks
  if (class(x)[1] == "SpatRaster") {
    n_cores <- ifelse(is.numeric(parallel_cores), parallel_cores, 1)
    x_blocks <- terra::blocks(x, n = 8*n_cores)
    n_blocks <- x_blocks$n
    terra::readStart(x)
    output_rast <- terra::rast(x[[1]])
    invisible(terra::writeStart(output_rast, filename = filename, ...))
  } else {
    n_blocks <- 1
  }
  if (debug) {
    message(paste("Range bagging predict writing output in", n_blocks,
                  "blocks\n"))
  }

  # Calculate the proportion of models that fit each location
  n_models <- length(object@ch_models)
  n_dim <- ncol(object@ch_models[[1]])
  if (n_dim > 1) {
    hulls <- lapply(object@ch_models, function(model) {
      geometry::convhulln(model, options='Pp')})
  }
  for (b in 1:n_blocks) {

    # Extract (block) data values for object variables from x
    if (class(x)[1] == "SpatRaster") {
      x_data <- terra::readValues(x,
                                  row = x_blocks$row[b],
                                  nrows = x_blocks$nrows[b],
                                  col = 1,
                                  ncols = terra::ncol(x),
                                  mat = TRUE)
      output_data <- rep(NA, nrow(x_data))
      x_data <- x_data[, object@variables, drop = FALSE]
      x_idx <- which(as.logical(rowSums(!is.na(x_data))))
      x_data <- x_data[x_idx,, drop = FALSE]
    } else if (is.data.frame(x) || is.matrix(x)) {
      x_data <- as.matrix(x[, object@variables, drop = FALSE])
    }

    # Count the proportion of convex hull model fits for each x data row/cell
    if (nrow(x_data) > 0) {

      # Parallel
      if (!is.null(parallel_cores) && parallel_cores > 1) {

        # Split the target data rows into groups
        n_groups <- min(nrow(x_data), parallel_cores*100)
        row_groups <- split(1:nrow(x_data),
                            cut(seq_along(1:nrow(x_data)),
                                breaks = n_groups, labels = FALSE))

        doParallel::registerDoParallel(cores = parallel_cores)
        ch_counts <- foreach(
          g = 1:n_groups,
          .combine = "c",
          .errorhandling = c("stop"),
          .export = c("x_data", "object", "hulls", "row_groups"),
          .noexport = c("x", "output_rast", "output_data", "x_idx")
        ) %dopar% {
          ch_counts_g <- numeric(length(row_groups[[g]]))
          idx <- row_groups[[g]]
          for (i in 1:n_models) {
            vars <- colnames(object@ch_models[[i]])
            if (n_dim == 1) {
              data_in_ch <- (x_data[idx, vars] >= object@ch_models[[i]][1, 1] &
                               x_data[idx, vars] <= object@ch_models[[i]][2, 1])
            } else {
              data_in_ch <- geometry::inhulln(hulls[[i]],
                                              x_data[idx, vars, drop = FALSE])
            }
            ch_counts_g <- ch_counts_g + data_in_ch
          }
          ch_counts_g
        }
        doParallel::stopImplicitCluster()

      } else { # serial
        ch_counts <- numeric(nrow(x_data))
        for (i in 1:n_models) {
          vars <- colnames(object@ch_models[[i]])
          if (n_dim == 1) {
            data_in_ch <- (x_data[, vars] >= object@ch_models[[i]][1, 1] &
                             x_data[, vars] <= object@ch_models[[i]][2, 1])
          } else {
            data_in_ch <- geometry::inhulln(hulls[[i]],
                                            x_data[, vars, drop = FALSE])
          }
          ch_counts <- ch_counts + data_in_ch
        }
      }

      if (class(x)[1] == "SpatRaster") {
        output_data[x_idx] <- ch_counts/n_models
      }
    }

    # Write output block
    if (class(x)[1] == "SpatRaster") {
      terra::writeValues(output_rast, v = output_data,
                         start = x_blocks$row[b], nrows = x_blocks$nrows[b])
    }
  }

  # Close block reading & writing
  if (class(x)[1] == "SpatRaster") {
    terra::readStop(x)
    invisible(terra::writeStop(output_rast))
  }

  # Return the count fraction as a raw vector, raster or data frame
  if (is.null(raw_output)) {
    if (is.data.frame(x) || is.matrix(x)) {
      raw_output <- TRUE
    } else {
      raw_output <- FALSE
    }
  }
  if (raw_output && n_blocks == 1) {
    return(ch_counts/n_models)
  } else if (x_is_raster) {
    return(raster::raster(output_rast))
  } else if (class(x)[1] == "SpatRaster") {
    return(output_rast)
  } else if (is.data.frame(x) || is.matrix(x)) {
    return(cbind(x[, which(!names(x) %in% object@variables)],
                 predicted = ch_counts/n_models))
  }
}
