#' Extrapolation Detection for Species Distribution Models
#'
#' Detect and quantify both univariate (Type 1) and multivariate (Type 2)
#' environmental novelty when projecting species distribution models.
#'
#' @param x Climate (or environmental) data as a `terra::SpatRaster`,
#'   \code{raster::Raster*}, `data.frame`, or `matrix` where each
#'   layer/column represents focal values of an environmental variable.
#' @param ref A `data.frame`, `matrix`, or `list` where each column/element
#'   represents reference values for an environmental variable (corresponding
#'   to those given in `x`).
#' @param mic Logical to indicate whether most influential covariates should be
#'   returned. If `TRUE`, the function returns which variables are most
#'   responsible for Type 1 and Type 2 novelty. Default is `FALSE`.
#' @param filename Optional filename for writing spatial raster output (i.e., 
#'   the `exdet` raster only). Default is "".
#' @param tol Tolerance value passed to `mahalanobis()` for matrix inversion.
#'   Default is `.Machine$double.eps`. See `?solve` for details.
#' @param ... Additional parameters.
#' @return If `x` is a `SpatRaster` or \code{Raster*} object, this function returns
#'   a list containing:
#'   - `exdet`: a `SpatRaster` layer giving the extrapolation detection scores
#'   where values < 0 indicate univariate novelty (Type 1), values between 0
#'   and 1 indicate analog conditions, and values > 1 indicate novel covariate
#'   combinations (Type 2);
#'   - `mic1`: a factor `SpatRaster` layer indicating which variable is most
#'   influential for Type 1 novelty. Value is "Not novel" where no Type 1
#'   novelty occurs (only included when `mic=TRUE`); and
#'   - `mic2`: a factor `SpatRaster` layer indicating which variable is most
#'   influential for Type 2 novelty. Value is "Not novel" where no Type 2
#'   novelty occurs (only included when `mic=TRUE`).
#'
#'   If `x` is a `data.frame` or `matrix`, the function will return a list as
#'   above, but with single layer `SpatRaster` objects replaced by vectors,
#'   and factor rasters replaced by factor vectors.
#' @details `exdet` implements the ExDet (Extrapolation Detection) method
#'   described in Mesgaran et al. (2014). It detects both novel univariate
#'   ranges (Type 1 novelty) and novel combinations of covariates (Type 2
#'   novelty) in the projection area compared to the reference area.
#'
#'   Type 1 novelty is assessed by comparing each environmental variable to its
#'   range in the reference data. Type 2 novelty is assessed using Mahalanobis
#'   distance to detect novel combinations of variables, even when individual
#'   variables are within their reference ranges.
#'
#'   The ExDet score prioritises univariate novelty. If one or more variables
#'   are outside their reference range, Type 2 novelty is not considered.
#'
#'   When `mic=TRUE`, the most influential covariates (MIC) are identified for
#'   both types of novelty. MIC1 identifies which variable contributes most to
#'   Type 1 novelty (most outside its reference range). MIC2 identifies which
#'   variable contributes most to Type 2 novelty (via leave-one-out influence
#'   analysis). Areas with no novelty are labeled "Not novel" in the MIC outputs.
#' @keywords extrapolation, novelty, exdet, sdm
#' @references Mesgaran, M. B., Cousens, R. D., and Webber, B. L. (2014). Here
#'   be dragons: a tool for quantifying novelty due to covariate range and
#'   correlation change when projecting species distribution models.
#'   \emph{Diversity and Distributions}, 20(10): 1147-1159.
#'   \doi{10.1111/ddi.12209}
#' @export
#' @examples
#' \dontrun{
#' library(geodata)
#' library(terra)
#' library(tmap)
#' bio <- worldclim_global("bio", res = 10, path = tempdir())
#' aus <- gadm("AUS", level = 0, resolution = 2, path = tempdir())
#' occ <- spatSample(aus, size = 100, method = "random")
#' ref <- terra::extract(bio, occ)
#' ex <- exdet(bio, ref, mic = TRUE)
#'
#' # Plot outputs
#' tmap::tm_shape(ex$exdet) +
#'   tmap::tm_raster(
#'     col.scale = tmap::tm_scale_continuous(values = "matplotlib.rd_bu", midpoint = 0),
#'     col.legend = tmap::tm_legend(
#'       position = tmap::tm_pos_out("right", "center"),
#'       title = "ExDet"
#'     )
#'   ) +
#'   tmap::tm_credits(
#'     paste(
#'       "Positive values indicate covariate correlations.",
#'       "Negative values indicate novel univariate ranges.",
#'       "Values around zero indicate analog conditions.",
#'       sep = "\n"
#'     ),
#'     position = tmap::tm_pos_out("center", "bottom")
#'   ) +
#'   tmap::tm_title("Extrapolation Detection")
#'
#' tmap::tm_shape(ex$mic1) +
#'   tmap::tm_raster(
#'     col.scale = tmap::tm_scale_categorical(
#'       values = "20",
#'       n.max = nrow(levels(ex$mic1)[[1]])
#'     ),
#'     col.legend = tmap::tm_legend(title = "Most influential\ncovariate (MIC1)"),
#'     col.chart = tmap::tm_chart_donut()
#'   ) +
#'   tmap::tm_title("Type 1 novelty: Most influential covariate")
#'
#' tmap::tm_shape(ex$mic2) +
#'   tmap::tm_raster(
#'     col.scale = tmap::tm_scale_categorical(
#'       values = "20",
#'       n.max = nrow(levels(ex$mic2)[[1]])
#'     ),
#'     col.legend = tmap::tm_legend(title = "Most influential\ncovariate (MIC2)"),
#'     col.chart = tmap::tm_chart_donut()
#'   ) +
#'   tmap::tm_title("Type 2 novelty: Most influential covariate")
#' }
exdet <- function(
  x,
  ref,
  mic = FALSE,
  filename = "",
  tol = .Machine$double.eps,
  ...
) {
  UseMethod("exdet")
}

#' @name exdet
#' @export
exdet.Raster <- function(
  x,
  ref,
  mic = FALSE,
  filename = "",
  tol = .Machine$double.eps,
  ...
) {
  # Convert to SpatRaster
  x <- terra::rast(x)

  # Call the terra version of the function
  exdet(x, ref, mic = mic, filename = filename, tol = tol, ...)
}

#' @name exdet
#' @export
exdet.SpatRaster <- function(
  x,
  ref,
  mic = FALSE,
  filename = "",
  tol = .Machine$double.eps,
  ...
) {
  # Prepare reference data and check variables
  ref_prep <- .prep_exdet_ref(ref, names(x), tol)
  cols <- ref_prep$cols
  n_vars <- length(cols)
  ref <- ref_prep$ref
  ref_min <- ref_prep$ref_min
  ref_max <- ref_prep$ref_max
  ref_mean <- ref_prep$ref_mean
  ref_cov <- ref_prep$ref_cov
  mah_max <- ref_prep$mah_max

  # Subset x to common columns
  x <- x[[cols]]

  # Set up block processing
  x_blocks <- terra::blocks(x, ...)
  n_blocks <- x_blocks$n
  terra::readStart(x)

  # Initialise output rasters and start writing
  out_exdet <- setNames(terra::init(x[[1]], NA), "exdet")
  terra::varnames(out_exdet) <- "exdet"
  invisible(terra::writeStart(out_exdet, filename = filename, ...))

  if (isTRUE(mic)) {
    out_mic1 <- terra::init(x[[1]], NA)
    terra::varnames(out_mic1) <- "mic1"
    invisible(terra::writeStart(out_mic1, filename = "", ...))
    out_mic2 <- terra::init(x[[1]], NA)
    terra::varnames(out_mic2) <- "mic2"
    invisible(terra::writeStart(out_mic2, filename = "", ...))
  }

  # Process each block
  for (b in 1:n_blocks) {
    # Read block
    x_data <- terra::readValues(
      x,
      row = x_blocks$row[b],
      nrows = x_blocks$nrows[b],
      mat = TRUE
    )

    # Find non-NA cells
    x_idx <- which(stats::complete.cases(x_data))

    if (length(x_idx) > 0) {
      x_vals <- x_data[x_idx, , drop = FALSE]

      # Calculate ExDet for this block
      exdet_result <- .calculate_exdet(
        x_vals,
        ref_min,
        ref_max,
        ref_mean,
        ref_cov,
        mah_max,
        cols,
        mic,
        tol
      )

      # Prepare output data for this block
      output_exdet <- rep(NA, nrow(x_data))
      output_exdet[x_idx] <- exdet_result$exdet

      if (isTRUE(mic)) {
        output_mic1 <- rep(NA, nrow(x_data))
        output_mic1[x_idx] <- exdet_result$mic1
        output_mic2 <- rep(NA, nrow(x_data))
        output_mic2[x_idx] <- exdet_result$mic2
      }
    } else {
      # All NA block
      output_exdet <- rep(NA, nrow(x_data))
      if (isTRUE(mic)) {
        output_mic1 <- rep(NA, nrow(x_data))
        output_mic2 <- rep(NA, nrow(x_data))
      }
    }

    # Write block outputs
    terra::writeValues(
      out_exdet,
      v = output_exdet,
      start = x_blocks$row[b],
      nrows = x_blocks$nrows[b]
    )
    if (isTRUE(mic)) {
      terra::writeValues(
        out_mic1,
        v = output_mic1,
        start = x_blocks$row[b],
        nrows = x_blocks$nrows[b]
      )
      terra::writeValues(
        out_mic2,
        v = output_mic2,
        start = x_blocks$row[b],
        nrows = x_blocks$nrows[b]
      )
    }
  }

  # Close block reading & writing
  terra::readStop(x)
  invisible(terra::writeStop(out_exdet))
  if (isTRUE(mic)) {
    invisible(terra::writeStop(out_mic1))
    invisible(terra::writeStop(out_mic2))
  }

  # Convert mic1 and mic2 to factors with variable names
  if (isTRUE(mic)) {
    out_mic1 <- suppressWarnings(terra::as.factor(out_mic1))
    levs <- data.frame(ID = 0:n_vars, variable = c("Not novel", cols))
    levels(out_mic1)[[1]] <- levs
    out_mic2 <- suppressWarnings(terra::as.factor(out_mic2))
    levels(out_mic2)[[1]] <- levs
  }

  # Return results
  if (isTRUE(mic)) {
    return(list(
      exdet = out_exdet,
      mic1 = out_mic1,
      mic2 = out_mic2
    ))
  } else {
    return(list(
      exdet = out_exdet
    ))
  }
}

#' @name exdet
#' @export
exdet.data.frame <- function(
  x,
  ref,
  mic = FALSE,
  tol = .Machine$double.eps,
  ...
) {
  # Prepare reference data and check variables
  ref_prep <- .prep_exdet_ref(ref, names(x), tol)
  cols <- ref_prep$cols
  ref <- ref_prep$ref
  ref_min <- ref_prep$ref_min
  ref_max <- ref_prep$ref_max
  ref_mean <- ref_prep$ref_mean
  ref_cov <- ref_prep$ref_cov
  mah_max <- ref_prep$mah_max

  # Subset x to common columns
  x <- x[, cols, drop = FALSE]

  # Calculate ExDet
  exdet_result <- .calculate_exdet(
    x,
    ref_min,
    ref_max,
    ref_mean,
    ref_cov,
    mah_max,
    cols,
    mic,
    tol
  )

  # Return results
  if (isTRUE(mic)) {
    # Convert MIC to factors with "Not novel" for 0 values
    mic1_factor <- factor(
      exdet_result$mic1,
      levels = 0:length(cols),
      labels = c("Not novel", cols)
    )
    mic2_factor <- factor(
      exdet_result$mic2,
      levels = 0:length(cols),
      labels = c("Not novel", cols)
    )
    return(list(
      exdet = exdet_result$exdet,
      mic1 = mic1_factor,
      mic2 = mic2_factor
    ))
  } else {
    return(list(
      exdet = exdet_result$exdet
    ))
  }
}

#' @name exdet
#' @export
exdet.matrix <- function(x, ref, mic = FALSE, tol = .Machine$double.eps, ...) {
  exdet(as.data.frame(x), ref, mic = mic, tol = tol, ...)
}

# Internal helper function for robust Mahalanobis distance calculation
.robust_mahalanobis <- function(x, center, cov, tol, warn = TRUE) {
  result <- tryCatch(
    mahalanobis(x = x, center = center, cov = cov, tol = tol),
    error = function(e) {
      if (grepl("singular", e$message)) {
        # Try increasingly SMALLER tolerances (more lenient)
        # Small tol = only singular if VERY singular (lenient)
        # Large tol = singular if moderately singular (strict)
        tol_values <- c(1e-10, 1e-12, 1e-14, 1e-16, 1e-18, 1e-20, 1e-25, 1e-30)
        tol_values <- tol_values[tol_values < tol]

        if (length(tol_values) == 0) {
          # User already provided very small tolerance, can't go smaller
          stop(
            sprintf(
              "Unable to compute Mahalanobis distance with tol=%.2e: covariance matrix is too singular. ",
              tol
            ),
            "Try reducing the number of variables or increasing sample size.",
            call. = FALSE
          )
        }

        for (new_tol in tol_values) {
          result <- tryCatch(
            mahalanobis(x = x, center = center, cov = cov, tol = new_tol),
            error = function(e2) NULL
          )
          if (!is.null(result)) {
            if (warn) {
              warning(
                sprintf(
                  "Covariance matrix is near-singular. Decreased tolerance from %.2e to %.2e. ",
                  tol,
                  new_tol
                ),
                "Consider reducing the number of variables or increasing sample size.",
                call. = FALSE
              )
            }
            return(result)
          }
        }
        # If all tolerances fail, give up
        stop(
          "Unable to compute Mahalanobis distance: covariance matrix is too singular. ",
          "Try reducing the number of variables or increasing sample size.",
          call. = FALSE
        )
      } else {
        stop(e)
      }
    }
  )
  return(result)
}

# Internal helper function to validate and preprocess reference data
.prep_exdet_ref <- function(ref, x_vars, tol) {
  # Ensure ref is a data frame
  if (!methods::is(ref, 'data.frame')) {
    ref <- as.data.frame(ref)
  }
  ref <- stats::na.omit(ref)

  # Check for overlapping variables
  ref_vars <- names(ref)
  x_vars_missing <- setdiff(ref_vars, x_vars)
  ref_vars_missing <- setdiff(x_vars, ref_vars)

  if (length(x_vars_missing) > 0) {
    warning(
      paste0(
        "The following variables are missing from x and will be ignored: ",
        paste(x_vars_missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (length(ref_vars_missing) > 0) {
    warning(
      paste0(
        "The following variables are missing from ref and will be ignored: ",
        paste(ref_vars_missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  cols <- intersect(x_vars, ref_vars)
  if (length(cols) == 0) {
    stop("No variables in common between x and ref.", call. = FALSE)
  }

  # Subset to common columns
  ref <- ref[, cols, drop = FALSE]

  # Pre-calculate reference statistics
  ref_min <- apply(ref, 2, min)
  ref_max <- apply(ref, 2, max)
  ref_mean <- colMeans(ref)
  ref_cov <- var(ref)

  # Calculate maximum Mahalanobis distance in reference data
  mah_ref <- .robust_mahalanobis(
    x = ref,
    center = ref_mean,
    cov = ref_cov,
    tol = tol,
    warn = TRUE
  )
  mah_max <- max(mah_ref[is.finite(mah_ref)])

  return(list(
    cols = cols,
    ref = ref,
    ref_min = ref_min,
    ref_max = ref_max,
    ref_mean = ref_mean,
    ref_cov = ref_cov,
    mah_max = mah_max
  ))
}

# Internal function to calculate ExDet statistics
.calculate_exdet <- function(
  x_vals,
  ref_min,
  ref_max,
  ref_mean,
  ref_cov,
  mah_max,
  cols,
  mic,
  tol
) {
  # Handle empty input
  if (nrow(x_vals) == 0) {
    return(list(
      exdet = numeric(0),
      mic1 = integer(0),
      mic2 = integer(0)
    ))
  }

  # Convert to data frame for mapply
  x_df <- as.data.frame(x_vals)
  names(x_df) <- cols

  # Calculate Type 1 novelty (univariate)
  ud <- mapply(
    function(p, ref_min, ref_max) {
      x <- findInterval(p, c(ref_min, ref_max))
      ifelse(
        x == 0,
        (p - ref_min) / (ref_max - ref_min),
        ifelse(x == 1, 0, (ref_max - p) / (ref_max - ref_min))
      )
    },
    x_df,
    ref_min,
    ref_max
  )

  # Ensure ud is a matrix even for single variable
  if (!is.matrix(ud)) {
    ud <- matrix(ud, ncol = 1)
  }

  nt1 <- rowSums(ud) # Type 1 novelty score

  # Calculate Type 2 novelty (multivariate)
  mah_p <- .robust_mahalanobis(
    x = x_vals,
    center = ref_mean,
    cov = ref_cov,
    tol = tol,
    warn = FALSE
  )
  nt2 <- mah_p / mah_max # Type 2 novelty score

  # Combine: prioritise Type 1 over Type 2
  exdet_score <- ifelse(nt1 == 0, nt2, nt1)

  # Calculate most influential covariates if requested
  if (isTRUE(mic)) {
    n_vars <- length(cols)

    # MIC1: Most influential covariate for Type 1 novelty
    if (ncol(ud) == 1) {
      mic1 <- rep(1L, nrow(x_vals))
    } else {
      mic1 <- integer(nrow(x_vals))
      type1_idx <- which(nt1 != 0)
      if (length(type1_idx) > 0) {
        mic1[type1_idx] <- apply(ud[type1_idx, , drop = FALSE], 1, which.min)
      }
    }

    # MIC2: Most influential covariate for Type 2 novelty
    mic2 <- integer(nrow(x_vals))
    type2_idx <- which(nt2 > 1)

    if (length(type2_idx) > 0 && n_vars > 1) {
      # For each variable, calculate influence on Type 2 novelty
      d_influence <- matrix(NA, nrow = length(type2_idx), ncol = n_vars)

      for (i in seq_len(n_vars)) {
        # Calculate Mahalanobis distance without variable i
        x_without_i <- x_vals[type2_idx, -i, drop = FALSE]
        mean_without_i <- ref_mean[-i]
        cov_without_i <- ref_cov[-i, -i, drop = FALSE]

        mah_without_i <- .robust_mahalanobis(
          x = x_without_i,
          center = mean_without_i,
          cov = cov_without_i,
          tol = tol,
          warn = FALSE
        )

        # Influence is the relative change in distance when variable is removed
        d_influence[, i] <- (mah_p[type2_idx] - mah_without_i) /
          mah_p[type2_idx]
      }

      # Most influential is the one that causes the largest decrease when removed
      mic2[type2_idx] <- apply(d_influence, 1, which.max)
    }

    return(list(
      exdet = exdet_score,
      mic1 = mic1,
      mic2 = mic2
    ))
  } else {
    return(list(
      exdet = exdet_score
    ))
  }
}
