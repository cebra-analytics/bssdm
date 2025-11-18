# This file contains code adapted from ecospat::ecospat.boyce
# Original code: Copyright (C) ecospat authors (Blaise Petitpierre,
#                Frank Breiner, Flavien Collart)
# Original licence: GPL-3.0
# Original source: https://github.com/ecospat/ecospat
#
# Modifications Copyright (C) 2025 bssdm authors (Sean Haythorne, James Camac,
#                                                 John Baumgartner)
#
# Changes made to the original code:
# - Restructured as S3 generic with methods for default, Rangebag, and Climatch
# - Added methods to work with Rangebag and Climatch model objects
# - Parameter names changed to snake_case (window_width, pe_plot, rm_duplicate)
# - Return value is a 'boyce' class object with additional 'indices' field
# - Return field names changed to snake_case (f_ratio, hs, cor)
# - Added comprehensive input validation with helpful error messages
# - Improved NA handling with explicit checks and validation
# - Changed error handling: errors instead of returning NA for insufficient data
# - Added print.boyce() and plot.boyce() methods
# - Changed default pe_plot to FALSE (was TRUE in original)
# - Updated documentation to roxygen2 format
# - Maintained core algorithm from original ecospat.boyce implementation
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public Licence as published by
# the Free Software Foundation, either version 3 of the Licence, or
# (at your option) any later version.

#' Calculate Boyce index for model validation
#'
#' Calculates the Boyce index as described in Hirzel et al. (2006) for
#' evaluating species distribution models using presence-only data.
#'
#' The Boyce index measures how much model predictions differ from a random
#' distribution of observed presences across prediction classes. It is a
#' presence-only metric that does not require absence data. Values range from
#' -1 to 1, where positive values indicate model predictions consistent with
#' the distribution of presences, values close to zero indicate predictions no
#' different than random, and negative values indicate counter predictions.
#'
#' This implementation is adapted from `ecospat.boyce` in the
#' [ecospat package](https://github.com/ecospat/ecospat), originally authored
#' by Blaise Petitpierre and Frank Breiner with updates by Flavien Collart
#' (GPL-3.0 licence).
#'
#' @param fit A vector or `terra::SpatRaster` containing predicted
#'   suitability values, or a fitted model object (`Rangebag` or
#'   `Climatch`).
#' @param ... Additional arguments passed to methods.
#' @return An object of class `boyce`, which is a list containing:
#'   - `cor`: Boyce index (correlation coefficient).
#'   - `f_ratio`: Predicted-to-expected ratio for each class interval.
#'   - `hs`: Mean habitat suitability values for each interval.
#'   - `indices`: Indices of values used in correlation (after removing
#'     successive duplicates if `rm_duplicate = TRUE`).
#' @references Hirzel, A. H., Le Lay, G., Helfer, V., Randin, C., & Guisan, A.
#'   (2006). Evaluating the ability of habitat suitability models to predict
#'   species presences. _Ecological Modelling_, 199(2), 142-152.
#'   \doi{10.1016/j.ecolmodel.2006.05.017}
#' @references Broennimann O, Di Cola V, Guisan A (2025). _ecospat: Spatial
#'   Ecology Miscellaneous Methods_. R package version 4.1.2,
#'   <https://CRAN.R-project.org/package=ecospat>.
#' @examples
#' \dontrun{
#' # With a fitted Rangebag model
#' model <- rangebag(presence_data, env_data)
#' validation_coords <- data.frame(x = c(...), y = c(...))
#' boyce_result <- boyce(model, env_data, validation_coords)
#'
#' # With prediction raster and presence values
#' predictions <- predict(model, env_data)
#' presence_values <- terra::extract(predictions, validation_coords)
#' boyce_result <- boyce(predictions, presence_values)
#' }
#' @export
boyce <- function(fit, ...) {
  UseMethod("boyce")
}

#' @describeIn boyce Default method for vectors and SpatRaster
#' @param obs A vector containing predicted suitability values at validation
#'   points (presence records), or a two-column matrix/data.frame of
#'   xy-coordinates (if `fit` is a `SpatRaster`).
#' @param nclass Number of classes (integer) or vector with class thresholds.
#'   If `nclass = 0` (default), the Boyce index is calculated with a
#'   moving window (see `window_width` and `res` parameters).
#' @param window_width Width of the moving window. Default is `"default"`,
#'   which uses 1/10 of the suitability range.
#' @param res Resolution of the moving window (number of bins). Default is 100.
#' @param pe_plot Logical. If `TRUE` (default), plots the
#'   predicted-to-expected ratio along the suitability classes.
#' @param rm_duplicate Logical. If `TRUE` (default), successive duplicated
#'   values are excluded from the correlation calculation.
#' @param method Correlation method used to compute the Boyce index. Default is
#'   `"spearman"`.
#' @export
boyce.default <- function(
  fit,
  obs,
  nclass = 0,
  window_width = "default",
  res = 100,
  pe_plot = FALSE,
  rm_duplicate = TRUE,
  method = "spearman",
  ...
) {
  # Validate fit
  if (!is.numeric(fit) && !inherits(fit, "SpatRaster")) {
    stop(
      "fit must be either a numeric vector or a SpatRaster object. ",
      "Received: ",
      class(fit)[1],
      call. = FALSE
    )
  }

  # Validate obs based on fit type
  if (inherits(fit, "SpatRaster")) {
    # When fit is a SpatRaster, obs can be coordinates or values
    if (!is.numeric(obs) && !is.data.frame(obs) && !is.matrix(obs)) {
      hint <- if (inherits(obs, "SpatVector")) {
        paste0(
          "\nTo extract coordinates from a SpatVector, use: ",
          "terra::crds(obs)"
        )
      } else {
        ""
      }
      stop(
        "When fit is a SpatRaster, obs must be one of:\n",
        "  - A numeric vector of predicted values at validation points\n",
        "  - A data.frame or matrix of xy-coordinates\n",
        "Received: ",
        class(obs)[1],
        hint,
        call. = FALSE
      )
    }
  } else {
    # When fit is numeric, obs must be numeric
    if (!is.numeric(obs)) {
      stop(
        "When fit is a numeric vector, obs must also be a numeric vector. ",
        "Received: ",
        class(obs)[1],
        call. = FALSE
      )
    }
  }

  # Internal function calculating predicted-to-expected ratio for each
  # class-interval
  boycei <- function(interval, obs, fit) {
    pi <- sum(as.numeric(obs >= interval[1] & obs <= interval[2])) / length(obs)
    ei <- sum(as.numeric(fit >= interval[1] & fit <= interval[2])) / length(fit)
    round(pi / ei, 10)
  }

  if (inherits(fit, "SpatRaster")) {
    if (is.data.frame(obs) || is.matrix(obs)) {
      obs <- terra::extract(fit, as.data.frame(obs), ID = FALSE)
      obs <- as.numeric(obs[, 1]) # need to be a vector
    }
    fit <- terra::values(fit, na.rm = TRUE)
    fit <- as.numeric(fit) # need to be a vector
  }

  # Remove NAs from both fit and obs
  fit <- fit[!is.na(fit)]
  obs <- obs[!is.na(obs)]

  # Check we have sufficient data after removing NAs
  if (length(fit) == 0) {
    stop("fit contains no valid (non-NA) values.", call. = FALSE)
  }
  if (length(obs) == 0) {
    stop("obs contains no valid (non-NA) values.", call. = FALSE)
  }
  if (length(obs) < 2) {
    stop(
      "obs must contain at least 2 validation points. Found: ",
      length(obs),
      call. = FALSE
    )
  }

  mini <- min(fit, obs)
  maxi <- max(fit, obs)

  if (length(nclass) == 1) {
    if (nclass == 0) {
      # moving window
      if (window_width == "default") {
        window_width <- (max(fit) - min(fit)) / 10
      }
      vec_mov <- seq(
        from = mini,
        to = maxi - window_width,
        by = (maxi - mini - window_width) / res
      )
      vec_mov[res + 1] <- vec_mov[res + 1] + 1
      # ^ Trick to avoid error with closed interval in R
      interval <- cbind(vec_mov, vec_mov + window_width)
    } else {
      # window based on nb of class
      vec_mov <- seq(from = mini, to = maxi, by = (maxi - mini) / nclass)
      interval <- cbind(vec_mov, c(vec_mov[-1], maxi))
    }
  } else {
    # user defined window
    vec_mov <- c(mini, sort(nclass[!nclass > maxi | nclass < mini]))
    interval <- cbind(vec_mov, c(vec_mov[-1], maxi))
  }

  f <- apply(interval, 1, boycei, obs, fit)
  to_keep <- which(f != "NaN") # index to keep no NaN data
  f <- f[to_keep]

  # Check we have enough bins to calculate correlation
  if (length(f) < 2) {
    stop(
      "Insufficient data to calculate Boyce index. ",
      "Only ",
      length(f),
      " bin(s) with data. ",
      "Need at least 2 bins. ",
      "Try reducing nclass or adjusting window.w/res parameters.",
      call. = FALSE
    )
  }

  # Identify which bins to use in correlation
  idx <- seq_along(f)
  if (rm_duplicate == TRUE) {
    idx <- c(seq_along(f))[f != c(f[-1], TRUE)]
    # ^ index to remove successive duplicates
  }

  # Calculate correlation (Boyce index)
  b <- cor(f[idx], vec_mov[to_keep][idx], method = method)
  hs <- apply(interval, 1, sum) / 2
  # ^ mean habitat suitability in the moving window
  if (length(nclass) == 1 && nclass == 0) {
    hs[length(hs)] <- hs[length(hs)] - 1
    # ^ Correction of the 'trick' to deal with closed interval
  }
  hs <- hs[to_keep] # exclude the NaN
  if (pe_plot == TRUE) {
    plot(b)
  }

  result <- list(
    cor = b,
    f_ratio = f,
    hs = hs,
    indices = idx
  )
  class(result) <- "boyce"
  return(result)
}

#' Print method for boyce objects
#'
#' @param x A boyce object
#' @param ... Additional arguments (not used)
#' @export
print.boyce <- function(x, ...) {
  cat("Boyce Index\n")
  cat("===========\n\n")
  cat("Correlation (Boyce index): ", round(x$cor, 3), "\n", sep = "")

  # Interpretation guide
  if (!is.na(x$cor)) {
    interpretation <- if (x$cor > 0.5) {
      "strong positive"
    } else if (x$cor > 0) {
      "positive"
    } else if (x$cor > -0.5) {
      "negative"
    } else {
      "strong negative"
    }
    cat("Interpretation:            ", interpretation, "\n", sep = "")
  }

  cat("\nDiagnostic components:\n")
  cat(
    "  - f_ratio: Predicted/expected ratios (",
    length(x$f_ratio),
    " bins)\n",
    sep = ""
  )
  cat(
    "  - hs:      Habitat suitability values (",
    length(x$hs),
    " bins)\n",
    sep = ""
  )
  cat(
    "  - indices: Indices used in correlation (",
    length(x$indices),
    " of ",
    length(x$f_ratio),
    " bins)\n",
    sep = ""
  )

  cat("\nUse plot(result) to visualise the predicted/expected ratio.\n")
  cat("Use str() or summary() to see all components.\n")
  invisible(x)
}

#' Plot method for boyce objects
#'
#' Plots the predicted-to-expected ratio along habitat suitability values.
#' Grey points show all bins, black points show bins used in the correlation
#' calculation (after removing successive duplicates if applicable).
#'
#' @param x A boyce object
#' @param ... Additional graphical parameters passed to plot
#' @export
plot.boyce <- function(x, ...) {
  # Plot all bins in grey
  plot(
    x$hs,
    x$f_ratio,
    xlab = "Habitat suitability",
    ylab = "Predicted/Expected ratio",
    col = "grey",
    cex = 0.75,
    las = 1,
    ...
  )

  # Overlay bins used in correlation in black
  points(x$hs[x$indices], x$f_ratio[x$indices], pch = 19, cex = 0.75)

  # Add reference line at 1
  abline(h = 1, lty = 2, col = "darkgrey")

  # Add title with Boyce index
  if (!is.na(x$cor)) {
    title(main = paste0("Boyce Index = ", round(x$cor, 3)))
  }
}

#' @describeIn boyce Method for Rangebag model objects
#' @param x Environmental/climate data with corresponding model variables as a
#'   `terra::SpatRaster`, `raster::Raster*`, `data.frame`, or
#'   `matrix`. Required when `fit` is a model object.
#' @export
boyce.Rangebag <- function(
  fit,
  x,
  obs,
  nclass = 0,
  window_width = "default",
  res = 100,
  pe_plot = TRUE,
  rm_duplicate = TRUE,
  method = "spearman",
  ...
) {
  # Generate predictions from the Rangebag model
  predictions <- predict(fit, x, raw_output = FALSE, ...)

  # Call the default method with the predictions
  boyce.default(
    fit = predictions,
    obs = obs,
    nclass = nclass,
    window_width = window_width,
    res = res,
    pe_plot = pe_plot,
    rm_duplicate = rm_duplicate,
    method = method
  )
}

#' @describeIn boyce Method for Climatch model objects
#' @export
boyce.Climatch <- function(
  fit,
  x,
  obs,
  nclass = 0,
  window_width = "default",
  res = 100,
  pe_plot = TRUE,
  rm_duplicate = TRUE,
  method = "spearman",
  ...
) {
  # Generate predictions from the Climatch model
  predictions <- predict(fit, x, raw_output = FALSE, ...)

  # Call the default method with the predictions
  boyce.default(
    fit = predictions,
    obs = obs,
    nclass = nclass,
    window_width = window_width,
    res = res,
    pe_plot = pe_plot,
    rm_duplicate = rm_duplicate,
    method = method
  )
}
