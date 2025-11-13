#' Calculate Multivariate Environmental Similarity
#'
#' Calculate Multivariate Environmental Similarity and most dissimilar/similar
#' variables with respect to a reference dataset, for a set of environmental
#' variables.
#'
#' @param x Climate (or environmental) data as a `terra::SpatRaster`,
#'   \code{raster::Raster*}, `list`, `data.frame`, or `matrix` where each
#'   layer/element/column represents focal values of an environmental 
#'   variable.
#' @param ref A `data.frame`, `matrix`, or `list` where each column/element
#'   represents reference values for an environmental variable (corresponding
#'   to those given in `x`).
#' @param full Logical to indicate whether similarity values should be returned
#'   for all variables. If `FALSE` (the default), then only the minimum
#'   similarity scores across variables will be returned.
#' @param filename Optional filename for writing spatial raster output (i.e., 
#'   the `mess` raster only). Default is "".
#' @param ... Additional parameters.
#' @return If `x` is a `SpatRaster` or \code{Raster*} object, this function returns
#'   a list containing:
#'   - `mess`: a `SpatRaster` layer giving the minimum similarity
#'   value across all variables for each location (i.e. the MESS);
#'   - `mess_by_variable`: a `SpatRaster` giving the environmental similarities 
#'   for each variable in `x` (only included when `full=TRUE`);
#'   - `mod`: a factor `SpatRaster` layer indicating which variable was most
#'   dissimilar to its reference range (i.e. the MoD map, Elith et al. 2010);
#'   and
#'   - `mos`: a factor `SpatRaster` layer indicating which variable was most
#'   similar to its reference range.
#'
#'   If `x` is a `list`, `matrix` or `data.frame`, the function will return a 
#'   list as above, but with multilayer and single layer `SpatRaster` objects
#'   replaced by matrix and vector objects, respectively.
#' @details `mess` uses the MESS algorithm described in Appendix S3 of Elith
#'   et al. 2010.
#' @keywords maxent, mess, similarity, environment
#' @references Elith, J., Kearney, M., and Phillips, S. (2010). The art of
#'   modelling range-shifting species. \emph{Methods in Ecology and Evolution},
#'   1: 330-342. \doi{10.1111/j.2041-210X.2010.00036.x}
#' @export
#' @examples
#' \dontrun{
#' library(geodata)
#' library(terra)
#' library(tmap)
#' bio <- worldclim_global("bio", res = 10, path = tempdir())
#' aus <- gadm("AUS", level = 0, resolution = 2, path = tempdir())
#' occ <- spatSample(aus, size = 100, method = "random")
#' ref <- terra::extract(bio, occ, ID = FALSE)
#' m <- mess(bio, ref, full = TRUE)
#'
#' # Plot outputs
#' tmap::tm_shape(m$mess) +
#'   tmap::tm_raster(
#'     col.scale = tmap::tm_scale_continuous(
#'       values = "matplotlib.rd_bu",
#'       midpoint = 0
#'     ),
#'     col.legend = tmap::tm_legend(title = "MESS")
#'   ) +
#'   tmap::tm_credits(
#'     "Positive values indicate similarity, negative values indicate dissimilarity", 
#'     position = tmap::tm_pos_out("center", "bottom")
#'   ) +
#'   tmap::tm_title("Multivariate Environmental Similarity")
#'
#' tmap::tm_shape(m$mod) +
#'   tmap::tm_raster(
#'     col.scale = tmap::tm_scale_categorical(
#'       values = "20",
#'       n.max = nrow(levels(m$mod)[[1]])
#'     ),
#'     col.legend = tmap::tm_legend(title = "Most dissimilar\ncovariate (MoD)"),
#'     col.chart = tmap::tm_chart_donut()
#'   ) +
#'   tmap::tm_title("Most dissimilar covariate (MoD)")
#'
#' tmap::tm_shape(m$mos) +
#'   tmap::tm_raster(
#'     col.scale = tmap::tm_scale_categorical(
#'       values = "20",
#'       n.max = nrow(levels(m$mos)[[1]])
#'     ),
#'     col.legend = tmap::tm_legend(title = "Most similar\ncovariate (MoS)"),
#'     col.chart = tmap::tm_chart_donut()
#'   ) +
#'   tmap::tm_title("Most similar covariate (MoS)")
#' }
mess <- function(x, ref, full = FALSE, filename = "", ...) {
  UseMethod("mess")
}

#' @name mess
#' @export
mess.Raster <- function(x, ref, full = FALSE, filename = "", ...) {
  # Convert to SpatRaster
  x <- terra::rast(x)
  
  # Call the terra version of the function
  mess(x, ref, full = full, filename = filename, ...)  
}

#' @name mess
#' @export
mess.SpatRaster <- function(x, ref, full = FALSE, filename = "", ...) {
  
  # Validate and preprocess reference data
  ref_prep <- .prep_mess_ref(ref, names(x))
  cols <- ref_prep$cols
  n_vars <- length(cols)
  ref_sorted <- ref_prep$ref_sorted
  ref_lengths <- ref_prep$ref_lengths
  rng <- ref_prep$rng
  
  # Subset x to common columns
  x <- x[[cols]]
  
  # Initialise output rasters
  if (isTRUE(full)) {
    out_sim <- terra::init(x, NA)
    names(out_sim) <- cols
  }
  out_min <- terra::init(x[[1]], NA)
  names(out_min) <- "mess"
  out_mod <- terra::init(x[[1]], NA)
  out_mos <- terra::init(x[[1]], NA)
  
  # Handle raster data in blocks
  x_blocks <- terra::blocks(x, n = 8)
  n_blocks <- x_blocks$n
  terra::readStart(x)
  
  if (isTRUE(full)) {
    invisible(terra::writeStart(out_sim, filename = "", ...))
  }
  invisible(terra::writeStart(out_min, filename = filename, ...))
  invisible(terra::writeStart(out_mod, filename = "", ...))
  invisible(terra::writeStart(out_mos, filename = "", ...))
  
  for (b in 1:n_blocks) {
    # Read block data
    x_data <- terra::readValues(x,
                                row = x_blocks$row[b],
                                nrows = x_blocks$nrows[b],
                                mat = TRUE)
    
    # Identify non-NA rows
    x_idx <- which(as.logical(rowSums(!is.na(x_data))))
    
    if (length(x_idx) > 0) {
      x_vals <- x_data[x_idx, , drop = FALSE]
      
      # Calculate MESS for this block
      mess_result <- .calculate_mess(x_vals, ref_sorted, ref_lengths, rng, cols)
      
      # Prepare output data for this block
      if (isTRUE(full)) {
        output_sim <- matrix(NA, nrow = nrow(x_data), ncol = n_vars)
        output_sim[x_idx, ] <- mess_result$sim
      }
      output_min <- rep(NA, nrow(x_data))
      output_min[x_idx] <- mess_result$min_sim
      output_mod <- rep(NA, nrow(x_data))
      output_mod[x_idx] <- mess_result$mod
      output_mos <- rep(NA, nrow(x_data))
      output_mos[x_idx] <- mess_result$mos
      
    } else {
      # All NA block
      if (isTRUE(full)) {
        output_sim <- matrix(NA, nrow = nrow(x_data), ncol = n_vars)
      }
      output_min <- rep(NA, nrow(x_data))
      output_mod <- rep(NA, nrow(x_data))
      output_mos <- rep(NA, nrow(x_data))
    }
    
    # Write block outputs
    if (isTRUE(full)) {
      terra::writeValues(out_sim, v = output_sim,
                        start = x_blocks$row[b], nrows = x_blocks$nrows[b])
    }
    terra::writeValues(out_min, v = output_min,
                      start = x_blocks$row[b], nrows = x_blocks$nrows[b])
    terra::writeValues(out_mod, v = output_mod,
                      start = x_blocks$row[b], nrows = x_blocks$nrows[b])
    terra::writeValues(out_mos, v = output_mos,
                      start = x_blocks$row[b], nrows = x_blocks$nrows[b])
  }
  
  # Close block reading & writing
  terra::readStop(x)
  if (isTRUE(full)) {
    invisible(terra::writeStop(out_sim))
  }
  invisible(terra::writeStop(out_min))
  invisible(terra::writeStop(out_mod))
  invisible(terra::writeStop(out_mos))
  
  # Convert mod and mos to factors with variable names
  # (terra::as.factor requires reading the raster, so do this after writeStop)
  # We suppressWarnings to avoid chatty warnings emitted by terra
  out_mod <- suppressWarnings(terra::as.factor(out_mod))
  levs <- data.frame(ID = seq_len(n_vars), variable = cols)
  levels(out_mod)[[1]] <- levs
  out_mos <- suppressWarnings(terra::as.factor(out_mos))
  levels(out_mos)[[1]] <- levs
  
  # Return results
  if (isTRUE(full)) {
    return(list(
      mess = out_min,
      mess_by_variable = out_sim,
      mod = out_mod,
      mos = out_mos
    ))
  } else {
    return(list(
      mess = out_min,
      mod = out_mod,
      mos = out_mos
    ))
  }
}

#' @name mess
#' @export
mess.data.frame <- function(x, ref, full = FALSE, ...) {
  
  # Validate and preprocess reference data
  ref_prep <- .prep_mess_ref(ref, names(x))
  cols <- ref_prep$cols
  ref_sorted <- ref_prep$ref_sorted
  ref_lengths <- ref_prep$ref_lengths
  rng <- ref_prep$rng
  
  # Subset x to common columns
  x <- x[, cols, drop = FALSE]
  
  # Calculate MESS
  mess_result <- .calculate_mess(x, ref_sorted, ref_lengths, rng, cols)
  
  # Return results
  if (isTRUE(full)) {
    return(list(
      mess = mess_result$min_sim,
      mess_by_variable = mess_result$sim,
      mod = mess_result$mod,
      mos = mess_result$mos
    ))
  } else {
    return(list(
      mess = mess_result$min_sim,
      mod = mess_result$mod,
      mos = mess_result$mos
    ))
  }
}

#' @name mess
#' @export
mess.matrix <- function(x, ref, full = FALSE, ...) {
  mess(as.data.frame(x), ref, full = full, ...)
}

# Internal helper function to validate and preprocess reference data
.prep_mess_ref <- function(ref, x_vars) {
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
    warning(paste0(
      "The following variables are missing from x and will be ignored: ",
      paste(x_vars_missing, collapse = ", ")
    ), call. = FALSE)
  }
  if (length(ref_vars_missing) > 0) {
    warning(paste0(
      "The following variables are missing from ref and will be ignored: ",
      paste(ref_vars_missing, collapse = ", ")
    ), call. = FALSE)
  }
  
  cols <- intersect(x_vars, ref_vars)
  if (length(cols) == 0) {
    stop("No variables in common between x and ref.", call. = FALSE)
  }
  
  # Subset to common columns and calculate reference statistics
  ref <- ref[, cols, drop = FALSE]
  rng <- as.data.frame(apply(ref, 2, range, na.rm = TRUE))
  ref_sorted <- lapply(ref, sort)
  ref_lengths <- lapply(ref, length)
  
  return(list(
    cols = cols,
    ref = ref,
    rng = rng,
    ref_sorted = ref_sorted,
    ref_lengths = ref_lengths
  ))
}

# Internal function to calculate MESS statistics
# This is called by both the SpatRaster and data.frame methods
.calculate_mess <- function(x_vals, ref_sorted, ref_lengths, rng, cols) {
  # Handle empty input
  if (nrow(x_vals) == 0) {
    n_vars <- length(cols)
    return(list(
      sim = matrix(numeric(0), nrow = 0, ncol = n_vars),
      min_sim = numeric(0),
      mod = integer(0),
      mos = integer(0)
    ))
  }
  
  # Convert to data frame for mapply
  x_df <- as.data.frame(x_vals)
  names(x_df) <- cols
  
  # Calculate percentile less than each value
  pct_less <- mapply(
    function(x, ref_sorted, ref_length) {
      findInterval(x, ref_sorted) / ref_length
    },
    x_df,
    ref_sorted,
    ref_lengths,
    SIMPLIFY = FALSE
  )
  
  # Calculate similarity for each variable
  sim <- mapply(
    function(f, rng, p) {
      ifelse(
        f == 0,
        (p - rng[1]) / diff(rng) * 100,
        ifelse(
          f > 0 & f <= 0.5,
          f * 200,
          ifelse(f > 0.5 & f < 1, (1 - f) * 200, (rng[2] - p) / diff(rng) * 100)
        )
      )
    },
    pct_less,
    rng,
    x_df
  )
  
  # Calculate minimum similarity and identify most/least similar variables
  if (ncol(sim) == 1) {
    min_sim <- sim[, 1]
    most_dissimilar_vec <- rep(1, nrow(sim))
    most_similar_vec <- rep(1, nrow(sim))
  } else {
    min_sim <- apply(sim, 1, min)
    mins <- apply(sim, 1, which.min)
    most_dissimilar_vec <- unlist(ifelse(lengths(mins) == 0, NA, mins))
    maxs <- apply(sim, 1, which.max)
    most_similar_vec <- unlist(ifelse(lengths(maxs) == 0, NA, maxs))
  }
  
  return(list(
    sim = sim,
    min_sim = min_sim,
    mod = most_dissimilar_vec,
    mos = most_similar_vec
  ))
}
