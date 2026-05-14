# Internal helper: build a categorical colour palette that scales with the
# number of levels in a factor SpatRaster. Uses the cols4all "20" palette for
# up to 20 categories; falls back to grDevices::hcl.colors() for larger sets.
.categorical_pal <- function(rast_layer) {
  levs <- terra::levels(rast_layer)[[1]]
  n <- if (!is.null(levs) && is.data.frame(levs) && nrow(levs) > 0) nrow(levs) else 20L
  if (n <= 20L) "20" else grDevices::hcl.colors(n, palette = "Dark 3")
}

#' Plot method for MessResult objects
#'
#' Produces a tmap plot for outputs of \code{\link{mess}}.
#'
#' @param x A \code{MessResult} object returned by \code{\link{mess}}.
#' @param which Character string specifying which output to plot. One of
#'   \code{"mess"} (default), \code{"mod"} (most dissimilar variable), or
#'   \code{"mos"} (most similar variable).
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns \code{x}. Called for its side-effect of rendering
#'   a tmap plot.
#' @note Requires the \pkg{tmap} package to be installed.
#' @export
plot.MessResult <- function(x, which = "mess", ...) {
  if (!requireNamespace("tmap", quietly = TRUE)) {
    stop(
      "The 'tmap' package is required for plotting MessResult objects. ",
      "Install it with: install.packages('tmap')",
      call. = FALSE
    )
  }

  which <- match.arg(which, c("mess", "mod", "mos"))

  p <- switch(which,
    mess = {
      tmap::tm_shape(x$mess) +
        tmap::tm_raster(
          col.scale = tmap::tm_scale_continuous(
            values = "-matplotlib.seismic",
            midpoint = 0
          ),
          col.legend = tmap::tm_legend(
            position = tmap::tm_pos_out("right", "center"),
            title = "MESS"
          )
        ) +
        tmap::tm_credits(
          paste(
            "Negative values (red) indicate at least one variable is beyond its reference range.",
            "Positive values (blue) indicate similarity to reference conditions.",
            sep = "\n"
          ),
          position = tmap::tm_pos_out("center", "bottom")
        ) +
        tmap::tm_title("Multivariate Environmental Similarity")
    },
    mod = {
      if (is.null(x$mod)) {
        stop(
          "MoD output not found. Re-run mess() to obtain this output.",
          call. = FALSE
        )
      }
      tmap::tm_shape(x$mod) +
        tmap::tm_raster(
          col.scale = tmap::tm_scale_categorical(
            values = .categorical_pal(x$mod)
          ),
          col.legend = tmap::tm_legend(
            position = tmap::tm_pos_out("right", "center"),
            group_id = "legend",
            title = "Most dissimilar\ncovariate (MoD)"
          ),
          col.chart = tmap::tm_chart_donut(
            position = tmap::tm_pos_out("right", "center"),
            group_id = "legend"
          )
        ) +
        tmap::tm_title("Most dissimilar covariate (MoD)") +
        tmap::tm_components(stack = "vertical")
    },
    mos = {
      if (is.null(x$mos)) {
        stop(
          "MoS output not found. Re-run mess() to obtain this output.",
          call. = FALSE
        )
      }
      tmap::tm_shape(x$mos) +
        tmap::tm_raster(
          col.scale = tmap::tm_scale_categorical(
            values = .categorical_pal(x$mos)
          ),
          col.legend = tmap::tm_legend(
            position = tmap::tm_pos_out("right", "center"),
            group_id = "legend",
            title = "Most similar\ncovariate (MoS)"
          ),
          col.chart = tmap::tm_chart_donut(
            position = tmap::tm_pos_out("right", "center"),
            group_id = "legend"
          )
        ) +
        tmap::tm_title("Most similar covariate (MoS)") +
        tmap::tm_components(stack = "vertical")
    }
  )

  print(p)
  invisible(x)
}

#' Plot method for ExdetResult objects
#'
#' Produces a tmap plot for outputs of \code{\link{exdet}}.
#'
#' @param x An \code{ExdetResult} object returned by \code{\link{exdet}}.
#' @param which Character string specifying which output to plot. One of
#'   \code{"exdet"} (default), \code{"mic1"} (most influential covariate for
#'   Type 1 novelty), or \code{"mic2"} (most influential covariate for Type 2
#'   novelty).
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns \code{x}. Called for its side-effect of rendering
#'   a tmap plot.
#' @note Requires the \pkg{tmap} package to be installed.
#' @export
plot.ExdetResult <- function(x, which = "exdet", ...) {
  if (!requireNamespace("tmap", quietly = TRUE)) {
    stop(
      "The 'tmap' package is required for plotting ExdetResult objects. ",
      "Install it with: install.packages('tmap')",
      call. = FALSE
    )
  }

  which <- match.arg(which, c("exdet", "mic1", "mic2"))

  p <- switch(which,
    exdet = {
      exdet_type1 <- terra::classify(x$exdet, cbind(0, Inf, NA))
      exdet_type2 <- terra::classify(x$exdet, cbind(-Inf, 1, NA))
      exdet_analog <- x$exdet >= 0 & x$exdet <= 1
      exdet_analog[!exdet_analog] <- NA

      tmap::tm_shape(exdet_type1) +
        tmap::tm_raster(
          col.scale = tmap::tm_scale_continuous(
            values = "reds3", values.range = c(0, 0.75)
          ),
          col.legend = tmap::tm_legend(
            position = tmap::tm_pos_out("right", "center"),
            title = "Type 1\nNovelty"
          )
        ) +
        tmap::tm_shape(exdet_type2) +
        tmap::tm_raster(
          col.scale = tmap::tm_scale_continuous(
            values = "-blues3", values.range = c(0, 0.75)
          ),
          col.legend = tmap::tm_legend(
            position = tmap::tm_pos_out("right", "center"),
            title = "Type 2\nNovelty"
          )
        ) +
        tmap::tm_shape(exdet_analog) +
        tmap::tm_raster(
          col.scale = tmap::tm_scale_categorical(
            values = "#FFEE88", labels = "Similar"
          ),
          col.legend = tmap::tm_legend(
            position = tmap::tm_pos_out("right", "center"),
            title = ""
          )
        ) +
        tmap::tm_credits(
          paste(
            "Values 0-1 indicate analog conditions.",
            "Values < 0 indicate novel univariate ranges.",
            "Values > 1 indicate novel covariate correlations.",
            sep = "\n"
          ),
          position = tmap::tm_pos_out("center", "bottom")
        ) +
        tmap::tm_title("Extrapolation Detection")
    },
    mic1 = {
      if (is.null(x$mic1)) {
        stop(
          "MIC1 output not found. Re-run exdet() with mic = TRUE.",
          call. = FALSE
        )
      }
      tmap::tm_shape(x$mic1) +
        tmap::tm_raster(
          col.scale = tmap::tm_scale_categorical(
            values = .categorical_pal(x$mic1)
          ),
          col.legend = tmap::tm_legend(
            position = tmap::tm_pos_out("right", "center"),
            group_id = "legend",
            title = "Most influential\ncovariate (MIC1)"
          ),
          col.chart = tmap::tm_chart_donut(
            position = tmap::tm_pos_out("right", "center"),
            group_id = "legend"
          )
        ) +
        tmap::tm_title("Type 1 novelty: Most influential covariate") +
        tmap::tm_components(stack = "vertical")
    },
    mic2 = {
      if (is.null(x$mic2)) {
        stop(
          "MIC2 output not found. Re-run exdet() with mic = TRUE.",
          call. = FALSE
        )
      }
      tmap::tm_shape(x$mic2) +
        tmap::tm_raster(
          col.scale = tmap::tm_scale_categorical(
            values = .categorical_pal(x$mic2)
          ),
          col.legend = tmap::tm_legend(
            position = tmap::tm_pos_out("right", "center"),
            group_id = "legend",
            title = "Most influential\ncovariate (MIC2)"
          ),
          col.chart = tmap::tm_chart_donut(
            position = tmap::tm_pos_out("right", "center"),
            group_id = "legend"
          )
        ) +
        tmap::tm_title("Type 2 novelty: Most influential covariate") +
        tmap::tm_components(stack = "vertical")
    }
  )

  print(p)
  invisible(x)
}
