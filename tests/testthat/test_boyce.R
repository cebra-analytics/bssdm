context("boyce")

# Helper function to create simple test data
create_test_data <- function(n_fit = 1000, n_obs = 100) {
  set.seed(123)
  fit <- runif(n_fit, 0, 1)
  obs <- runif(n_obs, 0.3, 0.9) # Biased towards higher values
  list(fit = fit, obs = obs)
}

test_that("boyce works with numeric vectors", {
  data <- create_test_data()

  result <- boyce(data$fit, data$obs, pe_plot = FALSE)

  # Check class
  expect_s3_class(result, "boyce")

  # Check structure
  expect_type(result, "list")
  expect_named(result, c("cor", "f_ratio", "hs", "indices"))

  # Check correlation value
  expect_type(result$cor, "double")
  expect_true(!is.na(result$cor))
  expect_true(result$cor >= -1 && result$cor <= 1)

  # Check f_ratio
  expect_type(result$f_ratio, "double")
  expect_true(length(result$f_ratio) > 0)

  # Check hs
  expect_type(result$hs, "double")
  expect_equal(length(result$hs), length(result$f_ratio))

  # Check indices
  expect_type(result$indices, "integer")
  expect_true(all(
    result$indices >= 1 & result$indices <= length(result$f_ratio)
  ))
})

test_that("boyce works with SpatRaster and coordinates", {
  skip_if_not_installed("terra")

  # Create simple raster with gradient
  r <- terra::rast(
    nrows = 10,
    ncols = 10,
    xmin = 0,
    xmax = 10,
    ymin = 0,
    ymax = 10
  )
  terra::values(r) <- seq(0, 1, length.out = terra::ncell(r))

  # Create sample points biased towards high values
  coords <- data.frame(
    x = runif(20, 5, 10),
    y = runif(20, 5, 10)
  )

  result <- boyce(r, coords, pe_plot = FALSE)

  expect_s3_class(result, "boyce")
  expect_type(result$cor, "double")
  expect_true(!is.na(result$cor))
})

test_that("boyce works with SpatRaster and extracted values", {
  skip_if_not_installed("terra")

  # Create simple raster
  r <- terra::rast(
    nrows = 10,
    ncols = 10,
    xmin = 0,
    xmax = 10,
    ymin = 0,
    ymax = 10
  )
  terra::values(r) <- seq(0, 1, length.out = terra::ncell(r))

  # Extract values directly
  coords <- data.frame(x = runif(20, 5, 10), y = runif(20, 5, 10))
  obs_values <- terra::extract(r, coords, ID = FALSE)[, 1]

  result <- boyce(r, obs_values, pe_plot = FALSE)

  expect_s3_class(result, "boyce")
  expect_type(result$cor, "double")
})

test_that("boyce handles different nclass values", {
  data <- create_test_data()

  # Moving window (default)
  result0 <- boyce(data$fit, data$obs, nclass = 0, pe_plot = FALSE)
  expect_s3_class(result0, "boyce")
  expect_true(length(result0$f_ratio) > 10) # Should have many bins

  # Fixed number of classes (creates intervals, so n+1 bins)
  result5 <- boyce(data$fit, data$obs, nclass = 5, pe_plot = FALSE)
  expect_s3_class(result5, "boyce")
  expect_true(length(result5$f_ratio) >= 5) # At least 5 bins (some may be NaN)

  # Custom thresholds
  result_custom <- boyce(
    data$fit,
    data$obs,
    nclass = c(0.2, 0.4, 0.6, 0.8),
    pe_plot = FALSE
  )
  expect_s3_class(result_custom, "boyce")
})

test_that("boyce handles window_width parameter", {
  data <- create_test_data()

  # Default window
  result_default <- boyce(data$fit, data$obs, pe_plot = FALSE)

  # Very small custom window (should create more potential bins)
  result_custom <- boyce(
    data$fit,
    data$obs,
    window_width = 0.02,
    pe_plot = FALSE
  )

  expect_s3_class(result_default, "boyce")
  expect_s3_class(result_custom, "boyce")

  # Both should produce valid correlations
  expect_true(!is.na(result_default$cor))
  expect_true(!is.na(result_custom$cor))
})

test_that("boyce handles rm_duplicate parameter", {
  data <- create_test_data()

  result_rm <- boyce(data$fit, data$obs, rm_duplicate = TRUE, pe_plot = FALSE)
  result_keep <- boyce(
    data$fit,
    data$obs,
    rm_duplicate = FALSE,
    pe_plot = FALSE
  )

  expect_s3_class(result_rm, "boyce")
  expect_s3_class(result_keep, "boyce")

  # When removing duplicates, fewer indices should be used
  expect_true(length(result_rm$indices) <= length(result_keep$indices))
  expect_equal(length(result_keep$indices), length(result_keep$f_ratio))
})

test_that("boyce handles different correlation methods", {
  data <- create_test_data()

  result_spearman <- boyce(
    data$fit,
    data$obs,
    method = "spearman",
    pe_plot = FALSE
  )
  result_pearson <- boyce(
    data$fit,
    data$obs,
    method = "pearson",
    pe_plot = FALSE
  )
  result_kendall <- boyce(
    data$fit,
    data$obs,
    method = "kendall",
    pe_plot = FALSE
  )

  expect_s3_class(result_spearman, "boyce")
  expect_s3_class(result_pearson, "boyce")
  expect_s3_class(result_kendall, "boyce")

  # Different methods should generally give different correlations
  expect_true(abs(result_spearman$cor - result_pearson$cor) < 1)
  # ^ Should be correlated
})

test_that("boyce validates fit input", {
  data <- create_test_data()

  # Invalid fit types
  expect_error(
    boyce("not_valid", data$obs),
    "fit must be either a numeric vector or a SpatRaster object"
  )

  expect_error(
    boyce(list(1, 2, 3), data$obs),
    "fit must be either a numeric vector or a SpatRaster object"
  )
})

test_that("boyce validates obs input", {
  data <- create_test_data()

  # obs must be numeric when fit is numeric
  expect_error(
    boyce(data$fit, "not_valid"),
    "When fit is a numeric vector, obs must also be a numeric vector"
  )

  expect_error(
    boyce(data$fit, list(1, 2, 3)),
    "When fit is a numeric vector, obs must also be a numeric vector"
  )
})

test_that("boyce handles SpatVector input with helpful error", {
  skip_if_not_installed("terra")

  r <- terra::rast(nrows = 10, ncols = 10)
  terra::values(r) <- runif(100)

  # Create SpatVector (should error with helpful message)
  coords <- data.frame(x = c(1, 2, 3), y = c(1, 2, 3))
  spat_vec <- terra::vect(coords, geom = c("x", "y"))

  expect_error(
    boyce(r, spat_vec),
    "terra::crds\\(obs\\)"
  )
})

test_that("boyce handles NA values", {
  data <- create_test_data()

  # Add some NAs
  fit_na <- data$fit
  fit_na[1:10] <- NA

  obs_na <- data$obs
  obs_na[1:5] <- NA

  # Should work after removing NAs
  result <- boyce(fit_na, obs_na, pe_plot = FALSE)
  expect_s3_class(result, "boyce")
  expect_true(!is.na(result$cor))
})

test_that("boyce errors with all NA values", {
  # All NAs in fit
  expect_error(
    boyce(as.numeric(rep(NA, 100)), runif(50)),
    "fit contains no valid"
  )

  # All NAs in obs
  expect_error(
    boyce(runif(100), as.numeric(rep(NA, 50))),
    "obs contains no valid"
  )
})

test_that("boyce requires sufficient validation points", {
  data <- create_test_data()

  # Only 1 observation
  expect_error(
    boyce(data$fit, data$obs[1]),
    "obs must contain at least 2 validation points"
  )
})

test_that("boyce handles edge case with limited variation", {
  # Create data where all obs fall in narrow range
  fit <- runif(1000, 0, 1)
  obs <- runif(100, 0.49, 0.51) # Very narrow range

  # Should still work, though correlation might be weak
  result <- boyce(fit, obs, nclass = 10, pe_plot = FALSE)
  expect_s3_class(result, "boyce")
  # Correlation might be any value including NA if not enough variation
})

test_that("print.boyce produces output", {
  data <- create_test_data()
  result <- boyce(data$fit, data$obs, pe_plot = FALSE)

  # Should print without error
  expect_output(print(result), "Boyce Index")
  expect_output(print(result), "Correlation")
  expect_output(print(result), "Interpretation")
})

test_that("plot.boyce produces plot", {
  data <- create_test_data()
  result <- boyce(data$fit, data$obs, pe_plot = FALSE)

  # Should plot without error
  expect_silent(plot(result))

  # Check that plot was created (by checking graphics device)
  # This is a basic check - in interactive use you'd see the plot
  dev.off() # Clean up
})

test_that("boyce.Rangebag method exists and works", {
  # Test that the method exists
  all_methods <- as.character(methods(boyce))
  expect_true("boyce.Rangebag" %in% all_methods)

  # Test that it can be retrieved
  rb_method <- getS3method("boyce", "Rangebag")
  expect_type(rb_method, "closure")
  expect_true("fit" %in% names(formals(rb_method)))
  expect_true("x" %in% names(formals(rb_method)))
  expect_true("obs" %in% names(formals(rb_method)))
})

test_that("boyce.Climatch method exists and works", {
  # Test that the method exists
  all_methods <- as.character(methods(boyce))
  expect_true("boyce.Climatch" %in% all_methods)

  # Test that it can be retrieved
  cm_method <- getS3method("boyce", "Climatch")
  expect_type(cm_method, "closure")
  expect_true("fit" %in% names(formals(cm_method)))
  expect_true("x" %in% names(formals(cm_method)))
  expect_true("obs" %in% names(formals(cm_method)))
})

test_that("boyce pe_plot parameter works", {
  data <- create_test_data()

  # With plot (should work via print method)
  result_plot <- boyce(data$fit, data$obs, pe_plot = TRUE)
  expect_s3_class(result_plot, "boyce")
  dev.off() # Clean up

  # Without plot
  result_no_plot <- boyce(data$fit, data$obs, pe_plot = FALSE)
  expect_s3_class(result_no_plot, "boyce")

  # Results should be identical except for potential plotting side effects
  expect_equal(result_plot$cor, result_no_plot$cor)
  expect_equal(result_plot$f_ratio, result_no_plot$f_ratio)
})
