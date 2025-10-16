context("exdet")

test_that("exdet runs with SpatRaster input", {
  TEST_DIRECTORY <- test_path("test_inputs")

  # Load climate rasters
  climate_rast <- terra::rast(file.path(
    TEST_DIRECTORY,
    sprintf("bioclim_10m/wc2.1bio%02d.tif", c(1, 4, 7, 12, 15))
  ))

  # Load coordinates and extract reference data
  sample_moths <- read.csv(file.path(TEST_DIRECTORY, "sample_moths.csv"))
  ref <- terra::extract(climate_rast, sample_moths[, c("lon", "lat")])
  ref$ID <- NULL # Remove ID column added by terra::extract

  # Test with mic = FALSE
  expect_silent(exdet_result <- exdet(climate_rast, ref, mic = FALSE))
  expect_is(exdet_result, "list")
  expect_equal(length(exdet_result), 1)
  expect_true("exdet" %in% names(exdet_result))

  # Check exdet output
  expect_is(exdet_result$exdet, "SpatRaster")
  expect_equal(terra::nlyr(exdet_result$exdet), 1)
  expect_equal(dim(exdet_result$exdet), dim(climate_rast[[1]]))
  expect_equal(names(exdet_result$exdet), "exdet")

  # Check that values are in reasonable range
  exdet_vals <- terra::values(exdet_result$exdet, na.rm = TRUE)
  expect_true(all(is.finite(exdet_vals)))

  # Test with mic = TRUE
  expect_silent(exdet_full <- exdet(climate_rast, ref, mic = TRUE))
  expect_is(exdet_full, "list")
  expect_equal(length(exdet_full), 3)
  expect_true(all(c("exdet", "mic1", "mic2") %in% names(exdet_full)))

  # Check mic1 output (most influential variable for Type 1)
  expect_is(exdet_full$mic1, "SpatRaster")
  expect_true(terra::is.factor(exdet_full$mic1))
  expect_equal(nrow(terra::levels(exdet_full$mic1)[[1]]), 6) # 5 vars + "Not novel"
  expect_true("Not novel" %in% terra::levels(exdet_full$mic1)[[1]]$variable)
  expect_true(all(
    names(climate_rast) %in%
      terra::levels(exdet_full$mic1)[[1]]$variable
  ))

  # Check mic2 output (most influential variable for Type 2)
  expect_is(exdet_full$mic2, "SpatRaster")
  expect_true(terra::is.factor(exdet_full$mic2))
  expect_equal(nrow(terra::levels(exdet_full$mic2)[[1]]), 6)
  expect_true("Not novel" %in% terra::levels(exdet_full$mic2)[[1]]$variable)

  # Check that exdet values are consistent between mic and non-mic
  expect_equal(
    terra::values(exdet_result$exdet),
    terra::values(exdet_full$exdet)
  )
})

test_that("exdet runs with data.frame input", {
  TEST_DIRECTORY <- test_path("test_inputs")

  # Load climate rasters and convert to data frame
  climate_rast <- terra::rast(file.path(
    TEST_DIRECTORY,
    sprintf("bioclim_10m/wc2.1bio%02d.tif", c(1, 4, 7, 12, 15))
  ))
  climate_df <- as.data.frame(climate_rast, xy = TRUE, na.rm = TRUE)

  # Load coordinates and extract reference data
  sample_moths <- read.csv(file.path(TEST_DIRECTORY, "sample_moths.csv"))
  ref <- terra::extract(climate_rast, sample_moths[, c("lon", "lat")])
  ref$ID <- NULL

  # Test with mic = FALSE
  expect_silent(exdet_result <- exdet(climate_df[, -c(1, 2)], ref, mic = FALSE))
  expect_is(exdet_result, "list")
  expect_equal(length(exdet_result), 1)
  expect_true("exdet" %in% names(exdet_result))

  # Check exdet output
  expect_is(exdet_result$exdet, "numeric")
  expect_equal(length(exdet_result$exdet), nrow(climate_df))

  # Check that values are reasonable
  expect_true(all(is.finite(exdet_result$exdet)))

  # Test with mic = TRUE
  expect_silent(exdet_full <- exdet(climate_df[, -c(1, 2)], ref, mic = TRUE))
  expect_is(exdet_full, "list")
  expect_equal(length(exdet_full), 3)
  expect_true(all(c("exdet", "mic1", "mic2") %in% names(exdet_full)))

  # Check mic1 output
  expect_is(exdet_full$mic1, "factor")
  expect_equal(length(exdet_full$mic1), nrow(climate_df))
  expect_equal(nlevels(exdet_full$mic1), 6) # 5 vars + "Not novel"
  expect_true("Not novel" %in% levels(exdet_full$mic1))

  # Check mic2 output
  expect_is(exdet_full$mic2, "factor")
  expect_equal(length(exdet_full$mic2), nrow(climate_df))
  expect_equal(nlevels(exdet_full$mic2), 6)
  expect_true("Not novel" %in% levels(exdet_full$mic2))

  # Check that exdet values are consistent
  expect_equal(exdet_result$exdet, exdet_full$exdet)
})

test_that("exdet runs with matrix input", {
  TEST_DIRECTORY <- test_path("test_inputs")

  # Load climate rasters and convert to matrix
  climate_rast <- terra::rast(file.path(
    TEST_DIRECTORY,
    sprintf("bioclim_10m/wc2.1bio%02d.tif", c(1, 4, 7, 12, 15))
  ))
  climate_mat <- as.matrix(as.data.frame(climate_rast, na.rm = TRUE))

  # Load coordinates and extract reference data
  sample_moths <- read.csv(file.path(TEST_DIRECTORY, "sample_moths.csv"))
  ref <- terra::extract(climate_rast, sample_moths[, c("lon", "lat")])
  ref$ID <- NULL

  # Test with matrix input
  expect_silent(exdet_result <- exdet(climate_mat, ref, mic = FALSE))
  expect_is(exdet_result, "list")
  expect_true("exdet" %in% names(exdet_result))
  expect_is(exdet_result$exdet, "numeric")
  expect_equal(length(exdet_result$exdet), nrow(climate_mat))
})

test_that("exdet handles missing variables correctly", {
  TEST_DIRECTORY <- test_path("test_inputs")

  # Load climate rasters
  climate_rast <- terra::rast(file.path(
    TEST_DIRECTORY,
    sprintf("bioclim_10m/wc2.1bio%02d.tif", c(1, 4, 7, 12, 15))
  ))

  # Load coordinates and extract reference data
  sample_moths <- read.csv(file.path(TEST_DIRECTORY, "sample_moths.csv"))
  ref <- terra::extract(climate_rast, sample_moths[, c("lon", "lat")])
  ref$ID <- NULL

  # Test with subset of variables in x
  expect_warning(
    exdet_result <- exdet(climate_rast[[1:3]], ref, mic = FALSE),
    "The following variables are missing from x"
  )
  expect_is(exdet_result, "list")

  # Test with subset of variables in ref
  expect_warning(
    exdet_result <- exdet(climate_rast, ref[, 1:3], mic = FALSE),
    "The following variables are missing from ref"
  )
  expect_is(exdet_result, "list")

  # Test with no overlapping variables
  names(ref) <- paste0("var_", 1:5)
  expect_error(
    expect_warning(
      exdet(climate_rast, ref, mic = FALSE),
      "missing from x"
    ),
    "No variables in common"
  )
})

test_that("exdet handles tolerance parameter correctly", {
  TEST_DIRECTORY <- test_path("test_inputs")

  # Load climate rasters
  climate_rast <- terra::rast(file.path(
    TEST_DIRECTORY,
    sprintf("bioclim_10m/wc2.1bio%02d.tif", c(1, 4, 7, 12, 15))
  ))

  # Load coordinates and extract reference data
  sample_moths <- read.csv(file.path(TEST_DIRECTORY, "sample_moths.csv"))
  ref <- terra::extract(climate_rast, sample_moths[, c("lon", "lat")])
  ref$ID <- NULL

  # Test with small tolerance (should work for 5 variables)
  expect_silent(
    exdet_result <- exdet(climate_rast, ref, mic = FALSE, tol = 1e-20)
  )
  expect_is(exdet_result, "list")

  # Test with larger tolerance (should also work)
  expect_silent(
    exdet_result <- exdet(climate_rast, ref, mic = FALSE, tol = 1e-12)
  )
  expect_is(exdet_result, "list")
})

test_that("exdet consistency between SpatRaster and data.frame", {
  TEST_DIRECTORY <- test_path("test_inputs")

  # Load climate rasters
  climate_rast <- terra::rast(file.path(
    TEST_DIRECTORY,
    sprintf("bioclim_10m/wc2.1bio%02d.tif", c(1, 4, 7, 12, 15))
  ))
  climate_df <- as.data.frame(climate_rast, na.rm = TRUE)

  # Load coordinates and extract reference data
  sample_moths <- read.csv(file.path(TEST_DIRECTORY, "sample_moths.csv"))
  ref <- terra::extract(climate_rast, sample_moths[, c("lon", "lat")])
  ref$ID <- NULL

  # Run ExDet on both input types
  expect_silent(exdet_rast <- exdet(climate_rast, ref, mic = TRUE))
  expect_silent(exdet_df <- exdet(climate_df, ref, mic = TRUE))

  # Extract values from rasters, ensuring we get non-NA cells
  exdet_vals <- terra::values(exdet_rast$exdet)
  non_na_idx <- which(!is.na(exdet_vals))

  # Compare exdet scores
  expect_equal(exdet_vals[non_na_idx], exdet_df$exdet)

  # Compare MIC factors
  # Extract factor IDs and convert to labels
  mic1_vals <- terra::values(exdet_rast$mic1)[non_na_idx]
  mic1_levels <- terra::levels(exdet_rast$mic1)[[1]]
  # Match by ID (levels start at 0, so we need to match, not index directly)
  mic1_rast_labels <- mic1_levels$variable[match(mic1_vals, mic1_levels$ID)]
  expect_equal(mic1_rast_labels, as.character(exdet_df$mic1))

  mic2_vals <- terra::values(exdet_rast$mic2)[non_na_idx]
  mic2_levels <- terra::levels(exdet_rast$mic2)[[1]]
  mic2_rast_labels <- mic2_levels$variable[match(mic2_vals, mic2_levels$ID)]
  expect_equal(mic2_rast_labels, as.character(exdet_df$mic2))
})

test_that("exdet with filename argument", {
  TEST_DIRECTORY <- test_path("test_inputs")

  # Load climate rasters
  climate_rast <- terra::rast(file.path(
    TEST_DIRECTORY,
    sprintf("bioclim_10m/wc2.1bio%02d.tif", c(1, 4, 7, 12, 15))
  ))

  # Load coordinates and extract reference data
  sample_moths <- read.csv(file.path(TEST_DIRECTORY, "sample_moths.csv"))
  ref <- terra::extract(climate_rast, sample_moths[, c("lon", "lat")])
  ref$ID <- NULL

  # Create temporary file
  temp_file <- tempfile(fileext = ".tif")

  # Test with filename
  expect_silent(
    exdet_result <- exdet(climate_rast, ref, mic = FALSE, filename = temp_file)
  )
  expect_is(exdet_result, "list")
  expect_true(file.exists(temp_file))

  # Read the file and compare
  exdet_from_file <- terra::rast(temp_file)
  expect_equal(
    terra::values(exdet_result$exdet),
    terra::values(exdet_from_file)
  )

  # Clean up
  unlink(temp_file)
})

test_that("exdet handles empty input correctly", {
  TEST_DIRECTORY <- test_path("test_inputs")

  # Load climate rasters to get variable names
  climate_rast <- terra::rast(file.path(
    TEST_DIRECTORY,
    sprintf("bioclim_10m/wc2.1bio%02d.tif", c(1, 4, 7, 12, 15))
  ))

  # Load coordinates and extract reference data
  sample_moths <- read.csv(file.path(TEST_DIRECTORY, "sample_moths.csv"))
  ref <- terra::extract(climate_rast, sample_moths[, c("lon", "lat")])
  ref$ID <- NULL

  # Create empty data frame with correct column names
  empty_df <- as.data.frame(matrix(numeric(0), ncol = 5))
  names(empty_df) <- names(climate_rast)

  # Test with empty data frame
  expect_silent(exdet_result <- exdet(empty_df, ref, mic = FALSE))
  expect_is(exdet_result, "list")
  expect_equal(length(exdet_result$exdet), 0)

  # Test with mic = TRUE
  expect_silent(exdet_full <- exdet(empty_df, ref, mic = TRUE))
  expect_is(exdet_full, "list")
  expect_equal(length(exdet_full$exdet), 0)
  expect_equal(length(exdet_full$mic1), 0)
  expect_equal(length(exdet_full$mic2), 0)
})

test_that("exdet score interpretation", {
  TEST_DIRECTORY <- test_path("test_inputs")

  # Load climate rasters
  climate_rast <- terra::rast(file.path(
    TEST_DIRECTORY,
    sprintf("bioclim_10m/wc2.1bio%02d.tif", c(1, 4, 7, 12, 15))
  ))
  climate_df <- as.data.frame(climate_rast, na.rm = TRUE)

  # Load coordinates and extract reference data
  sample_moths <- read.csv(file.path(TEST_DIRECTORY, "sample_moths.csv"))
  ref <- terra::extract(climate_rast, sample_moths[, c("lon", "lat")])
  ref$ID <- NULL

  # Run exdet
  exdet_result <- exdet(climate_df, ref, mic = TRUE)

  # Check score ranges make sense
  # Type 1 novelty: negative values
  type1_novel <- exdet_result$exdet < 0
  if (any(type1_novel)) {
    # MIC1 should not be "Not novel" for Type 1 novelty
    expect_true(all(exdet_result$mic1[type1_novel] != "Not novel"))
  }

  # Analog: 0 to 1
  analog <- exdet_result$exdet >= 0 & exdet_result$exdet <= 1
  if (any(analog)) {
    # Some analog locations should have "Not novel" for both MICs
    # (though not necessarily all, as nt2 could be between 0 and 1)
    expect_true("Not novel" %in% exdet_result$mic1[analog])
  }

  # Type 2 novelty: > 1
  type2_novel <- exdet_result$exdet > 1
  if (any(type2_novel)) {
    # MIC1 should be "Not novel" for Type 2 novelty (no univariate novelty)
    expect_true(all(exdet_result$mic1[type2_novel] == "Not novel"))
    # MIC2 should not be "Not novel" for Type 2 novelty
    expect_true(all(exdet_result$mic2[type2_novel] != "Not novel"))
  }
})
