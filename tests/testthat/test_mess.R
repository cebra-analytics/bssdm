context("mess")

test_that("mess runs with SpatRaster input", {
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

  # Test with full = FALSE
  expect_silent(mess_result <- mess(climate_rast, ref, full = FALSE))
  expect_is(mess_result, "list")
  expect_equal(length(mess_result), 3)
  expect_true(all(c("mess", "mod", "mos") %in% names(mess_result)))

  # Check mess output
  expect_is(mess_result$mess, "SpatRaster")
  expect_equal(terra::nlyr(mess_result$mess), 1)
  expect_equal(dim(mess_result$mess), dim(climate_rast[[1]]))

  # Check mod output (most dissimilar variable)
  expect_is(mess_result$mod, "SpatRaster")
  expect_true(terra::is.factor(mess_result$mod))
  expect_equal(nrow(terra::levels(mess_result$mod)[[1]]), 5)
  expect_true(all(
    terra::levels(mess_result$mod)[[1]]$variable == names(climate_rast)
  ))

  # Check mos output (most similar variable)
  expect_is(mess_result$mos, "SpatRaster")
  expect_true(terra::is.factor(mess_result$mos))
  expect_equal(nrow(terra::levels(mess_result$mos)[[1]]), 5)
  expect_true(all(
    terra::levels(mess_result$mos)[[1]]$variable == names(climate_rast)
  ))

  # Test with full = TRUE
  expect_silent(mess_full <- mess(climate_rast, ref, full = TRUE))
  expect_is(mess_full, "list")
  expect_equal(length(mess_full), 4)
  expect_true(all(
    c("mess", "mess_by_variable", "mod", "mos") %in%
      names(mess_full)
  ))

  # Check mess_by_variable output
  expect_is(mess_full$mess_by_variable, "SpatRaster")
  expect_equal(terra::nlyr(mess_full$mess_by_variable), 5)
  expect_equal(names(mess_full$mess_by_variable), names(climate_rast))

  # Check that mess values are consistent between full and non-full
  expect_equal(terra::values(mess_result$mess), terra::values(mess_full$mess))
  expect_equal(terra::values(mess_result$mod), terra::values(mess_full$mod))
  expect_equal(terra::values(mess_result$mos), terra::values(mess_full$mos))
})

test_that("mess runs with data.frame input", {
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
  ref$ID <- NULL # Remove ID column added by terra::extract

  # Test with full = FALSE
  expect_silent(mess_result <- mess(climate_df[, -c(1, 2)], ref, full = FALSE))
  expect_is(mess_result, "list")
  expect_equal(length(mess_result), 3)
  expect_true(all(c("mess", "mod", "mos") %in% names(mess_result)))

  # Check mess output
  expect_is(mess_result$mess, "numeric")
  expect_equal(length(mess_result$mess), nrow(climate_df))

  # Check mod output (most dissimilar variable)
  expect_is(mess_result$mod, "integer")
  expect_equal(length(mess_result$mod), nrow(climate_df))
  expect_true(all(mess_result$mod >= 1 & mess_result$mod <= 5))

  # Check mos output (most similar variable)
  expect_is(mess_result$mos, "integer")
  expect_equal(length(mess_result$mos), nrow(climate_df))
  expect_true(all(mess_result$mos >= 1 & mess_result$mos <= 5))

  # Test with full = TRUE
  expect_silent(mess_full <- mess(climate_df[, -c(1, 2)], ref, full = TRUE))
  expect_is(mess_full, "list")
  expect_equal(length(mess_full), 4)
  expect_true(all(
    c("mess", "mess_by_variable", "mod", "mos") %in%
      names(mess_full)
  ))

  # Check mess_by_variable output
  expect_is(mess_full$mess_by_variable, "matrix")
  expect_equal(ncol(mess_full$mess_by_variable), 5)
  expect_equal(nrow(mess_full$mess_by_variable), nrow(climate_df))
  expect_equal(colnames(mess_full$mess_by_variable), names(climate_rast))

  # Check that mess values are consistent between full and non-full
  expect_equal(mess_result$mess, mess_full$mess)
  expect_equal(mess_result$mod, mess_full$mod)
  expect_equal(mess_result$mos, mess_full$mos)
})

test_that("mess runs with matrix input", {
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
  ref$ID <- NULL # Remove ID column added by terra::extract

  # Test with matrix input
  expect_silent(mess_result <- mess(climate_mat, ref, full = FALSE))
  expect_is(mess_result, "list")
  expect_true(all(c("mess", "mod", "mos") %in% names(mess_result)))
  expect_is(mess_result$mess, "numeric")
  expect_equal(length(mess_result$mess), nrow(climate_mat))
})

test_that("mess handles missing variables correctly", {
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

  # Test with subset of variables in x
  expect_warning(
    mess_result <- mess(climate_rast[[1:3]], ref, full = FALSE),
    "The following variables are missing from x"
  )
  expect_is(mess_result, "list")
  expect_equal(nrow(terra::levels(mess_result$mod)[[1]]), 3)

  # Test with subset of variables in ref
  expect_warning(
    mess_result <- mess(climate_rast, ref[, 1:3], full = FALSE),
    "The following variables are missing from ref"
  )
  expect_is(mess_result, "list")
  expect_equal(nrow(terra::levels(mess_result$mod)[[1]]), 3)

  # Test with no overlapping variables
  names(ref) <- paste0("var_", 1:5)
  expect_error(
    expect_warning(
      mess(climate_rast, ref, full = FALSE),
      "missing from x"
    ),
    "No variables in common"
  )
})

test_that("mess consistency between SpatRaster and data.frame", {
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
  ref$ID <- NULL # Remove ID column added by terra::extract

  # Run MESS on both input types
  expect_silent(mess_rast <- mess(climate_rast, ref, full = TRUE))
  expect_silent(mess_df <- mess(climate_df, ref, full = TRUE))

  # Extract non-NA values from raster results
  mess_rast_vals <- as.vector(terra::values(mess_rast$mess, na.rm = TRUE))
  mod_rast_vals <- as.vector(terra::values(mess_rast$mod, na.rm = TRUE))
  mos_rast_vals <- as.vector(terra::values(mess_rast$mos, na.rm = TRUE))

  # Compare with data.frame results
  expect_equal(mess_rast_vals, mess_df$mess)
  expect_equal(mod_rast_vals, mess_df$mod)
  expect_equal(mos_rast_vals, mess_df$mos)

  # Compare mess_by_variable
  mess_by_var_rast <- terra::values(mess_rast$mess_by_variable, na.rm = TRUE)
  expect_equal(mess_by_var_rast, mess_df$mess_by_variable)
})

test_that("mess with filename argument", {
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

  # Create temporary file
  temp_file <- tempfile(fileext = ".tif")

  # Test with filename
  expect_silent(
    mess_result <- mess(climate_rast, ref, full = FALSE, filename = temp_file)
  )
  expect_is(mess_result, "list")
  expect_true(file.exists(temp_file))

  # Read the file and compare
  mess_from_file <- terra::rast(temp_file)
  expect_equal(terra::values(mess_result$mess), terra::values(mess_from_file))

  # Clean up
  unlink(temp_file)
})

test_that("mess handles empty input correctly", {
  TEST_DIRECTORY <- test_path("test_inputs")

  # Load climate rasters to get variable names
  climate_rast <- terra::rast(
    file.path(
      TEST_DIRECTORY,
      sprintf("bioclim_10m/wc2.1bio%02d.tif", c(1, 4, 7, 12, 15))
    )
  )

  # Load coordinates and extract reference data
  sample_moths <- read.csv(file.path(TEST_DIRECTORY, "sample_moths.csv"))
  ref <- terra::extract(climate_rast, sample_moths[, c("lon", "lat")])
  ref$ID <- NULL # Remove ID column added by terra::extract

  # Create empty data frame with correct column names
  empty_df <- as.data.frame(matrix(numeric(0), ncol = 5))
  names(empty_df) <- names(climate_rast)

  # Test with empty data frame
  expect_silent(mess_result <- mess(empty_df, ref, full = FALSE))
  expect_is(mess_result, "list")
  expect_equal(length(mess_result$mess), 0)
  expect_equal(length(mess_result$mod), 0)
  expect_equal(length(mess_result$mos), 0)

  # Test with full = TRUE
  expect_silent(mess_full <- mess(empty_df, ref, full = TRUE))
  expect_is(mess_full, "list")
  expect_equal(nrow(mess_full$mess_by_variable), 0)
  expect_equal(ncol(mess_full$mess_by_variable), 5)
  expect_equal(length(mess_full$mess), 0)
})
