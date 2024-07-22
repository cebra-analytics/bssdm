context("rangebag")

test_that("predicted output consistent with edmaps::range_bag", {
  TEST_DIRECTORY <- test_path("test_inputs")
  climate_rast <- terra::rast(file.path(TEST_DIRECTORY,
                                        sprintf("bioclim_10m/wc2.1bio%02d.tif",
                                                c(1,4,7,12,15))))
  climate_df <- as.data.frame(climate_rast, xy = TRUE, na.rm = FALSE)
  names(climate_df)[1:2] <- c("lon", "lat")
  sample_moths <- read.csv(file.path(TEST_DIRECTORY, "sample_moths.csv"))
  expected_rangebag <- terra::rast(file.path(TEST_DIRECTORY,
                                             "expected_rangebag.tif"))
  set.seed(1234)
  expect_silent(sdm.model <- rangebag(climate_rast, sample_moths))
  expect_is(sdm.model, "Rangebag")
  expect_silent(bsrb_output <- predict(sdm.model, climate_rast,
                                       raw_output = FALSE))
  expect_equal(round(bsrb_output[][,1], 6), round(expected_rangebag[][,1], 6))
  expect_warning(rangebag(climate_rast[[1]], sample_moths),
                 "Rangebag x data has fewer variables than n_dim.")
  expect_silent(bsrb_output <- predict(sdm.model, climate_rast,
                                       raw_output = FALSE))
  expect_true(all(is.finite(bsrb_output[][,1]) ==
                    is.finite(expected_rangebag[][,1])))
  set.seed(1234)
  expect_silent(sdm.model_df <- rangebag(climate_df, sample_moths))
  expect_is(sdm.model_df, "Rangebag")
  expect_equal(slot(sdm.model_df, "ch_models"), slot(sdm.model, "ch_models"))
  expect_silent(bsrb_output_df <- predict(sdm.model_df, climate_rast,
                                       raw_output = FALSE))
  expect_equal(round(bsrb_output_df[][,1], 6),
               round(expected_rangebag[][,1], 6))
  alt_ext <- c(1, 3, 40.5, 42.5)
  climate_rast_alt <- terra::extend(terra::crop(climate_rast, alt_ext),
                                    terra::ext(alt_ext))
  expected_rangebag_alt <- terra::extend(terra::crop(expected_rangebag,
                                                     alt_ext),
                                         terra::ext(alt_ext))
  expect_silent(bsrb_output_alt <- predict(sdm.model, climate_rast_alt,
                                           raw_output = FALSE))
  expect_true(terra::ext(bsrb_output_alt) == terra::ext(climate_rast_alt))
})
