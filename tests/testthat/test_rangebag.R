context("rangebag")

test_that("predicted output consistent with edmaps::range_bag", {
  TEST_DIRECTORY <- test_path("test_inputs")
  climate_rast <- terra::rast(file.path(TEST_DIRECTORY,
                                        sprintf("bioclim_10m/wc2.1bio%02d.tif",
                                                c(1,4,7,12,15))))
  sample_moths <- read.csv(file.path(TEST_DIRECTORY, "sample_moths.csv"))
  expected_rangebag <- terra::rast(file.path(TEST_DIRECTORY,
                                             "expected_rangebag.tif"))
  set.seed(1234)
  expect_silent(sdm.model <- rangebag(climate_rast, sample_moths))
  expect_is(sdm.model, "Rangebag")
  expect_silent(bsrb_output <- predict(sdm.model, climate_rast,
                                       raw_output = FALSE))
  expect_equal(round(bsrb_output[][,1], 6), round(expected_rangebag[][,1], 6))
  expect_warning(sdm.model <- rangebag(climate_rast[[1]], sample_moths),
                 "Rangebag x data has fewer variables than n_dim.")
  expect_silent(bsrb_output <- predict(sdm.model, climate_rast,
                                       raw_output = FALSE))
  expect_true(all(is.finite(bsrb_output[][,1]) ==
                    is.finite(expected_rangebag[][,1])))
})
