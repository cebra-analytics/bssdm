context("climatch")

test_that("predicted output consistent with Climatch site", {
  TEST_DIRECTORY <- test_path("test_inputs")
  world_climate_bc <- read.csv(file.path(TEST_DIRECTORY,
                                         "world_climate_bc.csv"))
  world_stations_sd <- read.csv(file.path(TEST_DIRECTORY,
                                          "world_stations_sd.csv"))
  sample_moths <- read.csv(file.path(TEST_DIRECTORY, "sample_moths.csv"))
  au_climate_sa <- read.csv(file.path(TEST_DIRECTORY, "au_climate_sa.csv"))
  expected_climatch <- read.csv(file.path(TEST_DIRECTORY,
                                          "expected_climatch.csv"))
  # 1. Euclidean method
  expect_silent(sdm.model <- climatch(world_climate_bc, sample_moths,
                                      algorithm = "euclidean", d_max = 100,
                                      sd_data = world_stations_sd[,2]))
  expect_is(sdm.model, "Climatch")
  expect_silent(bscm_output <- predict(sdm.model, au_climate_sa,
                                       raw_output = FALSE))
  expect_equal(bscm_output$predicted, expected_climatch$test_target_1)
  # 2. Closest standard score method
  expect_silent(sdm.model <- climatch(world_climate_bc, sample_moths,
                                      algorithm = "closest_standard_score",
                                      d_max = 100,
                                      sd_data = world_stations_sd[,2]))
  expect_silent(bscm_output <- predict(sdm.model, au_climate_sa,
                                       raw_output = FALSE))
  expect_equal(bscm_output$predicted, expected_climatch$test_target_2)
})
