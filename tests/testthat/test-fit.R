context("Fitting")

set.seed(1)

test_that("gaussian model - symmetric works", {
  fitted = fit(create_data(fishdist, asymmetric_model = FALSE))
    expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("gaussian model - asymmetric works", {
  fitted = fit(create_data(fishdist, asymmetric_model = TRUE))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("t model - symmetric works", {
  fitted = fit(create_data(fishdist, asymmetric_model = FALSE, est_t_model = TRUE))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("t model - asymmetric works", {
  fitted = fit(create_data(fishdist, asymmetric_model = TRUE, est_t_model = TRUE))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("gaussian model - symmetric no slope sigma", {
  fitted = fit(create_data(fishdist, asymmetric_model = FALSE, est_sigma_trend = FALSE))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("gaussian model - asymmetric no slope sigma", {
  fitted = fit(create_data(fishdist, asymmetric_model = TRUE, est_sigma_trend = FALSE))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("gaussian model - symmetric no slope mu", {
  fitted = fit(create_data(fishdist, asymmetric_model = FALSE, est_mu_trend = FALSE))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("gaussian model - asymmetric no slope mu", {
  fitted = fit(create_data(fishdist, asymmetric_model = TRUE, est_mu_trend = FALSE))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

