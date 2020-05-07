context("Fitting")

test_that("gaussian model - symmetric works", {
  set.seed(1)
  fitted = fit(create_data(fishdist, asymmetric_model = FALSE), silent=TRUE)
    expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("gaussian model - asymmetric works", {
  set.seed(1)
  fitted = fit(create_data(fishdist, asymmetric_model = TRUE), silent=TRUE)
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("t model - symmetric works", {
  set.seed(1)
  fitted = fit(create_data(fishdist, asymmetric_model = FALSE, tail_model = "student_t"), silent=TRUE)
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("t model - asymmetric works", {
  set.seed(2)
  fitted = fit(create_data(fishdist, asymmetric_model = TRUE, tail_model = "student_t"), silent=TRUE)
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("gaussian model - symmetric no slope sigma", {
  fitted = fit(create_data(fishdist, asymmetric_model = FALSE, est_sigma_trend = FALSE), silent=TRUE)
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("gaussian model - asymmetric no slope sigma", {
  set.seed(1)
  fitted = fit(create_data(fishdist, asymmetric_model = TRUE, est_sigma_trend = FALSE), silent=TRUE)
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("gaussian model - symmetric no slope mu", {
  set.seed(1)
  fitted = fit(create_data(fishdist, asymmetric_model = FALSE, est_mu_trend = FALSE), silent=TRUE)
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("gaussian model - asymmetric no slope mu", {
  set.seed(1)
  fitted = fit(create_data(fishdist, asymmetric_model = TRUE, est_mu_trend = FALSE), silent=TRUE)
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

