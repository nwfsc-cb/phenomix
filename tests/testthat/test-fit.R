context("Fitting")

set.seed(1)

test_that("gaussian model - symmetric works", {
  set.seed(1)
  fitted = fit(create_data(fishdist, asymmetric_model = FALSE))
    expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("gaussian model - asymmetric works", {
  set.seed(1)
  fitted = fit(create_data(fishdist, asymmetric_model = TRUE))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("t model - symmetric works", {
  set.seed(1)
  fitted = fit(create_data(fishdist, asymmetric_model = FALSE, est_t_model = TRUE))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("t model - asymmetric works", {
  set.seed(1)
  fitted = fit(create_data(fishdist, asymmetric_model = TRUE, est_t_model = TRUE))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})
