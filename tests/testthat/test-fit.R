context("Fitting")

test_that("gaussian model - symmetric works - 1 year", {
  set.seed(1)

  fitted <- fit(create_data(fishdist[which(fishdist$year == 1980), ],
                            asymmetric_model = FALSE
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year == 1980), ],
                            asymmetric_model = FALSE, family = "poisson"
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year == 1980), ],
                            asymmetric_model = FALSE,
                            family = "negbin"
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("student-t model - symmetric works - 1 year", {
  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year == 1980), ],
                            asymmetric_model = FALSE,
                            est_sigma_trend = FALSE,
                            est_mu_trend = FALSE,
                            tail_model = "student_t"
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year == 1980), ],
                            asymmetric_model = FALSE,
                            est_sigma_trend = FALSE,
                            est_mu_trend = FALSE,
                            tail_model = "student_t",
                            family = "poisson"
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year == 1980), ],
                            asymmetric_model = FALSE,
                            est_sigma_trend = FALSE,
                            est_mu_trend = FALSE,
                            tail_model = "student_t",
                            family = "negbin"
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("gnorm model - symmetric works - 1 year", {
  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year == 1980), ],
                            asymmetric_model = FALSE,
                            est_sigma_trend = FALSE,
                            est_mu_trend = FALSE,
                            tail_model = "gnorm"
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year == 1980), ],
                            asymmetric_model = FALSE,
                            est_sigma_trend = FALSE,
                            est_mu_trend = FALSE,
                            tail_model = "gnorm",
                            family = "poisson"
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year == 1980), ],
                            asymmetric_model = FALSE,
                            est_sigma_trend = FALSE,
                            est_mu_trend = FALSE,
                            tail_model = "gnorm",
                            family = "negbin"
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("gaussian model - symmetric works - 10 years - mu trend", {
  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year > 1970), ],
                            asymmetric_model = FALSE,
                            est_sigma_trend = FALSE,
                            est_mu_trend = TRUE,
                            est_sigma_re = FALSE,
                            est_mu_re = TRUE
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year > 1970), ],
                            asymmetric_model = FALSE,
                            est_sigma_trend = FALSE,
                            est_mu_trend = TRUE,
                            est_sigma_re = FALSE,
                            est_mu_re = TRUE, family =
                              "poisson"
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year > 1970), ],
                            asymmetric_model = FALSE,
                            est_sigma_trend = FALSE,
                            est_mu_trend = TRUE,
                            est_sigma_re = FALSE,
                            est_mu_re = TRUE,
                            family = "negbin"
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("gaussian model - symmetric works - 10 years - no trend but with re", {
  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year > 1970), ],
                            asymmetric_model = FALSE,
                            est_sigma_trend = FALSE,
                            est_mu_trend = FALSE,
                            est_sigma_re = FALSE,
                            est_mu_re = TRUE
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year > 1970), ],
                            asymmetric_model = FALSE,
                            est_sigma_trend = FALSE,
                            est_mu_trend = FALSE,
                            est_sigma_re = FALSE,
                            est_mu_re = TRUE, family =
                              "poisson"
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year > 1970), ],
                            asymmetric_model = FALSE,
                            est_sigma_trend = FALSE,
                            est_mu_trend = FALSE,
                            est_sigma_re = FALSE,
                            est_mu_re = TRUE,
                            family = "negbin"
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("gaussian model - symmetric works - 3 years - no trend but with re", {
  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year > 1977), ],
                            asymmetric_model = FALSE,
                            est_sigma_trend = FALSE,
                            est_mu_trend = FALSE,
                            est_sigma_re = TRUE,
                            est_mu_re = FALSE
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year > 1977), ],
                            asymmetric_model = FALSE,
                            est_sigma_trend = FALSE,
                            est_mu_trend = FALSE,
                            est_sigma_re = TRUE,
                            est_mu_re = FALSE,
                            family = "poisson"
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year > 1977), ],
                            asymmetric_model = FALSE,
                            est_sigma_trend = FALSE,
                            est_mu_trend = FALSE,
                            est_sigma_re = TRUE,
                            est_mu_re = FALSE,
                            family = "negbin"
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("gaussian model - symmetric works - 3 years - trend with re", {
  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year > 1977), ],
                            asymmetric_model = FALSE,
                            est_sigma_trend = TRUE,
                            est_mu_trend = FALSE,
                            est_sigma_re = TRUE,
                            est_mu_re = FALSE
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year > 1977), ],
                            asymmetric_model = FALSE,
                            est_sigma_trend = TRUE,
                            est_mu_trend = FALSE,
                            est_sigma_re = TRUE,
                            est_mu_re = FALSE,
                            family = "poisson"
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year > 1977), ],
                            asymmetric_model = FALSE,
                            est_sigma_trend = TRUE,
                            est_mu_trend = FALSE,
                            est_sigma_re = TRUE,
                            est_mu_re = FALSE,
                            family = "negbin"
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})


test_that("gaussian model - asymmetric works - 3 years - no trend but with re", {
  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year > 1977), ],
                            asymmetric_model = TRUE,
                            est_sigma_trend = FALSE,
                            est_mu_trend = FALSE,
                            est_sigma_re = TRUE,
                            est_mu_re = FALSE
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year > 1977), ],
                            asymmetric_model = TRUE,
                            est_sigma_trend = FALSE,
                            est_mu_trend = FALSE,
                            est_sigma_re = TRUE,
                            est_mu_re = FALSE,
                            family = "poisson"
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year > 1977), ],
                            asymmetric_model = TRUE,
                            est_sigma_trend = FALSE,
                            est_mu_trend = FALSE,
                            est_sigma_re = TRUE,
                            est_mu_re = FALSE,
                            family = "negbin"
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})


test_that("gaussian model - asymmetric works - 10 years - trend with re", {
  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year > 1970), ],
                            asymmetric_model = TRUE,
                            est_sigma_trend = TRUE,
                            est_mu_trend = FALSE,
                            est_sigma_re = TRUE,
                            est_mu_re = FALSE
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year > 1970), ],
                            asymmetric_model = TRUE,
                            est_sigma_trend = FALSE,
                            est_mu_trend = FALSE,
                            est_sigma_re = TRUE,
                            est_mu_re = FALSE,
                            family = "poisson"
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(fishdist[which(fishdist$year > 1970), ],
                            asymmetric_model = TRUE,
                            est_sigma_trend = FALSE,
                            est_mu_trend = FALSE,
                            est_sigma_re = TRUE,
                            est_mu_re = FALSE,
                            family = "negbin"
  ),
  silent = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})
