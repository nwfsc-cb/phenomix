context("Fitting")

set.seed(123)

# create 10 years of data
df = expand.grid("doy" = 100:200, "year" = 1:10, sig = 10)
df$mu = rnorm(10, 150, 5)[df$year]
df$pred = dnorm(df$doy, df$mu, sd = df$sig, log=TRUE)
df$pred = exp(df$pred + 10)
df$number = round(rnorm(nrow(df), df$pred, 0.1))

test_that("gaussian model - symmetric works - 1 year", {
  set.seed(1)

  d = df[which(df$year == 1), ]
  fitted <- fit(create_data(d, asymmetric_model = FALSE), silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(d, asymmetric_model = FALSE), silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  # expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  # set.seed(1)
  # fitted <- fit(create_data(d,
  #   asymmetric_model = FALSE,
  #   family = "negbin"
  # ),
  # silent = TRUE
  # )
  # expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("other obs error families work - 1 year", {
  set.seed(1)

  d = df[which(df$year == 1), ]
  fitted <- fit(create_data(d, asymmetric_model = FALSE, family="negbin"), silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(d, asymmetric_model = TRUE, family="negbin"), silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(d, asymmetric_model = FALSE, family="poisson"), silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(d, asymmetric_model = FALSE, family="poisson"), silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
  # expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  # set.seed(1)
  # fitted <- fit(create_data(d,
  #   asymmetric_model = FALSE,
  #   family = "negbin"
  # ),
  # silent = TRUE
  # )
  # expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})


test_that("gaussian model - symmetric works - multiple years", {
  set.seed(1)

  fitted <- fit(create_data(df, asymmetric_model = FALSE), silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(df, asymmetric_model = TRUE), silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(df, asymmetric_model = TRUE, est_sigma_re = FALSE), silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(df, asymmetric_model = TRUE, est_mu_re = FALSE), silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(df, asymmetric_model = FALSE, est_mu_re = FALSE), silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(df, asymmetric_model = FALSE, est_sigma_re = FALSE), silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(df, asymmetric_model = FALSE, est_sigma_re = FALSE, est_mu_re = FALSE), silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
  # set.seed(1)
  # fitted <- fit(create_data(d,
  #   asymmetric_model = FALSE,
  #   family = "negbin"
  # ),
  # silent = TRUE
  # )
  # expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})


test_that("gaussian model works - multiple years", {
  set.seed(1)

  fitted <- fit(create_data(df, asymmetric_model = FALSE), silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(df, asymmetric_model = TRUE), silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(df, asymmetric_model = TRUE, est_sigma_re = FALSE), silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(df, asymmetric_model = TRUE, est_mu_re = FALSE), silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(df, asymmetric_model = FALSE, est_mu_re = FALSE), silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(df, asymmetric_model = FALSE, est_sigma_re = FALSE), silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(df, asymmetric_model = FALSE, est_sigma_re = FALSE, est_mu_re = FALSE), silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  covar_data = data.frame(year = unique(df$year))
  covar_data$nyear = covar_data$year
  fitted <- fit(create_data(df, asymmetric_model = FALSE,
                        est_sigma_re = FALSE,
                        est_mu_re = FALSE,
                        mu = ~ nyear,
                        covar_data = covar_data),
                silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))

  covar_data$nyear =covar_data$year
  fitted <- fit(create_data(df, asymmetric_model = FALSE,
                            est_sigma_re = FALSE,
                            est_mu_re = FALSE,
                            sigma = ~ nyear,
                            covar_data = covar_data),
                silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))

  covar_data$nyear =covar_data$year
  fitted <- fit(create_data(df, asymmetric_model = FALSE,
                            est_sigma_re = TRUE,
                            est_mu_re = FALSE,
                            sigma = ~ nyear,
                            covar_data = covar_data),
                silent = TRUE,
                control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7))
  # set.seed(1)
  # fitted <- fit(create_data(d,
  #   asymmetric_model = FALSE,
  #   family = "negbin"
  # ),
  # silent = TRUE
  # )
  # expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})


test_that("student-t model - symmetric works - 1 year", {
  set.seed(1)
  d = df[which(df$year == 1), ]
  fitted <- fit(create_data(d,
    asymmetric_model = FALSE,
    tail_model = "student_t"
  ),
  silent = TRUE, control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  d = df[which(df$year == 1), ]
  fitted <- fit(create_data(d,
                            asymmetric_model = TRUE,
                            tail_model = "student_t"
  ),
  silent = TRUE, control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)


  fitted <- fit(create_data(df,
                            asymmetric_model = FALSE,
                            tail_model = "student_t"
  ),
  silent = TRUE, control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(d,
                            asymmetric_model = FALSE,
                            tail_model = "student_t"
  ),
  silent = TRUE, control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(df,
                            asymmetric_model = TRUE,
                            tail_model = "student_t"
  ),
  silent = TRUE, control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)


})

test_that("gnorm model - symmetric works - 1 year", {
  set.seed(1)
  d = df[which(df$year == 1), ]
  fitted <- fit(create_data(d,
                            asymmetric_model = FALSE,
                            tail_model = "gnorm"
  ),
  silent = TRUE, control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  # set.seed(1)
  # d = df[which(df$year == 1), ]
  # fitted <- fit(create_data(d,
  #                           asymmetric_model = TRUE,
  #                           tail_model = "gnorm"
  # ),
  # silent = TRUE, control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7)
  # )
  # expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(2)
  fitted <- fit(create_data(df,
                            asymmetric_model = FALSE,
                            tail_model = "gnorm"
  ),
  silent = TRUE, control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(df,
                            asymmetric_model = TRUE,
                            tail_model = "gnorm",
                            est_mu_re = FALSE, est_sigma_re = FALSE
  ),
  silent = TRUE, control = list(eval.max = 4000, iter.max = 5000, rel.tol = 1e-7)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})
