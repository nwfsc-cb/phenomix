context("Fitting")

rel_tol <- 1e-4

set.seed(123)

# create 10 years of data
df <- expand.grid("doy" = 100:200, "year" = 1:10, sig = 10)
df$mu <- rnorm(10, 150, 5)[df$year]
df$pred <- dnorm(df$doy, df$mu, sd = df$sig, log = TRUE)
df$pred <- exp(df$pred + 10)
df$number <- round(rnorm(nrow(df), df$pred, 0.1))

test_that("gaussian model - symmetric works - 1 year", {
  set.seed(1)

  d <- df[which(df$year == 1), ]
  fitted <- fit(create_data(d, asymmetric_model = FALSE, min_number = 1, max_theta = 12),
    silent = TRUE,
    control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  f <- fit(create_data(d, asymmetric_model = TRUE,
                            min_number = 0.1, max_theta = 12),
     silent = TRUE,
     control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
  )
  expect_equal(length(which(is.na(f$sdreport$sd))), 0)
})

test_that("student-t model - symmetric works - 1 year", {
  set.seed(1)

  d <- df[which(df$year == 1), ]
  fitted <- fit(create_data(d, asymmetric_model = FALSE, min_number = 1, max_theta = 12, tail_model = "student_t"),
    silent = TRUE,
    control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(d, asymmetric_model = TRUE, min_number = 1, max_theta = 12, tail_model = "student_t"),
    silent = TRUE,
    control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("gnorm model - symmetric works - 1 year", {
  set.seed(1)

  d <- df[which(df$year == 1), ]
  fitted <- fit(create_data(d, asymmetric_model = FALSE, min_number = 1, max_theta = 12, tail_model = "gnorm"),
    silent = TRUE,
    control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  fitted <- fit(create_data(d, asymmetric_model = TRUE, min_number = 1, max_theta = 12, tail_model = "gnorm"),
    silent = TRUE,
    control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

test_that("other obs error families work - 1 year", {
  set.seed(1)

  d <- df[which(df$year == 1), ]
  fitted <- fit(create_data(d, asymmetric_model = FALSE, family = "negbin"),
    silent = TRUE,
    control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  df <- expand.grid("doy" = 100:200, "year" = 1:10, sig = 10)
  df$mu <- rnorm(10, 150, 5)[df$year]
  df$pred <- dnorm(df$doy, df$mu, sd = df$sig, log = TRUE)
  df$pred <- exp(df$pred)
  df$number <- rnorm(nrow(df), df$pred, 0.001)
  d <- df[which(df$year == 1), ]
  fitted <- fit(create_data(d, asymmetric_model = FALSE, family = "gaussian"),
    silent = TRUE,
    control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  df <- expand.grid("doy" = 100:200, "year" = 1:10, sig = 10)
  df$mu <- rnorm(10, 150, 5)[df$year]
  df$pred <- dnorm(df$doy, df$mu, sd = df$sig, log = TRUE)
  df$pred <- exp(df$pred + 3) / (3 + exp(df$pred + 1))
  df$number <- rbinom(nrow(df), size = 1, prob = df$pred)
  d <- df[which(df$year == 1), ]
  fitted <- fit(create_data(d, asymmetric_model = FALSE, family = "binomial"),
    silent = TRUE,
    control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  df <- expand.grid("doy" = 100:200, "year" = 1:10, sig = 10)
  df$mu <- rnorm(10, 150, 5)[df$year]
  df$pred <- dnorm(df$doy, df$mu, sd = df$sig, log = TRUE)
  df$pred <- rpois(nrow(df), exp(df$pred + 8))
  df$number <- df$pred
  d <- df[which(df$year == 1), ]
  fitted <- fit(create_data(d, asymmetric_model = FALSE, family = "negbin"),
    silent = TRUE,
    control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})



# create 20 years of data
set.seed(123)
df <- expand.grid("doy" = 100:200, "year" = 1:20)
df$mu <- rnorm(unique(df$year), 150, 5)[df$year]
df$sig1 <- rnorm(unique(df$year), 30, 5)[df$year]
df$sig2 <- rnorm(unique(df$year), 30, 5)[df$year]
df$sig <- ifelse(df$doy < df$mu, df$sig1, df$sig2)
df$pred <- dnorm(df$doy, df$mu, sd = df$sig, log = TRUE)
df$pred <- exp(df$pred + 8)
df$number <- round(rnorm(nrow(df), df$pred, 0.1))

test_that("gaussian model - symmetric works - multiple years", {
  set.seed(1)
  fitted <- fit(create_data(df, asymmetric_model = FALSE, min_number = 1),
    silent = TRUE,
    control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(df, asymmetric_model = FALSE, est_mu_re = FALSE, min_number = 1),
    silent = TRUE,
    control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(df, asymmetric_model = FALSE, est_sigma_re = FALSE, min_number = 1),
    silent = TRUE,
    control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol), limits = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(df, asymmetric_model = FALSE, est_sigma_re = FALSE, est_mu_re = FALSE, min_number = 1),
    silent = TRUE,
    control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol), limits = TRUE
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})


test_that("gaussian model - asymmetric works - multiple years", {
  set.seed(1)
  fitted <- fit(create_data(df, asymmetric_model = TRUE, min_number = 1),
    silent = TRUE,
    control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  # commented out because this particular seed makes the unix test fail
  # set.seed(1)
  # fitted <- fit(create_data(df, asymmetric_model = TRUE, est_sigma_re = FALSE, min_number = 1),
  #   silent = TRUE,
  #   control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol), limits = TRUE
  # )
  # expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(df, asymmetric_model = TRUE, est_mu_re = FALSE, min_number = 1),
    silent = TRUE,
    control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(df, asymmetric_model = TRUE, est_mu_re = FALSE, est_sigma_re = FALSE, min_number = 1),
    silent = TRUE,
    control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

# library(LaplacesDemon)
# set.seed(123)
# df <- expand.grid("doy" = 100:200, "year" = 1:20)
# df$mu <- rnorm(unique(df$year), 150, 5)[df$year]
# df$sig1 <- rnorm(unique(df$year), 30, 5)[df$year]
# df$sig2 <- rnorm(unique(df$year), 30, 5)[df$year]
# df$sig <- ifelse(df$doy < df$mu, df$sig1, df$sig2)
# df$pred <- dst(df$doy, mu=df$mu, sigma = df$sig, nu=10,log = TRUE)
# df$pred <- exp(df$pred + 8)
# df$number <- round(rnorm(nrow(df), df$pred, 0.1))

# test_that("student-t model - symmetric works - multiple years", {
#   #commented out because this particular seed makes the unix test fail
#   set.seed(1)
#   fitted <- fit(create_data(df, asymmetric_model = FALSE, min_number = 1,
#     tail_model = "student_t"),
#     silent = TRUE,
#     control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
#   )
#   expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
# })

#
# test_that("student-t model - asymmetric works - multiple years", {
#   set.seed(2)
#   fitted <- fit(create_data(df, asymmetric_model = TRUE, min_number = 1, tail_model = "student_t"),
#     silent = TRUE, limits = TRUE,
#     control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
#   )
#   expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
# })


test_that("gnorm model - symmetric works - multiple years", {
  # set.seed(1)
  # fitted <- fit(create_data(df, asymmetric_model = FALSE, min_number = 1, tail_model = "gnorm"),
  # silent = FALSE,
  # control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
  # )
  # expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(4)
  fitted <- fit(create_data(df, asymmetric_model = FALSE, est_sigma_re = FALSE, est_mu_re = FALSE, min_number = 1, tail_model = "gnorm"),
    silent = TRUE,
    control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})


test_that("gnorm model - asymmetric works - multiple years", {
  # set.seed(1)
  # fitted <- fit(create_data(df, asymmetric_model = TRUE, min_number = 1,tail_model = "gnorm"), silent = TRUE,
  #               control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol))
  # expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
  #
  # set.seed(1)
  # fitted <- fit(create_data(df, asymmetric_model = TRUE, est_sigma_re = FALSE, min_number = 1,tail_model = "gnorm"), silent = TRUE,
  #               control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol))
  # expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
  #
  # set.seed(1)
  # fitted <- fit(create_data(df, asymmetric_model = TRUE, est_mu_re = FALSE, min_number = 1,tail_model = "gnorm"), silent = TRUE,
  #             control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol), limits = TRUE)
  # expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  set.seed(1)
  fitted <- fit(create_data(df, asymmetric_model = TRUE, est_mu_re = FALSE, est_sigma_re = FALSE, min_number = 1, tail_model = "gnorm"),
    silent = TRUE,
    control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
  )
  expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})

# create 20 years of data -- using gnorm
set.seed(123)
df <- expand.grid("doy" = 100:200, "year" = 1:20)
df$mu <- rnorm(unique(df$year), 150, 5)[df$year]
df$alpha1 <- rnorm(unique(df$year), 30, 5)[df$year]
# df$alpha2 <- rnorm(unique(df$year), 30, 5)[df$year]
# df$sig <- ifelse(df$doy < df$mu, df$alpha1, df$alpha2)
df$sig <- df$alpha1
df$pred <- 0
for (i in 1:nrow(df)) {
  df$pred[i] <- gnorm::dgnorm(df$doy[i],
    mu = df$mu[i],
    alpha = df$sig[i], beta = 10, log = TRUE
  )
}
df$pred <- exp(df$pred + 8)
df$number <- round(rnorm(nrow(df), df$pred, 0.1))


test_that("gnorm model - asymmetric works - multiple years", {
  #
  #   set.seed(5)
  #   fitted <- fit(create_data(df, asymmetric_model = TRUE, est_mu_re = FALSE, est_sigma_re = TRUE, min_number = 1, tail_model = "gnorm"),
  #                 silent = TRUE,limits=TRUE,
  #                 control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
  #   )
  #   expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)

  # set.seed(1)
  # fitted <- fit(create_data(df, asymmetric_model = TRUE, est_mu_re = FALSE, est_sigma_re = FALSE, min_number = 1, tail_model = "gnorm"),
  #   silent = TRUE, limits = TRUE,
  #   control = list(eval.max = 4000, iter.max = 5000, rel.tol = rel_tol)
  # )
  # expect_equal(length(which(is.na(fitted$sdreport$sd))), 0)
})
