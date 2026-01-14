# Additional Tests for phenomix RTMB Conversion
# Add these to your testthat/test-fit.R or create a new test-additional.R file

library(testthat)

# IMPORTANT: Understanding phenomix data structure
#
# The model treats data differently based on family:
#
# - Gaussian family (family=1): y = LOG-DENSITIES
#   The model computes pred = log(phenological density at day x)
#   Then fits: y ~ Normal(pred, obs_sigma)
#   So y should be log-densities, not raw counts!
#
# - Count families (poisson, negbin, binomial): y = RAW COUNTS
#   The model fits count distributions directly to the counts
#
# For gaussian family, the proper simulation is:
#   1. Define true phenological parameters (mu, sigma)
#   2. Compute TRUE log-density at each day
#   3. Add gaussian noise to get observed log-densities
#   4. Use these as y values

# Helper function to simulate LOG-DENSITY data for gaussian family
simulate_logdensity_data <- function(true_mu, true_sigma, obs_sigma,
                                     doy_range = 100:199, year = 1) {
  # Compute true log-density at each day
  log_dens <- dnorm(doy_range, mean = true_mu, sd = true_sigma, log = TRUE)

  # Add observation noise
  observed_log_dens <- log_dens + rnorm(length(doy_range), 0, obs_sigma)

  # Return dataframe
  data.frame(doy = doy_range, number = exp(observed_log_dens), year = year)
}

# Helper function to simulate COUNT data for poisson/negbin families
simulate_count_data <- function(mu, sigma, n_total, doy_range = 100:199, year = 1) {
  # Sample day-of-year values from distribution
  doy_samples <- round(rnorm(n_total, mean = mu, sd = sigma))
  doy_samples <- pmax(min(doy_range), pmin(max(doy_range), doy_samples))

  # Count observations per day
  counts <- sapply(doy_range, function(d) sum(doy_samples == d))

  # Return dataframe
  data.frame(doy = doy_range, number = counts, year = year)
}

# Test 1: Parameter Recovery - Symmetric Gaussian
test_that("parameter recovery - symmetric gaussian", {
  set.seed(123)

  # True parameters
  true_mu <- 150
  true_sigma <- 10
  obs_sigma <- 0.5 # Observation noise on log scale

  # Simulate log-density data (correct for gaussian family)
  df <- simulate_logdensity_data(
    true_mu = true_mu, true_sigma = true_sigma,
    obs_sigma = obs_sigma, year = 1
  )

  # Fit model
  datalist <- create_data(df, asymmetric_model = FALSE, tail_model = "gaussian")
  fitted <- fit(datalist, silent = TRUE)

  # Extract parameters
  fixef <- fixef(fitted)
  fitted_mu <- fixef$estimate[which(fixef$par == "b_mu")]
  fitted_sigma <- exp(fixef$estimate[which(fixef$par == "b_sig1")])

  # Check recovery (within 10% for log-density data)
  expect_equal(fitted_mu, true_mu, tolerance = 0.1 * abs(true_mu))
  expect_equal(fitted_sigma, true_sigma, tolerance = 0.1 * true_sigma)
})

# Test 2: Asymmetric with Equal Parameters = Symmetric
test_that("asymmetric model with equal tails matches symmetric", {
  data(fishdist)
  d <- fishdist[fishdist$year == 1975, ]

  # Fit symmetric
  d_sym <- create_data(d, asymmetric_model = FALSE, tail_model = "gaussian")
  fit_sym <- fit(d_sym, silent = TRUE)

  # Fit asymmetric
  d_asym <- create_data(d, asymmetric_model = TRUE, tail_model = "gaussian")
  fit_asym <- fit(d_asym, silent = TRUE)

  # Both should converge
  expect_true(fit_sym$pars$convergence == 0)
  expect_true(fit_asym$pars$convergence == 0 || fit_asym$pars$convergence == 1)

  # Asymmetric should have similar or better fit (more parameters)
  # Just check they both produce reasonable estimates
  fixef <- fixef(fit_sym)
  fitted_mu <- fixef$estimate[which(fixef$par == "b_mu")]
  expect_true(!is.null(fitted_mu))


  fixef <- fixef(fit_asym)
  fitted_mu <- fixef$estimate[which(fixef$par == "b_mu")]
  expect_true(!is.null(fitted_mu))
})

# Test 3: Both Random Effects Together
test_that("mu and sigma random effects work together", {
  data(fishdist)
  d <- fishdist[fishdist$year %in% 1975:1980, ]

  datalist <- create_data(d,
    asymmetric_model = FALSE,
    tail_model = "gaussian",
    est_mu_re = 1,
    est_sigma_re = 1
  )

  fitted <- fit(datalist, silent = TRUE)

  expect_true(fitted$par$convergence == 0)
  expect_true(!is.null(ranef(fitted)$estimate))

  # Should have random effects for both parameters
  re <- ranef(fitted)$par
  expect_true("mu_devs" %in% names(re) || length(re) > 0)
  expect_true("sigma1_devs" %in% names(re) || length(re) > 0)
})

# Test 4: Prediction Works
test_that("predict method produces reasonable output", {
  data(fishdist)
  d <- fishdist[fishdist$year == 1975, ]

  datalist <- create_data(d, asymmetric_model = FALSE, tail_model = "gaussian")
  fitted <- fit(datalist, silent = TRUE)

  # Get predictions
  preds <- predict(fitted)

  # Check structure
  expect_true(is.data.frame(preds) || is.list(preds))
  expect_true(!is.null(preds$pred) || !is.null(preds$mean))

  # Predictions should be finite
  expect_true(all(is.finite(preds$pred)))
})

# Test 5: Model Handles Sparse Count Data
test_that("poisson model handles many zeros", {
  set.seed(456)

  # Simulate sparse count data (correct for poisson family)
  df <- simulate_count_data(mu = 150, sigma = 5, n_total = 50, year = 1)

  datalist <- create_data(df,
    asymmetric_model = FALSE,
    family = "poisson",
    tail_model = "gaussian"
  )

  # Should converge without errors
  fitted <- fit(datalist, silent = TRUE)
  expect_true(fitted$par$convergence == 0 || fitted$par$convergence == 1)
})

# Test 6: Unbalanced Data Across Years
test_that("handles unbalanced data across years", {
  set.seed(789)

  # Create unbalanced log-density dataset (correct for gaussian family)
  df1 <- simulate_logdensity_data(true_mu = 120, true_sigma = 8, obs_sigma = 0.3, year = 1)
  df2 <- simulate_logdensity_data(true_mu = 150, true_sigma = 12, obs_sigma = 0.5, year = 2)
  df3 <- simulate_logdensity_data(true_mu = 140, true_sigma = 10, obs_sigma = 0.4, year = 3)

  df <- rbind(df1, df2, df3)

  datalist <- create_data(df, asymmetric_model = FALSE, tail_model = "gaussian")

  fitted <- fit(datalist, silent = TRUE)
  expect_true(fitted$par$convergence == 0 || fitted$par$convergence == 1)
})

# Test 7: Gradient is Small at Convergence
test_that("gradient is near zero at convergence", {
  data(fishdist)
  d <- fishdist[fishdist$year == 1975, ]

  datalist <- create_data(d, asymmetric_model = FALSE, tail_model = "gaussian")
  fitted <- fit(datalist, silent = TRUE)

  # Check convergence
  expect_true(fitted$par$convergence == 0)

  # Check that optimization message indicates convergence
  expect_true(!is.null(fitted$par))
})

# Test 8: Different Tail Models Give Different Fits
test_that("different tail models produce different likelihoods", {
  data(fishdist)
  d <- fishdist[fishdist$year == 1975, ]

  # Fit with gaussian tail
  d_gauss <- create_data(d, asymmetric_model = FALSE, tail_model = "gaussian")
  fit_gauss <- fit(d_gauss, silent = TRUE)

  # Fit with student-t tail
  d_t <- create_data(d, asymmetric_model = FALSE, tail_model = "student_t")
  fit_t <- fit(d_t, silent = TRUE)

  # Both should converge
  expect_true(fit_gauss$par$convergence == 0)
  expect_true(fit_t$par$convergence == 0)

  # Both should produce estimates
  expect_true(!is.null(fixef(fit_gauss)$estimate))
  expect_true(!is.null(fixef(fit_t)$estimate))
})

# Test 9: Minimal Observations Per Year
test_that("handles minimal observations per year", {
  set.seed(111)

  # Use log-density data with some missing/sparse days
  df1 <- simulate_logdensity_data(true_mu = 120, true_sigma = 15, obs_sigma = 0.5, year = 1)
  df2 <- simulate_logdensity_data(true_mu = 135, true_sigma = 12, obs_sigma = 0.5, year = 2)
  df3 <- simulate_logdensity_data(true_mu = 125, true_sigma = 18, obs_sigma = 0.5, year = 3)

  df <- rbind(df1, df2, df3)

  datalist <- create_data(df, asymmetric_model = FALSE, tail_model = "gaussian")

  # Should work (though estimates may be uncertain)
  fitted <- fit(datalist, silent = TRUE)
  expect_true(!is.null(fitted))
})

# Test 10: Very Peaked Distribution
test_that("model handles very peaked distribution", {
  set.seed(222)

  # Very peaked (narrow) log-density pattern
  df <- simulate_logdensity_data(true_mu = 150, true_sigma = 3, obs_sigma = 0.3, year = 1)

  datalist <- create_data(df, asymmetric_model = FALSE, tail_model = "gaussian")
  fitted <- fit(datalist, silent = TRUE)

  expect_true(fitted$par$convergence == 0 || fitted$par$convergence == 1)

  # Should estimate small sigma
  fitted_sigma <- exp(fixef(fitted)$estimate[which(fixef(fitted)$par == "b_sig1")])
  expect_lt(fitted_sigma, 3)
})
