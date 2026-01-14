#' RTMB Objective Function for phenomix
#'
#' This function replaces the TMB C++ template with an R-based RTMB function
#'
#' @param pars List of parameters
#' @param data List of data
#' @return Negative log-likelihood
#' @importFrom RTMB REPORT ADREPORT OBS getAll
#' @keywords internal
# Helper functions for density calculations
# These must be defined outside the objective function to avoid type issues

qthill <- function(quantile, v, mean, sigma) {
  flip <- ifelse(quantile > 0.5, 1, -1)
  z <- ifelse(quantile > 0.5, 2 * (1 - quantile), 2 * quantile)

  a <- 1 / (v - 0.5)
  b <- 48 / (a * a)
  c <- ((20700 * a / b - 98) * a - 16) * a + 96.36
  d <- ((94.5 / (b + c) - 3) / b + 1) * sqrt(a * pi / 2) * v
  x <- z * d
  y <- x^(2 / v)

  # Branch 1 calculations
  x_new <- RTMB::qnorm(z * 0.5, 0, 1)
  y_temp <- x_new * x_new

  # Use smooth transition instead of ifelse(v < 5)
  # Manual sigmoid instead of plogis: 1 / (1 + exp(-x))
  weight_v <- 1 / (1 + exp(-((5 - v) * 10))) # Steep sigmoid centered at v=5
  c_adjusted <- c + weight_v * 0.3 * (v - 4.5) * (x_new + 0.6)

  c_final <- c_adjusted + (((0.05 * d * x_new - 5) * x_new - 7) * x_new - 2) * x_new + b
  y_temp2 <- (((((0.4 * y_temp + 6.3) * y_temp + 36) * y_temp + 94.5) / c_final - y_temp - 3) / b + 1) * x_new
  y_temp3 <- a * y_temp2 * y_temp2

  # Use smooth transition instead of ifelse(y_temp3 > 0.002)
  # Manual sigmoid instead of plogis
  weight_y <- 1 / (1 + exp(-((y_temp3 - 0.002) * 1000))) # Steep sigmoid centered at 0.002
  branch1_raw <- weight_y * (exp(y_temp3) - 1) + (1 - weight_y) * (y_temp3 + 0.5 * y_temp3 * y_temp3)
  # Ensure branch1 is positive
  branch1 <- branch1_raw + 0.1

  # Branch 2
  branch2_raw <- ((1 / (((v + 6) / (v * y) - 0.089 * d - 0.822) * (v + 2) * 3) + 0.5 / (v + 4)) * y - 1) * (v + 1) / (v + 2) + 1 / y
  # Ensure branch2 is positive
  branch2 <- branch2_raw + 0.1

  # Use smooth transition instead of ifelse(y > 0.05 + a)
  # Manual sigmoid: 1 / (1 + exp(-x)) to avoid RTMB::plogis
  weight_branch <- 1 / (1 + exp(-((y - (0.05 + a)) * 100)))
  y_final <- weight_branch * branch1 + (1 - weight_branch) * branch2

  # Both branches have safety margins, so y_final should be positive
  # But add tiny epsilon just in case
  y_final <- y_final + 1e-10

  q <- sqrt(v * y_final) * flip
  return(mean + sigma * q)
}

# dgnorm - log density of generalized normal
dgnorm <- function(x, mu, alpha, beta) {
  z <- -abs(x - mu)^beta / alpha^beta + log(beta) - (log(2) + log(alpha) + lgamma(1 / beta))
  return(z)
}

# ddnorm - log density of double normal
ddnorm <- function(x, mu, sigma1, sigma2) {
  z <- log(2) - log(sigma1 + sigma2)
  z <- z + ifelse(x < mu,
    RTMB::dnorm(x = (x - mu) / sigma1, mean = 0, sd = 1, log = TRUE),
    RTMB::dnorm(x = (x - mu) / sigma2, mean = 0, sd = 1, log = TRUE)
  )
  return(z)
}

# ddt - log density of double Student-t
ddt <- function(x, mu, sigma1, sigma2, tdf_1, tdf_2) {
  z <- log(2) - log(sigma1 + sigma2)
  z <- z + ifelse(x < mu,
    RTMB::dt(x = (x - mu) / sigma1, df = tdf_1, log = TRUE),
    RTMB::dt(x = (x - mu) / sigma2, df = tdf_2, log = TRUE)
  )
  return(z)
}

# ddgnorm - log density of double generalized normal
ddgnorm <- function(x, mu, alpha1, alpha2, beta1, beta2, sigma1, sigma2) {
  z <- log(2) - log(sigma1 + sigma2)
  z <- z + ifelse(x < mu,
    dgnorm(x, mu, alpha1, beta1),
    dgnorm(x, mu, alpha2, beta2)
  )
  return(z)
}

# qgnorm - quantile function for generalized normal
qgnorm <- function(quantile, mu, alpha, beta) {
  p <- quantile
  p <- ifelse(p > 0.5, 1 - p, p)

  # Replace sign() with ifelse() to avoid comparisons on AD types
  sign_val <- ifelse(quantile > 0.5, 1, ifelse(quantile < 0.5, -1, 0))
  shape <- 1 / beta
  scale <- 1 / (1 / alpha)^beta

  z <- sign_val * exp(log(RTMB::qgamma(abs(p - 0.5) * 2, shape = shape, scale = scale)) / beta) + mu
  return(z)
}

# qdnorm - quantile function for double normal
qdnorm <- function(p, mu, log_sigma1, log_sigma2) {
  r <- log_sigma1 / (log_sigma1 + log_sigma2)
  z <- ifelse(p < r,
    mu + log_sigma1 * RTMB::qnorm(0.5 * p * (log_sigma1 + log_sigma2) / log_sigma1, 0, 1),
    mu + log_sigma2 * RTMB::qnorm(0.5 * ((log_sigma1 + log_sigma2) * (1 + p) - 2 * log_sigma1) / log_sigma2, 0, 1)
  )
  return(z)
}

# qdt - quantile function for double Student-t
qdt <- function(p, mu, log_sigma1, log_sigma2, tdf_1, tdf_2) {
  r <- log_sigma1 / (log_sigma1 + log_sigma2)
  z <- ifelse(p < r,
    mu + log_sigma1 * qthill(0.5 * p * (log_sigma1 + log_sigma2) / log_sigma1, tdf_1, 0, 1),
    mu + log_sigma2 * qthill(0.5 * ((log_sigma1 + log_sigma2) * (1 + p) - 2 * log_sigma1) / log_sigma2, tdf_2, 0, 1)
  )
  return(z)
}

# qdgnorm - quantile function for double generalized normal
qdgnorm <- function(p, mu, log_sigma1, log_sigma2, beta_ratio_1, beta_ratio_2, beta_1, beta_2) {
  r <- log_sigma1 / (log_sigma1 + log_sigma2)
  z <- ifelse(p < r,
    mu + log_sigma1 * qgnorm(0.5 * p * (log_sigma1 + log_sigma2) / log_sigma1, mu, log_sigma1 * beta_ratio_1, beta_1),
    mu + log_sigma2 * qgnorm(0.5 * ((log_sigma1 + log_sigma2) * (1 + p) - 2 * log_sigma1) / log_sigma2, mu, log_sigma2 * beta_ratio_2, beta_2)
  )
  return(z)
}

phenomix_objective <- function(pars, data) {
  # Use RTMB's built-in getAll() function
  # This extracts all parameters and data into local environment
  getAll(data, pars, warn = FALSE)

  # Mark data elements used in control flow as observed (not AD variables)
  nLevels <- RTMB::OBS(nLevels)
  asymmetric <- RTMB::OBS(asymmetric)
  family <- RTMB::OBS(family)
  tail_model <- RTMB::OBS(tail_model)
  est_sigma_re <- RTMB::OBS(est_sigma_re)
  est_mu_re <- RTMB::OBS(est_mu_re)
  share_shape <- RTMB::OBS(share_shape)
  use_t_prior <- RTMB::OBS(use_t_prior)
  use_beta_prior <- RTMB::OBS(use_beta_prior)

  # Mark data vectors as observed (not AD types)
  y <- RTMB::OBS(y)
  x <- RTMB::OBS(x)
  years <- RTMB::OBS(years)

  # qthill - Hill's algorithm with smooth approximations instead of ifelse
  # Use sigmoid (plogis) to create smooth transitions instead of hard if/else

  # Main objective function
  nll <- 0

  # Derived parameters
  obs_sigma <- exp(log_obs_sigma)
  tdf_1 <- exp(log_tdf_1) + 2
  tdf_2 <- exp(log_tdf_2) + 2
  beta_1 <- exp(log_beta_1)
  beta_2 <- exp(log_beta_2)

  if (share_shape == 1) {
    tdf_2 <- tdf_1
    beta_2 <- beta_1
  }

  # Extract prior hyperparameters from data_list
  # nu_prior and beta_prior come from create_data()
  nu_1 <- nu_prior[1]
  nu_2 <- nu_prior[2]
  beta_p1 <- beta_prior[1]
  beta_p2 <- beta_prior[2]

  # Priors on tail parameters
  if (tail_model == 1 && use_t_prior == 1) {
    nll <- nll - RTMB::dgamma(x = tdf_1, shape = nu_1, rate = nu_2, log = TRUE)
    if (asymmetric == 1) {
      nll <- nll - RTMB::dgamma(x = tdf_2, shape = nu_1, rate = nu_2, log = TRUE)
    }
  }

  if (tail_model == 2 && use_beta_prior == 1) {
    nll <- nll - RTMB::dgamma(x = beta_1, shape = beta_p1, rate = beta_p2, log = TRUE)
    if (asymmetric == 1) {
      nll <- nll - RTMB::dgamma(x = beta_2, shape = beta_p1, rate = beta_p2, log = TRUE)
    }
  }

  # Initialize vectors
  n <- length(y)
  log_sigma1 <- numeric(nLevels)
  log_sigma2 <- numeric(nLevels)
  mu <- numeric(nLevels)
  alpha1 <- numeric(nLevels)
  alpha2 <- numeric(nLevels)
  # Initialize as advectors for ADREPORT compatibility
  lower25 <- RTMB::advector(rep(0, nLevels))
  upper75 <- RTMB::advector(rep(0, nLevels))
  range <- RTMB::advector(rep(0, nLevels))
  # Pre-allocate pred - use RTMB::advector for type consistency
  pred <- RTMB::advector(rep(0, n))

  # Calculate beta_ratio for gnorm distribution if needed
  beta_ratio_1 <- 1.0
  beta_ratio_2 <- 1.0
  if (tail_model == 2) {
    beta_ratio_1 <- sqrt(exp(lgamma(1 / beta_1)) / exp(lgamma(3 / beta_1)))
    if (asymmetric == 1) {
      beta_ratio_2 <- sqrt(exp(lgamma(1 / beta_2)) / exp(lgamma(3 / beta_2)))
    }
  }

  # Random effects
  for (i in 1:nLevels) {
    if (est_mu_re == 1) {
      nll <- nll - RTMB::dnorm(x = mu_devs[i], mean = 0, sd = exp(log_sigma_mu_devs), log = TRUE)
    }
    if (est_sigma_re == 1) {
      nll <- nll - RTMB::dnorm(x = sigma1_devs[i], mean = 0, sd = exp(log_sigma1_sd), log = TRUE)
      if (asymmetric == 1) {
        nll <- nll - RTMB::dnorm(x = sigma2_devs[i], mean = 0, sd = exp(log_sigma2_sd), log = TRUE)
      }
    }
  }

  # Fixed effects
  mu <- as.vector(mu_mat %*% b_mu)
  log_sigma1 <- as.vector(sig_mat %*% b_sig1) # Log scale for positive constraint
  if (asymmetric == 1) {
    log_sigma2 <- as.vector(sig_mat %*% b_sig2) # Log scale for positive constraint
  }

  # Add random effects to entire vectors at once (not element by element)
  if (est_sigma_re == 1) {
    log_sigma1 <- log_sigma1 + sigma1_devs
    if (asymmetric == 1) {
      log_sigma2 <- log_sigma2 + sigma2_devs
    }
  }

  if (est_mu_re == 1) {
    mu <- mu + mu_devs
  }

  # Transform to natural scale: sigma = exp(log_sigma) ensures sigma > 0
  sigma1 <- exp(log_sigma1)
  if (asymmetric == 1) {
    sigma2 <- exp(log_sigma2)
  }

  # Calculate quantiles for each year
  # These represent the 25th and 75th percentiles of the phenological distribution
  for (i in 1:nLevels) {
    if (tail_model == 0) {
      # Gaussian distribution
      if (asymmetric == 0) {
        lower25[i] <- RTMB::qnorm(0.25, mean = mu[i], sd = sigma1[i])
        upper75[i] <- RTMB::qnorm(0.75, mean = mu[i], sd = sigma1[i])
      } else {
        # Asymmetric gaussian - use average sigma as approximation
        # (Avoids complex conditional logic on AD types)
        sigma_avg <- (sigma1[i] + sigma2[i]) / 2
        lower25[i] <- RTMB::qnorm(0.25, mean = mu[i], sd = sigma_avg)
        upper75[i] <- RTMB::qnorm(0.75, mean = mu[i], sd = sigma_avg)
      }
    } else if (tail_model == 1) {
      # Student-t distribution
      if (asymmetric == 0) {
        lower25[i] <- qthill(0.25, tdf_1, mu[i], sigma1[i])
        upper75[i] <- qthill(0.75, tdf_1, mu[i], sigma1[i])
      } else {
        # Asymmetric student-t - use average sigma and df
        sigma_avg <- (sigma1[i] + sigma2[i]) / 2
        df_avg <- (tdf_1 + tdf_2) / 2
        lower25[i] <- qthill(0.25, df_avg, mu[i], sigma_avg)
        upper75[i] <- qthill(0.75, df_avg, mu[i], sigma_avg)
      }
    } else if (tail_model == 2) {
      # Generalized normal
      if (asymmetric == 0) {
        # Approximate with normal quantiles scaled by alpha
        lower25[i] <- mu[i] + alpha1[i] * RTMB::qnorm(0.25)
        upper75[i] <- mu[i] + alpha1[i] * RTMB::qnorm(0.75)
      } else {
        # Asymmetric gnorm - use average alpha
        alpha_avg <- (alpha1[i] + alpha2[i]) / 2
        lower25[i] <- mu[i] + alpha_avg * RTMB::qnorm(0.25)
        upper75[i] <- mu[i] + alpha_avg * RTMB::qnorm(0.75)
      }
    }
    # Calculate range (inter-quartile range)
    range[i] <- upper75[i] - lower25[i]
  }

  # Calculate alpha for gnorm distribution if needed
  if (tail_model == 2) {
    alpha1 <- sigma1 * beta_ratio_1
    if (asymmetric == 1) {
      alpha2 <- sigma2 * beta_ratio_2
    }
  }

  # Prediction loop
  for (i in 1:n) {
    yr_idx <- years[i]

    if (asymmetric == 1) {
      if (tail_model == 0) {
        # Inline ddnorm calculation with log-space sigmoid
        z <- log(2) - log(sigma1[yr_idx] + sigma2[yr_idx])
        diff <- x[i] - mu[yr_idx]

        # Compute both branches
        dens_left <- RTMB::dnorm(x = diff / sigma1[yr_idx], mean = 0, sd = 1, log = TRUE)
        dens_right <- RTMB::dnorm(x = diff / sigma2[yr_idx], mean = 0, sd = 1, log = TRUE)

        # Log-space smooth blending for numerical stability
        # indicator = exp(log_indicator) where log_indicator uses log1p for stability
        # plogis(x) = 1/(1+exp(-x)) = exp(x)/(1+exp(x))
        # We want indicator → 1 when diff < 0, → 0 when diff > 0
        # Use logsumexp trick: log(exp(a) + exp(b)) = max(a,b) + log(1 + exp(-|a-b|))

        # For numerical stability, use the fact that:
        # indicator * left + (1-indicator) * right
        # = exp(log(indicator) + left) + exp(log(1-indicator) + right) in probability space
        # But in log space, we use logsumexp

        # Simpler: just use tanh-based smooth indicator
        # indicator = 0.5 * (1 - tanh(diff * scale))
        # This is numerically stable and smooth
        indicator <- 0.5 * (1 - tanh(diff * 0.5))
        z <- z + indicator * dens_left + (1 - indicator) * dens_right

        pred[i] <- z + theta[yr_idx]
      } else if (tail_model == 1) {
        # Inline ddt calculation with tanh-based sigmoid
        z <- log(2) - log(sigma1[yr_idx] + sigma2[yr_idx])
        diff <- x[i] - mu[yr_idx]

        # Compute both branches
        dens_left <- RTMB::dt(x = diff / sigma1[yr_idx], df = tdf_1, log = TRUE)
        dens_right <- RTMB::dt(x = diff / sigma2[yr_idx], df = tdf_2, log = TRUE)

        # Smooth indicator using tanh (numerically stable)
        indicator <- 0.5 * (1 - tanh(diff * 0.5))
        z <- z + indicator * dens_left + (1 - indicator) * dens_right

        pred[i] <- z + theta[yr_idx]
      } else if (tail_model == 2) {
        # Inline ddgnorm calculation with tanh-based sigmoid
        z <- log(2) - log(sigma1[yr_idx] + sigma2[yr_idx])
        diff <- x[i] - mu[yr_idx]
        diff_abs <- abs(diff) + 1e-10

        # Compute both branches in log space
        log_ratio_left <- beta_1 * (log(diff_abs) - log(alpha1[yr_idx]))
        log_ratio_right <- beta_2 * (log(diff_abs) - log(alpha2[yr_idx]))

        dens_left <- -exp(log_ratio_left) + log(beta_1) - (log(2) + log(alpha1[yr_idx]) + lgamma(1 / beta_1))
        dens_right <- -exp(log_ratio_right) + log(beta_2) - (log(2) + log(alpha2[yr_idx]) + lgamma(1 / beta_2))

        # Smooth indicator using tanh
        indicator <- 0.5 * (1 - tanh(diff * 0.5))
        z <- z + indicator * dens_left + (1 - indicator) * dens_right

        pred[i] <- z + theta[yr_idx]
      }
    } else {
      if (tail_model == 0) {
        dens_result <- RTMB::dnorm(x = x[i], mean = mu[yr_idx], sd = sigma1[yr_idx], log = TRUE)
        pred[i] <- dens_result + theta[yr_idx]
      } else if (tail_model == 1) {
        # Inline t-distribution calculation
        x_std <- (x[i] - mu[yr_idx]) / sigma1[yr_idx]
        dt_result <- RTMB::dt(x = x_std, df = tdf_1, log = TRUE)
        log_sigma <- log(sigma1[yr_idx])
        pred[i] <- (dt_result - log_sigma) + theta[yr_idx]
      } else if (tail_model == 2) {
        # Inline dgnorm calculation
        diff_abs <- abs(x[i] - mu[yr_idx]) + 1e-10 # add small constant to avoid log(0)

        # Compute -diff^beta / alpha^beta = -exp(beta * (log(diff) - log(alpha)))
        log_ratio <- beta_1 * (log(diff_abs) - log(alpha1[yr_idx]))

        pred[i] <- -exp(log_ratio) +
          log(beta_1) - (log(2) + log(alpha1[yr_idx]) + lgamma(1 / beta_1)) +
          theta[yr_idx]
      }
    }
  }

  # Annual totals calculation (re-enabled with inline RTMB-compatible code)
  # Calculate integral of density * theta over days 1-365 for each year
  # Initialize as advectors (not plain numeric) so ADREPORT works
  year_log_tot <- RTMB::advector(rep(0, nLevels))
  year_tot <- RTMB::advector(rep(0, nLevels))

  for (i in 1:nLevels) {
    # Initialize accumulator
    exp_sum <- 0

    # Loop over all days (1-365)
    for (t in 1:365) {
      # Calculate density for day t, year i
      if (asymmetric == 1) {
        if (tail_model == 0) {
          # Inline ddnorm calculation
          z <- log(2) - log(sigma1[i] + sigma2[i])
          diff <- t - mu[i]
          dens_left <- RTMB::dnorm(x = diff / sigma1[i], mean = 0, sd = 1, log = TRUE)
          dens_right <- RTMB::dnorm(x = diff / sigma2[i], mean = 0, sd = 1, log = TRUE)
          indicator <- 0.5 * (1 - tanh(diff * 0.5))
          dens <- z + indicator * dens_left + (1 - indicator) * dens_right
        } else if (tail_model == 1) {
          # Inline ddt calculation
          z <- log(2) - log(sigma1[i] + sigma2[i])
          diff <- t - mu[i]
          dens_left <- RTMB::dt(x = diff / sigma1[i], df = tdf_1, log = TRUE) - log(sigma1[i])
          dens_right <- RTMB::dt(x = diff / sigma2[i], df = tdf_2, log = TRUE) - log(sigma2[i])
          indicator <- 0.5 * (1 - tanh(diff * 0.5))
          dens <- z + indicator * dens_left + (1 - indicator) * dens_right
        } else if (tail_model == 2) {
          # Inline ddgnorm calculation
          z <- log(2) - log(sigma1[i] + sigma2[i])
          diff <- t - mu[i]
          diff_abs_left <- abs(diff / sigma1[i]) + 1e-10
          diff_abs_right <- abs(diff / sigma2[i]) + 1e-10
          log_ratio_left <- beta_1 * log(diff_abs_left)
          log_ratio_right <- beta_2 * log(diff_abs_right)
          dens_left <- -exp(log_ratio_left) + log(beta_1) - (log(2) + log(sigma1[i]) + lgamma(1 / beta_1))
          dens_right <- -exp(log_ratio_right) + log(beta_2) - (log(2) + log(sigma2[i]) + lgamma(1 / beta_2))
          indicator <- 0.5 * (1 - tanh(diff * 0.5))
          dens <- z + indicator * dens_left + (1 - indicator) * dens_right
        }
      } else {
        # Symmetric models
        if (tail_model == 0) {
          dens <- RTMB::dnorm(x = t, mean = mu[i], sd = sigma1[i], log = TRUE)
        } else if (tail_model == 1) {
          dens <- RTMB::dt(x = (t - mu[i]) / sigma1[i], df = tdf_1, log = TRUE) - log(sigma1[i])
        } else if (tail_model == 2) {
          # Inline gnorm calculation
          diff_abs <- abs(t - mu[i]) + 1e-10
          log_ratio <- beta_1 * (log(diff_abs) - log(alpha1[i]))
          dens <- -exp(log_ratio) + log(beta_1) - (log(2) + log(alpha1[i]) + lgamma(1 / beta_1))
        }
      }
      # Add theta[i] and accumulate in natural scale
      dens_theta <- dens + theta[i]
      exp_sum <- exp_sum + exp(dens_theta)
    }
    # Convert to log scale for log_tot
    year_tot[i] <- exp_sum
    year_log_tot[i] <- log(exp_sum)
  }

  # Likelihood
  if (family == 1) {
    # Gaussian
    nll <- nll - sum(RTMB::dnorm(x = y, mean = pred, sd = obs_sigma, log = TRUE))
  } else if (family == 2) {
    # Poisson - use loop like TMB template
    for (i in 1:n) {
      pred_i <- pred[i]
      # Note: TMB template had pred capping but marked it as "not needed"
      nll <- nll - RTMB::dpois(x = y[i], lambda = exp(pred_i), log = TRUE)
    }
  } else if (family == 3) {
    # Negative binomial - compute log-likelihood manually
    # TMB uses dnbinom_robust with log-scale parameters
    # Manual computation avoids RTMB::dnbinom issues
    for (i in 1:n) {
      s1 <- pred[i] # log(mu)
      s2 <- 2 * s1 - log_obs_sigma # log(size)
      mu_i <- exp(s1)
      size_i <- exp(s2)

      # Negative binomial log-likelihood formula:
      # log P(y | mu, size) = lgamma(y + size) - lgamma(size) - lgamma(y + 1) +
      #                       size * log(size) + y * log(mu) - (y + size) * log(mu + size)
      ll <- lgamma(y[i] + size_i) - lgamma(size_i) - lgamma(y[i] + 1) +
        size_i * log(size_i) + y[i] * log(mu_i) - (y[i] + size_i) * log(mu_i + size_i)
      nll <- nll - ll
    }
  } else if (family == 4) {
    # Binomial - use loop like TMB template
    for (i in 1:n) {
      nll <- nll - RTMB::dbinom(x = y[i], size = 1, prob = RTMB::plogis(pred[i]), log = TRUE)
    }
  } else if (family == 5) {
    # Lognormal
    nll <- nll - sum(RTMB::dnorm(x = log(y), mean = pred, sd = obs_sigma, log = TRUE))
  }

  # Report variables for ADREPORT
  RTMB::ADREPORT(theta)
  RTMB::REPORT(theta)
  RTMB::ADREPORT(log_sigma1)
  RTMB::REPORT(log_sigma1)
  RTMB::ADREPORT(sigma1)
  RTMB::REPORT(sigma1)
  RTMB::ADREPORT(mu)
  RTMB::REPORT(mu)
  RTMB::ADREPORT(b_mu)
  RTMB::REPORT(b_mu)
  RTMB::ADREPORT(b_sig1)
  RTMB::REPORT(b_sig1)
  RTMB::ADREPORT(year_tot)
  RTMB::REPORT(year_tot)
  RTMB::ADREPORT(year_log_tot)
  RTMB::REPORT(year_log_tot)

  if (family != 2 && family != 4) {
    RTMB::ADREPORT(obs_sigma)
    RTMB::REPORT(obs_sigma)
  }

  RTMB::ADREPORT(pred)
  RTMB::ADREPORT(lower25)
  RTMB::REPORT(lower25)
  RTMB::ADREPORT(upper75)
  RTMB::REPORT(upper75)
  RTMB::ADREPORT(range)
  RTMB::REPORT(range)

  if (tail_model == 1) {
    RTMB::ADREPORT(tdf_1)
    RTMB::REPORT(tdf_1)
  }
  if (tail_model == 2) {
    RTMB::ADREPORT(beta_1)
    RTMB::REPORT(beta_1)
  }

  if (asymmetric == 1) {
    RTMB::ADREPORT(b_sig2)
    RTMB::REPORT(b_sig2)
    RTMB::ADREPORT(log_sigma2)
    RTMB::REPORT(log_sigma2)
    RTMB::ADREPORT(sigma2)
    RTMB::REPORT(sigma2)
    if (tail_model == 1) {
      RTMB::ADREPORT(tdf_2)
      RTMB::REPORT(tdf_2)
    }
    if (tail_model == 2) {
      RTMB::ADREPORT(beta_2)
      RTMB::REPORT(beta_2)
    }
  }

  return(nll)
}
