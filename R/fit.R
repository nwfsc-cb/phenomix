#' Fitting function to be called by user
#'
#' This function creates a list of parameters, sets up RTMB object and attempts to
#' do fitting / estimation
#'
#' @param data_list A list of data, as output from create_data
#' @param silent Boolean passed to RTMB::MakeADFun, whether to be verbose or not (defaults to FALSE)
#' @param inits Optional named list of parameters for starting values, defaults to NULL
#' @param control Optional control list for stats::nlminb. For arguments see ?nlminb. Defaults to eval.max=2000, iter.max=1000, rel.tol=1e-10. For final model runs, the rel.tol should be even smaller
#' @param limits Whether to include limits for stats::nlminb. Can be a list of (lower, upper), or TRUE to use suggested hardcoded limits. Defaults to NULL,
#' where no limits used.
#' @param fit_model Whether to fit the model. If not, returns a list including the data, parameters, and initial values. Defaults to TRUE
#' @importFrom stats runif rnorm
#' @importFrom methods is
#' @importFrom RTMB MakeADFun sdreport
#' @export
#' @examples
#' \dontrun{
#' data(fishdist)
#'
#' # example of fitting fixed effects, no trends, no random effects
#' set.seed(1)
#' datalist <- create_data(fishdist[which(fishdist$year > 1970), ],
#'   asymmetric_model = FALSE,
#'   est_mu_re = FALSE, est_sigma_re = FALSE
#' )
#' fit <- fit(datalist)
#' }
fit <- function(data_list,
                silent = FALSE,
                inits = NULL,
                control = list(eval.max = 2000, iter.max = 1000, rel.tol = 1e-10),
                limits = NULL, # can also be a list, or TRUE
                fit_model = TRUE) {
  # create list of parameter starting values -- used in both
  # asymmetric and symmetric models
  parameters <- list(
    theta = rnorm(n = data_list$nLevels, log(mean(data_list$y[which(!is.na(data_list$y))]))),
    b_mu = rep(0, ncol(data_list$mu_mat)),
    log_sigma_mu_devs = 0,
    mu_devs = rep(0, data_list$nLevels),
    b_sig1 = rep(1, ncol(data_list$sig_mat)),
    b_sig2 = rep(1, ncol(data_list$sig_mat)),
    log_sigma1_sd = 0,
    sigma1_devs = rep(0, data_list$nLevels),
    log_sigma2_sd = 0,
    sigma2_devs = rep(0, data_list$nLevels),
    log_obs_sigma = 0.0,
    log_tdf_1 = 0,
    log_tdf_2 = 0,
    log_beta_1 = 0.1,
    log_beta_2 = 0.1
  )
  parameters$b_mu[1] <- mean(data_list$x[which(!is.na(data_list$y))])

  if (data_list$family == 1) {
    # for gaussian, don't log-transform theta
    parameters$theta <- rep(0, data_list$nLevels)
  }

  # Mapping off params as needed:
  rtmb_map <- list()
  if (data_list$family %in% c(2, 4)) {
    # don't include obs_sigma for poisson or binomial
    rtmb_map <- c(rtmb_map, list(log_obs_sigma = factor(NA)))
  }

  if (data_list$asymmetric == 0) {
    # map off pars not needed
    rtmb_map <- c(rtmb_map, list(
      b_sig2 = factor(rep(NA, ncol(data_list$sig_mat))),
      log_sigma2_sd = factor(NA),
      sigma2_devs = factor(rep(NA, data_list$nLevels))
    ))
  }

  if (data_list$tail_model == 0) {
    # then fit gaussian model and map off both parameters
    rtmb_map <- c(rtmb_map, list(
      log_tdf_1 = factor(NA),
      log_tdf_2 = factor(NA),
      log_beta_1 = factor(NA),
      log_beta_2 = factor(NA)
    ))
  }
  if (data_list$tail_model == 1) {
    # fit the t-model
    if (data_list$asymmetric == 0) {
      # then map off the tdf2, because model is symmetric
      rtmb_map <- c(rtmb_map, list(
        log_tdf_2 = factor(NA),
        log_beta_1 = factor(NA),
        log_beta_2 = factor(NA)
      ))
    } else {
      rtmb_map <- c(rtmb_map, list(
        log_beta_1 = factor(NA),
        log_beta_2 = factor(NA)
      ))
    }
  }
  if (data_list$tail_model == 2) {
    # then fit gnorm model and map off t parameters
    if (data_list$asymmetric == 0) {
      # then map off the beta2, because model is symmetric
      rtmb_map <- c(rtmb_map, list(
        log_beta_2 = factor(NA),
        log_tdf_1 = factor(NA),
        log_tdf_2 = factor(NA)
      ))
    } else {
      rtmb_map <- c(rtmb_map, list(
        log_tdf_1 = factor(NA),
        log_tdf_2 = factor(NA)
      ))
    }
  }
  if (data_list$asymmetric == 1) {
    if (data_list$share_shape == 1) {
      if (data_list$tail_model == 1) {
        # map off 2nd nu parameter
        rtmb_map <- c(rtmb_map, list(
          log_tdf_2 = factor(NA)
        ))
      }
      if (data_list$tail_model == 2) {
        rtmb_map <- c(rtmb_map, list(
          log_beta_2 = factor(NA)
        ))
      }
    }
  }

  if (data_list$nLevels == 1) {
    data_list$est_mu_re <- 0
    data_list$est_sigma_re <- 0
  }

  # if don't estimate mean random effects map them off
  if (data_list$est_mu_re == 0) {
    rtmb_map <- c(rtmb_map, list(
      log_sigma_mu_devs = factor(NA),
      mu_devs = factor(rep(NA, data_list$nLevels))
    ))
  }

  if (data_list$est_sigma_re == 0) {
    rtmb_map <- c(rtmb_map, list(
      log_sigma1_sd = factor(NA),
      log_sigma2_sd = factor(NA),
      sigma1_devs = factor(rep(NA, data_list$nLevels)),
      sigma2_devs = factor(rep(NA, data_list$nLevels))
    ))
  } else {
    if (data_list$asymmetric == 0) {
      rtmb_map <- c(rtmb_map, list(
        log_sigma2_sd = factor(NA),
        sigma2_devs = factor(rep(NA, data_list$nLevels))
      ))
    }
  }

  random <- ""
  if (data_list$est_mu_re == 1) {
    random <- c(random, "mu_devs")
  }
  if (data_list$est_sigma_re == 1) {
    random <- c(random, "sigma1_devs")
    if (data_list$asymmetric == TRUE) random <- c(random, "sigma2_devs")
  }
  random <- random[-1]
  if (length(random) == 0) random <- NULL

  # Pre-extract scalar values from data_list that will be used in the objective
  # This avoids indexing issues during RTMB tape construction
  nu_prior_shape <- data_list$nu_prior[1]
  nu_prior_rate <- data_list$nu_prior[2]
  beta_prior_shape <- data_list$beta_prior[1]
  beta_prior_rate <- data_list$beta_prior[2]

  # Add these to data_list
  data_list$nu_prior_shape <- nu_prior_shape
  data_list$nu_prior_rate <- nu_prior_rate
  data_list$beta_prior_shape <- beta_prior_shape
  data_list$beta_prior_rate <- beta_prior_rate

  # Capture data_list in a closure
  data_for_objective <- data_list
  obj_func <- function(pars) {
    phenomix_objective(pars, data_for_objective)
  }

  # Create RTMB objective function
  obj <- RTMB::MakeADFun(
    func = obj_func,
    parameters = parameters,
    map = rtmb_map,
    random = random,
    silent = silent
  )

  mod_list <- list(
    obj = obj,
    init_vals = obj$par,
    data_list = data_list
  )

  if (fit_model == TRUE) {
    # Attempt to do estimation
    if (!is.null(inits)) {
      init <- inits
    } else {
      init <- obj$par
    }
    if (is.null(limits)) {
      pars <- stats::nlminb(
        start = init,
        objective = obj$fn,
        gradient = obj$gr,
        control = control
      )
    } else {
      if (is(limits, "list")) {
        lower_limits <- limits$lower
        upper_limits <- limits$upper
      } else {
        lim <- limits(parnames = names(obj$par), max_theta = data_list$max_theta)
        lower_limits <- lim$lower
        upper_limits <- lim$upper
      }

      pars <- stats::nlminb(
        start = init, objective = obj$fn,
        gradient = obj$gr, control = control,
        lower = lower_limits,
        upper = upper_limits
      )
    }

    sdreport <- RTMB::sdreport(obj)
    mod_list <- c(mod_list, list(
      pars = pars, sdreport = sdreport,
      rtmb_map = rtmb_map, rtmb_random = random,
      rtmb_parameters = parameters
    ))
  }

  return(structure(mod_list, class = "phenomix"))
}
