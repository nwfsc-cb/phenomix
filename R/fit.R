#' @useDynLib salmix, .registration = TRUE
NULL

#' Fitting function to be called by user
#'
#' This function creates a list of parameters, sets up TMB object and attempts to
#' do fitting / estimation
#'
#' @param data_list A list of data, as output from create_data
#' @param silent Boolean passed to TMB::MakeADFun, whether to be verbose or not (defaults to FALSE)
#' @param inits Optional named list of parameters for starting values, defaults to NULL
#' @param control Optional control list for stats::nlminb. For arguments see ?nlminb. Defaults to eval.max=2000, iter.max=1000, rel.tol=1e-10. For final model runs, the rel.tol should be even smaller
#' @param limits Whether to include limits for stats::nlminb, defaults to FALSE
#' @param ... Additional arguments to be passed to
#' @importFrom stats runif rnorm
#' @export
#' @examples
#' data(fishdist)
#'
#' # example of fitting fixed effects, no trends
#' # set.seed(1)
#' # datalist <- create_data(fishdist[which(fishdist$year > 1970),], asymmetric_model = FALSE,
#' #    est_sigma_trend=FALSE, est_mu_trend=FALSE)
#' # fit <- fit(datalist)
#' #
#' # # example of model with random effects in means
#' # set.seed(1)
#' # datalist <- create_data(fishdist[which(fishdist$year > 1970),], asymmetric_model = FALSE,
#' #                         est_sigma_trend=FALSE, est_mu_trend=FALSE, est_mu_re = TRUE)
#' # fit <- fit(datalist)
#' # # example of model with random effects in variances
#' # set.seed(1)
#' # datalist <- create_data(fishdist[which(fishdist$year > 1970),], asymmetric_model = FALSE,
#' #                         est_sigma_trend=FALSE, est_mu_trend=FALSE, est_sigma_re = TRUE)
#' # fit <- fit(datalist)
#' #
#' # # example of model with poisson response
#' # set.seed(1)
#' # datalist <- create_data(fishdist[which(fishdist$year > 1970),], asymmetric_model = FALSE,
#' #                         est_sigma_trend=FALSE, est_mu_trend=FALSE, est_mu_re = TRUE,
#' #                         family="poisson")
#' # fit <- fit(datalist)
fit <- function(data_list,
                silent = FALSE,
                inits = NULL,
                control = list(eval.max = 2000, iter.max = 1000, rel.tol = 1e-10),
                limits = FALSE) {

  # create list of parameter starting values -- used in both
  # asymmetric and symmetric models
  parameters <- list(
    sigma1_devs = rep(1, data_list$nLevels),
    theta = rnorm(n = data_list$nLevels, log(mean(data_list$y))),
    log_mu_b0 = log(mean(data_list$x)),
    mu_b1 = 0.005,
    log_sigma_mu_devs = 0.0,
    sig1_b0 = 0.0,
    sig1_b1 = 0.0,
    log_sigma1 = 0.0,
    log_obs_sigma = 0.0,
    sigma2_devs = rep(0, data_list$nLevels),
    sig2_b0 = 0.0,
    sig2_b1 = 0.0,
    log_sigma2 = 0.0,
    log_tdf_1 = 0,
    log_tdf_2 = 0,
    log_beta_1 = 0,
    log_beta_2 = 0
  )

  if (data_list$est_mu_re == 0) {
    # parameters are fixed effects
    parameters <- append(
      parameters,
      list(mu_devs = rnorm(n = data_list$nLevels, mean(data_list$x)))
    )
  } else {
    # parameters are deviations in log-space from the mean
    parameters <- append(
      parameters,
      list(mu_devs = rep(0, data_list$nLevels))
    )
  }

  # If inits is included, use that instead of parameters
  if (!is.null(inits)) parameters <- inits

  # Mapping off params as needed:
  tmb_map <- list()
  if (data_list$family == 2) {
    tmb_map <- c(tmb_map, list(log_obs_sigma = as.factor(NA)))
  }

  if (data_list$asymmetric == 0) {
    # map off pars not needed
    tmb_map <- c(tmb_map, list(
      log_sigma2 = as.factor(NA),
      sigma2_devs = rep(as.factor(NA), data_list$nLevels)
    ))
  }

  if (data_list$sig_trend == 0) {
    tmb_map <- c(tmb_map, list(
      sig1_b1 = as.factor(NA)
    ))
    if (data_list$asymmetric == 0) {
      tmb_map <- c(tmb_map, list(
        sig2_b1 = as.factor(NA)
      ))
    }
  }
  if (data_list$mu_trend == 0) {
    tmb_map <- c(tmb_map, list(mu_b1 = as.factor(NA)))
  }

  if (data_list$tail_model == 0) {
    # then fit gaussian model and map off both parameters
    tmb_map <- c(tmb_map, list(
      log_tdf_1 = as.factor(NA),
      log_tdf_2 = as.factor(NA),
      log_beta_1 = as.factor(NA),
      log_beta_2 = as.factor(NA)
    ))
  }
  if (data_list$tail_model == 1) {
    # fit the t-model
    if (data_list$asymmetric == 0) {
      # then map off the tdf2, because model is symmetric
      tmb_map <- c(tmb_map, list(
        log_tdf_2 = as.factor(NA), log_beta_1 = as.factor(NA),
        log_beta_2 = as.factor(NA)
      ))
    } else {
      tmb_map <- c(tmb_map, list(
        log_beta_1 = as.factor(NA),
        log_beta_2 = as.factor(NA)
      ))
    }
  }
  if (data_list$tail_model == 2) {
    # then fit gaussian model and map off both parameters
    if (data_list$asymmetric == 0) {
      # then map off the beta2, because model is symmetric
      tmb_map <- c(tmb_map, list(
        log_beta_2 = as.factor(NA), log_tdf_1 = as.factor(NA),
        log_tdf_2 = as.factor(NA)
      ))
    } else {
      tmb_map <- c(tmb_map, list(
        log_tdf_1 = as.factor(NA),
        log_tdf_2 = as.factor(NA)
      ))
    }
  }

  if (data_list$nLevels == 1) {
    data_list$est_mu_re <- 0
    data_list$est_sigma_re <- 0
  }

  # if don't estimate mean random effects map them off
  if (data_list$est_mu_re == 0) {
    tmb_map <- c(tmb_map, list(
      log_sigma_mu_devs = as.factor(NA),
      log_mu_b0 = as.factor(NA)
    ))
    # mu_devs = rep(as.factor(NA), data_list$nLevels))
  }
  # if don't estimate sd random effects map them off
  if (data_list$est_sigma_re == 0 | data_list$asymmetric == 0) {
    tmb_map <- c(tmb_map, list(sig2_b0 = as.factor(NA)))
  }
  if (data_list$est_sigma_re == 0) {
    tmb_map <- c(tmb_map, list(
      sig1_b0 = as.factor(NA),
      log_sigma1 = as.factor(NA),
      log_sigma2 = as.factor(NA),
      sigma2_devs = rep(as.factor(NA), data_list$nLevels)
    ))
  } else {
    if (data_list$asymmetric == 0) {
      tmb_map <- c(tmb_map, list(
        sig2_b0 = as.factor(NA),
        log_sigma2 = as.factor(NA),
        sig2_b1 = as.factor(NA),
        sigma2_devs = rep(as.factor(NA), data_list$nLevels)
      ))
    } else {
      # estimate the random effects, just don't estimate the trend
      if (data_list$sig_trend == 0) {
        tmb_map <- c(tmb_map, list(sig2_b1 = as.factor(NA)))
      }
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

  obj <- TMB::MakeADFun(
    data = data_list,
    parameters = parameters,
    map = tmb_map,
    DLL = "salmix",
    random = random,
    silent = silent
  )

  # Attempt to do estimation
  if (!is.null(inits)) {
    init <- inits
  } else {
    init <- obj$par + runif(length(obj$par), -0.1, .1)
  }
  if (limits == FALSE) {
    pars <- stats::nlminb(
      start = init,
      objective = obj$fn,
      gradient = obj$gr,
      control = control
    )
  } else {
    pars <- stats::nlminb(
      start = init, objective = obj$fn,
      gradient = obj$gr, control = control,
      lower = limits(parnames = names(obj$par))$lower,
      upper = limits(parnames = names(obj$par))$upper,
      ...
    )
  }

  sdreport <- TMB::sdreport(obj)

  return(structure(list(
    obj = obj,
    pars = pars,
    sdreport = sdreport,
    init_values = parameters,
    data_list = data_list
  ), class = "salmix"))
}
