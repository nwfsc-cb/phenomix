#' @useDynLib phenomix, .registration = TRUE
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
#' @param limits Whether to include limits for stats::nlminb. Can be a list of (lower, upper), or TRUE to use suggested hardcoded limits. Defaults to NULL,
#' where no limits used.
#' @param fit_model Whether to fit the model. If not, returns a list including the data, parameters, and initial values. Defaults to TRUE
#' @importFrom stats runif rnorm
#' @export
#' @examples
#' data(fishdist)
#'
#' # example of fitting fixed effects, no trends, no random effects
#' set.seed(1)
#' datalist <- create_data(fishdist[which(fishdist$year > 1970), ],
#'   asymmetric_model = FALSE,
#'   est_mu_re = FALSE, est_sigma_re = FALSE
#' )
#' fit <- fit(datalist)
#' #
#' # # example of model with random effects in means only, and symmetric distribution
#' # set.seed(1)
#' # datalist <- create_data(fishdist[which(fishdist$year > 1970),], asymmetric_model = FALSE,
#' #                      est_sigma_re = FALSE)
#' # fit <- fit(datalist)
#' # # example of model with random effects in variances
#' # set.seed(1)
#' # datalist <- create_data(fishdist[which(fishdist$year > 1970),], asymmetric_model = TRUE,
#' #                          est_mu_re = TRUE)
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
                limits = NULL, # can also be a list, or TRUE
                fit_model = TRUE) {

  # create list of parameter starting values -- used in both
  # asymmetric and symmetric models
  parameters <- list(
    theta = rnorm(n = data_list$nLevels, log(mean(data_list$y[which(!is.na(data_list$y))]))),
    b_mu = rep(0, ncol(data_list$mu_mat)),#runif(ncol(data_list$mu_mat),0,10), #rep(0, ncol(data_list$mu_mat)),
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
  #CV <- 0.1
  #parameters$b_sig1[1] <- 1# CV * parameters$b_mu[1]
  #parameters$b_sig2[1] <- parameters$b_sig1[1]
  if (data_list$family == 1) {
    # for gaussian, don't log-transform theta
    parameters$theta <- rep(0, data_list$nLevels) # rnorm(n = data_list$nLevels, mean(data_list$y))
  }

  # If inits is included, use that instead of parameters
  # if (!is.null(inits)) parameters <- inits

  # Mapping off params as needed:
  tmb_map <- list()
  if (data_list$family %in% c(2, 4)) {
    # don't include obs_sigma for poisson or binomial
    tmb_map <- c(tmb_map, list(log_obs_sigma = as.factor(NA)))
  }

  if (data_list$asymmetric == 0) {
    # map off pars not needed
    tmb_map <- c(tmb_map, list(
      b_sig2 = rep(as.factor(NA), ncol(data_list$sig_mat)),
      log_sigma2_sd = as.factor(NA),
      sigma2_devs = rep(as.factor(NA), data_list$nLevels)
    ))
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
        log_tdf_2 = as.factor(NA),
        log_beta_1 = as.factor(NA),
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
        log_beta_2 = as.factor(NA),
        log_tdf_1 = as.factor(NA),
        log_tdf_2 = as.factor(NA)
      ))
    } else {
      tmb_map <- c(tmb_map, list(
        log_tdf_1 = as.factor(NA),
        log_tdf_2 = as.factor(NA)
      ))
    }
  }
  if (data_list$asymmetric == 1) {
    if (data_list$share_shape == 1) {
      if (data_list$tail_model == 1) {
        # map off 2nd nu parameter
        tmb_map <- c(tmb_map, list(
          log_tdf_2 = as.factor(NA)
        ))
      }
      if (data_list$tail_model == 2) {
        tmb_map <- c(tmb_map, list(
          log_beta_2 = as.factor(NA)
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
    tmb_map <- c(tmb_map, list(
      log_sigma_mu_devs = as.factor(NA),
      mu_devs = rep(as.factor(NA), data_list$nLevels)
    ))
  }

  if (data_list$est_sigma_re == 0) {
    tmb_map <- c(tmb_map, list(
      log_sigma1_sd = as.factor(NA),
      log_sigma2_sd = as.factor(NA),
      sigma1_devs = rep(as.factor(NA), data_list$nLevels),
      sigma2_devs = rep(as.factor(NA), data_list$nLevels)
    ))
  } else {
    if (data_list$asymmetric == 0) {
      tmb_map <- c(tmb_map, list(
        log_sigma2_sd = as.factor(NA),
        sigma2_devs = rep(as.factor(NA), data_list$nLevels)
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

  obj <- TMB::MakeADFun(
    data = data_list,
    parameters = parameters,
    map = tmb_map,
    DLL = "phenomix",
    random = random,
    silent = silent
  )

  mod_list <- list(
    obj = obj,
    init_vals = obj$par,
    # pars = pars,
    # sdreport = sdreport,
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
      if (class(limits) == "list") {
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

    sdreport <- TMB::sdreport(obj)
    mod_list <- c(mod_list, list(pars = pars, sdreport = sdreport))
  }

  return(structure(mod_list, class = "phenomix"))
}
