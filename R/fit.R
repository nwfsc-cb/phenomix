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
#' @importFrom stats runif
#' @export
#' @examples
#' data(fishdist)
#' datalist = create_data(fishdist)
fit <- function(data_list, silent=FALSE, inits = NULL, control=list(eval.max=2000, iter.max=1000,
  rel.tol=1e-10), limits=FALSE) {

  # create list of parameter starting values -- used in both
  # asymmetric and symmetric models
  parameters <- list(sigma1_devs = rep(0, data_list$nLevels),
    theta = runif(n=data_list$nLevels,5,7),
    mu_devs = rep(0, data_list$nLevels),
    log_mu_b0 = log(120),
    mu_b1 = 0,
    log_sigma_mu_devs = -1,
    sig1_b0 = 2.0,
    sig1_b1 = 0.0,
    log_sigma1 = -1,
    log_obs_sigma=0.005)

  # optional parameters to add for asymmetric model
  #if(data_list$asymmetric==1) {
    parameters <- append(parameters,
      list(sigma2_devs = rep(0, data_list$nLevels),
        sig2_b0 = 0.0,
        sig2_b1 = 0.0,
        log_sigma2 = 0.0,
        log_tdf_1 = 0,
        log_tdf_2 = 0))
  #}

  # If inits is included, use that instead of parameters
  if(!is.null(inits)) parameters = inits

  # Mapping off params as needed:
  tmb_map <- list()
  if (data_list$family == 2)
    tmb_map <- c(tmb_map, list(log_obs_sigma = as.factor(NA)))

  if(data_list$asymmetric==0) {
    # map off pars not needed
    tmb_map <- c(tmb_map, list(sig2_b0 = as.factor(NA),
      sig2_b1 = as.factor(NA),
      log_sigma2 = as.factor(NA),
      sigma2_devs = rep(as.factor(NA), data_list$nLevels)))
  }

  if(data_list$sig_trend == FALSE) {
    tmb_map = c(tmb_map, list(sig1_b1 = as.factor(NA),
      sig2_b1 = as.factor(NA)))
  }
  if(data_list$mu_trend == FALSE) {
    tmb_map = c(tmb_map, list(mu_b1 = as.factor(NA)))
  }

  if(data_list$t_model == 0) {
    # then fit gaussian model and map off both parameters
    tmb_map <- c(tmb_map, list(log_tdf_1 = as.factor(NA),
      log_tdf_2 = as.factor(NA)))
  } else {
    # fit the t-model
    if(data_list$asymmetric==0) {
      # then map off the tdf2, because model is symmetric
      tmb_map <- c(tmb_map, list(log_tdf_2 = as.factor(NA)))
    }
  }

  random = c("mu_devs","sigma1_devs")
  if(data_list$asymmetric==TRUE) random = c(random, "sigma2_devs")

  obj <- TMB::MakeADFun(data = data_list,
    parameters=parameters,
    map = tmb_map, DLL="salmix",
    random=random, silent=silent)

  # Attempt to do estimation
  if(!is.null(inits)) {
    init = inits
  } else {
    init = obj$par + runif(length(obj$par),-0.1,.1)
  }
  if(limits==FALSE) {
    pars = stats::nlminb(
      start = init, objective = obj$fn,
      gradient = obj$gr, control = control)
  } else {
    pars = stats::nlminb(
      start = init, objective = obj$fn,
      gradient = obj$gr, control = control,
      lower = limits(parnames = names(obj$par))$lower,
      upper = limits(parnames = names(obj$par))$upper)
  }

  sdreport = TMB::sdreport(obj)

  return(structure(list(obj = obj,
    pars=pars,
    sdreport=sdreport,
    init_values = parameters,
    data_list = data_list), class="salmix"))
}
