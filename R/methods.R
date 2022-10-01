#' Extract the number of observations of a phenomix model. Modified from sdmTMB implementation
#'
#' @param object The fitted sdmTMB model object
#' @importFrom stats nobs
#' @export
#' @noRd
nobs.phenomix <- function(object, ...) {
  length(object$data_list$y)
}


#' Extract the log likelihood of a phenomix model. Modified from sdmTMB implementation
#'
#' @param object The fitted phenomix model object
#' @importFrom stats logLik
#' @export
#' @noRd
logLik.phenomix <- function(object, ...) {
  val <- -object$par$objective

  nobs <- nobs.phenomix(object)
  df <- length(object$par) # fixed effects only
  structure(val,
    nobs = nobs, nall = nobs, df = df,
    class = "logLik"
  )
}

#' Extract the AIC of a phenomix model. Modified from sdmTMB implementation
#'
#' @param fit The fitted phenomix model
#' @param scale The scale (note used)
#' @param k Penalization parameter, defaults to 2
#' @param ... Anything else
#' @noRd
#'
#' @export
extractAIC.phenomix <- function(fit, scale, k = 2, ...) {
  L <- logLik(fit)
  edf <- attr(L, "df")
  return(data.frame("df" = edf, "AIC" = -2 * L + k * edf))
}

#' Get predicted values from model object, copying glmmTMB 'fast' implementation
#'
#' @param fit The fitted phenomix model
#' @param se.fit Boolean, whether to produce SEs defaults to FALSE
#' @param ... Extra parameters
#'
#' @export
predict.phenomix <- function (object, se.fit=FALSE, ...) {
  ee <- environment(object$obj$fn)
  lp <- ee$last.par.best # best maximum likelihood estimate, similar to what glmmTMB uses

  new_tmb_obj <- TMB::MakeADFun(
    data = object$data_list,
    parameters = ee$parList(lp),
    map = ee$map,
    random = ee$random,
    DLL = "phenomix",
    silent = TRUE
  )
  new_sdreport <- TMB::sdreport(new_tmb_obj)

  indx <- which(names(new_sdreport$value)=="pred")
  pred <- new_sdreport$value[indx]

  dat <- data.frame(y = ee$data$y,
                    yint = ee$data$yint,
                    x = ee$data$x,
                    years = ee$data$years,
                    pred = pred)

  if(se.fit==TRUE) {
    dat$`se.fit` <- new_sdreport$sd[indx]
  }

  return(dat)
}
