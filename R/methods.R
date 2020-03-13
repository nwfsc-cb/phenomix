#' Extract the number of observations of a salmix model. Modified from sdmTMB implementation
#'
#' @param object The fitted sdmTMB model object
#' @importFrom stats nobs
#' @export
#' @noRd
nobs.salmix <- function(object, ...) {
  length(fitted$data_list$y)
}


#' Extract the log likelihood of a salmix model. Modified from sdmTMB implementation
#'
#' @param object The fitted salmix model object
#' @importFrom stats logLik
#' @export
#' @noRd
logLik.salmix <- function(object, ...) {
  val <- -object$par$objective

  nobs <- nobs.salmix(object)
  df <- length(object$par) # fixed effects only
  structure(val,
    nobs = nobs, nall = nobs, df = df,
    class = "logLik"
  )
}

#' Extract the AIC of a salmix model. Modified from sdmTMB implementation
#'
#' @param fit The fitted salmix model
#' @param scale The scale (note used)
#' @param k Penalization parameter, defaults to 2
#' @param ... Anything else
#' @noRd
#'
#' @export
extractAIC.salmix <- function(fit, scale, k = 2, ...) {
  L <- logLik(fit)
  edf <- attr(L, "df")
  return(data.frame("df"=edf, "AIC"=-2 * L + k * edf))
}
