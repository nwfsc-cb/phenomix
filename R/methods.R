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

#' @importFrom stats predict
#' @export
fitted.phenomix <- function(object, ...) {
  predict(object)
}

#' Get predicted values from model object, copying glmmTMB 'fast' implementation
#'
#' @param object The fitted phenomix model
#' @param se.fit Boolean, whether to produce SEs defaults to FALSE
#' @param ... Extra parameters
#' @importFrom RTMB sdreport
#' @noRd
#'
#' @export
predict.phenomix <- function(object, se.fit = FALSE, ...) {
  new_sdreport <- get_sdreport(object)

  indx <- which(names(new_sdreport$value) == "pred")
  pred <- new_sdreport$value[indx]

  # For RTMB compatibility, get data from object$data_list instead of environment
  dat <- data.frame(
    y = object$data_list$y,
    yint = object$data_list$yint,
    x = object$data_list$x,
    years = object$data_list$years,
    pred = pred,
    stringsAsFactors = FALSE
  )

  if (se.fit == TRUE) {
    dat$`se.fit` <- new_sdreport$sd[indx]
  }

  return(dat)
}

#' Get parameters from model object, copying glmmTMB 'fast' implementation
#'
#' @param object The fitted phenomix model
#' @importFrom RTMB sdreport
#' @export
#' @examples
#' \donttest{
#' data(fishdist)
#'
#' # example of fitting fixed effects, no trends, no random effects
#' set.seed(1)
#' datalist <- create_data(fishdist[which(fishdist$year > 1970), ],
#'   asymmetric_model = FALSE,
#'   est_mu_re = FALSE, est_sigma_re = FALSE
#' )
#' fit <- fit(datalist)
#' p <- pars(fit)
#' names(p)
#' }
#' @export
#' @keywords internal
pars <- function(object) {
  new_sdreport <- get_sdreport(object)
  # new_sdreport[["pdHess"]] = NULL
  new_sdreport[["env"]] <- NULL
  # new_sdreport[["gradient.fixed"]] = NULL
  return(new_sdreport)
}


#' Get fixed effects parameters from model object, copying glmmTMB 'fast' implementation
#'
#' @param object The fitted phenomix model
#' @param ... Additional arguments
#' @importFrom RTMB sdreport
#' @importFrom nlme fixef
#' @export
fixef.phenomix <- function(object, ...) {
  new_sdreport <- get_sdreport(object)
  p <- data.frame(
    "par" = names(new_sdreport$par.fixed),
    "estimate" = new_sdreport$par.fixed,
    "se" = sqrt(diag(new_sdreport$cov.fixed))
  )
  return(p)
}


#' Get random effects parameters from model object, copying glmmTMB 'fast' implementation
#'
#' @param object The fitted phenomix model
#' @param ... Additional arguments
#' @importFrom RTMB sdreport
#' @importFrom nlme ranef
#' @export
ranef.phenomix <- function(object, ...) {
  new_sdreport <- get_sdreport(object)
  p <- data.frame(
    "par" = names(new_sdreport$par.random),
    "estimate" = new_sdreport$par.random,
    "se" = sqrt(new_sdreport$diag.cov.random)
  )
  return(p)
}


#' Get sdreport for predictions / coefficients, copying glmmTMB 'fast' implementation
#'
#' @param object The fitted phenomix model
#' @keywords internal
get_sdreport <- function(object) {
  # For RTMB, the sdreport is already computed and stored in the object
  # If it exists, return it; otherwise recompute it
  if (!is.null(object$sdreport)) {
    return(object$sdreport)
  } else {
    # Recompute if needed
    new_sdreport <- RTMB::sdreport(object$obj)
    return(new_sdreport)
  }
}
