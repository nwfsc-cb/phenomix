#' Output processing function to be called by user
#'
#' This function extracts the annual totals
#'
#' @param fit A fitted object returned from fit()
#' @param log Whether to return estimates in log space, defaults to TRUE
#' @export
#'
extract_annual <- function(fit, log = TRUE) {
  all_name <- names(fit$sdreport$value)
  if(log==TRUE) {
    idx <- grep("year_log_tot", all_name)
  } else {
    idx <- grep("year_tot", all_name)
  }

  par_name <- "year_tot"
  if(log==TRUE) par_name <- "year_log_tot"

  df <- data.frame("value" = fit$sdreport$value[idx],
                   "sd" = fit$sdreport$sd[idx],
                   "par" = par_name)
  return(df)
}


#' Output processing function to be called by user
#'
#' This function extracts the parameter means and respective sds
#'
#' @param fit A fitted object returned from fit()
#' @export
#'
extract_means <- function(fit) {
  all_name <- names(fit$sdreport$value)
  idx <- which(all_name == "mu")
  df <- data.frame("value" = fit$sdreport$value[idx],
                   "sd" = fit$sdreport$sd[idx],
                   "par" = "mu")
  return(df)
}

#' Output processing function to be called by user
#'
#' This function extracts the parameter sigma and respective sds
#'
#' @param fit A fitted object returned from fit()
#' @export
#'
extract_sigma <- function(fit) {
  all_name <- names(fit$sdreport$value)
  idx <- which(all_name == "sigma1")
  df <- data.frame("value" = fit$sdreport$value[idx],
                   "sd" = fit$sdreport$sd[idx],
                   "par" = "sigma1")
  # attempt to add in 2nd sigma
  idx <- which(all_name == "sigma2")
  if(length(idx) > 0) {
    df2 <- data.frame("value" = fit$sdreport$value[idx],
                     "sd" = fit$sdreport$sd[idx],
                     "par" = "sigma2")
    df <- rbind(df,df2)
  }
  return(df)
}

#' Output processing function to be called by user
#'
#' This function extracts the parameter theta and respective sds
#'
#' @param fit A fitted object returned from fit()
#' @export
#'
extract_theta <- function(fit) {
  all_name <- names(fit$sdreport$value)
  idx <- which(all_name == "theta")
  df <- data.frame("value" = fit$sdreport$value[idx],
                   "sd" = fit$sdreport$sd[idx],
                   "par" = "theta")
  return(df)
}

#' Output processing function to be called by user
#'
#' This function extracts the lower quartiles (25%) and respective sds
#'
#' @param fit A fitted object returned from fit()
#' @export
#'
extract_lower <- function(fit) {
  all_name <- names(fit$sdreport$value)
  idx <- which(all_name == "lower25")
  df <- data.frame("value" = fit$sdreport$value[idx],
                   "sd" = fit$sdreport$sd[idx],
                   "par" = "lower25")
  return(df)
}

#' Output processing function to be called by user
#'
#' This function extracts the upper quartiles (25%) and respective sds
#'
#' @param fit A fitted object returned from fit()
#' @export
#'
extract_upper <- function(fit) {
  all_name <- names(fit$sdreport$value)
  idx <- which(all_name == "upper75")
  df <- data.frame("value" = fit$sdreport$value[idx],
                   "sd" = fit$sdreport$sd[idx],
                   "par" = "upper75")
  return(df)
}

#' Output processing function to be called by user
#'
#' This function extracts the means, sigmas, thetas,
#' lower (25%) and upper (75%) quartiles, and respective sds
#'
#' @param fit A fitted object returned from fit()
#' @export
#'
extract_all <- function(fit) {
  lower <- extract_lower(fit)
  upper <- extract_upper(fit)
  m <- extract_means(fit)
  sig <- extract_sigma(fit)
  theta <- extract_theta(fit)
  df <- rbind(m, sig, theta,
              lower,
              upper)
  return(df)
}




