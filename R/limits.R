#' Internal function to assign upper and lower bounds to parameters, based on their names
#'
#' Arguments can be adjusted as needed
#'
#' @param parnames A vector of character strings or names used for estimation
#' @param max_theta A scalar (optional) giving the maximum value of theta, log(pred)
#'
limits <- function(parnames, max_theta) {
  df <- data.frame(name = parnames, lower = -1000, upper = 1000)

  for (i in 1:length(parnames)) {
    if (length(grep("log_sigma1", parnames[i])) > 0) df[i, c("lower", "upper")] <- c(-5, 6)
    if (length(grep("log_sigma2", parnames[i])) > 0) df[i, c("lower", "upper")] <- c(-5, 6)
    if (length(grep("theta", parnames[i])) > 0) df[i, c("lower", "upper")] <- c(0, max_theta)
    if (length(grep("log_sigma_mu_devs", parnames[i])) > 0) df[i, c("lower", "upper")] <- c(-5, 6)
    if (length(grep("log_obs_sigma", parnames[i])) > 0) df[i, c("lower", "upper")] <- c(-5, 5)
    if (length(grep("log_tdf", parnames[i])) > 0) df[i, c("lower", "upper")] <- c(-5, 5)
    if (length(grep("log_beta", parnames[i])) > 0) df[i, c("lower", "upper")] <- c(-5, 3)
  }
  return(df)
}
