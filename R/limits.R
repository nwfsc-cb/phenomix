#' Internal function to assign upper and lower bounds to parameters, based on their names
#'
#' Arguments can be adjusted as needed
#'
#' @param parnames A vector of character strings or names used for estimation
#'
limits <- function(parnames) {
  df = data.frame(name=parnames,lower=-10,upper=10)

  for(i in 1:length(parnames)) {
    if(length(grep("sigma1",parnames[i]))>0) df[i,c("lower","upper")] = c(-15,5)
    if(length(grep("sigma2",parnames[i]))>0) df[i,c("lower","upper")] = c(-15,5)
    if(length(grep("theta",parnames[i]))>0) df[i,c("lower","upper")] = c(2,25)
    if(length(grep("log_mu_b0",parnames[i]))>0) df[i,c("lower","upper")] = c(1,log(360))
    if(length(grep("mu_b1",parnames[i]))>0) df[i,c("lower","upper")] = c(-1,1)
    if(length(grep("log_sigma_mu_devs",parnames[i]))>0) df[i,c("lower","upper")] = c(-10,5)
    if(length(grep("sig1_b0",parnames[i]))>0) df[i,c("lower","upper")] = c(-1,5)
    if(length(grep("sig2_b0",parnames[i]))>0) df[i,c("lower","upper")] = c(-1,5)
    if(length(grep("sig1_b1",parnames[i]))>0) df[i,c("lower","upper")] = c(-0.2,0.2)
    if(length(grep("sig2_b1",parnames[i]))>0) df[i,c("lower","upper")] = c(-0.2,0.2)
    if(length(grep("log_obs_sigma",parnames[i]))>0) df[i,c("lower","upper")] = c(-5,3)
    if(length(grep("log_tdf",parnames[i]))>0) df[i,c("lower","upper")] = c(-10,5)
    if(length(grep("log_beta",parnames[i]))>0) df[i,c("lower","upper")] = c(-10,3)
  }
  return(df)
}
