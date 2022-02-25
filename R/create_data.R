#' Create data file for fitting time varying run timing distributions with TMB
#'
#' Does minimal processing of data to use as argument to fitting function
#'
#' @param data A data frame
#' @param min_number A minimum threshold to use, defaults to 0
#' @param variable A character string of the name of the variable in 'data' that contains the response (e.g. counts)
#' @param time A character string of the name of the variable in 'data' that contains the time variable (e.g. year)
#' @param date A character string of the name of the variable in 'data' that contains the response (e.g. day of year). The actual
#' #' column should contain a numeric response -- for example, the result from using lubridate::yday(x)
#' @param mu An optional formula allowing the mean to be a function of covariates. Random effects are not included in the formula
#' but specified with the `est_mu_re` argument
#' @param sigma An optional formula allowing the standard deviation to be a function of covariates. For asymmetric models,
#' each side of the distribution is allowed a different set of covariates. Random effects are not included in the formula
#' but specified with the `est_sigma_re` argument
#' @param covar_data a data frame containing covariates specific to each time step. These are used in the formulas `mu` and `sigma`
#' @param asymmetric_model Boolean, whether or not to let model be asymmetric (e.g. run timing before peak has a
#' different shape than run timing after peak)
#' @param est_sigma_re Whether to estimate random effects by year in sigma parameter controlling tail of distribution. Defaults to TRUE
#' @param est_mu_re Whether to estimate random effects by year in mu parameter controlling location of distribution. Defaults to TRUE
#' @param tail_model Whether to fit Gaussian ("gaussian" = default) or Student-t ("student_t") or generalized normal ("gnorm"). Defaults to "gaussian"
#' @param family Response for observation model, options are "gaussian", "poisson", "negbin"
#' @param max_theta Maximum value of log(pred) when `limits=TRUE`. Defaults to 10
#' @param share_shape Boolean argument for whether asymmetric student-t and generalized normal distributions should share the shape parameter (nu for the student-t;
#' beta for the generalized normal). Defaults to TRUE
#' @export
#' @importFrom stats model.matrix
#' @examples
#' data(fishdist)
#' datalist <- create_data(fishdist,
#'   min_number = 0, variable = "number", time = "year",
#'   date = "doy", asymmetric_model = TRUE, family = "gaussian"
#' )
create_data <- function(data,
                        min_number = 0,
                        variable = "number",
                        time = "year",
                        date = "doy",
                        asymmetric_model = TRUE,
                        mu = ~1,
                        sigma = ~1,
                        covar_data = NULL,
                        est_sigma_re = TRUE,
                        est_mu_re = TRUE,
                        tail_model = "gaussian",
                        family = "gaussian",
                        max_theta = 10,
                        share_shape = TRUE) {
  dist <- c("gaussian", "poisson", "negbin")
  fam <- match(family, dist)
  if (is.na(fam)) {
    stop("Make sure the entered family is in the list of accepted distributions")
  }

  tail <- c("gaussian", "student_t", "gnorm")
  tailmod <- match(tail_model, tail)
  if (is.na(tailmod)) {
    stop("Make sure the entered tail model is in the list of accepted distributions")
  }

  # check to make sure year and date are numeric
  if (!is.numeric(data[, time])) {
    stop("The time variable in the data frame (e.g. year) needs to be numeric")
  }
  if (is.numeric(data[, date])) {
    if (max(data[, date], na.rm = T) > 365) stop("The date variable in the data frame contains values greater than 365")
    if (min(data[, date], na.rm = T) < 1) stop("The date variable in the data frame contains values less than 1")
  } else {
    stop("The date variable in the data frame (e.g. day_of_year) needs to be numeric")
  }

  # if 1 level, turn off trend and random effect estimation
  if (length(unique(as.numeric(data[, time]))) == 1) {
    est_sigma_re <- FALSE
    est_mu_re <- FALSE
  }

  # drop rows below threshold or NAs
  drop_rows <- which(is.na(data[, variable]) | data[, variable] < min_number)
  if (length(drop_rows) > 0) data <- data[-drop_rows, ]

  # rescale year variable to start at 1 for indexing
  data$year <- data[, time] - min(data[, time]) + 1

  # parse formulas. covar_data contains covariates specific to each time step
  if (is.null(covar_data)) {
    covar_data <- data.frame(year = unique(data$year))
  }
  mu_mat <- model.matrix(mu, data = covar_data)
  sig_mat <- model.matrix(sigma, data = covar_data)

  data_list <- list(
    y = data[, variable],
    yint = round(data[, variable]),
    years = as.numeric(as.factor(data$year)),
    x = data[, date],
    year_levels = as.numeric(as.factor(unique(data$year))),
    unique_years = unique(data$year),
    nLevels = length(unique(data$year)),
    asymmetric = as.numeric(asymmetric_model),
    family = fam,
    mu_mat = mu_mat,
    sig_mat = sig_mat,
    tail_model = as.numeric(tailmod) - 1,
    est_sigma_re = as.numeric(est_sigma_re),
    est_mu_re = as.numeric(est_mu_re),
    max_theta = max_theta,
    share_shape = as.numeric(share_shape)
  )

  return(data_list)
}
