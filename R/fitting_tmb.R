#' Fit a trinomial mixture model with TMB
#'
#' Fit a trinomial mixture model that optionally includes covariates to estimate
#' effects of factor or continuous variables on proportions.
#'
#' @param formula The model formula for the design matrix. Does not need to have a response specified. If =NULL, then
#' the design matrix is ignored and all rows are treated as replicates
#' @param design_matrix A data frame, dimensioned as number of observations, and covariates in columns
#' @param data_matrix A matrix, with observations on rows and number of groups across columns
#' @param overdispersion Whether or not to include overdispersion parameter, defaults to FALSE
#' @param overdispersion_sd Prior standard deviation on 1/overdispersion parameter, Defaults to inv-Cauchy(0,5)
#' @param prior_sd Optional prior sd / penalty for fixed effects
#' @export
#' @import Rcpp
#' @importFrom stats model.frame
#'
#' @examples
#' \donttest{
#  #y <- matrix(c(3.77, 6.63, 2.60, 0.9, 1.44, 0.66, 2.10, 3.57, 1.33),
#  #            nrow = 3, byrow = TRUE
#  #)
#  #fit a model with no covariates
#  #fit <- fit_zoidTMB(data_matrix = y)
#'
#' # fit a model with 1 factor
#' #design <- data.frame("fac" = c("spring", "spring", "fall"))
#' #fit <- fit_zoidTMB(formula = ~fac, design_matrix = design, data_matrix = y)
#'
#' }
#'
fit_zoidTMB <- function(formula = NULL,
                        design_matrix,
                        data_matrix,
                        overdispersion = FALSE,
                        overdispersion_sd = 5,
                        prior_sd = NA) {

  # if a single observation
  if (class(data_matrix)[1] != "matrix") {
    data_matrix <- matrix(data_matrix, nrow = 1)
  }

  # fill with dummy values
  parsed_res <- list(design_matrix = matrix(0, nrow(data_matrix),ncol=1),
                     var_indx = 1,
                     n_re_by_group = 1,
                     tot_re = 1,
                     n_groups = 1)
  est_re <- FALSE
  re_group_names <- NA
  if (!is.null(formula)) {
    model_frame <- model.frame(formula, design_matrix)
    model_matrix <- model.matrix(formula, model_frame)
    # extract the random effects
    res <- parse_re_formula(formula, design_matrix)
    if(length(res$var_indx) > 0) {
      parsed_res <- res # only update if REs are in formula
      est_re <- TRUE
      model_matrix <- res$fixed_design_matrix
      re_group_names <- res$random_effect_group_names
    }
  } else {
    model_matrix <- matrix(1, nrow = nrow(data_matrix))
    colnames(model_matrix) <- "(Intercept)"
  }

  sd_prior <- 1 / ncol(data_matrix) # default if no covariates
  if (ncol(model_matrix) > 1) sd_prior <- 1
  use_prior_sd = 0L
  if (!is.na(prior_sd)) {
    sd_prior <- prior_sd
    use_prior_sd = 1L
  } else {
    sd_prior <- 0
  }

  par_names <- colnames(model_matrix)

  prod_idx <- matrix(0, ncol(data_matrix), ncol(data_matrix)-1)
  for(j in 1:ncol(data_matrix)){
    prod_idx[j,] <- seq(1,ncol(data_matrix),1)[-j]
  }

  N_bins = ncol(data_matrix)
  N_covar = ncol(model_matrix)
  n_groups = parsed_res$n_group
  tot_re = parsed_res$tot_re

  tmb_data <- list(
    X = data_matrix,
    prod_idx = prod_idx - 1L, # -1 for indexing to 0
    design_X = model_matrix,
    overdisp = ifelse(overdispersion == TRUE, 1, 0),
    overdispersion_sd = overdispersion_sd,
    use_prior_sd = use_prior_sd,
    prior_sd = sd_prior,
    #design_Z = parsed_res$design_matrix, # design matrix for Z (random int)
    re_var_indx = c(parsed_res$var_indx, 1) - 1L, # index of the group for each re
    #n_re_by_group = c(parsed_res$n_re_by_group, 1), # number of random ints per group
    #tot_re = tot_re, # total number of random ints, across all groups
    n_groups = n_groups,
    est_re = as.numeric(est_re)
  )

  # estimated parameters from TMB.
  # beta is included in all models
  # zeta_sds and zeta_vec are only in random effects models
  # phi_inv is only in overdispersion model
  tmb_pars <- list(beta_raw = matrix(0, N_bins-1, N_covar),
                   log_phi_inv = 0)
                   #log_zeta_sds = rep(0, n_groups),
                   #zeta_vec = rep(0, (N_bins - 1) * tot_re))

  tmb_map <- c()
  if(overdispersion == FALSE) {
    tmb_map <- c(tmb_map, list(log_phi_inv = as.factor(NA)))
  }

  #random <- "zeta_vec"
  #if(est_re == FALSE) {
  #  random <- NULL
    #tmb_map <- c(tmb_map,
    #             list(z = rep(as.factor(NA), (N_bins - 1) * tot_re)))
    #tmb_map <- c(tmb_map, list(log_zeta_sds = as.factor(rep(NA, n_groups)),#log_zeta_sds = as.factor(rep(NA, n_groups)),
    #                           zetavec = as.factor(rep(NA, (N_bins - 1) * tot_re))))
  #}

  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_pars,
    map = tmb_map,
    DLL = "zoidtmb",
    random = NULL,
    silent = TRUE
  )

  pars <- stats::nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr
  )
  sdreport <- TMB::sdreport(obj)

  mod_list <- list(pars = pars, sdreport = sdreport,
                   tmb_data = tmb_data, tmb_map = tmb_map,
                   #tmb_random = random,
                   tmb_parameters = tmb_pars)

  # tidy the mu and beta vectors
  idx <- grep("mu", names(sdreport$value))
  mod_list$mu <- matrix(sdreport$value[idx], nrow=nrow(data_matrix))
  mod_list$mu_se <- matrix(sdreport$sd[idx], nrow=nrow(data_matrix))

  return(mod_list)
}

#' Fit a trinomial mixture model that optionally includes covariates to estimate
#' effects of factor or continuous variables on proportions.
#'
#' @param formula The model formula for the design matrix.
#' @param data The data matrix used to construct RE design matrix
#' @importFrom stats model.matrix as.formula
parse_re_formula <- function(formula, data) {
  # Convert the formula to a character string
  formula_str <- as.character(formula)
  # Split the formula into parts based on '+' and '-' symbols
  formula_parts <- unlist(strsplit(formula_str, split = "[-+]", perl = TRUE))
  # Trim whitespace from each part
  formula_parts <- trimws(formula_parts)
  # Identify parts containing a bar '|'
  random_effects <- grep("\\|", formula_parts, value = TRUE)
  fixed_effects <- setdiff(formula_parts, random_effects)

  # Create design matrix for fixed effects. Catch the cases where no fixed
  # effects are included, or intercept-only models used
  if (length(fixed_effects) > 1 || (length(fixed_effects) == 1 && fixed_effects != "~")) {
    fixed_formula_str <- paste("~", paste(fixed_effects, collapse = "+"))
  } else {
    fixed_formula_str <- "~ 1" # Only intercept
  }
  fixed_design_matrix <- model.matrix(as.formula(fixed_formula_str), data)

  random_effect_group_names <- sapply(random_effects, function(part) {
    # Extract the part after the '|'
    split_part <- strsplit(part, "\\|", perl = TRUE)[[1]]
    # Remove the closing parenthesis and trim
    group_name <- gsub("\\)", "", split_part[2])
    trimws(group_name)
  })

  # create design matrices by group
  for(i in 1:length(random_effects)) {
    new_formula <- as.formula(paste("~", random_effect_group_names[i], "-1"))
    if(i ==1) {
      design_matrix <- model.matrix(new_formula, data)
      var_indx <- rep(1, ncol(design_matrix))
      n_re <- length(var_indx)
    } else {
      design_matrix <- cbind(design_matrix, model.matrix(new_formula, data))
      var_indx <- c(var_indx, rep(i, ncol(design_matrix)))
      n_re <- c(n_re, length(ncol(design_matrix)))
    }
  }
  n_groups <- 0
  if(length(var_indx) > 0) n_groups <- max(var_indx)
  return(list(design_matrix = design_matrix, var_indx = var_indx, n_re_by_group = n_re,
              tot_re = sum(n_re), n_groups = n_groups,
              fixed_design_matrix = fixed_design_matrix,
              random_effect_group_names = random_effect_group_names))
}

