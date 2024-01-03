#include <TMB.hpp>
// #include <omp.h>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;

  DATA_INTEGER(overdisp); // whether or not to include overdispersion term
  DATA_INTEGER(use_prior_sd); // whether to use penalty on fixed effects betas
  DATA_INTEGER(n_groups); // number of re groups
  DATA_INTEGER(est_re); // 0 or 1 indicator

  DATA_MATRIX(X); //[N_samples, N_bins] X; // proportions
  DATA_MATRIX(design_X); //[N_samples, N_covar] design_X;
  DATA_IMATRIX(prod_idx);//[N_bins,N_bins-1] int prod_idx;
  DATA_SCALAR(prior_sd);

  //DATA_MATRIX(design_Z);//[N_samples, tot_re] design_Z; // design matrix of random ints
  DATA_IVECTOR(re_var_indx); //[tot_re + 1]; // total index across acll groups

  // Parameters
  PARAMETER(log_phi_inv); // log-transform to ensure positivity
  PARAMETER_MATRIX(beta_raw);
  //PARAMETER_VECTOR(zeta_vec);
  //PARAMETER_VECTOR(log_zeta_sds); // log-transform for positivity

  // Transform parameters to original scale
  Type phi_inv = exp(log_phi_inv); // Exponentiating to ensure positivity
  //vector<Type> zeta_sds = exp(log_zeta_sds); // Same for standard deviations

  int N_samples = X.rows();
  int N_bins = X.cols();
  int N_covar = design_X.cols();
  //int tot_re = design_Z.cols();

  // Transformed data
  matrix<Type> is_zero(N_samples, N_bins);
  matrix<Type> is_proportion(N_samples, N_bins);
  matrix<Type> logX(N_samples, N_bins);
  matrix<Type> logNX(N_samples, N_bins);
  vector<Type> ESS(N_samples);
  vector<Type> ones(N_bins);

  // The negative log-likelihood
  Type nll = 0.0;

  // Calculate ones for Dirichlet
  ones.setZero();
  ones += ones;

  // Calculate ESS, summed across cols for each row
  ESS.setZero();
  for (int i = 0; i < N_samples; i++) {
    for (int j = 0; j < N_bins; j++) {
      ESS(i) += X(i, j);
    }
  }

  // Calculate is_zero and is_proportion indicators
  for (int i = 0; i < N_samples; i++) {
    for (int j = 0; j < N_bins; j++) {
      is_zero(i, j) = (X(i, j) == 0) ? Type(1) : Type(0);
      is_proportion(i, j) = (X(i, j) < ESS(i) && X(i, j) > 0) ? Type(1) : Type(0);
      if (is_proportion(i, j) == Type(1)) {
        logX(i, j) = log(X(i, j));
        logNX(i, j) = log(ESS(i) - X(i, j));
      }
    }
  }

  Type phi = 1.0; // Initialize phi
  if (overdisp == 1) {
    phi = Type(1) / phi_inv;
  }

  matrix<Type> p_zero(N_samples, N_bins); // Probability of 0 for each cell
  matrix<Type> p_one(N_samples, N_bins); // Probability of 1 for each cell
  matrix<Type> beta(N_bins, N_covar); // Coefficients
  //matrix<Type> zeta(N_bins, tot_re); // Coefficients for random effects
  matrix<Type> mu(N_samples, N_bins); // Estimates, in normal space

  // Fill beta and zeta matrices
  beta.row(N_bins - 1).setZero(); // For identifiability
  for (int k = 0; k < (N_bins - 1); k++) {
    beta.row(k) = beta_raw.row(k);
  }

  //if (est_re == 1) {
  //  int counter = 0;
    //zeta.row(N_bins - 1).setZero(); // For identifiability
  //  for (int k = 0; k < (N_bins - 1); k++) {
      //for (int j = 0; j < tot_re; j++) {
        //zeta(k,j) = zeta_vec(counter);
   //     counter += counter;
      //}
   // }
  //}

  // Calculate mu
  vector<Type> logits(N_bins);
  for (int n = 0; n < N_samples; n++) {
    logits.setZero();
    for (int m = 0; m < N_bins; m++) {
      logits(m) = (design_X.row(n) * beta.col(m)).sum();
      //if (est_re == 1) {
        //logits(m) += (design_Z.row(n) * zeta.col(m)).sum();
      //}
    }
    logits = logits.exp();
    mu.row(n) = logits / logits.sum(); // Softmax
  }

  // Calculate p_zero and p_one
  for (int i = 0; i < N_samples; i++) {
    for (int j = 0; j < N_bins; j++) {
      p_zero(i, j) = pow(1.0 - mu(i, j), ESS[i] * phi);

      // Calculation of p_one requires the product over a subset of p_zero
      Type prod_p_zero = 1.0;
      for (int idx = 0; idx < (N_bins - 1); idx++) {
        if (prod_idx(j, idx) != j) {
          prod_p_zero *= p_zero(i, prod_idx(j, idx));
        }
      }
      p_one(i, j) = (1.0 - p_zero(i, j)) * prod_p_zero;
    }
  }

  // Overdispersion prior
  //if (overdisp == 1) {
  //  nll -= dcauchy(phi_inv[1], Type(0), Type(5), true);
  //}

  // Priors for fixed effects for covariate factors
  if(use_prior_sd == 1) {
    for (int i = 0; i < N_covar; i++) {
      for (int j = 0; j < (N_bins - 1); j++) {
        nll -= dnorm(beta_raw(j, i), Type(0), prior_sd, true);
      }
    }
  }

  // Priors for random effects
  //if (est_re == 1) {
    // for (int i = 0; i < n_groups; i++) {
    //   nll -= dstudent_t(zeta_sds[i], Type(3), Type(0), Type(2), true);
    // }
    //for (int i = 0; i < (N_bins - 1); i++) {
      //for (int j = 0; j < tot_re; j++) {
        //nll -= dnorm(zeta(i, j), Type(0), zeta_sds(re_var_indx(j)), true);
      //}
    //}
  //}

  // Likelihood contributions for each sample and bin
  for (int i = 0; i < N_samples; i++) {
    for (int j = 0; j < N_bins; j++) {
      // Marginals of the trinomial are independent binomials
      nll -= dbinom(is_zero(i, j), Type(1), p_zero(i, j), true);

      if (is_proportion(i, j) == Type(1)) {
        Type alpha_temp = mu(i, j) * ESS(i) * phi;
        Type beta_temp = (Type(1) - mu(i, j)) * ESS(i) * phi;
        // Beta log probability mass function for 3-parameter model
        nll -= dbeta(X(i, j) / ESS(i), alpha_temp, beta_temp, true);
      }
    }
  }

  // report w derivatives
  ADREPORT(beta);
  ADREPORT(phi_inv);
  ADREPORT(log_phi_inv);
  //ADREPORT(zeta_sds);
  //ADREPORT(zeta);
  ADREPORT(mu);
  ADREPORT(p_zero);
  ADREPORT(p_one);

  REPORT(beta);
  REPORT(phi_inv);
  REPORT(log_phi_inv);
  //REPORT(zeta_sds);
  //REPORT(zeta);
  REPORT(mu);
  REPORT(p_zero);
  REPORT(p_one);
  return nll;
}
