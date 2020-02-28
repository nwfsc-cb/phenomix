#include <TMB.hpp>
// #include <omp.h>

template<class Type>
  Type objective_function<Type>::operator() ()
{
  using namespace density;

  DATA_VECTOR(y);
  DATA_VECTOR(x);
  DATA_IVECTOR(years);
  DATA_IVECTOR(year_levels);
  DATA_IVECTOR(unique_years);
  DATA_INTEGER(nLevels);
  PARAMETER_VECTOR(sigma1_devs);
  PARAMETER_VECTOR(theta);
  PARAMETER_VECTOR(mu_devs);
  PARAMETER(log_mu_b0);
  PARAMETER(mu_b1);
  PARAMETER(log_sigma_mu_devs);
  PARAMETER(sig1_b0);
  PARAMETER(sig1_b1);
  PARAMETER(log_sigma1);
  PARAMETER(log_obs_sigma);
  DATA_INTEGER(asymmetric); // 0 if false, 1 = true
  DATA_INTEGER(family); // 1 gaussian, 2 = poisson, 3 = neg bin
  PARAMETER_VECTOR(sigma2_devs);
  PARAMETER(sig2_b0);
  PARAMETER(sig2_b1);
  PARAMETER(log_sigma2);

  // derived parameters
  int n = x.size();
  Type obs_sigma=exp(log_obs_sigma);
  vector<Type> sigma1(nLevels), mu(nLevels);
  vector<Type> sigma2(nLevels),scalar(nLevels);
  vector<Type> lower25(nLevels), upper75(nLevels);
  int i = 0;

  Type nll=0;

  for(i = 0; i < nLevels; i++) {

    sigma1(i) = exp(sig1_b0 + sigma1_devs(i) + sig1_b1*Type(unique_years(i)));
    if(asymmetric == 1) {
      sigma2(i) = exp(sig2_b0 + sigma2_devs(i) + sig2_b1*Type(unique_years(i)));
      // scalar(i) is just log(sig2) - log(sig1)
      scalar(i) = sigma2_devs(i) + sigma2_devs(i) - (sigma1_devs(i) + sigma1_devs(i));
    }

    mu(i) = exp(log_mu_b0) + mu_devs(i) + mu_b1*Type(unique_years(i));

    // random effects contributions
    nll += dnorm(mu_devs(i), Type(0.0),exp(log_sigma_mu_devs),true);
    nll += dnorm(sigma1_devs(i),Type(0.0),exp(log_sigma1),true);
    if(asymmetric == 1) {
      nll += dnorm(sigma2_devs(i),Type(0.0),exp(log_sigma2),true);
    }

    lower25(i) = qnorm(Type(0.25), mu(i), sigma1(i));
    if(asymmetric == 1) {
      upper75(i) = qnorm(Type(0.75), mu(i), sigma2(i));
    } else {
      upper75(i) = qnorm(Type(0.75), mu(i), sigma1(i));
    }
  }

  vector<Type> log_dens(n), pred(n);
  for(i = 0; i < n; i++) {
    if(asymmetric == 1) {
      // model is asymmetric, left side smaller / right side bigger
      if(x(i) < mu(years(i)-1)) {
        log_dens(i) = dnorm(x(i), mu(years(i)-1), sigma1(years(i)-1), true);
        pred(i) = log_dens(i) + theta(years(i)-1);
      } else {
        log_dens(i) = dnorm(x(i), mu(years(i)-1), sigma2(years(i)-1), true);
        pred(i) = log_dens(i) + theta(years(i)-1) + scalar(years(i)-1);
      }
    } else {
      // model is symmetric around mu
      log_dens(i) = dnorm(x(i), mu(years(i)-1), sigma1(years(i)-1), true);
      pred(i) = log_dens(i) + theta(years(i)-1);
    }
    //nll += dnorm(y(i),pred(i),obs_sigma,true);
  }

  Type s1 = 0;
  Type s2 = 0;
  if(family==1) {
    // gaussian, both data and predictions in log space
    nll += sum(dnorm(log(y), pred, obs_sigma, true));
  }
  if(family==2) {
    nll += sum(dpois(y, exp(pred), true));
  }
  if(family==3) {
    for(i = 0; i < n; i++) {
      s1 = exp(pred(i));
      s2 = s1 + pow(s1, Type(2))*obs_sigma;
      nll += dnbinom2(y(i), s1, s2, true);
    }
  }

  // ADREPORT section
  //ADREPORT(theta); // nuisance parameter
  ADREPORT(sigma1);
  ADREPORT(mu);
  ADREPORT(obs_sigma);
  ADREPORT(pred);
  ADREPORT(log_mu_b0);
  ADREPORT(mu_b1);
  ADREPORT(mu_devs);
  ADREPORT(sig1_b0);
  ADREPORT(sig1_b1);
  ADREPORT(lower25);
  ADREPORT(upper75);

  if(asymmetric == 1) {
    // these are only reported for asymmetric model
    ADREPORT(sigma2);
    ADREPORT(sig2_b0);
    ADREPORT(sig2_b1);
  }

  return (-nll);
}
