#include <TMB.hpp>
// #include <omp.h>

template <class Type>
Type qthill(Type quantile, Type v, Type mean, Type sigma)
{
  // implements algorithm from Hill et al. 1970
  // this is source of qt in R, and elsewhere -- e.g.
  // https://repo.progsbase.com/repoviewer/no.inductive.libraries/BasicStatistics/0.1.14///HillsAlgorithm396/
  // I've extended this here to incldue mean and sigma
  Type z, flip;
  if(quantile > 0.5){
    flip = 1;
    z = 2*(1 - quantile);
  }else{
    flip = -1;
    z = 2*quantile;
  }

  Type a = 1/(v - 0.5);
  Type b = 48/(a*a);
  Type c = ((20700*a/b - 98)*a - 16)*a + 96.36;
  Type d = ((94.5/(b + c) - 3)/b + 1)*sqrt(a*3.14159265/2)*v;
  Type x = z*d;
  Type y = pow(x, 2/v);//x**(2/v);

  if(y > 0.05 + a){
    x = qnorm(z*0.5, Type(0), Type(1));
    y = x*x;
    if(v < 5){
      c = c + 0.3*(v - 4.5)*(x + 0.6);
    }
    c = c + (((0.05*d*x - 5)*x - 7)*x - 2)*x + b;
    y = (((((0.4*y + 6.3)*y + 36)*y + 94.5)/c - y - 3)/b + 1)*x;
    y = a*y*y;
    if(y > 0.002){
      y = exp(y) - 1;
    }else{
      y = y + 0.5*y*y;
    }
  }else{
    y = ((1/(((v + 6)/(v*y) - 0.089*d - 0.822)*(v + 2)*3) + 0.5/(v + 4))*y - 1)*(v + 1)/(v + 2) + 1/y;
  }

  Type q = sqrt(v*y);
  // flip sign if needed
  q = q * flip;

  return (mean + sigma*q);
}


template<class Type>
  Type objective_function<Type>::operator() ()
{
  using namespace density;

  DATA_VECTOR(y); // vector of counts
  DATA_VECTOR(x); // vector of calendar dates
  DATA_IVECTOR(years); // vector of years to assign each count to
  DATA_IVECTOR(year_levels);
  DATA_IVECTOR(unique_years); // vector containing unique years
  DATA_INTEGER(nLevels); // number of unique years
  DATA_INTEGER(asymmetric); // 0 if false, 1 = true. Whether to estimate same shape/scale parameters before/after mean
  DATA_INTEGER(family); // 1 gaussian, 2 = poisson, 3 = neg bin
  DATA_INTEGER(sig_trend); // 0 if false, 1 = true. Whether to estimate trend parameters with respect to sds
  DATA_INTEGER(mu_trend); // 0 if false, 1 = true. Whether to estimate trend parameters with respect to mean
  DATA_INTEGER(t_model); // 0 if false, 1 = true. Whether to estimate model with student-t tails or normal distribution

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
  PARAMETER(sig2_b0);
  PARAMETER(sig2_b1);
  PARAMETER(log_sigma2);
  PARAMETER_VECTOR(sigma2_devs);
  PARAMETER(log_tdf_1);
  PARAMETER(log_tdf_2);

  // derived parameters
  Type obs_sigma=exp(log_obs_sigma);
  Type tdf_1 = exp(log_tdf_1) + 2;
  Type tdf_2 = exp(log_tdf_2) + 2;
  vector<Type> sigma1(nLevels), mu(nLevels);
  vector<Type> sigma2(nLevels), scalar(nLevels);
  vector<Type> lower25(nLevels), upper75(nLevels);
  vector<Type> range(nLevels); // 75th - 25th percentile
  int i;
  int n = x.size();

  Type nll=0;

  for(i = 0; i < nLevels; i++) {

    if(sig_trend==1) {
      sigma1(i) = exp(sig1_b0 + sig1_b1*Type(unique_years(i)) + sigma1_devs(i));
    } else {
      sigma1(i) = exp(sig1_b0 + sigma1_devs(i));
    }

    if(asymmetric == 1) {
      if(sig_trend==1) {
      sigma2(i) = exp(sig2_b0 + sig2_b1*Type(unique_years(i)) + sigma2_devs(i));
      // scalar(i) is just log(sig2) - log(sig1)
      scalar(i) = sig2_b0 + sig2_b1*Type(unique_years(i)) + sigma2_devs(i) - (sig1_b0 + sig1_b1*Type(unique_years(i)) + sigma1_devs(i));
      } else {
        sigma2(i) = exp(sig2_b0 + sigma2_devs(i));
        scalar(i) = sig2_b0 + sigma2_devs(i) - (sig1_b0 + sigma1_devs(i));
      }
    }

    // trend in in normal space, e.g. not log-linear
    if(mu_trend == TRUE) {
    mu(i) = exp(log_mu_b0) + mu_devs(i) + mu_b1*Type(unique_years(i));
    } else {
      mu(i) = exp(log_mu_b0) + mu_devs(i);
    }

    // random effects contributions of mean and sigma1
    nll += dnorm(mu_devs(i), Type(0.0),exp(log_sigma_mu_devs),true);
    nll += dnorm(sigma1_devs(i),Type(0.0),exp(log_sigma1),true);
    if(asymmetric == 1) {
      nll += dnorm(sigma2_devs(i),Type(0.0),exp(log_sigma2),true);
    }

    if(t_model==0) {
      lower25(i) = qnorm(Type(0.25), mu(i), sigma1(i));
    } else {
      lower25(i) = qthill(Type(0.25),Type(tdf_1), mu(i), sigma1(i));
    }
    if(asymmetric == 1) {
      if(t_model == 0) {
        upper75(i) = qnorm(Type(0.75), mu(i), sigma2(i));
      } else {
        upper75(i) = qthill(Type(0.75),Type(tdf_2), mu(i), sigma2(i));
      }
    } else {
      if(t_model == 0) {
        upper75(i) = qnorm(Type(0.75), mu(i), sigma1(i));
      } else {
        upper75(i) = qthill(Type(0.75),Type(tdf_1), mu(i), sigma1(i));
      }
    }
    range(i) = upper75(i) - lower25(i);
  }

  vector<Type> log_dens(n), pred(n);
  for(i = 0; i < n; i++) {
    if(asymmetric == 1) {
      // model is asymmetric, left side smaller / right side bigger
      if(x(i) < mu(years(i)-1)) {
        if(t_model==0) {
          // model is asymmetric around mu, gaussian tails
          log_dens(i) = dnorm(x(i), mu(years(i)-1), sigma1(years(i)-1), true);
        } else {
          // model is asymmetric around mu, student-t tails
          log_dens(i) = dt((x(i) - mu(years(i)-1)) / sigma1(years(i)-1), Type(tdf_1), true) - log(sigma1(years(i)-1));
        }
        pred(i) = log_dens(i) + theta(years(i)-1);
      } else {
        if(t_model==0) {
          // model is asymmetric around mu, gaussian tails
          log_dens(i) = dnorm(x(i), mu(years(i)-1), sigma2(years(i)-1), true);
        } else {
          // model is asymmetric around mu, student-t tails
          log_dens(i) = dt((x(i) - mu(years(i)-1)) / sigma2(years(i)-1), tdf_2, true) - log(sigma2(years(i)-1));
        }
        pred(i) = log_dens(i) + theta(years(i)-1) + scalar(years(i)-1);
      }
    } else {
      if(t_model==0) {
        // model is symmetric around mu, gaussian tails
        log_dens(i) = dnorm(x(i), mu(years(i)-1), sigma1(years(i)-1), true);
      } else {
        // model is symmetric around mu, student-t tails
        log_dens(i) = dt((x(i) - mu(years(i)-1)) / sigma1(years(i)-1), tdf_1, true) - log(sigma1(years(i)-1));
      }
      pred(i) = log_dens(i) + theta(years(i)-1);
    }
  }

  Type s1 = 0;
  Type s2 = 0;
  if(family==1) {
    // gaussian, both data and predictions in log space
    nll += sum(dnorm(log(y), pred, obs_sigma, true));
  }
  if(family==2) {
    REPORT(nll);
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
  ADREPORT(theta); // nuisance parameter
  ADREPORT(sigma1); // sigma, LHS
  ADREPORT(mu); // mean of curves by year
  ADREPORT(obs_sigma); // obs sd (or phi, NB)
  ADREPORT(pred); // predictions in link space (log)
  ADREPORT(log_mu_b0); // hypermean, log space
  if(mu_trend==1) {
    ADREPORT(mu_b1); // trend in mean
  }
  ADREPORT(mu_devs); // deviations year to year from hypermean
  ADREPORT(sig1_b0); // mean sigma for LHS
  if(sig_trend==1) {
    ADREPORT(sig1_b1); // optional trend parameter for LHS sigmas
  }
  ADREPORT(lower25); // lower quartile
  ADREPORT(upper75); // upper quartile
  ADREPORT(range); // diff between upper and lower quartiles
  if(t_model==1) {
    ADREPORT(tdf_1); // tdf for LHS
  }
  if(asymmetric == 1) {
    // these are only reported for asymmetric model
    ADREPORT(sigma2); // same as above, but RHS optionally
    ADREPORT(sig2_b0);
    if(t_model==1) {
      ADREPORT(tdf_2);
    }
    if(sig_trend==1) {
      ADREPORT(sig2_b1);
    }
    //ADREPORT(scalar);
  }

  return (-nll);
}
