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

template <class Type>
Type dgnorm(Type x, Type mu, Type alpha, Type beta)
{
  // implements dgnorm
  // copied from https://github.com/maryclare/gnorm/blob/master/R/gnorm.R
  Type z = -pow(fabs(x - mu)/alpha,beta) + log(beta) - (log(2.0) +
    log(alpha) + lgamma(1.0/beta));
  return(z);
}

template <class Type>
Type qgnorm(Type quantile, Type mu, Type alpha, Type beta)
{
  // implements qgnorm
  // copied from https://github.com/maryclare/gnorm/blob/master/R/gnorm.R
  Type p = quantile;
  if(p > 0.5) {
    p = 1 - p;
  }
  Type sign = 0;
  if(p - 0.5 > 0.0) {
    sign = 1;
  }
  if(p - 0.5 < 0.0){
    sign = -1;
  }
  Type shape =1.0/beta;
  Type scale = 1.0/pow(1.0/alpha, beta);
  //return(sign(p - 0.5) * qgamma(abs(p - 0.5) * 2, shape = 1/beta, scale = 1/lambda)^(1/beta) + mu)
  //return (sign*pow(qgamma(fabs(p - 0.5)*2, shape = shape, scale = scale), 1.0/beta) + mu);
  //Type z = sign*exp(log(qgamma(fabs(p - 0.5)*2, shape = shape, scale = scale))/beta) + mu;
  Type z = sign*exp(log(qgamma(fabs(p - 0.5)*2, shape, scale))/beta) + mu;
  return (z);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;

  DATA_VECTOR(y); // vector of counts
  //DATA_IVECTOR(yint); // counts for poisson/negbin models
  DATA_VECTOR(x); // vector of calendar dates
  DATA_IVECTOR(years); // vector of years to assign each count to
  DATA_IVECTOR(year_levels);
  DATA_IVECTOR(unique_years); // vector containing unique years
  DATA_INTEGER(nLevels); // number of unique years
  DATA_INTEGER(asymmetric); // 0 if false, 1 = true. Whether to estimate same shape/scale parameters before/after mean
  DATA_INTEGER(family); // 1 gaussian, 2 = poisson, 3 = neg bin
  //DATA_INTEGER(sig_trend); // 0 if false, 1 = true. Whether to estimate trend parameters with respect to sds
  //DATA_INTEGER(mu_trend); // 0 if false, 1 = true. Whether to estimate trend parameters with respect to mean
  DATA_INTEGER(tail_model); // 0 if gaussian, 1 = student_t, 2 = generalized normal for tails
  DATA_INTEGER(est_sigma_re); // 0 if FALSE, 1 = TRUE. Whether to estimate deviations as random effects
  DATA_INTEGER(est_mu_re); // 0 if FALSE, 1 = TRUE. Whether to estimate deviations as random effects
  //DATA_VECTOR(mu_cov); // vector of covariates for mean trend
  //DATA_VECTOR(sigma_cov); // vector of covariates for sigma trend
  DATA_MATRIX(mu_mat); // matrix of covariates for mean trend
  DATA_MATRIX(sig_mat);  // matrix of covariates for sigma trend

  PARAMETER(log_sigma1);// log sd of random effects on sigma1
  PARAMETER_VECTOR(sigma1_devs);// random effects on sigma1
  PARAMETER(log_sigma2);// log sd of random effects on sigma2
  PARAMETER_VECTOR(sigma2_devs);// random effects on sigma2
  PARAMETER_VECTOR(theta); // scaling parameter
  PARAMETER_VECTOR(mu_devs);// random effects on mean
  PARAMETER(log_sigma_mu_devs); // log sd of random effects on mean
  PARAMETER(log_tdf_1);// log of Student t df 1
  PARAMETER(log_tdf_2);// log of Student t df 2
  PARAMETER(log_beta_1);
  PARAMETER(log_beta_2);
  PARAMETER(log_obs_sigma);// log sd of observation error

  PARAMETER_VECTOR(b_mu); // coefficients for covariates on mean
  PARAMETER_VECTOR(b_sig1); // coefficients for covariates on sigma1
  PARAMETER_VECTOR(b_sig2); // coefficients for covariates on sigma2

  // derived parameters
  Type obs_sigma=exp(log_obs_sigma);
  Type tdf_1 = exp(log_tdf_1) + 2;
  Type tdf_2 = exp(log_tdf_2) + 2;
  Type beta_1 = exp(log_beta_1);///(1+exp(log_beta_1))*20;
  Type beta_2 = exp(log_beta_2);///(1+exp(log_beta_2))*20;
  vector<Type> sigma1(nLevels), mu(nLevels);
  vector<Type> sigma2(nLevels), scalar(nLevels);
  vector<Type> alpha1(nLevels), alpha2(nLevels);
  vector<Type> logalpha1(nLevels), logalpha2(nLevels);
  vector<Type> lower25(nLevels), upper75(nLevels);
  vector<Type> range(nLevels); // 75th - 25th percentile
  int i;
  int n = y.size();
  vector<Type> log_dens(n);
  vector<Type> pred(n);

  // calculations for beta for gnorm dist if implemented
  vector<Type> beta_ratio(2);
  if(tail_model == 2) {
    beta_ratio(1) = sqrt(exp(lgamma(1.0/Type(beta_1))) / exp(lgamma(3.0/Type(beta_1))));
    if(asymmetric == 1) {
      beta_ratio(2) = sqrt(exp(lgamma(1.0/Type(beta_2))) / exp(lgamma(3.0/Type(beta_2))));
    }
  }

  Type nll=0;
  //Type temp_calc=0;

  // random effects components
  for(i = 0; i < nLevels; i++) {
    // random effects contributions of mean and sigma1
    if(est_mu_re==1) {
      nll += dnorm(mu_devs(i), Type(0.0), exp(log_sigma_mu_devs), true);
    }
    if(est_sigma_re==1) {
      nll += dnorm(sigma1_devs(i),Type(0.0),exp(log_sigma1),true);
      if(asymmetric == 1) {
        nll += dnorm(sigma2_devs(i),Type(0.0),exp(log_sigma2),true);
      }
    }
  }

  // fixed effects for mu and both sigmas
  mu = mu_mat * b_mu;
  sigma1 = sig_mat * b_sig1;
  if(asymmetric == 1) sigma2 = sig_mat * b_sig2;

  for(i = 0; i < nLevels; i++) {

    //temp_calc = sigma1_devs(i);
    //if(sig_trend==1) {
    //  temp_calc += temp_calc + sig1_b1*Type(unique_years(i));
    //}
    //sigma1(i) = exp(temp_calc);
    if(est_sigma_re == 1) {
      sigma1(i) += sigma1_devs(i);
      if(asymmetric == 1) sigma2(i) += sigma2_devs(i);
    }

    if(asymmetric == 1) scalar(i) = sigma2(i) - sigma1(i);//{
      //temp_calc = sigma2_devs(i);
      //if(sig_trend==1) {
        //temp_calc += sigma2_devs(0) + sig2_b1*Type(unique_years(i));
        //temp_calc += sigma2_devs(0) + sig2_b1*Type(sigma_cov(i));
        // scalar(i) is just log(sig2) - log(sig1)
        //scalar(i) = sig2_b1*Type(unique_years(i)) + sigma2_devs(i) - (sig1_b1*Type(unique_years(i)) + sigma1_devs(i));
        //scalar(i) = sig2_b1*Type(sigma_cov(i)) + sigma2_devs(i) - (sig1_b1*Type(sigma_cov(i)) + sigma1_devs(i));
      //} else {
        //scalar(i) = sigma2_devs(i) - sigma1_devs(i);
      //}
      //sigma2(i) = exp(temp_calc);
    //}

    // calculate alphas if the gnorm model is used
    if(tail_model == 2) {
      alpha1(i) = sigma1(years(i)-1)*beta_ratio(1);
      if(asymmetric==1) {
        alpha2(i) = sigma2(years(i)-1)*beta_ratio(2);
      }
    }

    // trend in in normal space, e.g. not log-linear
    //mu(i) = mu_devs(i);
    if(est_mu_re == 1) {
      mu(i) += mu_devs(i);
    }
    //if(mu_trend == 1) {
      //mu(i) += mu_b1*Type(unique_years(i));
    //  mu(i) += mu_b1*Type(mu_cov(i));
    //}

    // this is all for calculating quantiles on LHS
    if(tail_model==0) {
      lower25(i) = qnorm(Type(0.25), mu(i), sigma1(i));
    } else {
      if(tail_model == 1) {
        lower25(i) = qthill(Type(0.25),Type(tdf_1), mu(i), sigma1(i));
      } else {
        // gnorm
        lower25(i) = qgnorm(Type(0.25), mu(i), sigma1(i)*beta_ratio(1), beta_1);
      }
    }
    // this is all for calculating quantiles on RHS
    if(asymmetric == 1) {
      if(tail_model == 0) {
        upper75(i) = qnorm(Type(0.75), mu(i), sigma2(i));
      } else {
        if(tail_model == 1) {
          upper75(i) = qthill(Type(0.75),Type(tdf_2), mu(i), sigma2(i));
        } else {
          // gnorm
          upper75(i) = qgnorm(Type(0.75), mu(i), sigma2(i)*beta_ratio(2), beta_2);
        }
      }
    } else {
      if(tail_model == 0) {
        upper75(i) = qnorm(Type(0.75), mu(i), sigma1(i));
      } else {
        if(tail_model==1) {
          upper75(i) = qthill(Type(0.75),Type(tdf_1), mu(i), sigma1(i));
        } else {
          // gnorm
          upper75(i) = qgnorm(Type(0.75), mu(i), sigma1(i)*beta_ratio(1), beta_1);
        }
      }
    }
    range(i) = upper75(i) - lower25(i);
  }

  // this is for the predictions
  for(i = 0; i < y.size(); i++) {
    if(asymmetric == 1) {
      // model is asymmetric, left side smaller / right side bigger
      if(x(i) < mu(years(i)-1)) {
        if(tail_model==0) {
          // model is asymmetric around mu, gaussian tails
          log_dens(i) = dnorm(x(i), mu(years(i)-1), sigma1(years(i)-1), true);
        } else {
          if(tail_model==1) {
            // model is asymmetric around mu, student-t tails
            log_dens(i) = dt((x(i) - mu(years(i)-1)) / sigma1(years(i)-1), Type(tdf_1), true) - log(sigma1(years(i)-1));
          } else {
            // gnorm, copied from maryclare/gnorm
            // alpha = sqrt( var * gamma(1/beta) / gamma(3/beta) ), alpha = sigma(1)*beta_ratio(1)
            log_dens(i) = dgnorm(x(i), mu(years(i)-1), alpha1(years(i)-1), beta_1);
          }
        }
        pred(i) = log_dens(i) + theta(years(i)-1);
      } else {
        if(tail_model==0) {
          // model is asymmetric around mu, gaussian tails
          log_dens(i) = dnorm(x(i), mu(years(i)-1), sigma2(years(i)-1), true);
        } else {
          if(tail_model==1) {
            // model is asymmetric around mu, student-t tails
            log_dens(i) = dt((x(i) - mu(years(i)-1)) / sigma2(years(i)-1), tdf_2, true) - log(sigma2(years(i)-1));
          } else {
            // gnorm, copied from maryclare/gnorm
            // alpha = sqrt( var * gamma(1/beta) / gamma(3/beta) ), alpha = sigma(1)*beta_ratio(1)
            log_dens(i) = dgnorm(x(i), mu(years(i)-1), alpha2(years(i)-1), beta_2);
          }
        }
        pred(i) = log_dens(i) + theta(years(i)-1) + scalar(years(i)-1);
      }
    } else {
      if(tail_model==0) {
        // model is symmetric around mu, gaussian tails
        log_dens(i) = dnorm(x(i), mu(years(i)-1), sigma1(years(i)-1), true);
      } else {
        if(tail_model==1) {
          // model is symmetric around mu, student-t tails
          log_dens(i) = dt((x(i) - mu(years(i)-1)) / sigma1(years(i)-1), tdf_1, true) - log(sigma1(years(i)-1));
        } else {
          // gnorm, copied from maryclare/gnorm
          // alpha = sqrt( var * gamma(1/beta) / gamma(3/beta) ), alpha = sigma(1)*beta_ratio(1)
          log_dens(i) = dgnorm(x(i), mu(years(i)-1), alpha1(years(i)-1), beta_1);
        }
      }
      pred(i) = log_dens(i) + theta(years(i)-1);
    }
  }

  // this is the likelihood
  Type s1 = 0;
  Type s2 = 0;

  if(family==1) {
    // gaussian, both data and predictions in log space
    nll += sum(dnorm(log(y), pred, obs_sigma, true));
  }
  if(family==2) {
    for(i = 0; i < n; i++) {
      //std::cout << pred(i) << std::endl;
      nll += dpois(y(i), exp(pred(i)), true);
    }
  }
  if(family==3) {
    for(i = 0; i < n; i++) {
      s1 = exp(pred(i));
      //s2 = s1 + pow(s1, Type(2))*obs_sigma;
      s2 = 2. * s1 - log_obs_sigma; // log(var - mu)
      nll += dnbinom_robust(y(i), s1, s2, true);
      //std::cout << dnbinom2(y(i), s1, s2, true) << std::endl;
      //std::cout << s2 << std::endl;
      //nll += dnbinom2(y(i), s1, s2, true);
    }
  }

  // ADREPORT section
  ADREPORT(theta); // nuisance parameter
  ADREPORT(sigma1); // sigma, LHS
  ADREPORT(mu); // mean of curves by year
  ADREPORT(b_mu); // sigma, LHS
  ADREPORT(b_sig1); // mean of curves by year
  if(family != 2) {
    ADREPORT(obs_sigma); // obs sd (or phi, NB)
  }
  ADREPORT(pred); // predictions in link space (log)

  ADREPORT(lower25); // lower quartile
  ADREPORT(upper75); // upper quartile
  ADREPORT(range); // diff between upper and lower quartiles
  if(tail_model==1) {
    ADREPORT(tdf_1); // tdf for LHS
  }
  if(tail_model==2) {
    ADREPORT(beta_1); // tdf for LHS
  }
  if(asymmetric == 1) {
    // these are only reported for asymmetric model
    ADREPORT(b_sig2); // mean of curves by year
    if(est_sigma_re==1) {
      ADREPORT(sigma2); // same as above, but RHS optionally
    }
    if(tail_model==1) {
      ADREPORT(tdf_2);
    }

  }
  return (-nll);
}
