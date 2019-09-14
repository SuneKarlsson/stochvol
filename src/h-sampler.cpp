#include <RcppArmadillo.h>
#include "aug-kalman-filter.h"
#include "h-sampler.h"
#include "h-utils.h"
#include "mixture-state-sampler.h"
#include "simulation-smoother.h"
#include "auxmix.h"
#include "parameterization.h"
#include <chrono>
#include <cmath>

using namespace Rcpp;

namespace stochvol {
  extern std::chrono::steady_clock::duration measure;
}

arma::vec draw_latent(
    const arma::vec& y,
    const arma::vec& y_star,
    const arma::ivec& d,
    const arma::vec& h,
    const arma::vec& ht,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const double priormu_mu,
    const double priormu_sigma,
    const bool correct) {
  // Draw h from AUX
  const arma::vec s = draw_s_auxiliary(y_star, d, h, ht, phi, rho, sigma2, mu, Parameterization::CENTERED);
  const auto time_start = std::chrono::steady_clock::now();
  const arma::vec proposed = draw_h_auxiliary(y_star, d, s, phi, rho, sigma2, mu, priormu_mu, priormu_sigma, Parameterization::CENTERED);
  const auto time_end = std::chrono::steady_clock::now();
  stochvol::measure += time_end-time_start;

  if (correct) {
    return correct_latent_auxiliaryMH(y, y_star, d, h, ht, proposed, phi, rho, sigma2, mu);
  } else {
    return proposed;
  }
}
 
arma::vec draw_h_auxiliary(
    const arma::vec& y_star,
    const arma::ivec& d,
    const arma::vec& s,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const double priormu_mu,
    const double priormu_sigma,
    const Parameterization centering) {
  static arma::vec mixing_a, mixing_b, mixing_m, mixing_v;
  mixing_a.resize(s.size()), mixing_b.resize(s.size()), mixing_m.resize(s.size()), mixing_v.resize(s.size());
  std::transform(s.cbegin(), s.cend(), mixing_a.begin(), [](const int selem) -> double {return mix_a[selem];});
  std::transform(s.cbegin(), s.cend(), mixing_b.begin(), [](const int selem) -> double {return mix_b[selem];});
  std::transform(s.cbegin(), s.cend(), mixing_m.begin(), [](const int selem) -> double {return mix_mean[selem];});
  std::transform(s.cbegin(), s.cend(), mixing_v.begin(), [](const int selem) -> double {return sqrt(mix_var[selem]);});
  
  const List filter_result = aug_kalman_filter(phi, rho, sigma2, mixing_a, mixing_b, mixing_m, mixing_v, d, y_star, priormu_mu, pow(priormu_sigma, 2), centering);
  
  const List smoothing_result = simulation_smoother(mu, filter_result, centering);
  const arma::vec eta = smoothing_result["eta"];
  const double eta0 = as<NumericVector>(smoothing_result["eta0"])[0];

  const int n = as<NumericVector>(filter_result["D"]).size();
  arma::vec h(n, arma::fill::zeros);
  arma::vec dt;
  switch (centering) {
    case Parameterization::CENTERED:
    h[0] = mu + eta0;
    dt = mu*(1-phi) + rho*sqrt(sigma2)*(d%mixing_a%exp(mixing_m/2));
    break;
    case Parameterization::NONCENTERED:
    h[0] = eta0;
    dt = rho*(d%mixing_a%exp(mixing_m/2));
    break;
  }

  for (int i = 0; i < n-1; i++) {
    h[i+1] = dt[i] + phi*h[i] + eta[i];
  }

  return h;
}
 
arma::vec draw_h_auxiliary_awol(
    const arma::vec& y_star,
    const arma::ivec& d,
    const arma::vec& s,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const double priormu_mu,
    const double priormu_sigma,
    const Parameterization centering) {
  const int n = s.size();
  const int mix_count = sizeof(mix_prob)/sizeof(mix_prob[0]);

  static arma::vec mixing_a, mixing_b, mixing_m, mixing_v;
  mixing_a.resize(n), mixing_b.resize(n), mixing_m.resize(n), mixing_v.resize(n);
  std::transform(s.cbegin(), s.cend(), mixing_a.begin(), [](const int selem) -> double {return mix_a[selem];});
  std::transform(s.cbegin(), s.cend(), mixing_b.begin(), [](const int selem) -> double {return mix_b[selem];});
  std::transform(s.cbegin(), s.cend(), mixing_m.begin(), [](const int selem) -> double {return mix_mean[selem];});
  std::transform(s.cbegin(), s.cend(), mixing_v.begin(), [](const int selem) -> double {return sqrt(mix_var[selem]);});
  
  arma::vec h(n, arma::fill::zeros);
  arma::vec omega_diag(n, arma::fill::zeros);
  arma::vec omega_offdiag(n-1, arma::fill::zeros);
  arma::vec e_12m(mix_count, arma::fill::zeros);
  for (int i = 0; i < mix_count; i++) {
    e_12m[i] = std::exp(.5*mix_mean[i]);
    // TODO implement awol based on McCausland and the Mathematica notebook on the wu cloud
    // run devtools::load_all() and then svlsample(exrates[, "CHF"]) to do measurements using the manual timing here
  }

  return h;
}

arma::vec correct_latent_auxiliaryMH(
    const arma::vec& y,
    const arma::vec& y_star,
    const arma::ivec& d,
    const arma::vec& h,
    const arma::vec& ht,
    const arma::vec& proposed,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu) {
    //const CharacterVector centering,

  // Calculate MH acceptance ratio
  const double hlp1 = h_log_posterior(proposed, y, phi, rho, sigma2, mu);
  const double hlp2 = h_log_posterior(h, y, phi, rho, sigma2, mu);
  const double halp1 = h_aux_log_posterior(proposed, y_star, d, phi, rho, sigma2, mu);
  const double halp2 = h_aux_log_posterior(h, y_star, d, phi, rho, sigma2, mu);
  const double log_acceptance = hlp1-hlp2-(halp1-halp2);
  arma::vec result;
  if (log_acceptance > 0 || exp(log_acceptance) > R::runif(0, 1)) {
    result = proposed;
  } else {
    result = h;
  }

  return result;
}

