#include <RcppArmadillo.h>
// #include <expert.hpp>
// #include <type_definitions.hpp>
// #include "utils_parameters.h"


namespace stochvol {

namespace fast_sv {

// utils_parameters.h. line 57

// Encapsulate the four parameters mu, phi, and sigma,
// as one returned object, including whether
// they were updated/accepted.
// ADDED: full conditional posterior for sigma evaluated at zero, signed sigma
struct SampledTheta4SD{
  double mu,
  phi,
  sigma;
  bool mu_accepted,
  phi_accepted,
  sigma_accepted;
  double sigma_density,
  sigma_signed;
};

namespace noncentered {

// utils_parameters.h. line 107

// Two-block update of mu, phi, and sigma -- noncentered.
// For more details, see Kastner and Fr?hwirth-Schnatter (2014)
// MODIFIED to also return full conditional posterior for sigma evaluated at zero, signed sigma
// SampledTheta4SD draw_theta_2block4SD(
//     const arma::vec& log_data2,
//     const double mu,
//     const double phi,
//     const double sigma,
//     const double ht0,
//     const arma::vec& ht,
//     const arma::uvec& r,
//     const PriorSpec& prior_spec,
//     const ExpertSpec_FastSV& expert);



  
  
} // END namespace noncentered

} // END namespace fast_sv

// single_update.h, line 74

// void update_fast_sv4SD(
//     const arma::vec& log_data2,
//     double& mu,
//     double& phi,
//     double& sigma,
//     double& h0,
//     // SK
//     double& sigma_density,
//     double& sigma_signed,
//     //SK
//     arma::vec& h,
//     arma::uvec& r,
//     const PriorSpec& prior_spec,  // old parameters: C0, cT, Bsigma, a0, b0, bmu, Bmu, Gammaprior, truncnormal, dontupdatemu, priorlatent0 feed into this
//     const ExpertSpec_FastSV& expert);  // old parameters: parameterization, centered_baseline, B011inv, B022inv, Mhcontrol, MHsteps feed into this

Rcpp::List svsample_fast_cpp4SD(
    const arma::vec& data,
    const int draws,
    const int burnin,
    const Rcpp::List& priorspec_in,
    const Rcpp::List& startpara,
    const arma::vec& startlatent,
    const Rcpp::List& expert_in);
  
} // END namespace stochvol
