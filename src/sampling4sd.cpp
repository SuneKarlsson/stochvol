#include <RcppArmadillo.h>
#include <expert.hpp>
#include <type_definitions.hpp>
#include "sampling4sd.h"
#include "utils.h"
#include "utils_main.h"
#include "utils_latent_states.h"
#include "utils_parameters.h"
#include "sampling_latent_states.h"

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

namespace stochvol {

namespace fast_sv {

// Borrow from single_update.cc

std::vector<Parameterization> expert_to_strategy(
    const ExpertSpec_FastSV& expert);


namespace noncentered {

// Borrow definition from utils_parameters.cc

std::pair<bool, double> sample_phi(
    const double phi,
    const double ht0,
    const arma::vec& ht,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert);


// utils_parameters.cc, line 455

// MODIFIED to also return full conditional posterior for sigma evaluated at zero, signed sigma
SampledTheta4SD draw_theta_2block4SD(
    const arma::vec& log_data2,
    const double mu,
    const double phi,
    const double sigma,
    const double ht0,
    const arma::vec& ht,
    const arma::uvec& r,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert) {
  R_assert(prior_spec.mu.distribution == PriorSpec::Mu::NORMAL, "The non-centered 2-block sampler is only implemented with a normal prior for mu");
  
  // Draw mu and sigma
  const int T = ht.n_elem;
  double BT11 = std::pow(prior_spec.mu.normal.sd, -2),
    BT12 = 0,
    BT22 = 2*prior_spec.sigma2.gamma.rate,
    bT1 = 0,
    bT2 = prior_spec.mu.normal.mean * BT11,
    // SK
    bTsigma = 0;
  
  for (int j = 0; j < T; j++) {
    const double tmp1 = mix_varinv[r[j]],
                                  tmp2 = (log_data2[j]-mix_mean[r[j]])*tmp1,
                                  tmp3 = ht[j]*tmp1;
    BT11 += tmp1;
    BT12 -= tmp3;
    BT22 += ht[j]*tmp3;
    bT1 += ht[j]*tmp2;
    bT2 += tmp2;
    // SK
    bTsigma += (log_data2[j]-mix_mean[r[j]]-mu)*tmp1;
  }
  
  {
    const double det = BT11*BT22-BT12*BT12;
    BT11 /= det;
    BT12 /= det;
    BT22 /= det;
  }
  
  {
    const double bT1_old = bT1;
    bT1 = BT11*bT1_old + BT12*bT2;
    bT2 = BT12*bT1_old + BT22*bT2;
  }
  
  const double chol11 = std::sqrt(BT11),
    chol12 = (BT12/chol11),
    chol22 = std::sqrt(BT22-chol12*chol12);
  
  const double innov = R::norm_rand(),
    sigma_new = bT1 + chol11*innov,
    mu_new = bT2 + chol12*innov + chol22*R::norm_rand();
  
  // Draw phi
  const auto phi_draw = sample_phi(phi, ht0, ht, prior_spec, expert);
  
  // SK: evaluate full conditional for sigma at zero
  const double BTsigma = 1.0/BT22;
  const double sigden = R::dnorm(0.0,0.0,bTsigma,std::sqrt(BTsigma));
  
  return {mu_new, phi_draw.second, std::fabs(sigma_new), true, phi_draw.first, true, sigden, sigma_new};
}

} // END namespace noncentered

} // END namspace fast_sv

} // END namespace stochvol



namespace stochvol {

namespace fast_sv {

namespace centered {

// sampling_parameters.cc, line 49

// Fast SV step (b): sample mu, phi, sigma - __CENTERED__ version
SampledTheta4SD regression4SD(
    const double mu,
    const double phi,
    const double sigma,
    const double h0,
    const arma::vec& h,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert) {
  SampledTheta ret;
  switch (expert.mh_blocking_steps) {
  case 1:
    // return centered::draw_theta_1block(mu, phi, sigma, h0, h, prior_spec, expert);
    ::Rf_error("Parameter fast_sv$mh_blocking_steps can only be 2 with the modified stochvol.");
  case 2:
    ret = centered::draw_theta_2block(mu, phi, sigma, h0, h, prior_spec, expert);
    return {ret.mu, ret.phi, ret.sigma, ret.mu_accepted, ret.phi_accepted, ret.sigma_accepted, 0.0, 0.0};
  case 3:
    // return centered::draw_theta_3block(mu, phi, sigma, h0, h, prior_spec, expert);
    ::Rf_error("Parameter fast_sv$mh_blocking_steps can only be 2 with the modified stochvol.");
  default:
    ::Rf_error("Parameter fast_sv$mh_blocking_steps should an integer between 1 and 3.");
  }
}

}  // END namespace centered

namespace noncentered {

// sampling_parameters.cc, line 74

// Fast SV step (b): sample mu, phi, sigma - __NONCENTERED__ version
SampledTheta4SD regression4SD(
    const arma::vec& log_data2,
    const double mu,
    const double phi,
    const double sigma,
    const double h0,
    const arma::vec& h,
    const arma::uvec& r,
    const PriorSpec& prior_spec,
    const ExpertSpec_FastSV& expert) {
  switch (expert.mh_blocking_steps) {
  case 1: //[[fallthrough]];
    ::Rf_error("Parameter fast_sv$mh_blocking_steps can only be 2 with the modified stochvol.");
  case 2:
    return noncentered::draw_theta_2block4SD(log_data2, mu, phi, sigma, h0, h, r, prior_spec, expert);
  case 3:
    // SK This will not work without modifying this function as well
    // return noncentered::draw_theta_3block(log_data2, mu, phi, sigma, h0, h, r, prior_spec, expert);
    ::Rf_error("Parameter fast_sv$mh_blocking_steps can only be 2 with the modified stochvol.");
  default:
    ::Rf_error("Parameter fast_sv$mh_blocking_steps should an integer between 1 and 3.");
  }
}

} // END namespace noncentered


// sampling_parameters.cc, line 99

SampledTheta4SD draw_theta4SD(
    const arma::vec& log_data2,
   const double mu,
   const double phi,
   const double sigma,
   const double h0,
   const double ht0,
   const arma::vec& h,
   const arma::vec& ht,
   const arma::uvec& r,
   const PriorSpec& prior_spec,
   const ExpertSpec_FastSV& expert,
   const Parameterization parameterization) {
 switch (parameterization) {
   case Parameterization::CENTERED:
     return centered::regression4SD(mu, phi, sigma, h0, h, prior_spec, expert);
   case Parameterization::NONCENTERED:
     return noncentered::regression4SD(log_data2, mu, phi, sigma, ht0, ht, r, prior_spec, expert);
   default:
     ::Rf_error("draw_theta: Mistake in the switch-case");
 }
}

} // END namspace fast_sv

// single_update.cc, line 72

void update_fast_sv4SD(
    const arma::vec& log_data2,
    double& mu,
    double& phi,
    double& sigma,
    double& h0,
    // SK
    double& sigma_density,
    double& sigma_signed,
    //SK
    arma::vec& h,
    arma::uvec& r,
    const PriorSpec& prior_spec,  // C0, cT, Bsigma, a0, b0, bmu, Bmu, Gammaprior, dontupdatemu, priorlatent0 feed into this
    const ExpertSpec_FastSV& expert) {  // parameterization, centered_baseline, B011inv, B022inv, truncnormal, Mhcontrol, MHsteps feed into this
  // TODO setup validation (put in a "debug" environment in the end)
  // phi is beta
  // sigma2 is either inverse_gamma with shape == 2.5 or gamma with shape == 0.5
  // inverse_gamma prior and mh_blocking_steps != 2 is nyi
  // constant mu implies mh_blocking_steps == 3 (mh_blocking_steps == 1 is nyi and 2 does not make sense because mu is drawn jointly with either phi or sigma2 in the 2-block)
  
  double ht0 = centered_to_noncentered(mu, sigma, h0);
  arma::vec ht = centered_to_noncentered(mu, sigma, h);
  
  if (prior_spec.mu.distribution == PriorSpec::Mu::CONSTANT) {
    mu = 0;
  }
  
  if (expert.update.mixture_indicators) {  // Step (c): sample indicators
    r = fast_sv::draw_mixture_indicators(log_data2, mu, phi, sigma, h);
  }
  
  if (expert.update.latent_vector) {  // Step (a): sample the latent volatilities h:
    const auto latent_new = fast_sv::draw_latent(log_data2, mu, phi, sigma, r, prior_spec, expert);
    
    switch (expert.baseline) {
    case Parameterization::CENTERED:
      h = latent_new.h;
      h0 = latent_new.h0;
      ht = centered_to_noncentered(mu, sigma, h);
      ht0 = centered_to_noncentered(mu, sigma, h0);
      break;
    case Parameterization::NONCENTERED:
      ht = latent_new.h;
      ht0 = latent_new.h0;
      h = noncentered_to_centered(mu, sigma, ht);
      h0 = noncentered_to_centered(mu, sigma, ht0);
      break;
    }
  }
  
  if (expert.update.parameters) {  // Step (b): sample mu, phi, sigma
    const auto strategy = fast_sv::expert_to_strategy(expert);
    for (const auto parameterization : strategy) {
      // Draw theta
      const auto parameter_draw = fast_sv::draw_theta4SD(log_data2, mu, phi, sigma, h0, ht0, h, ht, r, prior_spec, expert, parameterization);
      mu = parameter_draw.mu;
      phi = parameter_draw.phi;
      sigma = parameter_draw.sigma;
      sigma_density = parameter_draw.sigma_density;
      sigma_signed = parameter_draw.sigma_signed;
      // Update latent vectors
      switch (parameterization) {
        case Parameterization::CENTERED:
          ht = centered_to_noncentered(mu, sigma, h);
          ht0 = centered_to_noncentered(mu, sigma, h0);
          break;
        case Parameterization::NONCENTERED:
          h = noncentered_to_centered(mu, sigma, ht);
          h0 = noncentered_to_centered(mu, sigma, ht0);
          break;
      }
    }
  }
}

// utils_main.cc, line 65

void transpose_and_rename4SD(
    const int T,
    NumericMatrix& para,
    NumericMatrix& latent0,
    NumericMatrix& latent,
    NumericMatrix& savagedickey) {
  para = Rcpp::transpose(para);
  latent = Rcpp::transpose(latent);
  savagedickey = Rcpp::transpose(savagedickey);
  
  {  // colnames in para
    const Rcpp::CharacterVector col_names {"mu", "phi", "sigma"};
    colnames(para) = col_names;
  }
  {  // colnames in latent0
    colnames(latent0) = Rcpp::CharacterVector({"h_0"});
  }
  {  // colnames in latent
    const unsigned int ncol = latent.ncol();
    Rcpp::CharacterVector col_names(ncol);
    for (unsigned int c = 1; c <= ncol; c++) {
      std::string name = "h_";
      name += std::to_string(T-ncol+c);
      col_names[c-1] = name;
    }
    colnames(latent) = col_names;
  }
  {  // colnames in savagedickey 
    colnames(savagedickey) = Rcpp::CharacterVector({"sigma_density","sigma_signed"});
  }
}

// utils_main.cc, line 116

List cleanup4SD(
    const int T,
    NumericMatrix& para,
    NumericMatrix& latent0,
    NumericMatrix& latent,
    IntegerMatrix& mixture_indicators,
    NumericMatrix& savagedickey) {
  transpose_and_rename4SD(T, para, latent0, latent, savagedickey);
  
  mixture_indicators = Rcpp::transpose(mixture_indicators);
  
  {  // colnames in mixture_indicators
    const unsigned int ncol = mixture_indicators.ncol();
    Rcpp::CharacterVector col_names(ncol);
    for (unsigned int c = 1; c <= ncol; c++) {
      std::string name = "r_";
      name += std::to_string(T-ncol+c);
      col_names[c-1] = name;
    }
    colnames(mixture_indicators) = col_names;
  }
  
  List val = List::create(
    _["para"] = para,
    _["latent"] = latent,
    _["latent0"] = latent0,
    _["indicators"] = mixture_indicators + 1u,
    _["savagedickey"] = savagedickey);
  
  return val;
}


// sampling_man.cc, line 48

// Wrapper function around fast SV.
// See documentation above the declaration


List svsample_fast_cpp4SD(
    const arma::vec& data,
    const int draws,
    const int burnin,
    const Rcpp::List& priorspec_in,
    const Rcpp::List& startpara,
    const arma::vec& startlatent,
    const Rcpp::List& expert_in) {
  ::Rprintf("A\n");
  const unsigned int T = data.size();
  const int T_int = T;
  ::Rprintf("B\n");
  //  const unsigned int p = X.n_cols;
  
  // const int thinpara = 1,
  //   thinlatent = 1,
  //   thintime = 1;
  ::Rprintf("BA\n");
  const SEXP priorlatent0_sexp = priorspec_in["latent0_variance"];
  ::Rprintf("BB\n");
  const List priormu = priorspec_in["mu"],
                           priorphi = priorspec_in["phi"],
                                          priorsigma2 = priorspec_in["sigma2"],
                                                            priornu = priorspec_in["nu"],
                                                                          priorrho = priorspec_in["rho"],
                                                                                         priorbeta = priorspec_in["beta"];
  
  ::Rprintf("BBB\n");
  
  
  
  const PriorSpec prior_spec = list_to_priorspec(priorspec_in);
  ::Rprintf("C\n");
  const ExpertSpec_FastSV expert = list_to_fast_sv(expert_in, true);
  ::Rprintf("D\n");
  
  // const bool is_regression = !::ISNA(X.at(0,0)),
  const bool is_heavy_tail = prior_spec.nu.distribution != PriorSpec::Nu::INFINITE;
  if ( is_heavy_tail ) {
    ::Rf_error("Heavy tails not implemented with svsample4SD");
  }
  // if (prior_spec.mu.distribution == PriorSpec::Mu::CONSTANT && expert.mh_blocking_steps == 1) { // not implemented (would be easy, though)
  //   ::Rf_error("Single block update leaving mu constant is not yet implemented");
  // }
  if ( expert.mh_blocking_steps !=2 ) {
    ::Rf_error("Only 2-block updating with svsample4SD");
  }
  
  // shortcuts / precomputations
//  const int thintime = determine_thintime(T, keeptime_in);
  
  // number of MCMC draws
  // const int N = burnin + draws;
  
  // verbosity control
  // const int chain = print_settings["chain"],
  //                                 n_chains = print_settings["n_chains"];
  // const bool quiet = print_settings["quiet"],
  //                                  verbose = !quiet,
  //                                  single_chain = n_chains == 1;
  
  // initialize the variables:
  double mu = startpara["mu"],
                       phi = startpara["phi"],
                                      sigma = startpara["sigma"]; //,
//                                                       nu = startpara["nu"];
  // extra variables for Savage-Dickey
  double sigma_density, sigma_signed;
  
  arma::vec h = startlatent;  // contains h1 to hT, but not h0!
  double h0 = startpara["latent0"];
  // arma::vec beta = startpara["beta"];
  // arma::vec tau = expert_in["init_tau"];
  // if (tau.n_elem == 1) {
  //   const double tau_elem = tau[0];
  //   tau.set_size(T);
  //   tau.fill(tau_elem);
  // } else if (tau.n_elem != T) {
  //   ::Rf_error("Bad initialization for the vector tau. Should have length %d, received length %d, first element %f", T, tau.n_elem, tau[0]);
  // }
  arma::uvec r = expert_in["init_indicators"];  // mixture indicators
  if (r.n_elem == 1) {
    const double r_elem = r[0];
    r.set_size(T);
    r.fill(r_elem);
  } else if (r.n_elem != T) {
    ::Rf_error("Bad initialization for the vector of mixture indicators. Should have length %d, received length %d, first element %u", T, r.n_elem, r[0]);
  }
  r -= 1u;
  R_assert(arma::all(r <= 9u), "Initial values of the mixture indicators need to be between 1 and 10 inclusive.");
  
  // keep some arrays cached
  // arma::vec data_demean = is_regression ? data - X * beta : data,
  //   log_data2_normal = arma::log(arma::square(data_demean) / tau + offset),  // standardized "data" (different for t-errors, "normal data")
  //   exp_h_half_inv;
  arma::vec log_data2_normal = arma::log(arma::square(data)),
    exp_h_half_inv;
  clamp_log_data2(log_data2_normal);
  const arma::vec conditional_mean(data.n_elem, arma::fill::zeros),
  conditional_sd(data.n_elem, arma::fill::ones);
  
  // storage
  const bool keep_r = expert_in["store_indicators"];
  // const int para_draws = draws / thinpara;
  // Rcpp::NumericMatrix para_store(5, para_draws);
  Rcpp::NumericMatrix para_store(3, draws);
//  Rcpp::NumericMatrix beta_store(p, is_regression * para_draws);
  // const int latent_length = T / thintime;  // thintime must be either 1 or T
  // const int latent_draws = draws / thinlatent;
  // Rcpp::NumericVector latent0_store(latent_draws);
  // Rcpp::NumericMatrix latent_store(latent_length, latent_draws);
  Rcpp::NumericVector latent0_store(draws);
  Rcpp::NumericMatrix latent_store(T, draws);
  //  Rcpp::NumericMatrix tau_store(latent_length, keep_tau * latent_draws);
  // Rcpp::IntegerMatrix r_store(latent_length, keep_r * latent_draws);
  Rcpp::IntegerMatrix r_store(T, keep_r);
//  Rcpp::NumericVector correction_weight_latent(correct_model_specification * latent_draws);
//  Rcpp::NumericVector correction_weight_para(correct_model_specification * para_draws);
  Rcpp::NumericMatrix sd_store(2,draws);
  // // a warning about intended use
  // if (correct_model_specification and (is_regression or is_heavy_tail)) {
  //   Rcpp::warning("Function 'svsample_fast_cpp' is not meant to do correction for model misspecification along with regression/heavy tails. Function 'svsample_general_cpp' can do this correction.");
  // }
  
  // initializes the progress bar
  // "show" holds the number of iterations per progress sign
  // int show = 0;
  // if (verbose) {
  //   if (single_chain) {
  //     show = progressbar_init(N);
  //   } else {
  //     show = chain_print_init(chain, burnin, draws);
  //   }
  // }
  
  for (int i = -burnin + 1; i < draws + 1; i++) {  // BEGIN main MCMC loop
    if (i % 20 == 0) {
      ::R_CheckUserInterrupt();
    }
    
    // const bool parasave_round = (i - 1) % thinpara == thinpara - 1,  // is this a parameter saving round?
    //   latentsave_round = (i - 1) % thinlatent == thinlatent - 1;  // is this a latent saving round?
    
    // // print a progress sign every "show" iterations
    // if (verbose and i % show == 0) {
    //   if (single_chain) {
    //     progressbar_print();
    //   } else {
    //     show = chain_print(chain, i, burnin, draws);
    //   }
    // }
    
    // a single MCMC update: update indicators, latent volatilities,
    // and parameters ONCE
    update_fast_sv4SD(log_data2_normal, mu, phi, sigma, h0, sigma_density, sigma_signed, h, r, prior_spec, expert);
    // if (correct_model_specification or is_regression or is_heavy_tail) {
    //   exp_h_half_inv = arma::exp(-.5 * h);
    // }
    
    // // update tau and nu
    // if (is_heavy_tail) {
    //   update_t_error(data_demean % exp_h_half_inv, tau, conditional_mean, conditional_sd, nu, prior_spec, false);
    //   // update cached data arrays after the regression
    // }
    // 
    // // update beta
    // if (is_regression) {
    //   const arma::vec normalizer = exp_h_half_inv / arma::sqrt(tau);
    //   update_regressors(
    //     data % normalizer,
    //     X.each_col() % normalizer,
    //     beta, prior_spec);
    //   // update cached data arrays
    //   data_demean = data - X * beta;
    //   if (is_heavy_tail) {
    //     log_data2_normal = arma::log(arma::square(data_demean) / tau + offset);
    //   } else {
    //     log_data2_normal = arma::log(arma::square(data_demean) + offset);
    //   }
    //   clamp_log_data2(log_data2_normal);
    // } else if (is_heavy_tail) {
    //   log_data2_normal = arma::log(arma::square(data_demean) / tau + offset);
    //   clamp_log_data2(log_data2_normal);
    // } else {
    //   ;  // no update needed
    // }
    
    
    // store draws
    // double correction_weight_i = 1;
    // if ((parasave_round or latentsave_round) and correct_model_specification) {
    //   correction_weight_i = fast_sv::compute_correction_weight(data_demean, log_data2_normal, h, 1/exp_h_half_inv);
    // }
    if (i >= 1 ) { // and parasave_round) {
      // const unsigned int index = (i - 1) / thinpara;
      // save_para_sample(index, mu, phi, sigma, nu, beta, para_store, beta_store, is_regression);
      // instead of inlined function
      para_store(0, i) = mu;
      para_store(1, i) = phi;
      para_store(2, i) = sigma;
      
      // if (correct_model_specification) {
      //   correction_weight_para[index] = correction_weight_i;
      // }
    // }
    // if (i >= 1 ) { // and latentsave_round) {
      // const unsigned int index = (i - 1) / thinlatent;
      // save_latent_sample(index, h0, h, tau, r, thintime, latent_length, latent0_store, latent_store, tau_store, r_store, keep_tau and is_heavy_tail, keep_r);
      // instead of inlined function
      latent0_store[i] = h0;
      // for (int volind = 0, thincol = 0; thincol < latent_length; volind++, thincol++) {
      //   latent_store(volind, 0) = h[thintime * (thincol + 1) - 1];
      // }
      for (int volind = 0; volind < T_int; volind++) {
        latent_store(volind,i) = h[volind];
      }
      if ( keep_r) {
        for (int volind = 0; volind < T_int; volind++) {
          r_store.at(volind, 0) = r[volind];
        }
      }
      
      // if (correct_model_specification) {
      //   correction_weight_latent[index] = correction_weight_i;
      // }
      sd_store(0,i) = sigma_density;
      sd_store(1,i) = sigma_signed;
    }
  }  // END main MCMC loop
  
  // if (verbose) {
  //   if (single_chain) {
  //     progressbar_finish(N);  // finalize progress bar
  //   } else {
  //     chain_print_finish(chain);
  //   }
  // }
  
  Rcpp::NumericMatrix latent0_store_mat(latent0_store.length(), 1);
  latent0_store_mat(_, 0) = latent0_store;
  return cleanup4SD(T, para_store, latent0_store_mat, latent_store, r_store,sd_store);
}


} // END namespace stochvol

