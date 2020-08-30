#ifndef _UTILS_LATENT_STATES_H_
#define _UTILS_LATENT_STATES_H_

// Constants and functions related to the sampling of the latent states.
// This includes the all-without-a-loop (AWOL, McCausland et al., 2011)
// sampling of the latent vector in the auxiliary model from Omori et al. (2007)
// and the corresponding correction step thereafter

#include <RcppArmadillo.h>

namespace stochvol {

const arma::vec::fixed<10> mix_prob {.00609, .04775, .13057, .20674, .22715, .18842, .12047, .05591, .01575, .00115};
const arma::vec::fixed<10> mix_mean {1.92677, 1.34744, .73504, .02266, -.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000};
const arma::vec::fixed<10> mix_var {.11265, .17788, .26768, .40611, .62699, .98583, 1.57469, 2.54498, 4.16591, 7.33342};
const arma::vec::fixed<10> mix_a {1.01418, 1.02248, 1.03403, 1.05207, 1.08153, 1.13114, 1.21754, 1.37454, 1.68327, 2.50097};
const arma::vec::fixed<10> mix_b {0.50710, 0.51124, 0.51701, 0.52604, 0.54076, 0.56557, 0.60877, 0.68728, 0.84163, 1.25049};
const arma::vec::fixed<10> mix_varinv {
8.8770528184642696345463264151476323604583740234375000000,
5.6217674836968738460996064532082527875900268554687500000,
3.7358039450089663979781562375137582421302795410156250000,
2.4623870379946319886244054941926151514053344726562500000,
1.5949217690872261599110970564652234315872192382812500000,
1.0143736749743870184659044753061607480049133300781250000,
0.6350456280283738319525355109362863004207611083984375000,
0.3929303963095977514363710270117735490202903747558593750,
0.2400435919162919873315331642515957355499267578125000000,
0.1363620248124340350592831327958265319466590881347656250};
const arma::vec::fixed<10> mix_2varinv {
4.4385264092321348172731632075738161802291870117187500000,
2.8108837418484369230498032266041263937950134277343750000,
1.8679019725044831989890781187568791210651397705078125000,
1.2311935189973159943122027470963075757026672363281250000,
0.7974608845436130799555485282326117157936096191406250000,
0.5071868374871935092329522376530803740024566650390625000,
0.3175228140141869159762677554681431502103805541992187500,
0.1964651981547988757181855135058867745101451873779296875,
0.1200217959581459936657665821257978677749633789062500000,
0.0681810124062170175296415663979132659733295440673828125};
const arma::vec::fixed<10> mix_pre {
-4.0093723912083900628999799664597958326339721679687500000,
-2.1784531553855770447114537091692909598350524902343750000,
-1.3768642766903782526100030736415646970272064208984375000,
-1.1257277037836319610875079888501204550266265869140625000,
-1.2487323430568648685579091761610470712184906005859375000,
-1.6619460888428292388852014482836239039897918701171875000,
-2.3433837334574310062862423365004360675811767578125000000,
-3.3510734196563021214387845247983932495117187500000000000,
-4.8643822832849297199686589010525494813919067382812500000,
-7.7642143280080739842219372803810983896255493164062500000};

inline
Rcpp::List get_omori_constants () {
  return Rcpp::List::create(
      Rcpp::_["prob"] = mix_prob,
      Rcpp::_["mean"] = mix_mean,
      Rcpp::_["var"] = mix_var,
      Rcpp::_["a"] = mix_a,
      Rcpp::_["b"] = mix_b);
}

// Cholesky factor for a tridiagonal matrix with constant off-diagonal
void cholesky_tridiagonal(
    const arma::vec& omega_diag,
    const double omega_offdiag,
    arma::vec& chol_diag,
    arma::vec& chol_offdiag);

// Solves Chol*x = covector ("forward algorithm")
void forward_algorithm(
    const arma::vec& chol_diag,
    const arma::vec& chol_offdiag,
    const arma::vec& covector,
    arma::vec& htmp);

// Solves (Chol')*x = htmp ("backward algorithm")
void backward_algorithm(
    const arma::vec& chol_diag,
    const arma::vec& chol_offdiag,
    const arma::vec& htmp,
    arma::vec& h);

// Draws length(r) RVs, expects the non-normalized CDF mixprob
void inverse_transform_sampling(
    const arma::vec& mixprob,
    arma::ivec& r,
    const int T);

// Computes the CDF of the mixture indicators
void find_mixture_indicator_cdf(
    arma::vec& mixprob,
    const arma::vec& datanorm);

double h_log_posterior(
    const arma::vec& h,
    const arma::vec& y,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const double h0);

double h_aux_log_posterior(
    const arma::vec& h,
    const arma::vec& y_star,
    const arma::ivec& d,
    const double phi,
    const double rho,
    const double sigma2,
    const double mu,
    const double h0);

}

#endif  // H_UTILS_H
