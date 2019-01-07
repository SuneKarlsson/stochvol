#ifndef _AUXMIX_H_
#define _AUXMIX_H_

#include <Rcpp.h>

// Some constants relating to the approximation of log(chisq) trough
// normal mixture (from Omori et al., 2007), and
// corresponding functions related to sampling the indicators

const double mix_prob[10] = {.00609, .04775, .13057, .20674, .22715, .18842, .12047, .05591, .01575, .00115};

const double mix_mean[10] = {1.92677, 1.34744, .73504, .02266, -.85173, -1.97278, -3.46788, -5.55246, -8.68384, -14.65000};

const double mix_var[10] = {.11265, .17788, .26768, .40611, .62699, .98583, 1.57469, 2.54498, 4.16591, 7.33342};

const double mix_a[10] = {1.01418, 1.02248, 1.03403, 1.05207, 1.08153, 1.13114, 1.21754, 1.37454, 1.68327, 2.50097};

const double mix_b[10] = {0.50710, 0.51124, 0.51701, 0.52604, 0.54076, 0.56557, 0.60877, 0.68728, 0.84163, 1.25049};

const double mix_varinv[10] = {
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

const double mix_2varinv[10] = {
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

const double mix_pre[10] = {
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

// Non-normalized posterior probabilities
void findMixprobs(double * mixprob, const Rcpp::NumericVector & datanorm);

// Cumulative sum over columns of a matrix
void colCumsums(double * x, int const nrow, int const ncol);

// Combines findMixprobs() and colCumsums() (see above) into one function
void findMixCDF(double * mixprob, const Rcpp::NumericVector & datanorm);

#endif
