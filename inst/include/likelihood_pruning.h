// Function to compute the likelihood of the mvBM model.

#ifndef likelihood_pruning_h
#define likelihood_pruning_h

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

double logLikNode_C(arma::vec ss, arma::mat sigma_len, arma::mat sigma_len_inv, int k) {
  double val;
  double signal;
  arma::log_det(val, signal, sigma_len); // val is logdeterminant and signal should be 0.
  return -0.5 * ( k * log(2 * arma::datum::pi) + val + as_scalar(trans(ss) * sigma_len_inv * ss));
}

#endif
