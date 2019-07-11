// Define some proposal functions.

#ifndef proposals_h
#define proposals_h

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// #######################################################
// ###### The proposal functions
// #######################################################

arma::vec multiplierProposal_C(int size, arma::vec w_sd){
  // A proposal that scales with the absolute value of the parameter. (Always positive.)
  // Get this return factor and multiply by the current value for the proposal.
  // Return is a vector of length 'size', so it will work for multiple parameters.
  // The proposal ratio is the sum of the output vector.
  return exp( (randu(size) - 0.5) % w_sd );
}

arma::vec slideWindowLogSpace_C(arma::vec mu, arma::vec w_mu){
  // This will take the vector of root values mu and a vector of widths and make a draw.
  // Difference is that here we do the sliding window in log-space.
  // This proposal strategy is symmetrical in log-space but assymetrical in normal space.
  // Not a problem because the symmetry need to happen in the proposal step.
  // in R this is doing:
  // y <- sapply(1:length(mu), function(x) exp( log(mu) - (log(w_mu)/2) + (runif(1) * log(w_mu)) ) )
  return exp( log(mu) - (log(w_mu)/2) + (randu(mu.n_elem) % log(w_mu)) );
}

arma::vec slideWindow_C(arma::vec mu, arma::vec w_mu){
  // This will take the vector of root values mu and a vector of widths and make a draw.
  // in R this is doing:
  // y <- sapply(1:length(mu), function(x) runif(1, min = mu[i] - (w_mu[i]/2), max = mu[i] + (w_mu[i]/2) ) )
  return (mu - (w_mu/2)) + ( randu( mu.n_elem ) % w_mu );
}

// arma::vec slideWindowPositive_C(arma::vec sd, arma::vec w_sd){
//   // This will take the vector of root values mu and a vector of widths and make a draw.
//   arma::vec prop_sd = (sd - (w_sd/2)) + ( randu(sd.n_elem) % w_sd );
//   for(int i=0; i < sd.n_elem; i++){
//     prop_sd[i] = std::abs(prop_sd[i]);
//   }
//   return prop_sd;
// }

arma::mat makePropIWish_C(arma::mat vcv, double k, double v){
  // Function to make a proposal for the correlation structure.
  // Proposal for a single R matrix, so vcv is not a cube.
  // Defining the parameters here as double because of the multiplication.
  // v here is the size if the step. This will usually be a large number.
  arma::mat center = (v-k-1) * vcv;
  // Need the riwish function here.
  return riwish_C(v, center);
}

#endif
