// Define multiple functions to log posterior samples to files.

#ifndef log_fn_h
#define log_fn_h

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

void writeToMultFile_C(std::ostream& mcmc_stream, arma::uword p, arma::uword k, arma::cube R, arma::vec mu){
  // Note the 'std::ostream&' argument here is the use of a reference.
  for( arma::uword i=0; i < p; i++ ){
    for( arma::uword j=0; j < k; j++ ){
      for( arma::uword z=0; z < k; z++ ){
	mcmc_stream << R.slice(i)(j,z);
	mcmc_stream << "; ";
      }
    }
  }
  for( arma::uword i=0; i < k-1; i++ ){
    mcmc_stream << mu[i];
    mcmc_stream << "; ";
  }
  mcmc_stream << mu.tail(1);
  // mcmc_stream << "\n";
}

void writeQToFile(std::ostream& Q_mcmc_stream, arma::vec vec_Q, arma::uword k, std::string model_Q){
  // Note the 'std::ostream&' argument here is the use of a reference.
  if( model_Q == "ER" ){
    Q_mcmc_stream << vec_Q;
  } else{
    arma::uword print_size = vec_Q.n_rows;
    for( arma::uword i=0; i < (print_size-1); i++ ){
      Q_mcmc_stream << vec_Q[i];
      Q_mcmc_stream << "; ";
    }
    Q_mcmc_stream << vec_Q[print_size-1];
    Q_mcmc_stream << "\n";
  }
}

#endif
