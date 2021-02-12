// Define some custom distributions and some helping functions.

#ifndef dist_help_h
#define dist_help_h

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// #######################################################
// ###### The distribution functions
// #######################################################
// These are some distribution functions that are not available on Rcpp or on RcppArmadillo.

int rMultinom(arma::vec p) {
  // This is a function to make a single draw from a multinominal distribution.
  // p is a vector of probabilities and need to sum to 1.
  p = p / sum(p);

  if( any(p < 0.0) ){
    // If probabilities < 0 then bounce to 0.
    for( arma::uword w = 0; w < p.n_rows; w++ ){
      if( p[w] < 0.0 ){
	p[w] = 0.0;
      }
    }
  }
  
  double unif_draw = as_scalar(randu(1));

  // Use this random draw to select one of the outcomes. Note that this is a simple map computation, the continuous draw is mapped to an integer depending on the value.

  arma::uword i = 0;
  double map_ref = p[i];
  while( unif_draw >= map_ref ) {
    // The loop will not break because p sums to 1.
    i++;
    map_ref = map_ref + p[i];
  }

  return i;  
}

double logDensityIWish_C(arma::mat W, double v, arma::mat S){
  // Function for the density of a inverse Wishart distribution.
  // W the covariance matrix.
  // v the degrees of freedom parameter.
  // S the standard matrix for the distribution.
  double valS;
  double valW;
  double sign_sink;
  double lgammapart = 0;
    
  double k = S.n_cols;
  for(arma::uword i=0; i < S.n_cols; i++) {
    lgammapart = lgammapart + lgamma((v-i)/2);
  }
  log_det(valS, sign_sink, S);
  log_det(valW, sign_sink, W);
  
  double ldenom = lgammapart + ( ( (v*k)/2.0 ) * log( 2.0 ) ) + ( ( (k*(k-1.0))/4.0 ) * log( arma::datum::pi ) );
  // Need to make sure that we are doing 'double' operations here!
  double lnum = ( ( v/2.0 ) * valS ) + ( ( -(v + k + 1.0)/2.0 ) * valW ) + ( -0.5 * trace( S * inv(W) ) );
  return lnum - ldenom;
}

arma::mat riwish_C(int v, arma::mat S){
  // Generates a random draw from a inverse-Wishart distribution.
  // arma::mat CC = chol( inv(S) );
  arma::mat CC = chol( inv(S) );
  int p = S.n_cols;
  // Make a diagonal matrix with the elements:
  // R::rchisq( df ) // with df in a sequence v:(v - p + 1)
  arma::vec chi_sample(p);
  // regspace will make a sequence.
  arma::vec df = regspace(v, v-p+1); // The degrees of freedom when sampling.
  for( int i=0; i < p; i++ ) {
    chi_sample[i] = sqrt( R::rchisq( df[i] ) );
  }
  arma::mat Z = diagmat( chi_sample );
  // Need to fill the upper triangular elements with samples from a standard normal.
  for( int i=0; i < p-1; i++) {
    // randn uses a normal distribution to sample.
    Z(i,span((i+1), (p-1))) = trans(randn(p-(i+1)));
  }
  arma::mat out = Z * CC;
  return inv( trans(out) * out );
}

double hastingsDensity_C(arma::cube R, arma::cube R_prop, int k, arma::vec v, int Rp){
  // The hasting is only computed for the regime that is updated (Rp).
  arma::mat center_curr = (v[Rp]-k-1) * R.slice(Rp);
  arma::mat center_prop = (v[Rp]-k-1) * R_prop.slice(Rp);
  return logDensityIWish_C(R.slice(Rp), v[Rp], center_prop) - logDensityIWish_C(R_prop.slice(Rp), v[Rp], center_curr);
}  

arma::vec extractQ(arma::mat Q, arma::uword size, std::string model_Q){
  // Function to extract a column vector from the Q matrix.
  // Length of the vector will depend on the type of the model for the Q matrix.
  // Need to use the same pattern to extract and rebuild the matrix.
  arma::vec vec_Q;
  
  if( model_Q == "ER" ){
    vec_Q = vec(1, fill::zeros);
    vec_Q[0] = Q(0,1); // All off-diagonals are the same.
  } else if( model_Q == "SYM" ){
    arma::uword count = 0;
    int size_vec = ( ( size * size ) - size ) / 2;
    vec_Q = vec(size_vec, fill::zeros);
    for( arma::uword i=0; i < size; i++ ){
      for( arma::uword j=0; j < size; j++ ){
	if( i >= j ) continue;
	vec_Q[count] = Q(i,j);
	count++;
      }
    }
  } else{ // model_Q == "ARD"
    arma::uword count = 0;
    int size_vec = ( size * size ) - size;
    vec_Q = vec(size_vec, fill::zeros);
    for( arma::uword i=0; i < size; i++ ){
      for( arma::uword j=0; j < size; j++ ){
	if( i == j ) continue;
	vec_Q[count] = Q(i,j);
	count++;
      }
    }
  }

  return vec_Q;
}

arma::mat cov2cor_C(arma::mat V){
  // This is a **brute force** function for the correlation matrix.
  arma::mat Vdiag = inv( sqrt( diagmat(V) ) );
  return Vdiag * V * Vdiag;
}

#endif
