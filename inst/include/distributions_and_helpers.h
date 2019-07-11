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
  // arma::mat CC = chol( inv_sympd(S) );
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

arma::mat buildQ(arma::vec vec_Q, arma::uword size, std::string model_Q){
  // Function to re-build the Q matrix.
  // Need to follow the same pattern used to extract the vector.
  arma::mat Q = mat(size, size, fill::zeros);
  
  if( model_Q == "ER" ){
    Q.fill(vec_Q[0]);
    // Now fill the diagonal.
    for( arma::uword i=0; i < size; i++ ){
      Q(i,i) = -1.0 * ( sum( Q.row(i) ) - vec_Q[0] );
    }
  } else if( model_Q == "SYM" ){
    Q.fill(0); // Fill the matrix with 0.
    arma::uword count = 0;
    // Go over the matrix and fill the upper and lower-tri.
    for( arma::uword i=0; i < size; i++ ){
      for( arma::uword j=0; j < size; j++ ){
	if( i >= j ) continue;
	Q(i,j) = vec_Q[count];
	Q(j,i) = vec_Q[count]; // The trick to fill the lower-tri.
	count++;
      }
    }
    // Now fill the diagonal.
    for( arma::uword i=0; i < size; i++ ){
      Q(i,i) = -1.0 * sum( Q.row(i) );
    }
  } else{ // model_Q == "ARD"
    Q.fill(0); // Fill the matrix with 0.
    arma::uword count = 0;
    // Go over the matrix and fill the upper and lower-tri.
    for( arma::uword i=0; i < size; i++ ){
      for( arma::uword j=0; j < size; j++ ){
	if( i == j ) continue;
	Q(i,j) = vec_Q[count];
	count++;
      }
    }
    // Now fill the diagonal.
    for( arma::uword i=0; i < size; i++ ){
      Q(i,i) = -1.0 * sum( Q.row(i) );
    }
  }

  return Q;
}

arma::mat cov2cor_C(arma::mat V){
  // This is a **brute force** function for the correlation matrix.
  arma::mat Vdiag = inv( sqrt( diagmat(V) ) );
  return Vdiag * V * Vdiag;
}

// Some functions to work with polytopes:

arma::mat samplePolytope(arma::mat edge) {
  // Function to sample a single point from the polytope given the edges with min and max for the traits.
  // This assumes that the space is a n-dimensional cube, with straight lines making the limits of the space.
  // We sample from the volume delimited by the max and min values for each of the dimentions.

  arma::uword dim_poly = edge.n_cols / 2; // number of dimensions.
  arma::uword n_samples = edge.n_rows; // number of samples to be taken.
  // Retuning container need to have samples as rows and traits as columns (same format as other functions).
  arma::mat samples = mat(n_samples, dim_poly);
  arma::vec rand_sample = vec(dim_poly);
  arma::vec range_axis = vec(dim_poly);
  arma::vec min_axis = vec(dim_poly);
  int axis_count;

  for(uword i=0; i < n_samples; i++) {
    // Need to sample from a uniform distribution for each of the dimensions.
    rand_sample = randu(dim_poly);
    axis_count = 0;

    for(uword j=0; j < edge.n_cols; ) { // No advancement statement here.
      // Create the vector with the range.
      range_axis(axis_count) = edge(i,j+1) - edge(i,j);
      min_axis(axis_count) = edge(i,j);
      // Make the advancement:
      j += 2;
      axis_count++;
    }

    // Here need to transpose the col vector to a row vector.
    if( all(range_axis) ){
      // If all elements of the range are non-zero:
    samples.row(i) = trans( min_axis + ( rand_sample % range_axis ) );
    } else{
      // At least one of the elements is zero. Do not compute the difference.
      // Keep the min values, since the max and min are equal.
      samples.row(i) = trans( min_axis );
    }
  }

  arma::mat samples_trans = trans(samples);
  return samples_trans;
}

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  // Define a function to sample from a multivariate normal.
  // Source: http://gallery.rcpp.org/articles/simulate-multivariate-normal/
  // Author: Ahmadou Dicko - Mar 12, 2013
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

arma::vec Gibbs_get_trio(arma::vec a, arma::vec d1, arma::vec d2, arma::vec s, arma::vec t1, arma::vec t2, arma::cube R, arma::uword p, int is_root){
  // Sample the value for a node given the descendant values and the ancestor value.
  // This follows algorithm described in Quintero et al. (2015).
  // The branch lengths s, t1, t2 are just (transposed) rows from the mapped.edge matrix.
  if( p == 1){
    // Single regime only.
    double muA;
    arma::vec muB;
    arma::vec mu; // Final mean of the multivariate normal.
    double scaler;
    arma::mat V;
    arma::mat node_sample; // But this should be a column!
    if( is_root == 1 ){
      mu = ( (d1*t2[0])+(d2*t1[0]) ) / (t1[0]+t2[0]);
      scaler = (t1[0] * t2[0]) / (t1[0] + t2[0]);
      V = scaler * R.slice(0);      
    } else{
      // Compute the mean:
      muA = 1.0 / ( (t1[0]*t2[0]) + ((t1[0]+t2[0])*s[0]) );
      muB = (t1[0]*t2[0]*a) + (s[0]*t2[0]*d1) + (s[0]*t1[0]*d2);
      mu = muA * muB;
      // Compute the variance:
      scaler = (t1[0]*t2[0]*s[0]) / ( (t1[0]*t2[0]) + ( (t1[0]+t2[0])*s[0] ) );
      V = R.slice(0) * scaler;
    }
    // We can draw from the posterior distribution:
    node_sample = mvrnormArma(1, mu, V);
    arma::vec output_fn = trans( node_sample.row(0) );
    return trans( node_sample.row(0) ); // Return need to be a row vector.
  } else{
    // Make a draw taking into account multiple regimes.
    // In this case we need to use similar steps to the pruning algorithm.
    // When computing the variance we need to scale each of the rate matrices by their respective regimes.
    // DUMMY RETURN. THIS ELSE BRACHET DOES NOT WORK IN THE CURRENT VERSION.
    // Using object of correct format and scope for the dummy return.
    arma::vec dummy;
    return dummy; // Return need to be a row vector.
  }
}


#endif
