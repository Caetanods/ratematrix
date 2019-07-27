// Define the prior functions.

#ifndef priors_h
#define priors_h

// #######################################################
// ###### The prior functions
// #######################################################

double priorRoot_C(arma::vec mu, arma::mat par_prior_mu, std::string den_mu){
  // Here the 'par_prior_mu' need to be a matrix constructed before.
  // The value for 'den_mu' will also be given before.
  // These objects will be constructed by 'makePrior' function.
  // If you have a problem with the vector type you can use:
  // NumericVector(a.begin(),a.end()) to transform into a Rcpp vector.
  double pp = 0.0;
  if( den_mu == "unif" ){
    for( arma::uword i=0; i < mu.n_elem; i++ ){
      pp = pp + R::dunif(mu[i], par_prior_mu(i,0), par_prior_mu(i,1), true);
    }
  } else{
    for( arma::uword i=0; i < mu.n_elem; i++ ){
      pp = pp + R::dnorm(mu[i], par_prior_mu(i,0), par_prior_mu(i,1), true);
    }
  }  
  return pp;
}

double priorSD_C(arma::mat sd, arma::mat par_prior_sd, std::string den_sd, int p){
  // FUNCTION FOR 2 OR MORE RATE REGIMES.
  // 'sd' is a matrix with number of columns equal to 'p', number of regimes.
  // 'par_prior_sd' is a matrix with parameters, each row for each regime.
  // p is the number of regimes in the model.

  // Here we are checking the number of regimes instead of the dimensions of the matrix.

  double pp = 0.0;
  if( den_sd == "unif" ){
    for( arma::uword i=0; i < sd.n_rows; i++ ){
      for( arma::uword j=0; j < sd.n_cols; j++){
	// Each line of the 'par_prior_sd' correspond to each of the sd vectors stored as the columns of the 'sd' matrix.
	pp = pp + R::dunif(sd(i,j), par_prior_sd(j,0), par_prior_sd(j,1), true);
      }
    }
  } else{
    for( arma::uword i=0; i < sd.n_rows; i++ ){
      for( arma::uword j=0; j < sd.n_cols; j++){
	// Each line of the 'par_prior_sd' correspond to each of the sd vectors stored as the columns of the 'sd' matrix.
	pp = pp + R::dlnorm(sd(i,j), par_prior_sd(j,0), par_prior_sd(j,1), true);
      }
    }
  }  
  return pp;
}

double priorSD_vec(arma::rowvec sd, arma::rowvec par_prior_sd, std::string den_sd, int p){
  // SAME FUNCTION AS priorSD_C BUT WORKS FOR THE CASE OF A SINGLE REGIME.
  // 'sd' is a matrix with number of columns equal to 'p', number of regimes.
  // 'par_prior_sd' is a matrix with parameters, each row for each regime.
  // p is the number of regimes in the model.

  double pp = 0.0;
  if( den_sd == "unif" ){
    for( arma::uword i=0; i < sd.n_cols; i++ ){
	// Each line of the 'par_prior_sd' correspond to each of the sd vectors stored as the columns of the 'sd' matrix.
	pp = pp + R::dunif(sd(i), par_prior_sd(0), par_prior_sd(1), true);
    }
  } else{
    for( arma::uword i=0; i < sd.n_cols; i++ ){
      // Each line of the 'par_prior_sd' correspond to each of the sd vectors stored as the columns of the 'sd' matrix.
      pp = pp + R::dlnorm(sd(i), par_prior_sd(0), par_prior_sd(1), true);
    }
  }
  return pp;
  
}

double priorCorr_C(arma::cube corr, arma::vec nu, arma::cube sigma){
  // FUNCTION FOR 2 OR MORE RATE REGIMES.
  // nu is a vector with the correspondent degree of freedom for the rate regime.
  // This is the more complicated one. Need to deal with the Wishart distributions.
  // If the correlation prior was set to "uniform'. Then we just need to set sigma and v to the standard values when doing the 'makePrior' step. No need for a if test here.
  // Will treat the correlation and sigma as arrays (cube). This will work even if they are matrices. I think.
  arma::uword p = corr.n_slices;
  double pp = 0.0;
  for( arma::uword i=0; i < p; i++ ) {
    pp = pp + logDensityIWish_C(corr.slice(i), nu[i], sigma.slice(i)); // Need to define this one.
  }
  return pp;
}

double priorQ(arma::vec vec_Q, arma::vec par_prior_Q, std::string den_Q){
  // Compute the prior probability for the Q matrix rates.
  // Can be uniform or exponential prior.
  // The 'par_prior_Q' can be a vector with 1 or 2 elements, depeding on the density.
  double pp = 0.0;
  if( den_Q == "uniform" ){
    for( arma::uword i=0; i < vec_Q.n_rows; i++ ){ // vec_Q is a column.
      // If "unif" then 'par_prior_Q' is a vector with 2 elements.
      pp = pp + R::dunif(vec_Q[i], par_prior_Q[0], par_prior_Q[1], true);
    }
  } else{ // Then it is an exponential prior.
    for( arma::uword i=0; i < vec_Q.n_rows; i++ ){
      // In Rcpp the exponential is '1/rate' whereas in R it is simply 'rate' .
      // Here only one parameter is needed.
      pp = pp + R::dexp(vec_Q[i], 1/par_prior_Q[0], true);
    }
  }
  
  return pp;
}

#endif

// Local Variables:
// mode: c++
// End:
