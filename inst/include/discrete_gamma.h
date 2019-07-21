// Define the functions for the discrete-Gamma distribution.

#ifndef discrete_Gamma_h
#define discrete_Gamma_h

// Functions to implement a relaxed random walk using the discrete-Gamma distributed rates.

// Sets the first integration.
double f_int1(double x, void *params){
  double beta = *(double *) params;
  double den = x * gsl_ran_gamma_pdf(x, beta, 1.0/beta);
  return den;
}
  
// Sets the second integration.
double f_int2(double x, void *params){
  double beta = *(double *) params;
  double den = gsl_ran_gamma_pdf(x, beta, 1.0/beta);
  return den;
}

// Define the integration function when a and b are finite.
double integrateGammaAB(double a, double b, double beta) {
  // This function will return the mean value for a given category for a discrete-Gamma distribution.

  // Define the integration parameters (hard coded at the moment).
  double abserr = 0.0;
  double relerr = 1.e-7;
  double result_int1; // the integral value
  double error_int1; // the error estimate
  double result_int2; // the integral value
  double error_int2; // the error estimate
  size_t np = 1000; // work area size
  
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(np);

  gsl_function F_int1;  
  gsl_function F_int2;

  F_int1.function = &f_int1;
  F_int1.params = &beta;
  F_int2.function = &f_int2;
  F_int2.params = &beta;  
  // The function now only have a single parameter. So this is not necessary.
  // F.params = &n;

  gsl_integration_qag(&F_int1, a, b, abserr, relerr, np, GSL_INTEG_GAUSS41, w, &result_int1, &error_int1);
  gsl_integration_qag(&F_int2, a, b, abserr, relerr, np, GSL_INTEG_GAUSS41, w, &result_int2, &error_int2);

  double final_int = result_int1 / result_int2;

  return final_int;  
}

// Defines the integration function from a to Inf.
// This is a different case and need a special integrator.
double integrateGammaAInf(double a, double beta) {
  // This function will return the mean value for a given category for a discrete-Gamma distribution.

  // Define the integration parameters (hard coded at the moment).
  double abserr = 0.0;
  double relerr = 1.e-7;
  double result_int1; // the integral value
  double error_int1; // the error estimate
  double result_int2; // the integral value
  double error_int2; // the error estimate
  size_t np = 1000; // work area size
  
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(np);

  gsl_function F_int1;
  gsl_function F_int2;

  F_int1.function = &f_int1;
  F_int1.params = &beta;
  F_int2.function = &f_int2;
  F_int2.params = &beta;  
  // The function now only have a single parameter. So this is not necessary.
  // F.params = &n;

  // These integrations work from a to +Inf.
  gsl_integration_qagiu(&F_int1, a, abserr, relerr, np, w, &result_int1, &error_int1);
  gsl_integration_qagiu(&F_int2, a, abserr, relerr, np, w, &result_int2, &error_int2);

  double final_int = result_int1 / result_int2;

  return final_int;  
}

// Define function to generate the vector of rates.
arma::vec getGammaRates(double beta, arma::uword ncat) {

  // Get the quantiles to set the chunks:
  arma::vec quantiles = vec(ncat);
  // NOTE: The upper bound is ALWAYS +Inf! It is not shown in this vector.
  quantiles(0) = 0.0;
  
  for( arma::uword i=1; i < ncat; i++ ){
    // Need to cast the integers as doubles to get the correct output for the division below.
    quantiles(i) = gsl_cdf_gamma_Pinv((double)i/(double)ncat, beta, 1.0/beta);
  }
 
  // Find the mean of the rates for each of the chunks:
  arma::vec rates = vec(ncat);
  for( arma::uword j=0; j < ncat-1; j++){
    rates(j) = integrateGammaAB(quantiles(j), quantiles(j+1), beta);
  }

  // Last interval always is bounded by +Inf.
  // Need to use a different functions for this case.
  rates(ncat-1) = integrateGammaAInf(quantiles(ncat-1), beta);
  
  return rates;
  
}

#endif

// Local Variables:
// mode: c++
// End:
