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

// Define the integration functions.
double integrateGamma(double a, double b, double beta) {
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

  gsl_integration_qag(&F_int1, a, b, abserr, relerr, np, GSL_INTEG_GAUSS15, w, &result_int1, &error_int1);
  gsl_integration_qag(&F_int2, a, b, abserr, relerr, np, GSL_INTEG_GAUSS15, w, &result_int2, &error_int2);

  double final_int = result_int1 / result_int2;

  return final_int;  
}

#endif
