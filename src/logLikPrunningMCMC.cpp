#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

/* Function to calculate the likelihood of the model using the prunning algorithm.
*/

/* This code will depend on lists of matrices. To do that:

std::vector<MatrixXd> my_list;

This will define the vector of matrices based on the Eigen library.
But need to initialize with values.
*/

/*
COMMENTS OF IMPLEMENTATION
cache$ll is an incremental sum. So keep it as a double.

All these need to have dynamic lengths in the analysis:
cache$key is a vector that will grow in the analysis.
cache$X0 is a two dimension array, do not need to be a matrix.
cache$V0 is a vector of matrices.

Not using the 'cache' scheme here. These are independent objects.
And copying the lists over and over again might be part of the problem
   with the speed of the code anyways.
*/

// [[Rcpp::export]]
double logLikPrunningMCMC(NumericMatrix X, unsigned int k, IntegerVector nodes, IntegerVector anc, NumericMatrix mappedEdge, List R, NumericVector mu){
  
  // Declare objects.
  int nd_size = nodes.size();
  double X0[k][nd_size+1]; // Empty 2D array.
  std::vector<MatrixXd> V0; // Empty vector of matrices.
  NumericVector key(); // Empty numeric vector.
  double ll; //Empty double.

  // Traverse the tree:
  IntegerVector::iterator it; // The iterator for the loop.
  int ndend = nodes.end();
  for(it = nodes.begin(); it != 

  // return total / n;
}
