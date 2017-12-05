#include <Rcpp.h>

using namespace Rcpp;
//' Multiply a number by three
//'
//' @param x A single integer.
//' @export
// [[Rcpp::export]]
int timesThree(int x){
   return x * 3;
}
