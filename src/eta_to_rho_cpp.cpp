#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]] 
double eta_to_rho_cpp(double eta, double a) {
  return(a * (exp(eta) - 1) / (1 + exp(eta)));
}