#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double rho_to_eta_cpp(double rho, double a) {
  return(log((rho + a) / (a - rho)));
}