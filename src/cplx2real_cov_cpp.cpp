#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat cplx2real_cov_cpp(arma::mat var, arma::cx_double rho) {
  arma::cx_mat rel = var * rho;
  arma::mat Vxx = arma::real(var + rel);
  arma::mat Vyx = arma::imag(var + rel);
  arma::mat Vxy = arma::imag(-var + rel);
  arma::mat Vyy = arma::real(var - rel);
  
  return arma::join_rows(arma::join_cols(Vxx, Vyx), arma::join_cols(Vxy, Vyy)) / 2;
}