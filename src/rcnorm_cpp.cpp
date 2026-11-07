#include <RcppArmadillo.h>
using namespace Rcpp;

// Forward declarations
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);
arma::mat cplx2real_cov_cpp(arma::mat var, arma::cx_double rho);

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat rcnorm_cpp(int n, arma::vec mu, arma::mat var, arma::cx_double rho,
                     bool is_real = true, bool is_circular = false) {
  // n: number of samples
  // mu: mean vector
  // var: variance matrix
  // rho: complex correlation
  // is_real: flag for real output
  // is_circular: flag for circular complex normal distribution
  
  int n2 = 2 * n;
  arma::mat I = arma::eye(n2, n2);
  int p;
  arma::mat out;
  
  // Determine the dimension
  p = (var.n_elem == 1) ? 1 : var.n_cols;
  
  if (is_circular) {
    if (real(rho) == 0 && imag(rho) == 0) {
      // Special case for circular with zero correlation
      if (p == 1) {
        out = mvrnormArma(1, arma::zeros(n2), (var/2) * I).t();
      } else {
        out = mvrnormArma(n2, arma::zeros(p), var/2);
      }
      // Redundant return, but keeping original logic
      return out;
    }
  }
  
  // Generate samples for non-circular case
  arma::mat V = cplx2real_cov_cpp(var, rho);
  arma::vec mu_vec(2 * p);
  
  // Set mean vector
  mu_vec = (mu.n_elem == 1) ? arma::ones(2 * p) * mu : mu;
  
  out = mvrnormArma(n, mu_vec, V);
  
  // Return result based on is_real flag
  // Note: Both branches are identical in original code
  return join_cols(out.cols(0, p - 1), out.cols(p, 2 * p - 1));
}