#include <RcppArmadillo.h>
using namespace Rcpp;

// Forward declaration
arma::vec rrinvgauss(int n, arma::vec mu, double lambda);

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double invtau2j_bl_cpp(arma::mat beta, int j, double sig2, double lambda2) {
  arma::vec bb = arma::vectorise(beta);
  
  // Create a single-element arma::vec for mu
  arma::vec mu = arma::vec{sqrt(lambda2 * sig2) / std::abs(bb[j - 1])};
  
  // Use the single-element vector in rrinvgauss
  return arma::as_scalar(rrinvgauss(1, mu, lambda2));
}