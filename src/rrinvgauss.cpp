#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec rrinvgauss(int n, arma::vec mu, double lambda) {
  arma::vec random_vector(n);
  
  for(int i = 0; i < n; ++i) {
    double z = R::rnorm(0, 1);
    double y = z * z;
    double x = mu(i) + (mu(i) * mu(i) * y) / (2 * lambda) - (mu(i) / (2 * lambda)) * 
      sqrt(4 * mu(i) * lambda * y + mu(i) * mu(i) * y * y);
    double u = R::runif(0, 1);
    if(u <= mu(i) / (mu(i) + x)) {
      random_vector(i) = x;
    } else {
      random_vector(i) = mu(i) * mu(i) / x;
    }
  }
  return random_vector;
}