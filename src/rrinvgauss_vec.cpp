#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec rrinvgauss_vec(int n, arma::vec mu, arma::vec lambda){
  arma::vec random_vector(n);
  double z, y, x, u;
  
  for(int i = 0; i < n; ++i) {
    z = R::rnorm(0, 1);
    y = z * z;
    x = mu(i) + (mu(i) * mu(i) * y) / (2 * lambda(i)) - (mu(i) / (2 * lambda(i))) * 
      sqrt(4 * mu(i) * lambda(i) * y + mu(i) * mu(i) * y * y);
    u = R::runif(0, 1);
    if(u <= mu(i) / (mu(i) + x)) {
      random_vector(i) = x;
    }else{
      random_vector(i) = mu(i) * mu(i) / x;
    };
  }
  return(random_vector);
}