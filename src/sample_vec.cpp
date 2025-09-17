#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec sample_vec(const arma::vec& x, int size, bool replace, 
                     const arma::vec& prob) {
  return Rcpp::RcppArmadillo::sample(x, size, replace, prob);
}