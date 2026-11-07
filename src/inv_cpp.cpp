#include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;
using namespace arma;

 // [[Rcpp::depends(RcppArmadillo)]]
 // [[Rcpp::export]]
 arma::mat inv_cpp(arma::mat m1) {
   // Check if the matrix is square
   if (m1.n_rows != m1.n_cols) {
     stop("The input matrix must be square.");
   }
   
   // Try to compute the inverse using inv_sympd
   try {
     return inv_sympd(m1);
   } catch (std::runtime_error &e) {
     stop("Matrix is not symmetric positive definite or not invertible.");
   }
 }
 