#include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;
using namespace arma;

// Forward declarations
double eta_to_rho_cpp(double eta, double s);
arma::vec dmvn_cpp(arma::mat x,  
                   arma::rowvec mean,  
                   arma::mat sigma, 
                   bool logd = false);


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double mh_log_ratio_im_cpp(double rho_re, double eta_im, double eta_im_star, 
                           arma::rowvec mut, arma::mat ymat, arma::mat sig2I) {
  
  double a = sqrt(1 - pow(rho_re, 2));
  double rho_im = eta_to_rho_cpp(eta_im, a);
  double rho_im_star = eta_to_rho_cpp(eta_im_star, a);
  arma::mat rho_mat(2, 2);
  rho_mat(0, 0) = (1 + rho_re);
  rho_mat(0, 1) = rho_im;
  rho_mat(1, 0) = rho_im;
  rho_mat(1, 1) = (1 - rho_re);
  arma::mat rho_mat_star = rho_mat;
  rho_mat_star(0, 1) = rho_im_star;
  rho_mat_star(1, 0) = rho_im_star;
  arma::mat Sig = kron(rho_mat, sig2I);
  arma::mat Sig_star = kron(rho_mat_star, sig2I);
  arma::vec A_star = dmvn_cpp(ymat, mut, Sig_star, true) + 
    eta_im_star - 2 * log(1 + exp(eta_im_star));
  arma::vec A = dmvn_cpp(ymat, mut, Sig, true) + 
    eta_im - 2 * log(1 + exp(eta_im));
  return(as_scalar(A_star - A));
}