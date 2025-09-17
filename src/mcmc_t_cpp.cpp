#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat beta_t_cpp(arma::mat XtX, arma::mat XtY, double sig2,
                     arma::vec tau2) {
  arma::mat A = XtX + diagmat(1/tau2);
  arma::mat Ainv = inv_sympd(A);
  arma::mat Sig_beta = sig2 * Ainv;
  arma::vec mu_beta = vectorise(Ainv * XtY);
  return mvrnormArma(1, mu_beta, Sig_beta);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec tau2_t_cpp(arma::mat beta, double atau, double btau, double sig2, int p) {
  arma::vec bb = arma::vectorise(beta);
  double shape = atau + 0.5;
  arma::vec rate = btau + 0.5 * (pow(bb, 2) / sig2);
  arma::vec tau2_vec(p);
  
  for(int i = 0; i < p; ++i) {
    tau2_vec(i) = 1 / R::rgamma(shape, 1/rate(i));
  }
  return(tau2_vec);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double sig2_t_cpp(arma::vec YtY, arma::mat XtY, arma::mat XtX, arma::mat beta, 
                  arma::vec tau2, int n, int p, double asig, double bsig) {
  double a_sig = (0.5) * (n + p) + asig;
  double b_sig;
  b_sig = arma::as_scalar(YtY - 2 * beta * XtY + beta * XtX * beta.t());
  b_sig = (0.5) * (b_sig + sum(pow(beta.t(), 2)/tau2)) + bsig;
  return 1 / R::rgamma(a_sig, 1/b_sig);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List mcmc_t_cpp(int n_mcmc, int warm, int thin, 
                arma::mat beta, arma::vec tau2, double sig2,
                arma::vec& y, arma::mat& X, double a_tau, double b_tau, 
                double a_sig, double b_sig,
                arma::mat draws) {
  
  // #-------------------------------
  // # Here create some objects that will be used later.
  // #-------------------------------
  
  const int n = X.n_rows;
  const int p = X.n_cols;
  arma::mat XtX = X.t() * X;
  arma::mat Xty = X.t() * y.reshape(n, 1);
  arma::mat yty = y.t() * y.reshape(n, 1);
  
  
  // #-------------------------------
  // # Gibbs sampling algorithm
  // #------------------------------- 
  
  for (int iter = 0; iter < n_mcmc; ++iter) {
    
    // # -----------------------------------------------------
    // #  Gibbs sampler functions
    // # -----------------------------------------------------
    
    beta = beta_t_cpp(XtX, Xty, sig2, tau2);
    
    tau2 = tau2_t_cpp(beta, a_tau, b_tau, sig2, p);
    
    sig2 = sig2_t_cpp(yty, Xty, XtX, beta, tau2, n, p, a_sig, b_sig);
    
    // # -----------------------------------------------------
    // #  Save samples
    // # -----------------------------------------------------
    if ((iter + 1) > warm) {
      if (fmod(iter + 1 , thin) == 0) {
        int d = (iter + 1 - warm) / thin;
        draws(d - 1, span(0, p - 1)) = beta;
        draws(d - 1, span(p, 2 * p - 1)) = tau2.t();
        draws(d - 1, draws.n_cols - 1) = sig2;
      }
    }
    
    // # -----------------------------------------------------
    // #  Print iterations
    // # -----------------------------------------------------
    if ((iter + 1) % 100 == 0) {
      Rcout << "\r" << "Student t MCMC Iter: " << iter + 1 << std::flush;
    }
  }
  
  // #--------------
  // # Write output
  // #--------------
  // Object draws does not have colnames. Need to assign.
  return List::create(
    _["draws"] = draws,
    _["warm"] = warm, 
    _["thin"] = thin, 
    _["nmcmc"] = n_mcmc
  // _["n"] = n,
  // _["p"] = p,
  // _["XtX"] = XtX,
  // _["Xty"] = Xty,
  // _["yty"] = yty
  );
}

