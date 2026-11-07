#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);
arma::vec rrinvgauss_vec(int n, arma::vec mu, arma::vec lambda);
arma::vec sample_vec(const arma::vec& x, int size, bool replace, const arma::vec& prob);

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat beta_gdp_cpp(arma::mat XtX, arma::mat XtY, double sig2, 
                       arma::vec invtau) {
  arma::mat A = XtX + diagmat(invtau);
  arma::mat Ainv = inv_sympd(A);
  arma::vec mu_beta = Ainv * XtY;
  arma::mat Sig_beta = sig2 * Ainv;
  // arma::mat S = mvrnormArma(1, mu_beta, Sig_beta);
  // this is a 1 by p row vector
  return (mvrnormArma(1, mu_beta, Sig_beta));
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec invtau_gdp_cpp(arma::mat beta, double sig2, arma::vec lam) {
  arma::vec bb = arma::vectorise(beta);
  arma::vec lam2 = pow(lam, 2);
  return (rrinvgauss_vec(beta.n_elem, sqrt(lam2 * sig2) / abs(bb), lam2));
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double sig2_gdp_cpp(arma::mat YtY, arma::mat XtY, arma::mat XtX, arma::mat beta, arma::vec invtau,
                    int n, int p, double asig, double bsig) {
  double a_sig = (0.5) * (n + p) + asig;
  double b_sig;
  b_sig = arma::as_scalar(YtY - 2 * beta * XtY + beta * XtX * beta.t());
  b_sig = (0.5) * (b_sig + sum(invtau % pow(beta.t(), 2))) + bsig;
  return 1 / R::rgamma(a_sig, 1/b_sig);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double lamj_gdp_cpp(double betaj, double sig2, double alpha, double eta) {
  return (R::rgamma(alpha + 1, 1 / ((abs(betaj) / sqrt(sig2)) + eta)));
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec a_cond_log_cpp(arma::vec a, int p, arma::mat beta, double sig2, double eta) {
  arma::vec bb = arma::vectorise(beta);
  return(p * (log(1 - a) - log(a)) - 
         (1 / a) * sum(log(1 + (abs(bb) / (sqrt(sig2) * eta)))));
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec e_cond_log_cpp(arma::vec e, int p, arma::mat beta, double sig2, double alpha) {
  arma::vec bb = arma::vectorise(beta);
  arma::vec d = e / (1 - e);
  double sig = sqrt(sig2);
  arma::vec log_sum = zeros(d.n_elem);
  for(int i = 0; i < d.n_elem; ++i) {
    log_sum(i) = sum(log(1 + d(i)/sig * abs(bb)));
  }
  return(p * log(d) + (-alpha - 1) * log_sum);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List mcmc_gdp_cpp(int n_mcmc, int warm, int thin, 
                  arma::mat beta, arma::vec invtau, arma::vec lam, double sig2,
                  arma::vec& y, arma::mat& X,
                  double a_sig, double b_sig, double alpha, double eta, 
                  bool griddy,
                  arma::mat draws) {
  
  // #-------------------------------
  // # Here create some objects that will be used later.
  // #-------------------------------
  
  int n = X.n_rows;
  int p = X.n_cols;
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
    
    beta = beta_gdp_cpp(XtX, Xty, sig2, invtau);
    
    invtau = invtau_gdp_cpp(beta, sig2, lam);
    
    sig2 = sig2_gdp_cpp(yty, Xty, XtX, beta, invtau, n, p, a_sig, b_sig);
    
    for (int j = 0; j < p; ++j) {
      lam(j) = lamj_gdp_cpp(beta(0, j), sig2, alpha, eta);
    }
    
    if (griddy == true) {
      arma::vec aa = linspace(0.01, 0.99, 100);
      arma::vec ee = linspace(0.01, 0.99, 100);
      arma::vec weight_a = exp(a_cond_log_cpp(aa, p, beta, sig2, eta));
      arma::vec a = sample_vec(aa, 1, false, weight_a/sum(weight_a));
      alpha = (1 / a[0]) - 1;
      
      arma::vec weight_e = exp(e_cond_log_cpp(ee, p, beta, sig2, alpha));
      arma::vec e = sample_vec(ee, 1, false, weight_e/sum(weight_e));
      eta = (1 / e[0]) - 1;
    }
    
    
    // # -----------------------------------------------------
    // #  Save samples
    // # -----------------------------------------------------
    if (griddy == true) {
      if ((iter + 1) > warm) {
        if (fmod(iter + 1 , thin) == 0) {
          int d = (iter + 1 - warm) / thin;
          draws(d - 1, span(0, p - 1)) = beta;
          draws(d - 1, span(p, 2 * p - 1)) = trans(invtau);
          draws(d - 1, span(2 * p, 3 * p - 1)) = trans(lam);
          draws(d - 1, 3 * p) = sig2;
          draws(d - 1, draws.n_cols - 2) = alpha;
          draws(d - 1, draws.n_cols - 1) = eta;
        }
      }
    } else {
      if ((iter + 1) > warm) {
        if (fmod(iter + 1 , thin) == 0) {
          int d = (iter + 1 - warm) / thin;
          draws(d - 1, span(0, p - 1)) = beta;
          draws(d - 1, span(p, 2 * p - 1)) = trans(invtau);
          draws(d - 1, span(2 * p, 3 * p - 1)) = trans(lam);
          draws(d - 1, 3 * p) = sig2;
        }
      }
    }
    
    // # -----------------------------------------------------
    // #  Print iterations
    // # -----------------------------------------------------
    if ((iter + 1) % 100 == 0) {
      Rcout << "\r" << "GDP MCMC Iter: " << iter + 1 << std::flush;
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

