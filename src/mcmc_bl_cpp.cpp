#include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;
using namespace arma;

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);
arma::vec rrinvgauss(int n, arma::vec mu, double lambda);

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat beta_bl_cpp(arma::mat& XtX, arma::mat& XtY, double sig2, arma::vec& invtau2) {
  arma::mat A = XtX + diagmat(invtau2);
  arma::mat AA = (A + A.t()) / 2;
  int count = 0;
  while (AA.is_sympd() == 0) {
    count = count + 1;
    AA = AA + eye(A.n_cols, A.n_rows) * count * 0.001;
  }
  arma::mat Ainv = inv_sympd(AA);
  arma::mat Sig_beta = sig2 * Ainv;
  arma::vec mu_beta = Ainv * XtY;
  return mvrnormArma(1, mu_beta, Sig_beta);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double lam2_bl_cpp(arma::vec& invtau2, int p, double r, double delta) {
  return R::rgamma(p + r, 1 / ((0.5) * sum(1 / invtau2) + delta));
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec invtau2_bl_cpp(arma::mat& beta, double sig2, double lam2) {
  arma::vec bb = arma::vectorise(beta);
  return (rrinvgauss(beta.n_elem, sqrt(lam2 * sig2) / abs(bb), lam2));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double sig2_bl_cpp(arma::mat& YtY, arma::mat& XtY, arma::mat& XtX, 
                   arma::mat& beta, arma::vec& invtau2, 
                   int n, int p, double asig, double bsig) {
  double a_sig = (0.5) * (n + p) + asig;
  double b_sig;
  b_sig = arma::as_scalar(YtY - 2 * beta * XtY + beta * XtX * beta.t());
  b_sig = (0.5) * (b_sig + sum(invtau2 % pow(beta.t(), 2))) + bsig;
  return (1 / R::rgamma(a_sig, 1 / b_sig));
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List mcmc_bl_cpp(int n_mcmc, int warm, int thin, 
                 arma::mat beta, arma::vec tau2, double lam2, double sig2,
                 arma::vec& y, arma::mat& X, double r, double delta,
                 double a_sig, double b_sig,
                 arma::mat draws) {
  
  int n = X.n_rows;
  int p = X.n_cols;
  arma::mat XtX = X.t() * X;
  arma::mat Xty = X.t() * y.reshape(n, 1);
  arma::mat yty = y.t() * y.reshape(n, 1);
  arma::vec invtau2 = 1 / tau2;
  
  for (int iter = 0; iter < n_mcmc; ++iter) {
    beta = beta_bl_cpp(XtX, Xty, sig2, invtau2);
    invtau2 = invtau2_bl_cpp(beta, sig2, lam2);
    lam2 = lam2_bl_cpp(invtau2, p, r, delta);
    sig2 = sig2_bl_cpp(yty, Xty, XtX, beta, invtau2, n, p, a_sig, b_sig);
    
    if ((iter + 1) > warm) {
      if (fmod(iter + 1, thin) == 0) {
        int d = (iter + 1 - warm) / thin;
        draws(d - 1, span(0, p - 1)) = beta;
        draws(d - 1, span(p, 2 * p - 1)) = trans(1 / invtau2);
        draws(d - 1, draws.n_cols - 2) = lam2;
        draws(d - 1, draws.n_cols - 1) = sig2;
      }
    }
    
    if ((iter + 1) % 100 == 0) {
      Rcout << "\r" << "Blasso MCMC Iter: " << iter + 1 << std::flush;
    }
  }
  
  return List::create(
    _["draws"] = draws,
    _["warm"] = warm, 
    _["thin"] = thin, 
    _["nmcmc"] = n_mcmc
  );
}
