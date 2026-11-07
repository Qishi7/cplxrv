#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

arma::sp_mat cplx2real_cov_inv_cpp(arma::mat var, arma::cx_double rho);
double eta_to_rho_cpp(double eta, double a);
double rho_to_eta_cpp(double rho, double a);
arma::uvec keep_elem(const arma::vec& x, const arma::ivec& to_remove);
arma::vec rrinvgauss(int n, arma::vec mu, double lambda);
arma::vec dmvn_cpp(arma::mat x, arma::rowvec mean, arma::mat sigma, bool logd = false);
double get_tune_cpp(double tune, double keep, int iter, double step, double exp_term, double target);
double mh_log_ratio_re_cpp(double eta_re, double eta_re_star, double rho_im, arma::rowvec mut, arma::mat ymat, arma::mat sig2I);
double mh_log_ratio_im_cpp(double rho_re, double eta_im, double eta_im_star, arma::rowvec mut, arma::mat ymat, arma::mat sig2I);

// [[Rcpp::depends("RcppArmadillo")]]
arma::sp_mat cplx2real_cov_inv_cpp(arma::mat var, arma::cx_double rho) {
  double rho_re = real(rho);
  double rho_im = imag(rho);
  double dd = 1/ ((1 - rho_re) - pow(rho_im, 2) / (1 + rho_re));
  double aa = 1 + rho_re;
  arma::mat rho_mat = mat(2, 2, fill::ones);
  rho_mat(0, 0) = (aa + pow(rho_im, 2)) / pow(aa, 2);
  rho_mat(0, 1) = -rho_im / aa;
  rho_mat(1, 0) = -rho_im / aa;
  
  arma::sp_mat Vinv = sp_mat(kron(rho_mat, dd * var));
  
  return(Vinv);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mvrnormArma1(arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(1, ncols);
  return mu.t() + Y * arma::chol(sigma);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat mvrnormArma_inv(arma::vec mu, arma::mat sigma_inv) {
  int ncols = sigma_inv.n_cols;
  arma::mat Y = arma::randn(1, ncols);
  arma::mat L = arma::chol(sigma_inv, "lower");
  arma::mat B = arma::solve(arma::trimatl(L), Y.t()); 
  arma::mat C = arma::solve(arma::trimatu(L.t()), B);
  return (mu + C).t();
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat beta_cbl_cpp(arma::vec y, arma::mat X, arma::mat beta,
                       List Sigj_inv_lst, arma::sp_mat& Vinv, double sig2, int p,
                       arma::vec all_ind, List select_idx_lst, List Xj_lst) {
  arma::ivec remove_idx;
  arma::mat Sigj;
  arma::uvec idx;
  arma::mat Aj;
  arma::vec mu;
  
  for (int j = 0; j < p; ++j) {
    arma::mat Sigj_inv = Sigj_inv_lst[j];
    Sigj = inv_sympd(Sigj_inv);
    remove_idx = {j, j + p};
    idx = keep_elem(all_ind, remove_idx);
    Aj = y - X.cols(idx) * beta.cols(idx).t();
    arma::mat Xj = Xj_lst[j];
    mu = Sigj * Xj.t() * Vinv * Aj;
    beta.cols(as<uvec>(select_idx_lst[j])) = mvrnormArma1(mu, sig2 * Sigj);
  }
  
  return(beta);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec ga_cbl_cpp(double lam2, double sig2, arma::vec btb, int p) {
  return(rrinvgauss(p, sqrt(lam2 * sig2 / btb), lam2));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double sig2_cbl_cpp(double asig, double bsig, int n, int p,
                    arma::vec y, arma::vec mu, arma::vec tau2, arma::vec btb,
                    arma::sp_mat Vinv) {
  arma::vec e = y - mu;
  double b_sig = as_scalar(bsig + 0.5 * (e.t() * Vinv * e + sum(btb / tau2)));
  return (1 / R::rgamma(n + p + asig, 1/b_sig));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double lam2_cbl_cpp(arma::vec tau2, double r, double delta, int p) {
  return(R::rgamma(3 * p / 2 + r, 1/(delta + (sum(tau2) / 2))));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double mh_log_re_cpp(double eta_re, double rho_im, 
                     arma::rowvec mut, arma::mat ymat, arma::mat sig2I) {
  double rho_re = eta_to_rho_cpp(eta_re, 1);
  arma::mat rho_mat(2, 2);
  rho_mat(0, 0) = (1 + rho_re);
  rho_mat(0, 1) = rho_im;
  rho_mat(1, 0) = rho_im;
  rho_mat(1, 1) = (1 - rho_re);
  
  arma::mat Sig = kron(rho_mat, sig2I);
  
  arma::vec A = dmvn_cpp(ymat, mut, Sig, true) + 
    eta_re - 2 * log(1 + exp(eta_re));
  return(as_scalar(A));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double mh_log_im_cpp(double rho_re, double eta_im,
                     arma::rowvec mut, arma::mat ymat, arma::mat sig2I) {
  double a = sqrt(1 - pow(rho_re, 2));
  double rho_im = eta_to_rho_cpp(eta_im, a);
  
  arma::mat rho_mat(2, 2);
  rho_mat(0, 0) = (1 + rho_re);
  rho_mat(0, 1) = rho_im;
  rho_mat(1, 0) = rho_im;
  rho_mat(1, 1) = (1 - rho_re);
  
  arma::mat Sig = kron(rho_mat, sig2I);
  
  arma::vec A = dmvn_cpp(ymat, mut, Sig, true) +
    eta_im - 2 * log(1 + exp(eta_im));
  return(as_scalar(A));
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List mcmc_cbl_cpp(const int& n_mcmc, const int& warm, const int& thin, 
                  arma::mat beta, arma::vec tau2, double lam2, double sig2,
                  arma::vec y, arma::mat& X, arma::cx_double rho_b, 
                  arma::cx_double rho_ep, const double& r, const double& delta, 
                  double a_sig, double b_sig,
                  const bool& learn_rho, const int& n_chains,       
                  double keep_re, double keep_im,
                  double tune_re, double tune_im,
                  int& tune_len, bool& adapt, const double& adapt_step, 
                  const double& target_accept_rate,
                  arma::mat draws) {
  int Tb = tune_len;
  double keep_tmp_re = keep_re;
  double keep_tmp_im = keep_im;
  double rho_re, rho_im, eta_re, eta_im;
  double dd, aa;
  double eta_re_star, rho_re_star, logpost_re, eta_im_star, logpost_im;
  const double eps = 0.0001;
  int n = X.n_rows / 2;
  int p = X.n_cols / 2;
  
  arma::mat I = arma::eye(n, n);
  arma::sp_mat Vinv = cplx2real_cov_inv_cpp(I, rho_ep);
  
  arma::vec all_ind = linspace(0, 2 * p - 1, 2 * p);
  arma::vec ga;
  arma::vec mu;
  arma::ivec remove_idx;
  arma::uvec select_idx;
  arma::mat select_beta;
  arma::vec btVbinvb(p);
  arma::mat Vbinv = eye(2, 2);
  int j;
  
  List select_idx_lst(p), Xj_lst(p), Sigj_inv_lst(p);
  
  for (j = 0; j < p; ++j) {
    select_idx = {(unsigned int) j, (unsigned int) j + p};
    select_idx_lst[j] = select_idx;
    arma::mat Xj = X.cols(select_idx);
    Xj_lst[j] = Xj;
    Sigj_inv_lst[j] = Xj.t() * Vinv * Xj + Vbinv / tau2(j);
  }
  
  if (learn_rho) {
    rho_re = 0.5;
    rho_im = 0.1;
    rho_ep = arma::cx_double(rho_re, rho_im);
    eta_re = rho_to_eta_cpp(rho_re, 1);
    eta_im = rho_to_eta_cpp(rho_re, sqrt(1 - pow(rho_re, 2)));
  }
  
  int iter;
  
  for (iter = 0; iter < n_mcmc; ++iter) {
    beta = beta_cbl_cpp(y, X, beta, Sigj_inv_lst, Vinv, sig2, p, all_ind,
                        select_idx_lst, Xj_lst);
    
    for (j = 0; j < p; ++j) {
      select_beta = beta.cols(as<uvec>(select_idx_lst[j]));
      btVbinvb(j) = dot(select_beta, select_beta);
    }
    
    ga = ga_cbl_cpp(lam2, sig2, btVbinvb, p);
    tau2 = 1 / ga;
    lam2 = lam2_cbl_cpp(tau2, r, delta, p);
    mu = X * beta.t();
    sig2 = sig2_cbl_cpp(a_sig, b_sig, n, p, y,
                        mu, tau2, btVbinvb, Vinv);
    
    if (learn_rho) {
      mat rho_mat = mat(2, 2, fill::ones);
      arma::mat ymat = reshape(y, 1, 2 * n);
      arma::mat sig2I = sig2 * I;
      arma::rowvec mut = mu.t();
      
      if (adapt && fmod(iter, Tb) == 0) {
        keep_tmp_re = keep_tmp_re / Tb;
        keep_tmp_im = keep_tmp_im / Tb;
        tune_re = get_tune_cpp(tune_re, keep_tmp_re, iter, adapt_step,
                               1/2, target_accept_rate);
        tune_im = get_tune_cpp(tune_im, keep_tmp_im, iter, adapt_step,
                               1/2, target_accept_rate);
        keep_tmp_re = 0;
        keep_tmp_im = 0;
      }
      
      eta_re_star = R::rnorm(eta_re, tune_re);
      rho_re_star = eta_to_rho_cpp(eta_re_star, 1);
      while(pow(rho_re_star, 2) + pow(rho_im, 2) > 1) {
        eta_re_star = R::rnorm(eta_re, tune_re);
        rho_re_star = eta_to_rho_cpp(eta_re_star, 1);
      }
      
      logpost_re = mh_log_ratio_re_cpp(eta_re, eta_re_star, rho_im,
                                       mut, ymat, sig2I);
      
      if(log(R::runif(0, 1)) < logpost_re) {
        eta_re = eta_re_star;
        ++keep_re;
        ++keep_tmp_re;
      }
      
      rho_re = eta_to_rho_cpp(eta_re, 1);
      eta_im_star = R::rnorm(eta_im, tune_im);
      logpost_im = mh_log_ratio_im_cpp(rho_re, eta_im, eta_im_star,
                                       mut, ymat, sig2I);
      
      if(log(R::runif(0, 1)) < logpost_im) {
        eta_im = eta_im_star;
        ++keep_im;
        ++keep_tmp_im;
      }
      
      rho_im = eta_to_rho_cpp(eta_im, sqrt(1 - pow(rho_re, 2)));
      dd = 1/ ((1 - rho_re) - pow(rho_im, 2) / (1 + rho_re));
      aa = 1 + rho_re;
      rho_mat(0, 0) = (aa + pow(rho_im, 2)) / pow(aa, 2);
      rho_mat(0, 1) = -rho_im / aa;
      rho_mat(1, 0) = -rho_im / aa;
      Vinv = sp_mat(kron(rho_mat, dd * I));
    }
    
    for (int j = 0; j < p; ++j) {
      arma::mat Xj = Xj_lst[j];
      arma::mat S = Xj.t() * Vinv * Xj + Vbinv / tau2(j);
      Sigj_inv_lst[j] = (S + S.t()) / 2 + eps * eye(2, 2);
    }
    
    if ((iter + 1) > warm) {
      if (fmod(iter + 1 , thin) == 0) {
        int d = (iter + 1 - warm) / thin;
        draws(d - 1, span(0, 2 * p - 1)) = beta;
        draws(d - 1, span(2 * p, 3 * p - 1)) = tau2.t();
        draws(d - 1, 3 * p) = lam2;
        draws(d - 1, 3 * p + 1) = sig2;
        
        if (learn_rho) {
          draws(d - 1, draws.n_cols - 2) = rho_re;
          draws(d - 1, draws.n_cols - 1) = rho_im;
        }
      }
    }
  }
  
  keep_im = keep_im / n_mcmc;
  keep_re = keep_re / n_mcmc;
  
  return List::create(
    _["draws"] = draws, 
    _["warm"] = warm, 
    _["thin"] = thin, 
    _["nmcmc"] = n_mcmc,
    _["accept_re"] = keep_re, 
    _["accept_im"] = keep_im, 
    _["tune_re"] = tune_re, 
    _["tune_im"] = tune_im
  );
}