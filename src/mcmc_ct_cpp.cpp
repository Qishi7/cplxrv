#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

double get_tune_cpp(double tune, double keep, int iter, double step, 
                    double exp_term, double target);
arma::mat cplx2real_cov_cpp(arma::mat var, arma::cx_double rho);
double rho_to_eta_cpp(double rho, double a);
double eta_to_rho_cpp(double eta, double a);
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);
arma::vec dmvn_cpp(arma::mat x,  arma::rowvec mean,  arma::mat sigma, bool logd = false);
double mh_log_ratio_re_cpp(double eta_re, double eta_re_star, double rho_im, 
                           arma::rowvec mut, arma::mat ymat, arma::mat sig2I);
double mh_log_ratio_im_cpp(double rho_re, double eta_im, double eta_im_star, 
                           arma::rowvec mut, arma::mat ymat, arma::mat sig2I);

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat beta_ct_cpp(arma::mat& XtVinvy, arma::mat& Siginv, double sig2) {
  // arma::mat A = XtX + diagmat(1/tau2);
  arma::mat Sig = inv_sympd(Siginv);
  arma::vec mu_beta = Sig * XtVinvy;
  arma::mat Sig_beta = sig2 * Sig;
  // arma::mat S = mvrnormArma(1, mu_beta, Sig_beta);
  // this is a 1 by p row vector
  return mvrnormArma(1, mu_beta, Sig_beta);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec tau2_ct_cpp(arma::mat& beta, arma::cx_double rho_b, double df, int p, 
                      double sig2) {
  arma::vec bb = arma::vectorise(beta);
  arma::vec bb_re = bb.subvec(0, p - 1);
  arma::vec bb_im = bb.subvec(p, 2*p - 1);
  arma::cx_vec bb_cplx = cx_vec(bb_re, bb_im);
  arma::vec bb_mod_sq = pow(bb_re, 2) + pow(bb_im, 2);
  arma::vec bbb;
  if (real(rho_b) == 0 && imag(rho_b) == 0) {
    bbb = bb_mod_sq / sig2;
  } else {
    bbb = (bb_mod_sq - real(rho_b * pow(conj(bb_cplx), 2))) / 
      (sig2 * (1 - (pow(real(rho_b), 2) + pow(imag(rho_b), 2))));
  }
  
  arma::vec tau2_vec(p);
  arma::vec rate = df + bbb;
  
  for(int i = 0; i < p; ++i) {
    tau2_vec(i) = 1 / R::rgamma(df + 1, 1/rate(i));
  }
  return(tau2_vec);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double sig2_ct_cpp(double asig, double bsig, int n, int p, 
                   arma::vec y, arma::mat beta, arma::mat X, arma::mat Vinv, 
                   arma::mat Vbinv) {
  double a_sig = n + p + asig;
  arma::vec e = vectorise(y - X * beta.t());
  double b_sig = as_scalar(bsig + 0.5 * (e.t() * Vinv * e + beta * Vbinv * beta.t()));
  return 1 / R::rgamma(a_sig, 1/b_sig);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List mcmc_ct_cpp(int n_mcmc, int warm, int thin, 
                 arma::mat beta, arma::vec tau2, double sig2,
                 arma::vec& y, arma::mat& X, arma::cx_double rho_b, 
                 arma::cx_double rho_ep, double df, double a_sig, double b_sig,
                 bool learn_rho, int n_chains,       
                 double keep_re, double keep_im,
                 double tune_re, double tune_im,
                 int tune_len, bool adapt, double& adapt_step, 
                 double target_accept_rate,
                 arma::mat draws) {
  
  // #-------------------------------
  // # Here create some objects that will be used later.
  // #-------------------------------
  
  int n = X.n_rows / 2;
  int p = X.n_cols / 2;
  arma::mat I = eye(n, n);
  arma::mat V = cplx2real_cov_cpp(2 * I, rho_ep);
  arma::mat Vinv = inv_sympd(V);
  arma::mat Xt = X.t();
  arma::mat XtVinvX = Xt * Vinv * X;
  arma::mat XtVinvy = Xt * Vinv * y;
  
  // So far only consider rho_b = 0 case
  arma::mat Vbinv = diagmat(1/join_cols(tau2, tau2));
  arma::mat Siginv = XtVinvX + Vbinv;
  
  // Adaptive tuning
  int Tb = tune_len;
  double keep_tmp_re = keep_re;
  double keep_tmp_im = keep_im;
  double rho_re;
  double rho_im;
  double eta_re;
  double eta_im;
  
  if (learn_rho) {
    rho_re = 0.5;
    rho_im = 0.1;
    rho_ep = arma::cx_double(rho_re, rho_im);
    eta_re = rho_to_eta_cpp(rho_re, 1);
    eta_im = rho_to_eta_cpp(rho_re, sqrt(1 - pow(rho_re, 2)));
  }
  
  // #-------------------------------
  // # Gibbs sampling algorithm
  // #------------------------------- 
  
  for (int iter = 0; iter < n_mcmc; ++iter) {
    
    beta = beta_ct_cpp(XtVinvy, Siginv, sig2);
    
    tau2 = tau2_ct_cpp(beta, rho_b, df, p, sig2);
    
    sig2 = sig2_ct_cpp(a_sig, b_sig, n, p, y, beta, X, Vinv, Vbinv);
    
    // update ---
    Vbinv = diagmat(1/join_cols(tau2, tau2));
    
    
    // #-------------------------------
    // # HERE put Metropolis steps if needed
    // #-------------------------------
    
    if (learn_rho) {
      arma::vec mu = X * beta.t();
      arma::mat sig2I = sig2 * I;  
      
      // Rcout << "mu: " << mu << endl;
      // # Update tuning parameter
      // #------------------------
      
      if (adapt && fmod(iter, Tb) == 0) {
        
        // # Adaptive tuning
        keep_tmp_re = keep_tmp_re / Tb;
        keep_tmp_im = keep_tmp_im / Tb;
        tune_re = get_tune_cpp(tune_re, keep_tmp_re, iter, adapt_step,
                               1/2, target_accept_rate);
        tune_im = get_tune_cpp(tune_im, keep_tmp_im, iter, adapt_step,
                               1/2, target_accept_rate);
        keep_tmp_re = 0;
        keep_tmp_im = 0;
      }
      
      // #----------------------------
      // # Random walk Normal proposal
      // #----------------------------
      
      double eta_re_star = R::rnorm(eta_re, tune_re);
      double rho_re_star = eta_to_rho_cpp(eta_re_star, 1);
      while(pow(rho_re_star, 2) + pow(rho_im, 2) > 1) {
        eta_re_star = R::rnorm(eta_re, tune_re);
        rho_re_star = eta_to_rho_cpp(eta_re_star, 1);
      }
      
      
      // # Acceptance Probability
      // #----------------------------
      double logpost_re = mh_log_ratio_re_cpp(eta_re, eta_re_star, rho_im,
                                              mu, y, sig2I);
      
      if(log(R::runif(0, 1)) < logpost_re) {
        eta_re = eta_re_star;
        // keep_re = keep_re + 1;
        // keep_tmp_re = keep_tmp_re + 1;
        
        ++keep_re;
        ++keep_tmp_re;
      }
      
      rho_re = eta_to_rho_cpp(eta_re, 1);
      
      double eta_im_star = R::rnorm(eta_im, tune_im);
      double logpost_im = mh_log_ratio_im_cpp(rho_re, eta_im, eta_im_star,
                                              mu, y, sig2I);
      
      if(log(R::runif(0, 1)) < logpost_im) {
        eta_im = eta_im_star;
        // keep_im = keep_im + 1;
        // keep_tmp_im = keep_tmp_im + 1;
        ++keep_im;
        ++keep_tmp_im;
      }
      
      rho_im = eta_to_rho_cpp(eta_im, sqrt(1 - pow(rho_re, 2)));
      
      V = join_cols(join_rows((1 + rho_re) * I, rho_im * I),
                    join_rows(rho_im * I, (1 - rho_re) * I));
      
      Vinv = inv_sympd(V);
      XtVinvX = Xt * Vinv * X;
      XtVinvy = Xt * Vinv * y;
      Siginv = XtVinvX + Vbinv;
    }
    
    // # -----------------------------------------------------
    // #  Save samples
    // # -----------------------------------------------------
    
    if ((iter + 1) > warm) {
      if (fmod(iter + 1 , thin) == 0) {
        int d = (iter + 1 - warm) / thin;
        draws(d - 1, span(0, 2 * p - 1)) = beta;
        draws(d - 1, span(2 * p, 3 * p - 1)) = tau2.t();
        draws(d - 1, 3 * p) = sig2;
        
        if (learn_rho) {
          draws(d - 1, draws.n_cols - 2) = rho_re;
          draws(d - 1, draws.n_cols - 1) = rho_im;
        }
      }
    }
    
    // # -----------------------------------------------------
    // #  Print iterations
    // # -----------------------------------------------------
    if ((iter + 1) % 500 == 0) {
      Rcout << "\r" << "Cplx t MCMC Iter: " << iter + 1 << std::flush;
    }
  }
  
  // # Acceptance Probability
  // #----------------------------
  keep_im = keep_im / n_mcmc;
  keep_re = keep_re / n_mcmc;
  
  
  // # Write output
  // #--------------
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

