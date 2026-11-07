#include <stdio.h>
#include <math.h> 
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

arma::vec dmvn_cpp(arma::mat x,arma::rowvec mean,arma::mat sigma,bool logd = false);
double get_tune_cpp(double tune, double keep, int iter, double step, 
                      double exp_term, double target);
arma::mat cplx2real_cov_cpp(arma::mat var, arma::cx_double rho);  
double rho_to_eta_cpp(double rho, double a);
double eta_to_rho_cpp(double eta, double a);
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);
arma::vec rrinvgauss_vec(int n, arma::vec mu, arma::vec lambda);
arma::uvec keep_elem(const arma::vec& x, const arma::ivec& to_remove);
arma::vec sample_vec(const arma::vec& x, int size, bool replace, const arma::vec& prob);

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat beta_cgdp_cpp(arma::vec& y, arma::mat& X, arma::mat beta,
                        List& Sigj_inv_lst, arma::mat& Vinv, double sig2, int& p,
                        arma::vec& all_ind, List& select_idx_lst, List& Xj_lst) {
  // arma::mat Sigj_inv;
  arma::ivec remove_idx;
  arma::mat Sigj;
  arma::uvec idx;
  arma::mat Aj;
  // arma::uvec select_idx;
  arma::vec mu;
  arma::mat sigma;
  
  for (int j = 0; j < p; ++j) {
    arma::mat Sigj_inv = Sigj_inv_lst[j];
    Sigj = inv_sympd(Sigj_inv);
    remove_idx = {j, j + p};
    // arma::uvec idx = {(unsigned int) j, (unsigned int) j + p};
    // arma:: mat XX = X.cols(idx);
    idx = keep_elem(all_ind, remove_idx);
    Aj = y - X.cols(idx) * beta.cols(idx).t();
    // select_idx = {(unsigned int) j, (unsigned int) j + p};
    // mu = Sigj * X.cols(select_idx).t() * Vinv * Aj;
    arma::mat Xj = Xj_lst[j];
    mu = Sigj * Xj.t() * Vinv * Aj;
    sigma = sig2 * Sigj;
    // select_idx = as<uvec>(select_idx_lst[j]);
    beta.cols(as<uvec>(select_idx_lst[j])) = mvrnormArma(1, mu, sigma);
  }
  return(beta);
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec ga_cgdp_cpp(arma::vec lam, double sig2, arma::vec& btb, int& p) {
  arma::vec lam2 = pow(lam, 2);
  return(rrinvgauss_vec(p, sqrt(lam2 * sig2 / btb), lam2));
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double sig2_cgdp_cpp(double asig, double bsig, int n, int p, 
                     arma::vec& y, arma::vec& mu, arma::vec tau2, arma::vec& btb, 
                     arma::mat& Vinv) {
  arma::vec e = y - mu;
  double a_sig = n + p + asig;
  // arma::vec e = vectorise(y - X * beta.t());
  double b_sig = as_scalar(bsig + 0.5 * (e.t() * Vinv * e + sum(btb / tau2)));
  return (1 / R::rgamma(a_sig, 1/b_sig));
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec lam_cgdp_cpp(arma::vec& btb, double alpha, double eta, double sig2, 
                       int& p) {
  
  arma::vec lam_vec(p);
  arma::vec rate = sqrt(btb / sig2) + 2 * eta;
  double shape = 2 * alpha + 2;
  
  for(int j = 0; j < p; ++j) {
    lam_vec(j) = R::rgamma(shape, 1/rate(j));
  }
  return(lam_vec);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double mh_log_ratio_re_cpp_gdp(double eta_re, double eta_re_star, double rho_im, 
                           arma::vec& mu, arma::vec& yr, double sig2,
                           int& n, arma::mat& I) {
  
  double rho_re = eta_to_rho_cpp(eta_re, 1);
  double rho_re_star = eta_to_rho_cpp(eta_re_star, 1);
  // int p = 0.5 * yr.n_elem;
  // arma::mat I = eye(p, p);
  arma::mat Vrr = (1 + rho_re) * I;
  arma::mat Vrr_star = (1 + rho_re_star) * I;
  arma::mat Vii = (1 - rho_re) * I;
  arma::mat Vii_star = (1 - rho_re_star) * I;
  arma::mat Vri = rho_im * I;
  arma::mat Sig = sig2 * 
    join_cols(join_rows(Vrr, Vri), join_rows(Vri, Vii));
  arma::mat Sig_star = sig2 * 
    join_cols(join_rows(Vrr_star, Vri), join_rows(Vri, Vii_star));
  
  arma::mat ymat = reshape(yr, 1, 2 * n);
  arma::vec A_star = dmvn_cpp(ymat, mu.t(), Sig_star, true) + 
    eta_re_star - 2 * log(1 + exp(eta_re_star));
  arma::vec A = dmvn_cpp(ymat, mu.t(), Sig, true) + 
    eta_re - 2 * log(1 + exp(eta_re));
  return(as_scalar(A_star - A));
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double mh_log_ratio_im_cpp_gdp(double rho_re, double eta_im, double eta_im_star, 
                           arma::vec& mu, arma::vec& yr, double sig2,
                           int& n, arma::mat& I) {
  
  double a = sqrt(1 - pow(rho_re, 2));
  double rho_im = eta_to_rho_cpp(eta_im, a);
  double rho_im_star = eta_to_rho_cpp(eta_im_star, a);
  // int p = 0.5 * yr.n_elem;
  // arma::mat I = eye(p, p);
  arma::mat Vrr = (1 + rho_re) * I;
  arma::mat Vii = (1 - rho_re) * I;
  arma::mat Vri = rho_im * I;
  arma::mat Vri_star = rho_im_star * I;
  arma::mat Sig = sig2 * 
    join_cols(join_rows(Vrr, Vri), join_rows(Vri, Vii));
  arma::mat Sig_star = sig2 * 
    join_cols(join_rows(Vrr, Vri_star), join_rows(Vri_star, Vii));
  
  arma::mat ymat = reshape(yr, 1, 2 * n);

  arma::vec A_star = dmvn_cpp(ymat, mu.t(), Sig_star, true) + 
    eta_im_star - 2 * log(1 + exp(eta_im_star));
  arma::vec A = dmvn_cpp(ymat, mu.t(), Sig, true) + 
    eta_im - 2 * log(1 + exp(eta_im));
  return(as_scalar(A_star - A));
  // return(Sig);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec a_cond_log_cgdp_cpp(arma::vec a, int& p, arma::vec& btb, double sig2, 
                              double eta) {
  arma::vec a_value = p * log((1 - a) / a) - p * log(a) -
    ((1 / a) + 1) * sum(log(1 + (sqrt(btb) / (sqrt(sig2) * eta))));
  
  return(a_value);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec e_cond_log_cgdp_cpp(arma::vec e, int& p, arma::vec& btb, double sig2, 
                              double alpha) {
  arma::vec d = e / (1 - e);
  arma::vec e_value(d.n_elem);
  for (int k = 0; k < d.n_elem; ++k) {
    e_value(k) = 2 * p * log(d(k)) - 
      (alpha + 2) * sum(log(1 + (sqrt(btb) * d(k) / sqrt(sig2))));
  }
  return(e_value);
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
List mcmc_cgdp_cpp(const int& n_mcmc, const int& warm, const int& thin, 
                   arma::mat beta, arma::vec tau2, arma::vec lam, double sig2,
                   arma::vec& y, arma::mat& X, arma::cx_double rho_b, 
                   arma::cx_double rho_ep, double alpha, double eta,
                   double a_sig, double b_sig,
                   const bool& learn_rho, const bool& griddy, const int& n_chains,       
                   double keep_re, double keep_im,
                   double tune_re, double tune_im,
                   int& tune_len, bool& adapt, const double& adapt_step, 
                   const double& target_accept_rate,
                   arma::mat draws) {
  
  // #-------------------------------
  // # Here create some objects that will be used later.
  // #-------------------------------
  
  // Adaptive tuning
  // Adaptive tuning
  double rho_re, rho_im, eta_re, eta_im;
  const int Tb = tune_len;
  double keep_tmp_re = keep_re;
  double keep_tmp_im = keep_im;
  
  int n = X.n_rows / 2;
  int p = X.n_cols / 2;
  arma::mat I = eye(n, n);
  arma::mat V = cplx2real_cov_cpp(2 * I, rho_ep);
  arma::mat Vinv = inv_sympd(V);
  arma::mat Xt = X.t();
  arma::mat XtVinvX = Xt * Vinv * X;
  arma::mat XtVinvy = Xt * Vinv * y;
  arma::vec all_ind = linspace(0, 2 * p - 1, 2 * p);
  arma::vec ga;
  arma::vec mu;
  arma::ivec remove_idx;
  arma::uvec select_idx;
  arma::mat select_beta;
  arma::vec btVbinvb = zeros(p);
  
  arma::mat Vbinv = eye(2, 2);
  
  List select_idx_lst(p);
  for (int j = 0; j < p; ++j) {
    select_idx = {(unsigned int) j, (unsigned int) j + p};
    // arma::mat Xj = X.cols(select_idx);
    select_idx_lst[j] = select_idx;
    
  }
  
  List Xj_lst(p);
  for (int j = 0; j < p; ++j) {
    // select_idx = {(unsigned int) j, (unsigned int) j + p};
    // arma::mat Xj = X.cols(select_idx);
    arma::mat Xj = X.cols(as<uvec>(select_idx_lst[j]));
    Xj_lst[j] = Xj;
    
  }
  
  
  List Sigj_inv_lst(p);
  for (int j = 0; j < p; ++j) {
    // select_idx = {(unsigned int) j, (unsigned int) j + p};
    // arma::uvec select_idx = select_idx_lst[j];
    // arma::mat Xj = X.cols(as<uvec>(select_idx_lst[j]));
    arma::mat Xj = Xj_lst[j];
    Sigj_inv_lst[j] = Xj.t() * Vinv * Xj + Vbinv / tau2(j);
  }
  
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
    
    beta = beta_cgdp_cpp(y, X, beta, Sigj_inv_lst, Vinv, sig2, p, all_ind,
                         select_idx_lst, Xj_lst);
    
    //  only consider rho_b == 0 case
    for (int j = 0; j < p; ++j) {
      select_beta = beta.cols(as<uvec>(select_idx_lst[j]));
      btVbinvb(j) = as_scalar(select_beta * select_beta.t());
    }
    
    mu = X * beta.t();
    
    ga = ga_cgdp_cpp(lam, sig2, btVbinvb, p);
    
    tau2 = 1 / ga;
    
    lam = lam_cgdp_cpp(btVbinvb, alpha, eta, sig2, p);
    
    sig2 = sig2_cgdp_cpp(a_sig, b_sig, n, p, y, mu, tau2, btVbinvb, Vinv);
    
    if (griddy == true) {
      arma::vec aa = linspace(0.1, 0.9, 100);
      arma::vec ee = linspace(0.1, 0.9, 100);
      arma::vec weight_a = exp(a_cond_log_cgdp_cpp(aa, p, btVbinvb, sig2,
                                                   eta));
      arma::vec a = sample_vec(aa, 1, false,
                                                weight_a/sum(weight_a));
      alpha = (1 / a[0]) - 1;
      
      arma::vec weight_e = exp(e_cond_log_cgdp_cpp(ee, p, btVbinvb, sig2,
                                                   alpha));
      arma::vec e = sample_vec(ee, 1, false,
                                                weight_e/sum(weight_e));
      eta = (1 / e[0]) - 1;
    }
    
    // #-------------------------------
    // # HERE put Metropolis steps if needed
    // #-------------------------------
    
    if (learn_rho) {
      
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
      double logpost_re = mh_log_ratio_re_cpp_gdp(eta_re, eta_re_star, rho_im,
                                              mu, y, sig2, n, I);
      if(log(R::runif(0, 1)) < logpost_re) {
        eta_re = eta_re_star;
        ++keep_re;
        ++keep_tmp_re;
      }
      
      rho_re = eta_to_rho_cpp(eta_re, 1);
      
      double eta_im_star = R::rnorm(eta_im, tune_im);
      double logpost_im = mh_log_ratio_im_cpp_gdp(rho_re, eta_im, eta_im_star,
                                              mu, y, sig2, n, I);
      
      if(log(R::runif(0, 1)) < logpost_im) {
        eta_im = eta_im_star;
        ++keep_im;
        ++keep_tmp_im;
      }
      
      rho_im = eta_to_rho_cpp(eta_im, sqrt(1 - pow(rho_re, 2)));
      
      V = join_cols(join_rows((1 + rho_re) * I, rho_im * I),
                    join_rows(rho_im * I, (1 - rho_re) * I));
      
      Vinv = inv_sympd(V);
    }
    
    // # ----------
    // ## update
    // # ---------
    
    for (int j = 0; j < p; ++j) {
      arma::mat Xj = Xj_lst[j];
      Sigj_inv_lst[j] = Xj.t() * Vinv * Xj + Vbinv / tau2(j);
    }
    
    // # -----------------------------------------------------
    // #  Save samples
    // # ----------------------------------------------------
    
    if ((iter + 1) > warm) {
      if (fmod(iter + 1 , thin) == 0) {
        int d = (iter + 1 - warm) / thin;
        draws(d - 1, span(0, 2 * p - 1)) = beta;
        draws(d - 1, span(2 * p, 3 * p - 1)) = tau2.t();
        draws(d - 1, span(3 * p, 4 * p - 1)) = lam.t();
        draws(d - 1, 4 * p) = sig2;
        
        if (griddy) {
          draws(d - 1, 4 * p + 1) = alpha;
          draws(d - 1, 4 * p + 2) = eta;
        }
        
        if (learn_rho) {
          draws(d - 1, draws.n_cols - 2) = rho_re;
          draws(d - 1, draws.n_cols - 1) = rho_im;
        }
      }
    }
    
    
    // # -----------------------------------------------------
    // #  Print iterations
    // # -----------------------------------------------------
    if ((iter + 1) % 100 == 0) {
      Rcout << "\r" << "Cplx GDP MCMC Iter: " << iter + 1 << std::flush;
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

