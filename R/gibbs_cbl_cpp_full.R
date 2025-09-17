#' Gibbs Sampling for Complex Bayesian Lasso Regression
#'
#' Implements Gibbs sampling for a complex Bayesian Lasso regression model
#'
#' @param n_mcmc Number of MCMC iterations
#' @param start_lst List of initial values for the MCMC chain
#' @param warm Number of warm-up iterations to discard (default: 0)
#' @param thin Thinning interval for MCMC samples (default: 1)
#' @param number_par Total number of parameters (default: length of name_par)
#' @param name_par Vector of parameter names
#' @param rho_b Correlation parameter for beta (default: 0)
#' @param rho_ep Correlation parameter for errors (default: 0)
#' @param y Response vector
#' @param X Design matrix
#' @param r Shape parameter for lam2 prior (default: 0.5)
#' @param delta Scale parameter for lam2 prior (default: 0.5)
#' @param a_sig Shape parameter for sigma2 prior (default: 0.5)
#' @param b_sig Scale parameter for sigma2 prior (default: 0.5)
#' @param learn_rho Logical, whether to learn correlation parameters (default: FALSE)
#' @param n_chains Number of MCMC chains (default: 1)
#' @param keep.re Acceptance rate for real part of correlation (default: 0)
#' @param keep.im Acceptance rate for imaginary part of correlation (default: 0)
#' @param tune.re Tuning parameter for real part of correlation (default: 0.2)
#' @param tune.im Tuning parameter for imaginary part of correlation (default: 0.2)
#' @param tune.len Length of tuning period (default: 50)
#' @param target.accept.rate Target acceptance rate (default: 0.3)
#' @param adapt Logical, whether to use adaptive tuning (default: TRUE)
#' @param adapt_step Adaptive tuning step size (default: 0.2)
#'
#' @return List containing MCMC draws and metadata
#' @export
gibbs_cbl_cpp_full <- function(n_mcmc, start_lst, warm = 0, thin = 1, 
                               number_par = length(name_par), name_par, rho_b = 0,
                               rho_ep = 0,
                               y, X, r = 0.5, delta = 0.5, a_sig = 0.5, b_sig = 0.5,
                               learn_rho = FALSE, 
                               n_chains = 1,
                               keep.re = 0, keep.im = 0,
                               tune.re = 0.2, tune.im = 0.2,
                               tune.len = 50, 
                               target.accept.rate = 0.3,
                               adapt = TRUE,
                               adapt_step = 0.2) {
  
  #-------------------------------
  # Starting values
  #-------------------------------
  beta <- matrix(start_lst$beta, nrow = 1)
  tau2 <- start_lst$tau2
  lam2 <- start_lst$lam2
  sig2 <- start_lst$sig2
  
  
  #-------------------------------
  # Storage
  #-------------------------------
  # idx recording which iterations are sampled
  sampleidx <- seq(from = (warm + thin), to = n_mcmc, by = thin)
  # sample storage
  draws <- matrix(NA, nrow = length(sampleidx), ncol = number_par)
  
  mcmc_cbl_cpp(n_mcmc, warm, thin, beta, tau2, lam2, sig2,
               y, X, rho_b, rho_ep, r, delta, a_sig, b_sig,
               learn_rho, n_chains, keep.re, keep.im,
               tune.re, tune.im, tune.len, adapt, adapt_step, 
               target.accept.rate, draws)
  
}



