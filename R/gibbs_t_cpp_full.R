#' @title Gibbs Sampling for Bayesian Linear Regression with Student-t
#' @description Implements Gibbs sampling for Bayesian linear regression with 
#'             Student-t using Rcpp for computational efficiency.
#'
#' @import Rcpp
#' @import stats
#' @import emulator
#' @import LaplacesDemon
#'
#' @param n_mcmc Integer. Number of MCMC iterations
#' @param start_lst List. Initial values for MCMC chain containing:
#'                 - beta: Initial regression coefficients
#'                 - tau2: Initial scale parameter
#'                 - sig2: Initial variance parameter
#' @param warm Integer. Number of warm-up draws to discard (default: 0)
#' @param thin Integer. Keep every thin-th draw for thinning (default: 1)
#' @param number_par Integer. Total number of parameters (default: length(name_par))
#' @param name_par Character vector. Names of parameters
#' @param y Numeric vector/matrix. Response variable (n x 1)
#' @param X Numeric matrix. Design matrix (n x p)
#' @param a_tau Numeric. Shape parameter for tau2 prior (default: 0.5)
#' @param b_tau Numeric. Scale parameter for tau2 prior (default: 0.5)
#' @param a_sig Numeric. Shape parameter for sig2 prior (default: 0.5)
#' @param b_sig Numeric. Scale parameter for sig2 prior (default: 0.5)
#' @param metro Logical. Whether to use Metropolis-Hastings (default: FALSE)
#' @param metro.tune Numeric. Tuning parameter for Metropolis-Hastings (default: NULL)
#' @param n_chains Integer. Number of MCMC chains to run (default: 1)
#'
#' @return List containing:
#' \itemize{
#'   \item draws - Matrix of posterior draws (rows: samples, cols: parameters)
#'   \item n_mcmc - Number of MCMC iterations
#'   \item warm - Number of warm-up iterations
#'   \item thin - Thinning interval
#'   \item sampleidx - Vector indicating which iterations were kept
#' }
#'
#' @examples
#' \dontrun{
#' # Generate sample data
#' n <- 100
#' p <- 3
#' X <- matrix(rnorm(n*p), n, p)
#' true_beta <- c(1, -0.5, 2)
#' y <- X %*% true_beta + rt(n, df=3)
#'
#' # Set initial values
#' start <- list(
#'   beta = rep(0, p),
#'   tau2 = 1,
#'   sig2 = 1
#' )
#'
#' # Run Gibbs sampler
#' res <- gibbs_t_cpp_full(
#'   n_mcmc = 5000,
#'   start_lst = start,
#'   warm = 1000,
#'   thin = 2,
#'   name_par = c("beta1", "beta2", "beta3", "tau2", "sig2"),
#'   y = y,
#'   X = X
#' )
#' }
#'
#' @export
gibbs_t_cpp_full <- function(n_mcmc, 
                             start_lst, 
                             warm = 0, 
                             thin = 1,
                             number_par = length(name_par), 
                             name_par,
                             y, 
                             X, 
                             a_tau = 0.5, 
                             b_tau = 0.5,
                             a_sig = 0.5, 
                             b_sig = 0.5,
                             metro = FALSE, 
                             metro.tune = NULL,
                             n_chains = 1) {
  
  #-------------------------------
  # Input validation
  #-------------------------------
  if (!is.list(start_lst) || !all(c("beta", "tau2", "sig2") %in% names(start_lst))) {
    stop("start_lst must be a list containing beta, tau2, and sig2")
  }
  
  if (nrow(X) != length(y)) {
    stop("Number of rows in X must match length of y")
  }
  
  #-------------------------------
  # Starting values
  #-------------------------------
  beta <- matrix(start_lst$beta, nrow = 1)
  tau2 <- start_lst$tau2
  sig2 <- start_lst$sig2
  
  #-------------------------------
  # Storage initialization
  #-------------------------------
  sampleidx <- seq(from = (warm + thin), to = n_mcmc, by = thin)
  draws <- matrix(NA, nrow = length(sampleidx), ncol = number_par)
  colnames(draws) <- name_par
  
  #-------------------------------
  # Main MCMC sampling
  #-------------------------------
  mcmc_result <- mcmc_t_cpp(
    n_mcmc = n_mcmc,
    warm = warm,
    thin = thin,
    beta = beta,
    tau2 = tau2,
    sig2 = sig2,
    y = y,
    X = X,
    a_tau = a_tau,
    b_tau = b_tau,
    a_sig = a_sig,
    b_sig = b_sig,
    draws = draws
  )
  
  #-------------------------------
  # Return results
  #-------------------------------
  return(list(
    draws = mcmc_result,
    n_mcmc = n_mcmc,
    warm = warm,
    thin = thin,
    sampleidx = sampleidx
  ))
}