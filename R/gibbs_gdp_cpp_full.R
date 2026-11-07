#' Gibbs Sampling for GDP Prior Model
#' 
#' @description Implements Gibbs sampling for Bayesian linear regression with 
#' Generalized Double Pareto (GDP) prior using Rcpp for computational efficiency.
#'
#' @param n_mcmc Integer. Number of MCMC iterations
#' @param start_lst List. Initial values containing:
#'        - beta: Initial regression coefficients
#'        - invtau: Initial inverse tau parameter
#'        - lam: Initial lambda parameter
#'        - sig2: Initial variance parameter
#' @param warm Integer. Number of warm-up draws to discard (default: 0)
#' @param thin Integer. Keep every thin-th draw for thinning (default: 1)
#' @param number_par Integer. Total number of parameters (default: length(name_par))
#' @param name_par Character vector. Names of parameters
#' @param y Numeric vector/matrix. Response variable (n x 1)
#' @param X Numeric matrix. Design matrix (n x p)
#' @param alpha Numeric. Shape parameter for GDP prior (default: 1)
#' @param eta Numeric. Scale parameter for GDP prior (default: 1)
#' @param griddy Logical. Whether to use griddy Gibbs sampling (default: FALSE)
#' @param a_sig,b_sig Numeric. Shape and scale parameters for sig2 prior (default: 0.5)
#' @param metro Logical. Whether to use Metropolis-Hastings (default: FALSE)
#' @param metro.tune Numeric. Tuning parameter for Metropolis-Hastings (default: NULL)
#' @param n_chains Integer. Number of MCMC chains to run (default: 1)
#'
#' @return List containing:
#' \itemize{
#'   \item draws - Matrix of posterior draws
#'   \item n_mcmc - Number of MCMC iterations
#'   \item warm - Number of warm-up iterations
#'   \item thin - Thinning interval
#'   \item sampleidx - Vector indicating which iterations were kept
#' }
#'
#' @import Rcpp
#' @importFrom stats rgamma
#'
#' @export
gibbs_gdp_cpp_full <- function(n_mcmc, 
                               start_lst, 
                               warm = 0, 
                               thin = 1,
                               number_par = length(name_par), 
                               name_par, 
                               y, 
                               X, 
                               alpha = 1, 
                               eta = 1, 
                               griddy = FALSE,
                               a_sig = 0.5, 
                               b_sig = 0.5,
                               metro = FALSE, 
                               metro.tune = NULL, 
                               n_chains = 1) {
  
  #-------------------------------
  # Input validation
  #-------------------------------
  if (!is.list(start_lst) || !all(c("beta", "invtau", "lam", "sig2") %in% names(start_lst))) {
    stop("start_lst must be a list containing beta, invtau, lam, and sig2")
  }
  
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(y)) y <- as.matrix(y)
  
  if (nrow(X) != length(y)) {
    stop("Number of rows in X must match length of y")
  }
  
  #-------------------------------
  # Initialize parameters
  #-------------------------------
  beta <- matrix(start_lst$beta, nrow = 1)
  invtau <- start_lst$invtau
  lam <- start_lst$lam
  sig2 <- start_lst$sig2
  
  #-------------------------------
  # Setup storage
  #-------------------------------
  sampleidx <- seq(from = (warm + thin), to = n_mcmc, by = thin)
  draws <- matrix(NA, nrow = length(sampleidx), ncol = number_par)
  colnames(draws) <- name_par
  
  #-------------------------------
  # Main MCMC sampling
  #-------------------------------
  mcmc_result <- mcmc_gdp_cpp(
    n_mcmc = n_mcmc,
    warm = warm,
    thin = thin,
    beta = beta,
    invtau = invtau,
    lam = lam,
    sig2 = sig2,
    y = y,
    X = X,
    a_sig = a_sig,
    b_sig = b_sig,
    alpha = alpha,
    eta = eta,
    griddy = griddy,
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