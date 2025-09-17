#' Full Bayesian Lasso with Gibbs Sampling via Rcpp
#' @description Implements Gibbs sampling for Bayesian Lasso regression using Rcpp
#' @param n_mcmc Number of MCMC iterations
#' @param start_lst List of starting values containing beta, tau2, lam2, and sig2
#' @param warm Number of warmup iterations to discard
#' @param thin Thinning interval for MCMC chain
#' @param number_par Total number of parameters
#' @param name_par Vector of parameter names 
#' @param y Response vector
#' @param X Design matrix
#' @param r Parameter for Bayesian Lasso prior (default = 0.5)
#' @param delta Parameter for Bayesian Lasso prior (default = 0.5)
#' @param a_sig Shape parameter for inverse gamma prior on sig2 (default = 0.5)
#' @param b_sig Scale parameter for inverse gamma prior on sig2 (default = 0.5)
#' @param metro Logical indicating whether to use Metropolis step (default = FALSE)
#' @param metro.tune Tuning parameter for Metropolis step (default = NULL)
#' @param n_chains Number of MCMC chains to run (default = 1)
#' @return List containing:
#'   \item{draws}{Matrix of posterior draws}
#'   \item{n_mcmc}{Number of MCMC iterations}
#'   \item{warm}{Number of warmup iterations}
#'   \item{thin}{Thinning interval}
#'   \item{sampleidx}{Indices of kept samples}
#' @examples
#' # Generate example data
#' n <- 100
#' p <- 10
#' X <- matrix(rnorm(n*p), n, p)
#' beta_true <- c(rep(1, 3), rep(0, p-3))
#' y <- X %*% beta_true + rnorm(n, 0, 1)
#' 
#' # Set starting values
#' start <- list(
#'   beta = rep(0, p),
#'   tau2 = rep(1, p), 
#'   lam2 = 1,
#'   sig2 = 1
#' )
#' 
#' # Run MCMC
#' result <- gibbs_bl_cpp_full(
#'   n_mcmc = 1000,
#'   start_lst = start,
#'   warm = 100,
#'   thin = 2,
#'   name_par = paste0("beta", 1:p),
#'   y = y,
#'   X = X
#' )
#' @export
gibbs_bl_cpp_full <- function(n_mcmc, start_lst, warm = 0, thin = 1, 
                              number_par = length(name_par), name_par, 
                              y, X, r = 0.5, delta = 0.5,
                              a_sig = 0.5, b_sig = 0.5,
                              metro = FALSE, metro.tune = NULL, 
                              n_chains = 1) {
  
  beta <- matrix(start_lst$beta, nrow = 1)
  tau2 <- start_lst$tau2
  lam2 <- start_lst$lam2
  sig2 <- start_lst$sig2
  invtau2 <- 1 / tau2
  
  sampleidx <- seq(from = (warm + thin), to = n_mcmc, by = thin)
  draws <- matrix(NA, nrow = length(sampleidx), ncol = number_par)
  
  mcmc_bl_cpp(n_mcmc, warm, thin, beta, tau2, lam2, sig2, y, X, r, delta,
              a_sig, b_sig, draws)
}