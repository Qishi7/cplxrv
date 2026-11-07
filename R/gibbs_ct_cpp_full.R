#' Complex-valued Student-t Gibbs Sampler
#' 
#' @description Implements Gibbs sampling for complex-valued data with Student-t
#' using Rcpp for computational efficiency.
#'
#' @param n_mcmc Integer. Number of MCMC iterations
#' @param start_lst List. Initial values containing:
#'        - beta: Initial complex-valued regression coefficients
#'        - tau2: Initial scale parameter
#'        - sig2: Initial variance parameter
#' @param warm Integer. Number of warm-up draws to discard (default: 0)
#' @param thin Integer. Keep every thin-th draw for thinning (default: 1)
#' @param number_par Integer. Total number of parameters (default: length(name_par))
#' @param name_par Character vector. Names of parameters
#' @param y Complex vector/matrix. Response variable (n x 1)
#' @param X Numeric matrix. Design matrix (n x p)
#' @param rho_b Numeric. Prior correlation for beta (default: 0)
#' @param rho_ep Numeric. Prior correlation for error term (default: 0)
#' @param df Numeric. Degrees of freedom for Student-t distribution (default: 1)
#' @param a_sig,b_sig Numeric. Shape and scale parameters for sig2 prior (default: 0.5)
#' @param learn_rho Logical. Whether to learn correlation parameters (default: FALSE)
#' @param n_chains Integer. Number of MCMC chains (default: 1)
#' @param keep.re,keep.im Numeric. Initial acceptance rates for real and imaginary parts (default: 0)
#' @param tune.re,tune.im Numeric. Initial tuning parameters for MH steps (default: 0.2)
#' @param tune.len Integer. Frequency of adaptive tuning (default: 50)
#' @param target.accept.rate Numeric. Target acceptance rate for MH steps (default: 0.3)
#' @param adapt Logical. Whether to use adaptive tuning (default: TRUE)
#' @param adapt_step Numeric. Step size for adaptive tuning (default: 0.2)
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
#' @export
gibbs_ct_cpp_full <- function(n_mcmc, 
                              start_lst, 
                              warm = 0, 
                              thin = 1,
                              number_par = length(name_par), 
                              name_par,
                              y, 
                              X, 
                              rho_b = 0, 
                              rho_ep = 0, 
                              df = 1, 
                              a_sig = 0.5, 
                              b_sig = 0.5,
                              learn_rho = FALSE,
                              n_chains = 1,       
                              keep.re = 0, 
                              keep.im = 0,
                              tune.re = 0.2, 
                              tune.im = 0.2,
                              tune.len = 50, 
                              target.accept.rate = 0.3,
                              adapt = TRUE,
                              adapt_step = 0.2) {
  
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
  # Adaptive tuning setup
  #-------------------------------
  keep.tmp.re <- keep.re 
  keep.tmp.im <- keep.im 
  Tb <- tune.len
  
  #-------------------------------
  # Initialize parameters
  #-------------------------------
  beta <- matrix(start_lst$beta, nrow = 1)
  tau2 <- start_lst$tau2
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
  mcmc_result <- mcmc_ct_cpp(
    n_mcmc = n_mcmc,
    warm = warm,
    thin = thin,
    beta = beta,
    tau2 = tau2,
    sig2 = sig2,
    y = y,
    X = X,
    rho_b = rho_b,
    rho_ep = rho_ep,
    df = df,
    a_sig = a_sig,
    b_sig = b_sig,
    learn_rho = learn_rho,
    n_chains = n_chains,
    keep_re = keep.re,       
    keep_im = keep.im,        
    tune_re = tune.re,        
    tune_im = tune.im,        
    tune_len = tune.len,
    adapt = adapt,
    adapt_step = adapt_step,
    target_accept_rate = target.accept.rate,  
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