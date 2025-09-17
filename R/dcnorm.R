#' Complex-valued Density Representation
#'
#' @import emulator
#' @import mvtnorm
#'
#' @param z Complex-valued input (vector or matrix)
#' @param mean Mean of the distribution
#' @param var Variance-covariance matrix
#' @param rel Relation matrix (only used when is.circular = FALSE)
#' @param is.circular Logical; if TRUE, uses circular complex normal distribution
#' @param log Logical; if TRUE, returns log-density
#'
#' @return Density or log-density value(s)
#'
#' @export
dcnorm <- function(z, mean, var, rel, is.circular = FALSE, log = FALSE) {
  
  # Error checking and default values
  if (missing(z) || missing(mean) || missing(var)) {
    stop("Arguments 'z', 'mean', and 'var' are required.")
  }
  
  if (!is.circular && missing(rel)) {
    stop("Argument 'rel' is required when is.circular = FALSE.")
  }
  
  if (!is.matrix(var)) {
    stop("Argument 'var' must be a matrix.")
  }
  
  # Convert vector input to matrix
  if (is.vector(z)) {
    z <- matrix(z, nrow = length(z))
  }
  
  # Get the dimension of the variance matrix
  p <- ncol(var)
  
  if (is.circular) {
    # Circular complex normal distribution
    log_den <- - p * log(pi) - emulator::quad.form.inv(M = var, x = z - mean) - 
      2*sum(log(diag(chol(var))))
    # Note: 2*sum(log(diag(chol(var)))) = sum(log(svd(var)$d))
  } else {
    # Non-circular complex normal distribution
    zz <- sweep(z, 1, mean)  # Center the data
    var_inv <- solve(var)    # Inverse of variance matrix
    R <- emulator::cprod(rel, var_inv)
    P <- Conj(var) - R %*% rel
    p_conj_inv <- solve(Conj(P))
    qqbar <- emulator::quad.form(p_conj_inv, zz) - Re(t(R %*% zz) %*% p_conj_inv %*% zz)
    
    log_den <- -(p * log(pi) + sum(log(diag(chol(Re(var %*% P))))) + qqbar)
  }
  
  # Return result based on log argument
  if (log) {
    return(log_den)
  } else {
    return(exp(log_den))
  }
}