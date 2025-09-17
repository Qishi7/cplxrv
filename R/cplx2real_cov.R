#' Convert Complex Covariance Matrix to Real Covariance Matrix
#' @description Convert a complex covariance matrix to its real representation
#' @param var Complex covariance matrix
#' @param rel Complex correlation matrix (optional)
#' @param rho Complex correlation coefficient (optional)
#' @param is.circular Logical indicating if using circular covariance structure
#' @return Real covariance matrix
#' @examples
#' var <- matrix(c(1+0i, 0.5+0.3i, 0.5-0.3i, 1+0i), 2, 2) 
#' cplx2real_cov(var, is.circular = TRUE)
#' cplx2real_cov(var, rho = 0.5+0.3i)
#' @export
cplx2real_cov <- function(var, rel = NULL, rho = NULL, 
                          is.circular = TRUE) {
  
  ## rho is const across variables. 
  ## For each variable, its real-imag correlation is assumed the same
  
  if (is.circular) {
    return(as.matrix(Matrix::bdiag(var, var))/2)
  }
  if (!is.null(rho)) rel <- var * rho
  
  Vxx <- Re(var + rel)
  Vyx <- Im(var + rel)
  Vxy <- Im(-var + rel)
  Vyy <- Re(var - rel)
  
  return(cbind(rbind(Vxx, Vyx), rbind(Vxy, Vyy))/2)
  
}