#' Convert Complex Variance to Real Variance
#' @description Convert a complex variance to its real representation
#' @param var Complex variance matrix
#' @param rel Complex correlation matrix (optional)
#' @param rho Complex correlation coefficient (optional)
#' @param is.circular Logical indicating if using circular variance structure
#' @return Real variance matrix
#' @details 
#' For circular structure, returns a block diagonal matrix scaled by 1/2 using Matrix::bdiag.
#' For non-circular structure, constructs variance matrix using real and imaginary parts.
#' If rho is provided but rel is NULL, rel is computed as var * rho.
#' @examples
#' # Create a 2x2 complex variance matrix
#' var <- matrix(c(1, 0.5, 0.5, 2), 2, 2) + 
#'       1i * matrix(c(0, 0.3, -0.3, 0), 2, 2)
#' 
#' # Circular case
#' result1 <- cplx2real_variance(var, is.circular = TRUE)
#' 
#' # Non-circular case with rho
#' result2 <- cplx2real_variance(var, rho = 0.5 + 0.3i)
#' @export

cplx2real_variance <- function(var, rel = NULL, rho = NULL, 
                               is.circular = TRUE) {
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
