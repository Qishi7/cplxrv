#' Convert Real Matrix to Complex Variance Matrix
#' @description Convert a real matrix to complex variance matrix with optional block structure
#' @param V Full real matrix (optional)
#' @param Vx Top-left block matrix (optional)
#' @param Vy Bottom-right block matrix (optional)
#' @param Vxy Top-right block matrix (optional)
#' @param Vyx Bottom-left block matrix (optional)
#' @param is.circular Logical indicating if using circular structure
#' @return List containing complex variance matrix (var) and relation matrix (rel)
#' @examples
#' # Example 1: Using full matrix with circular structure
#' V <- matrix(c(1, 0.5, 0, 0,
#'               0.5, 1, 0, 0,
#'               0, 0, 1, 0.5,
#'               0, 0, 0.5, 1), 4, 4)
#' result1 <- real2cplx_var(V, is.circular = TRUE)
#'
#' # Example 2: Using block matrices
#' Vx <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' Vy <- matrix(c(1, -0.3, -0.3, 1), 2, 2)
#' Vxy <- matrix(c(0.2, 0.1, 0.1, 0.2), 2, 2)
#' Vyx <- matrix(c(0.2, -0.1, -0.1, 0.2), 2, 2)
#' V <- rbind(cbind(Vx, Vxy),cbind(Vyx, Vy))
#' result2 <- real2cplx_var(V, Vx, Vy, Vxy, Vyx)
#' @export
real2cplx_var <- function(V, Vx, Vy, Vxy, Vyx, is.circular = FALSE) {
  p <- ncol(V)/2
  idx <- 1:p
  
  if(missing(Vx)) {
    Vx <- V[idx, idx]       
  }
  
  if (is.circular) {
    var <- 2 * Vx
    rel <- 0
    return(list(var = var, rel = rel))
  }
  
  if(missing(Vy)) {
    Vy <- V[-idx, -idx]
  }
  
  if(missing(Vxy)) {
    Vy <- V[idx, -idx]
  }
  
  if(missing(Vyx)) {
    Vy <- V[-idx, idx]
  }
  
  var <- Vx + Vy + (Vyx - Vxy) * 1i
  rel <- Vx - Vy + (Vyx + Vxy) * 1i
  
  return(list(var = var, rel = rel))
  
}