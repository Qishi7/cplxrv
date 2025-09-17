get_V_from_eps <- function (eps, sig2) {
  n <- length(eps) / 2
  Vrr <- var(eps[1:n])/sig2
  Vii <- var(eps[-c(1:n)])/sig2
  Vri <- cov(eps[1:n], eps[-c(1:n)])/sig2
  # rho_re <- (Vrr + Vii)/2
  V <- rbind(cbind(Vrr * diag(n), Vri * diag(n)),
             cbind(Vri * diag(n), Vii * diag(n)))
  return(list(V = V, rho_ep = (Vrr - Vii)/2 + Vri * 1i))
  
}