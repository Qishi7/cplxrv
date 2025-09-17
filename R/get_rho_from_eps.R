get_rho_from_eps <- function (eps) {
  n <- length(eps) / 2
  Vrr <- var(eps[1:n])
  Vii <- var(eps[-c(1:n)])
  rho_im <- cov(eps[1:n], eps[-c(1:n)])
  rho_re <- (Vrr + Vii)/2
  return(rho_re + rho_im * 1i)
}