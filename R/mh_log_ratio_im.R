mh_log_ratio_im <- function(rho_re, eta_im, eta_im_star, mu, yr, sig2) {
  # rho_re <- 1 / (1 + exp(-eta_re))
  a <- sqrt(1 - rho_re ^ 2)
  rho_im <- eta_to_rho(eta_im, a = a)
  rho_im_star <- eta_to_rho(eta_im_star, a = a)
  
  # mu <- Xr %*% betar
  I <- diag(length(yr) / 2)
  # n <- length(yr) / 2
  Vrr <- (1 + rho_re) * I
  Vii <- (1 - rho_re) * I
  Vri <- rho_im * I
  Vri_star <- rho_im_star * I
  Sig <- sig2 * rbind(cbind(Vrr, Vri), cbind(Vri, Vii))
  Sig_star <- sig2 * rbind(cbind(Vrr, Vri_star), cbind(Vri_star, Vii))
  
  A_star <- mvnfast::dmvn(yr, mu = mu, sigma = Sig_star, log = TRUE) + 
    eta_im_star + log(2 * a) - 2 * log(1 + exp(eta_im_star))
  
  A <- mvnfast::dmvn(yr, mu = mu, sigma = Sig, log = TRUE) + 
    eta_im + log(2 * a) - 2 * log(1 + exp(eta_im))
  
  return(A_star - A)
}