set_par <- function (p, prior = "t", is_cplx = FALSE, learn_rho = FALSE,
                     griddy = FALSE) {
  # prior = "t", "bl", "gdp"
  
  ## add stop() message
  
  b_re_vec <- paste("beta_re", 1:p, sep = "_")
  b_im_vec <- paste("beta_im", 1:p, sep = "_")
  b_vec <- paste("beta", 1:p, sep = "_")
  tau2_vec <- paste("tau2", 1:p, sep = "_")
  
  if (prior == "t") {
    if (is_cplx) {
      name_par <- c(b_re_vec, b_im_vec, tau2_vec, "sig2")
      start_lst <- list("beta" = rep(0.5, 2 * p),
                        "tau2" = rep(1, p),
                        "sig2" = 0.5)
      
      if (learn_rho) {
        name_par <- c(name_par, "rho_re", "rho_im")
        start_lst <- c(start_lst, list("rho_re" = 0.5, "rho_im" = 0.1))
      }
      
      
    } else {
      name_par <- c(b_vec, tau2_vec, "sig2") 
      start_lst <- list("beta" = rep(0.5, p),
                        "tau2" = rep(1, p),
                        "sig2" = 1)
    }
    
    
  }
  
  ################
  if (prior == "bl") {
    if (is_cplx) {
      name_par <- c(b_re_vec, b_im_vec, tau2_vec, "lam2", "sig2")
      
      start_lst <- list("beta" = rep(0.5, 2*p),
                        "tau2" = rep(1, p),
                        "lam2" = 1,
                        "sig2" = 0.5)
      
      if (learn_rho) {
        name_par <- c(name_par, "rho_re", "rho_im")
        start_lst <- c(start_lst, list("rho_re" = 0.5, "rho_im" = 0.1))
      }
      
    } else {
      name_par <- c(b_vec, tau2_vec, "lam2", "sig2") 
      start_lst <- list("beta" = rep(0.5, p),
                        "tau2" = rep(1, p),
                        "lam2" = 1,
                        "sig2" = 1)
    }
  }
  
  
  ################
  if (prior == "gdp") {
    if (is_cplx) {
      start_lst <- list("beta" = rep(0.5, 2*p),
                        "tau2" = rep(1, p),
                        "lam" = rep(1, p),
                        "sig2" = 0.5)
      name_par <- c(b_re_vec, b_im_vec, tau2_vec,
                    paste("lam", 1:p, sep = "_"),
                    "sig2")
      
      # if (!griddy) {
      # 
      # } else {
      #     name_par <- c(name_par, "alpha", "eta")
      # }
      
      if (griddy) {
        name_par <- c(name_par, "alpha", "eta")
      }
      
      if (learn_rho) {
        name_par <- c(name_par, "rho_re", "rho_im")
        start_lst <- c(start_lst, list("rho_re" = 0.5, "rho_im" = 0.1))
      }
      
      
    } else {
      start_lst <- list("beta" = rep(0.5, p),
                        "invtau" = rep(1, p),
                        "lam" = rep(1, p),
                        "sig2" = 1)
      name_par <- c(b_vec,
                    paste("invtau", 1:p, sep = "_"),
                    paste("lam", 1:p, sep = "_"), 
                    "sig2")
      # if (!griddy) {
      # 
      # } else {
      #     name_par <- c(name_par, "alpha", "eta")
      # }
      
      if (griddy) {
        name_par <- c(name_par, "alpha", "eta")
      }
      
    }
  }
  
  return(list(name_par = name_par, start_lst = start_lst,
              number_par = length(name_par)))
}
