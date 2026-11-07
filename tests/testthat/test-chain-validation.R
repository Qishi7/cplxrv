test_that("multi-coefficient chains expose exact variables and reject inconsistent structure", {
  X <- matrix(seq(-1,1,length.out=24),8,3)
  for(model in c("real","complex")) for(prior in c("student","lasso","gdp")) {
    z <- cplx_shrink_chains(1:8,X,prior=prior,model=model,chains=2,iterations=14,warm=2)
    a <- as_draws_array(z)
    beta <- if(model=="complex") c(paste0("beta_re",1:3),paste0("beta_im",1:3)) else paste0("beta",1:3)
    h <- switch(prior,student=character(),lasso="L",gdp=paste0("lambda",1:3))
    expect_identical(dimnames(a)[[3]],c(beta,"s",paste0("w",1:3),h,paste0("magnitude",1:3)))
    expect_equal(dim(a),c(12,2,length(beta)+7+length(h)))
    expect_equal(mcmc_diagnostics(z)$variable,dimnames(a)[[3]])
    bad <- z; bad$fits[[2]]$sampleidx[1] <- 0L
    expect_error(as_draws_array(bad),"inconsistent")
    bad <- z; bad$fits[[2]]$thin <- 2L
    expect_error(as_draws_array(bad),"inconsistent")
    bad <- z; bad$fits[[2]]$s[1] <- NaN
    expect_error(as_draws_array(bad),"scale")
    bad <- z; colnames(bad$fits[[2]]$beta)[1] <- "changed"
    expect_error(as_draws_array(bad),"inconsistent")
    if(prior!="student") {
      bad <- z; bad$fits[[2]]$hyper <- matrix(numeric(),12,0)
      expect_error(as_draws_array(bad),"hyperparameter")
    }
  }
})
