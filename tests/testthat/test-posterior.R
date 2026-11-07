make_fit <- function(beta, d=2L, rho=0, s=rep(1,nrow(beta)), family="student") {
  p <- ncol(beta)/d; n <- nrow(beta)
  colnames(beta) <- if(d==2L) c(paste0("beta_re",seq_len(p)),paste0("beta_im",seq_len(p))) else paste0("beta",seq_len(p))
  w <- matrix(1,n,p,dimnames=list(NULL,paste0("w",seq_len(p))))
  h <- matrix(numeric(),n,0)
  structure(list(beta=beta,s=s,w=w,hyper=h,model=list(d=d,rho=rho,family=family),
    sampleidx=seq_len(n), iterations=n, warm=0L,thin=1L),class="cplxshrink_fit")
}

test_that("magnitude uses draw-wise modulus with independent hand oracles", {
  f <- make_fit(cbind(c(1,-1),c(0,0)))
  set.seed(19); rng <- .Random.seed; saved <- f
  m <- magnitude_summary(f)
  expect_equal(m$mean,1); expect_equal(Mod(mean(c(1+0i,-1+0i))),0)
  expect_equal(m$sd,0); expect_equal(m$lower,1); expect_equal(m$width,0)
  expect_identical(f,saved); expect_identical(.Random.seed,rng)
  # Magnitudes are exactly (5, 0, 13) and (0, 2, 4).
  f <- make_fit(rbind(c(3,0,4,0), c(0,0,0,2), c(5,4,12,0)))
  f$design <- list(columns=c("a","b"))
  m <- magnitude_summary(f,level=.5)
  expect_equal(m$mean,c(6,2)); expect_equal(m$sd,c(sqrt(43),2))
  expect_equal(m$lower,c(2.5,1)); expect_equal(m$upper,c(9,3))
  expect_equal(m$width,c(6.5,2)); expect_equal(m$level,c(.5,.5))
  expect_equal(m$n_draws,c(3L,3L)); expect_equal(m$coefficient,c("a","b"))
  f <- make_fit(matrix(0,2,2)); expect_equal(magnitude_summary(f)$mean,0)
})

test_that("magnitude rejects unsupported and malformed inputs", {
  f <- make_fit(matrix(1,2,2))
  for(l in list(0,1,NA,Inf,1i,"a",c(.5,.9))) expect_error(magnitude_summary(f,l),"level")
  expect_error(magnitude_summary(list()),"fit")
  expect_error(magnitude_summary(make_fit(matrix(1,2,1),d=1L)),"complex model")
  expect_error(magnitude_summary(make_fit(matrix(1,1,2))),"at least 2")
  f$beta[1,1] <- Inf; expect_error(magnitude_summary(f),"finite")
  f$beta <- matrix(numeric(),2,0); expect_error(magnitude_summary(f),"nonempty")
  f$beta <- matrix(1,2,3); expect_error(magnitude_summary(f),"nonempty")
})

test_that("streams reproduce fits and preserve caller RNG and iteration ordering", {
  X <- matrix(seq(-1,1,length.out=12),12,1); y <- X[,1]+.1
  for(model in c("real","complex")) for(prior in c("student","lasso","gdp")) {
    yy <- if(model=="complex") y+1i*rev(y) else y
    set.seed(33); saved <- .Random.seed; kind <- RNGkind()
    a <- cplx_shrink_chains(yy,X,model=model,prior=prior,chains=2,seed=81,
      starts=list(list(s=.3),list(s=2)),iterations=30,warm=3,thin=2)
    b <- cplx_shrink_chains(yy,X,model=model,prior=prior,chains=2,seed=81,
      starts=list(list(s=.3),list(s=2)),iterations=30,warm=3,thin=2)
    expect_identical(a,b); expect_identical(.Random.seed,saved); expect_identical(RNGkind(),kind)
    expect_false(identical(a$streams[[1]],a$streams[[2]]))
    expect_identical(a$streams[[2]],parallel::nextRNGStream(a$streams[[1]]))
    expect_false(identical(a$fits[[1]]$beta,a$fits[[2]]$beta))
    expect_equal(a$sampleidx,seq(5,29,2))
    ar <- as_draws_array(a)
    d <- if(model=="complex") 2 else 1
    h <- switch(prior, student=0,lasso=1,gdp=1)
    expect_equal(dim(ar),c(13,2,d+3+h))
    expect_equal(as.numeric(ar[,1,1]), a$fits[[1]]$beta[,1])
    if(model=="complex") expect_true(all(is.finite(magnitude_summary(a$fits[[1]])$mean)))
    # A stream replay calls the original interface and reproduces a chain exactly.
    replay <- cplxshrink:::.with_rng(function() {
      do.call(RNGkind,as.list(a$rng_kind)); assign(".Random.seed",a$streams[[1]],.GlobalEnv)
      cplx_shrink(yy,X,model=model,prior=prior,start=list(s=.3),iterations=30,warm=3,thin=2)
    })
    expect_identical(replay,a$fits[[1]])
  }
})

test_that("diagnostics match independently constructed posterior arrays", {
  set.seed(123)
  raw <- array(rnorm(200*4*2),c(200,4,2))
  fits <- lapply(1:4,function(k) make_fit(raw[,k,]))
  x <- structure(list(fits=fits),class="cplxshrink_chains")
  # Assemble all variables without calling any new summary/conversion helper.
  ref <- array(NA_real_, c(200,4,5), dimnames=list(NULL,NULL,
    c("beta_re1","beta_im1","s","w1","magnitude1")))
  for(k in 1:4) ref[,k,] <- cbind(raw[,k,1],raw[,k,2],1,1,
                                  sqrt(raw[,k,1]^2+raw[,k,2]^2))
  oracle <- posterior::as_draws_array(ref)
  expect_equal(unname(as.array(as_draws_array(x))),unname(as.array(oracle)),ignore_attr=TRUE)
  dg <- mcmc_diagnostics(x)
  expected <- posterior::summarise_draws(oracle, rhat=posterior::rhat,
    ess_bulk=posterior::ess_bulk, ess_tail=posterior::ess_tail, mcse_mean=posterior::mcse_mean)
  for(j in c(1,2,5)) {
    expect_equal(dg$rhat[j],expected$rhat[j])
    expect_equal(dg$ess_bulk[j],expected$ess_bulk[j])
    expect_equal(dg$ess_tail[j],expected$ess_tail[j])
    expect_equal(dg$mcse_mean[j],expected$mcse_mean[j])
  }
  expect_true(all(is.na(dg$rhat[3:4]))); expect_true(all(nzchar(dg$note[3:4])))
  for(k in 1:4) x$fits[[k]]$beta[,1] <- raw[,k,1]+k*10
  expect_gt(mcmc_diagnostics(x)$rhat[1],1.5)
  one <- mcmc_diagnostics(fits[[1]]); expect_true(all(is.na(one$rhat)))
  expect_true(all(grepl("two chains",one$note)))
  short <- structure(list(fits=lapply(1:2,function(k) make_fit(matrix(1,3,2)))),class="cplxshrink_chains")
  expect_true(all(is.na(mcmc_diagnostics(short)$ess_bulk)))
  x$fits[[2]]$model$rho <- .2; expect_error(as_draws_array(x),"inconsistent")
})

test_that("chain validation errors restore RNG including previously absent state", {
  X <- matrix(1:6,6,1)
  for(n in list(0,NA,1.5,1i,"four",Inf)) expect_error(cplx_shrink_chains(1:6,X,chains=n),"chains")
  for(s in list(-1,NA,Inf,1.5,1i)) expect_error(cplx_shrink_chains(1:6,X,seed=s),"seed")
  expect_error(cplx_shrink_chains(1:6,X,starts=list(NULL)),"starts")
  set.seed(19); old <- .Random.seed
  expect_error(cplx_shrink_chains(1:6,X,relation=1),"rho")
  expect_identical(.Random.seed,old)
  rm(".Random.seed",envir=.GlobalEnv)
  on.exit(assign(".Random.seed",old,.GlobalEnv))
  z <- cplx_shrink_chains(1:6,X,model="real",iterations=8,warm=2)
  expect_false(exists(".Random.seed",.GlobalEnv,inherits=FALSE))
  z$fits[[2]]$w <- matrix(1,6,2); expect_error(as_draws_array(z),"scale")
  expect_error(as_draws_array(list()),"chains")
})
