test_that("complex mean prediction matches independent complex multiplication", {
  # Construct directly rather than rely on helpers in other test files.
  b <- rbind(c(1,2,3,4),c(-1,0,2,-2))
  f <- structure(list(beta=b,model=list(d=2L)),class="cplxshrink_fit")
  X <- rbind(c(1+2i,-1i),c(-1,2+1i))
  expected <- rbind(drop(X%*%c(1+3i,2+4i)),drop(X%*%c(-1+2i,-2i)))
  set.seed(92); rng <- .Random.seed; old <- f
  actual <- predict(f,X)
  expect_equal(unname(actual),unname(expected))
  expect_identical(.Random.seed,rng); expect_identical(f,old)
  sm <- predict(f,X,summary=TRUE,level=.5)
  expect_equal(sm$mean[1:2],colMeans(Re(expected)))
  expect_equal(sm$mean[3:4],colMeans(Im(expected)))
  expect_equal(sm$mean[5:6],colMeans(Mod(expected)))
  f$beta <- matrix(c(1,-1,0,0),2,2)
  m <- predict(f,matrix(1,1,1),summary=TRUE)
  expect_equal(m$mean,c(0,0,1)); expect_equal(m$lower,c(-.95,0,1))
  expect_equal(dim(predict(f,matrix(1,1,1))),c(2L,1L))
  f$beta <- matrix(c(1,2),1,2)
  expect_equal(dim(predict(f,matrix(1,1,1))),c(1L,1L))
  expect_error(predict(f,matrix(1,1,1),summary=TRUE),"at least 2")
})

test_that("predictive covariance matches nonzero relation and real variance", {
  n <- 100000L; s <- 1.7; rho <- .3-.4i
  f <- structure(list(beta=matrix(0,n,2),s=rep(s,n),
    model=list(d=2L,rho=rho)),class="cplxshrink_fit")
  set.seed(4); rng <- .Random.seed; old <- f
  a <- predict(f,matrix(0,1,1),type="predictive",seed=714)
  expect_identical(a,predict(f,matrix(0,1,1),type="predictive",seed=714))
  expect_identical(.Random.seed,rng); expect_identical(f,old)
  target <- s*matrix(c(1.3,-.4,-.4,.7),2,2)
  emp <- stats::cov(cbind(Re(a[,1]),Im(a[,1])))
  # Gaussian covariance standard error sqrt((Sigma_ij^2+Sigma_ii*Sigma_jj)/(N-1)).
  se <- sqrt((target^2+outer(diag(target),diag(target)))/(n-1))
  expect_true(all(abs(emp-target)<6*se))
  expect_true(all(abs(c(mean(Re(a)),mean(Im(a))))<6*sqrt(diag(target)/n)))
  # Separate new observations have independent errors given fixed parameters.
  two <- predict(f,matrix(0,2,1),type="predictive",seed=58)
  expect_lt(abs(cor(Re(two[,1]),Re(two[,2]))),.02)
  f$model <- list(d=1L,rho=0); f$beta <- matrix(0,n,1)
  real <- predict(f,matrix(0,1,1),type="predictive",seed=8)
  expect_false(is.complex(real))
  expect_lt(abs(var(real[,1])-s),6*s*sqrt(2/(n-1)))
  expect_true(all(predict(f,matrix(0,1,1),type="mean")==0))
  # Null seed consumes RNG, seeded calls restore it; variance varies by draw.
  f$s <- seq(.5,2.5,length.out=n)
  v <- predict(f,matrix(0,1,1),type="predictive",seed=71)[,1]/sqrt(f$s)
  expect_lt(abs(var(v)-1),.03)
  predict(f,matrix(0,1,1),type="predictive")
  expect_false(identical(.Random.seed,rng))
})

test_that("prediction rejects invalid inputs and mismatched columns", {
  X <- matrix(1:12,6,2,dimnames=list(NULL,c("a","b")))
  f <- cplx_shrink(1:6,X,model="real",iterations=20,warm=5)
  for(xx in list(1:2,matrix(1,1,1),matrix(NA_real_,1,2),matrix(numeric(),0,2),matrix("a",1,2)))
    expect_error(predict(f,xx),"newdata")
  expect_error(predict(f,X+1i),"real newdata")
  expect_error(predict(f,X[,2:1]),"column names")
  expect_error(predict(f,X,summary=1),"logical")
  expect_error(predict(f,X,level=1),"level")
  expect_error(predict(f,X,seed=-1),"seed")
  expect_error(predict(f,X,type="bad"),"arg")
  f$s[1] <- -1; expect_error(predict(f,X,type="predictive"),"noise scale")
  f$s[] <- 1; f$model$rho <- 1
  expect_error(predict(f,X,type="predictive"),"rho")
  f$beta[1,1] <- Inf; expect_error(predict(f,X),"finite")
})

test_that("summaries and plots work for single and multiple chains headlessly", {
  file <- tempfile(fileext=".pdf"); grDevices::pdf(file)
  on.exit({grDevices::dev.off(); unlink(file)},add=TRUE)
  X <- matrix(seq(-1,1,length.out=12),12,1)
  for(model in c("real","complex")) for(prior in c("student","lasso","gdp")) {
    z <- cplx_shrink_chains(X[,1],X,model=model,prior=prior,iterations=25,warm=5)
    f <- z$fits[[1]]
    ss <- summary(f); ms <- summary(z)
    expect_s3_class(ss,"summary.cplxshrink_fit")
    expect_s3_class(ms,"summary.cplxshrink_chains")
    expect_equal(ms$parameters$mean[1],mean(vapply(z$fits,function(k) mean(k$beta[,1]),numeric(1))))
    expect_output(print(ss),"Not assessed")
    expect_output(print(ms),"no convergence classification")
    saved <- graphics::par("mfrow"); set.seed(21); rng <- .Random.seed
    for(type in c("trace","density")) {
      expect_silent(plot(f,type=type)); expect_silent(plot(z,type=type))
    }
    expect_identical(graphics::par("mfrow"),saved); expect_identical(.Random.seed,rng)
    expect_error(plot(z,parameters="missing"),"parameters")
    expect_error(plot(f,parameters=character()),"parameters")
  }
  expect_gt(file.info(file)$size,0)
  f <- cplx_shrink(1:6,matrix(1:6,6,1),model="real",iterations=2,warm=1)
  expect_error(summary(f),"at least 2")
  expect_error(coefficient_summary(f),"at least 2")
  expect_error(plot(f),"at least 2")
})
