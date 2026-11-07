# Deterministic software-use validation; these short chains do not establish convergence.
library(cplxshrink)
set.seed(104)
X <- matrix(rnorm(60),20,3)
colnames(X) <- c("signal","noise1","noise2")
y <- drop(X %*% c(1+1i,0,0)) + rnorm(20) + 1i*rnorm(20)
fit <- cplx_shrink(y,X,prior="gdp",relation=.2+.3i,iterations=80,warm=20,thin=2)
mag <- magnitude_summary(fit)
stopifnot(nrow(mag)==3L, all(mag$mean>=0),all(mag$width>=0),all(mag$n_draws==30))
z <- cplx_shrink_chains(y,X,prior="gdp",relation=.2+.3i,
  chains=4,seed=201,iterations=80,warm=20,thin=2,
  starts=lapply(c(.2,.5,1,2),function(s) list(s=s)))
a <- as_draws_array(z); dg <- mcmc_diagnostics(z)
stopifnot(identical(dim(a),c(30L,4L,16L)), nrow(dg)==16L,
          all(c("rhat","ess_bulk","ess_tail","mcse_mean") %in% names(dg)))
new <- X[1:2,,drop=FALSE]
mu <- predict(fit,new,type="mean")
yp <- predict(fit,new,type="predictive",seed=77)
b <- fit$beta[,1:3,drop=FALSE]+1i*fit$beta[,4:6,drop=FALSE]
stopifnot(isTRUE(all.equal(unname(mu),unname(b%*%t(new)))),
  all(is.finite(yp)),!identical(mu,yp),
  identical(yp,predict(fit,new,type="predictive",seed=77)))
print(mag)
print(head(dg))
print(predict(fit,new,type="predictive",seed=77,summary=TRUE))
print(summary(fit))
cat("Validation demo completed. Short chains do not establish convergence.\n")
