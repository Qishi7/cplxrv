# Short code-use demonstration; assess convergence separately.
library(cplxshrink)
set.seed(42)
n <- 100
p <- 5
X <- matrix(rnorm(n * p), n, p)
beta <- c(1 + 1i, -0.5 + 0.5i, 0, 0, 0)
y <- drop(X %*% beta) + rnorm(n) + 1i * rnorm(n)
fit <- cplx_shrink(y, X, prior = "gdp")
print(coef(fit))
print(coefficient_summary(fit))
