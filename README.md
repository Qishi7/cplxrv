# cplxshrink

Bayesian shrinkage regression for real and complex-valued data in R.
Complex coefficients share a local shrinkage scale across their real and
imaginary components.

- Student-t, Bayesian lasso and generalized double Pareto priors.
- Posterior coefficient and magnitude summaries.
- Reproducible multiple chains with R-hat, effective sample size and Monte Carlo error diagnostics.
- Conditional-mean and posterior-predictive draws.
- Summary, trace and density plotting methods.

## Installation

Requires R >= 4.1 and a C++ toolchain.

```r
install.packages("remotes")
remotes::install_github("Qishi7/cplxrv")
library(cplxshrink)
```

## Quick start

```r
set.seed(42)
n <- 100
p <- 5
X <- matrix(rnorm(n * p), n, p)
beta <- c(1 + 1i, -0.5 + 0.5i, 0, 0, 0)
y <- drop(X %*% beta) + rnorm(n) + 1i * rnorm(n)

fit <- cplx_shrink(y, X, prior = "gdp")
coef(fit)
coefficient_summary(fit)
magnitude_summary(fit)
summary(fit)
```

Choose `prior = "student"`, `"lasso"` or `"gdp"`. For scalar real regression,
use `model = "real"` with a real response and design matrix.

`magnitude_summary()` summarizes the magnitude of each posterior coefficient
draw. This differs from taking the magnitude of the posterior mean coefficient.

## Multiple-chain diagnostics

```r
chains <- cplx_shrink_chains(
  y, X, prior = "gdp", chains = 4, seed = 123,
  iterations = 2000, warm = 500,
  starts = lapply(c(0.2, 0.5, 1, 2), function(s) list(s = s))
)

mcmc_diagnostics(chains)
summary(chains)
plot(chains, parameters = c("beta_re1", "s"), type = "trace")
plot(chains, parameters = "beta_re1", type = "density")
```

Diagnostics include rank-normalized R-hat, bulk and tail ESS, and MCSE of
posterior means, computed with the `posterior` package. Inspect multiple chains
and Monte Carlo error before drawing substantive conclusions.

## Prediction

```r
newX <- X[1:2, , drop = FALSE]

# Conditional means: coefficient uncertainty.
mean_draws <- predict(fit, newX, type = "mean")

# New responses: coefficient uncertainty and observation noise.
response_draws <- predict(fit, newX, type = "predictive", seed = 45)
predict(fit, newX, type = "predictive", seed = 45, summary = TRUE)
```

Supply an unstacked design matrix in the original column order. Complex
prediction summaries include real, imaginary and draw-wise magnitude summaries.

## Model conventions

- No intercept or preprocessing is added automatically. An explicitly included intercept is shrunk.
- Complex likelihood relation `rho` is fixed, with `Mod(rho) < 1`.
- Noise scale `s` gives total complex error variance `2*s`; real-model error variance is `s`.
- Coefficient and prediction intervals are marginal intervals.
- Package defaults differ from the manuscript's matched-prior calibration.

See the [model conventions](docs/model-conventions.md),
[manuscript calibration](docs/manuscript-calibration.md) and `?cplx_shrink`
for details on covariance, priors, draw indexing and prediction.

## Examples and local checks

Run the installed feature demonstration:

```r
source(system.file("examples", "validate_new_features.R",
                   package = "cplxshrink", mustWork = TRUE))
```

From the repository root, run the tests with:

```r
testthat::test_local(".")
```
