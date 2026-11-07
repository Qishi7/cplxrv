# cplxshrink

Bayesian shrinkage regression for complex data with a shared local scale for
each real/imaginary coefficient pair. Scalar real regression is also available.
No intercept is added automatically; include a column of ones if needed
(the intercept coefficient is then subject to the same shrinkage prior).

Development version **0.1.0.9000**, for review; this is not an official v0.2 release.
The frozen manuscript baseline is upstream commit
`d985957c4254593ae2cc5c449d6dc3d108dd1b5b` (0.1.0).
All four certified sampler/binding files and reference fixtures are unchanged.
This development package is hosted in Qishi7/cplxrv at the owner's request.
The repository previously held the older cplxrv package. This replacement installs
as **cplxshrink**: use `library(cplxshrink)`, not the former `library(cplxrv)`.
The previous repository contents remain recoverable from commit
`7ef7f9e3decb6294e6774d08d9444768f6ad9f58`. This is an existing repository,
not a GitHub fork of the upstream manuscript package.
Historical local validation reports describe the earlier, pre-upload handoff;
see PUBLISHING_REPORT.md for this repository replacement.

## Install

A working R installation (>= 4.1) and C++ toolchain are required.
In R, install from GitHub:

```r
install.packages("remotes")  # if needed
remotes::install_github("Qishi7/cplxrv")
library(cplxshrink)
```

For source installation, first install `Rcpp`, `RcppArmadillo` and `posterior`.

From a terminal, install the supplied source archive:

```sh
R CMD INSTALL cplxshrink_0.1.0.9000.tar.gz
```

Alternatively, run `R CMD INSTALL cplxshrink` from the directory containing
this source folder. To run its tests, also install `testthat` and `jsonlite`. Install `posterior` for diagnostics.

## Fit and summarize

```r
library(cplxshrink)
set.seed(42)
n <- 100
p <- 5
X <- matrix(rnorm(n * p), n, p)
beta <- c(1 + 1i, -0.5 + 0.5i, 0, 0, 0)
y <- drop(X %*% beta) + rnorm(100) + 1i * rnorm(100)
fit <- cplx_shrink(y, X, prior = "gdp")
coef(fit)
coefficient_summary(fit)
```

`coef(fit)` returns posterior mean coefficients, complex for a complex fit
and real for a real fit. `coefficient_summary()` gives means, standard
deviations, and marginal central intervals for each real and imaginary
coordinate. Draws remain available in `fit$beta` (real columns followed by
imaginary columns), `fit$w`, `fit$s`, and `fit$hyper`.

Choose `prior = "student"` (Student-t, the default), `"lasso"` (Bayesian
lasso), or `"gdp"` (generalized double Pareto). Use `model = "real"` for
real data. Optional `relation`, `hyper`, `iterations`, `warm`, `thin`, and
`start` arguments control fixed error relation, prior hyperparameters, total
iterations, warmup, thinning, and initial values; see `?cplx_shrink` for details. Set a seed immediately before fitting for repeatability.
After installation, run the bundled demonstration:

```r
source(system.file("examples", "demo.R", package = "cplxshrink", mustWork = TRUE))
```

## Manuscript calibration (separate from the generic demo)

The primary simulation comparison uses real `X`, proper errors (`relation = 0`),
and `a_s = b_s = 0.5` for each error variance. Specify the following `hyper`
lists explicitly; package defaults are unchanged.

| Prior | Paired complex fit (`model = "complex"`) | Each scalar fit (`model = "real"`) |
| --- | --- | --- |
| Student-t | `list(a_s=0.5, b_s=0.5, nu=1)` | `list(a_s=0.5, b_s=0.5, a_tau=1, b_tau=1)` |
| Lasso | `list(a_s=0.5, b_s=0.5, r=0.5, delta=0.5)` | `list(a_s=0.5, b_s=0.5, r=0.5, delta=0.8966996574)` |
| GDP | `list(a_s=0.5, b_s=0.5, alpha=1, eta=1)` | `list(a_s=0.5, b_s=0.5, alpha=1, eta=1.3403844863)` |

For example, using `y` and real `X` from the generic demo, the calibrated
independent scalar GDP fits use:

```r
paper_gdp <- list(a_s=0.5, b_s=0.5, alpha=1, eta=1.3403844863)
fit_R <- cplx_shrink(Re(y), X, prior="gdp", model="real", hyper=paper_gdp)
fit_I <- cplx_shrink(Im(y), X, prior="gdp", model="real", hyper=paper_gdp)
```

M-ReIm requires two independent fits. M-Scalar-Pooled instead requires an
explicit block-diagonal real design and stacked response, sharing one error
variance and, for lasso, one global rate. A single ordinary fitting call does
not construct these simulation strategies. These parameter examples do not
supply the paper's datasets, seed mapping, four-chain schedules, diagnostics,
scoring, or full ablation and fMRI workflows. For the H1 sensitivity only,
the scalar lasso uses `delta=0.5` and scalar GDP uses `eta=1`.

## Parameter notation

The package's noise variance `s` is the manuscript's sigma-squared (\(\sigma^2\)),
the variance of each real error component, not total complex variance.
`a_s` and `b_s` correspond to \(a_\sigma\) and \(b_\sigma\), with
`s ~ IG(a_s, b_s)` and density proportional to `s^(-a_s-1)*exp(-b_s/s)`.
`w` denotes the local scale. Complex Student `nu` is \(\nu\); scalar Student
`a_tau` and `b_tau` are the shape and scale of the inverse-gamma local scale
(`a_tau=b_tau=1` gives the paper's scalar t with two degrees of freedom).
Lasso `delta` is \(\delta_L\), the rate in `L ~ Gamma(r, delta)`;
GDP `alpha` and `eta` are the shape and rate in `lambda[j] ~ Gamma(alpha, eta)`.

## Limits

The likelihood relation is fixed, with modulus below one; coefficient priors
are proper. The noise scale `s` is real-component variance, so complex noise
covariance is `2*s*I`. Gamma distributions use shape and rate. Hyperparameter
defaults depend on the model; they do not automatically implement the paper's
matched-prior comparisons. The model conventions are detailed in the function help.

A short run demonstrates code use; it does not establish MCMC convergence or
reproduce every manuscript simulation or the complete fMRI analysis. The
package does not implement likelihood-relation learning or the full structural
ablation workflows. Assess convergence and Monte
Carlo error with suitable multiple-chain diagnostics for substantive inference.
Marginal coefficient intervals are not joint selection regions. This minimal
package contains no human data, empirical outputs, or full replication workflows.

## Draw-wise magnitudes and four-chain diagnostics

```r
magnitude_summary(fit, level = 0.9)
# mean(Mod(beta draw)), not Mod(mean(beta draw)):
# draws 1+0i and -1+0i have mean magnitude 1, but mean coefficient 0.
chains <- cplx_shrink_chains(y, X, prior="gdp", chains=4, seed=123,
  iterations=1000, warm=500, thin=1,
  starts=lapply(c(0.2, 0.5, 1, 2), function(s) list(s=s)))
mcmc_diagnostics(chains)
dim(as_draws_array(chains))  # retained iterations x chains x variables
summary(chains)
plot(chains, parameters=c("beta_re1", "s"), type="trace")
plot(fit, parameters="beta_re1", type="density")
```

Chains run serially using separate L'Ecuyer-CMRG streams, whose initial states
are stored in `chains$streams`. The wrapper restores caller RNG kind and state.
`posterior` supplies rank-normalized split/folded R-hat, bulk and tail ESS, and
MCSE of means. Variables include coefficient coordinates, s, w, prior-specific
rates and draw-wise magnitudes (absolute coefficients for real fits).
Undefined statistics return NA with an explanation. Fewer than two chains or
four draws per chain return NA for all diagnostics. Summaries pool draws only
for descriptive marginal summaries. A short run or favorable diagnostics alone
does not establish convergence. Single-chain summaries report it as unassessed.
Memory scales with draws x chains x variables; array conversion makes a copy.
Plots default to four parameters and allow at most sixteen per call.

## Conditional means and predictive responses

```r
newX <- X[1:2, , drop=FALSE]
mean_draws <- predict(fit, newX, type="mean")
response_draws <- predict(fit, newX, type="predictive", seed=45)
predict(fit, newX, type="predictive", seed=45, summary=TRUE, level=0.9)
summary(fit)
```

`newdata` is an original unstacked real/complex matrix in the original design
column order; names must match when both matrices are named. There is no automatic
intercept or preprocessing. Raw results are draws x observations. Complex
summaries report real, imaginary and draw-wise magnitude means, SDs and marginal
central intervals. These are not joint credible sets or activation classifications.
Mean prediction reflects coefficient uncertainty; predictive responses additionally
include independent observation errors conditional on each retained draw.
For fixed `rho=a+iu`, their real-coordinate covariance is
`s * matrix(c(1+a,u,u,1-a),2,2)`. Thus `s` is the noise component scale, marginal
variances are `s*(1+a)` and `s*(1-a)`, and total complex variance is `2*s`.
Real errors have variance `s`. Seeded prediction restores caller RNG state;
`seed=NULL` advances it for predictive draws. Mean prediction uses no randomness.
No user transformations are inverted and no fMRI original-scale prediction is claimed.

Summaries require at least two finite retained draws; raw predictions support one.
`magnitude_summary` requires a complex fit. Invalid levels, dimensions, nonfinite
values and malformed scales produce errors. `fit$design` stores n, p and column
names; original data are not retained. Other fitted fields and coefficient indexing
are described above and in `?cplx_shrink`.

After installing the package, run from the source folder:

```sh
Rscript examples/validate_new_features.R
```

The deterministic synthetic demonstration checks all required interfaces. It is
software validation, not full manuscript reproduction or a convergence study.
CI is configured for release R on Linux, macOS and Windows plus R 4.1 on Linux;
only the local platform results in the validation reports have been executed here.
