# Model conventions

## Design and coefficient draws

Supply the original, unstacked response and design matrix. Both real and complex
design matrices are supported for complex regression; real regression requires
real data. No intercept, centering, scaling or whitening is added automatically.
An intercept column supplied by the user is subject to the same shrinkage prior.

For a complex fit with p coefficients, `fit$beta` is a retained-draws by 2p matrix:
real coefficient columns first, imaginary coefficient columns second. For real
fits it has p columns. `fit$w` stores local scales, `fit$s` stores the noise scale,
and `fit$hyper` stores prior-specific rates. `fit$sampleidx` gives original retained
iteration indices. `fit$design` records the design size and column names.

## Error covariance

The relation parameter `rho = a + i*u` is fixed and must satisfy `Mod(rho) < 1`.
For each observation, real and imaginary error coordinates have covariance

```r
s * matrix(c(1 + a, u, u, 1 - a), 2, 2)
```

The marginal variances are `s*(1+a)` and `s*(1-a)`, and total complex variance
is `2*s`. For a real model, the error variance is `s` and relation must be zero.
The package does not estimate relation or spatial/temporal error dependence.

## Prior parameterization

Gamma distributions use shape and rate. An inverse-gamma distribution with
shape a and scale b has density proportional to `x^(-a-1)*exp(-b/x)`.

- `a_s`, `b_s`: inverse-gamma parameters for the noise scale s.
- Complex Student-t: `nu` controls the inverse-gamma local scales.
- Real Student-t: `a_tau`, `b_tau` are inverse-gamma local-scale parameters.
- Lasso: `L ~ Gamma(r, delta)`; `delta` is a rate.
- GDP: `lambda[j] ~ Gamma(alpha, eta)`; `eta` is a rate.

See `?cplx_shrink` for the full hierarchy and defaults. Package defaults do not
automatically implement the [manuscript's matched-prior calibration](manuscript-calibration.md).

## Summaries and diagnostics

Magnitude summaries use the modulus of each paired complex coefficient draw,
then summarize those moduli. They differ from the modulus of the posterior mean:
draws +1 and -1 have mean magnitude 1 and mean coefficient 0.
Central intervals use R type-7 quantiles. All real, imaginary and magnitude
intervals are marginal; they are not joint credible regions or activation rules.
Summaries require at least two finite retained draws.

Multiple chains run serially with reproducible L'Ecuyer-CMRG streams and restore
the caller's RNG state. Diagnostics use posterior's rank-normalized split/folded
R-hat, bulk/tail ESS and MCSE without pooling the chain axis. Undefined statistics
and posterior warnings are explained in the returned note column. At least two
chains and four retained draws per chain are required by the diagnostics wrapper.
Short examples and passing software tests do not establish MCMC convergence.

## Prediction

`predict(fit, newdata, type="mean")` returns conditional linear-predictor draws.
`type="predictive"` additionally samples independent observation errors using
each draw's s and fixed rho. Raw results are draws by new observations.
Complex summaries report separate real, imaginary and draw-wise magnitude summaries.

Supply new data in the original design-column order. Seeded predictive calls
restore caller RNG state; unseeded predictive calls advance it. No transformations
are automatically inverted. Predictions on an original measurement scale require
the user to undo any preprocessing.

## Storage

Fits retain draws in memory. Multiple-chain storage and conversion scale with
retained draws, chains and variables; conversion creates a draws-array copy.
Prediction stores retained draws by new observations. Plot methods default to
four variables and accept at most sixteen per call.
