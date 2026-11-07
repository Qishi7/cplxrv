# Checkpoint 0: frozen baseline

Recorded before feature edits: upstream main/HEAD is
`d985957c4254593ae2cc5c449d6dc3d108dd1b5b`, authored/committed
2026-11-04 22:58:25 -0600 (2026-11-05 UTC), package 0.1.0.
The expected commit has not advanced. Local branch: feature/posterior-interface-v02.
Upstream fetch URL: https://github.com/chenghanyustats/cplxshrink-package;
push URL: DISABLED. No upstream writes made.

## Access deviation
GitHub connector authenticated as Qishi7 but exposes no fork operation.
GitHub browser is signed out; gh is unavailable. No fork was created and
no push is claimed. Following the supplied fallback, work is on a local branch;
final handoff includes source ZIP, patch and exact fork/push instructions.
No AI-owned GitHub account is available.

## Reproduction and results
From the parent directory, a pristine `git archive HEAD` was extracted to baseline.
- `Rscript -e 'testthat::test_local(".", reporter="summary")'` in unchanged checkout: passed.
- `LC_ALL=en_US.UTF-8 Rscript -e 'testthat::test_local("baseline", reporter="summary")'`:
  8 test blocks, 151 expectations, 0 failures/errors/warnings/skips; counted log and CSV saved.
- `R CMD build baseline`: succeeded, cplxshrink_0.1.0.tar.gz.
- `R CMD check --as-cran cplxshrink_0.1.0.tar.gz`: initial environment failure,
  1 error, 1 warning, 2 notes. Invalid inherited C.UTF-8 locale caused R startup
  warnings in DESCRIPTION checks; incoming-submission and current-time notes also present.
- `LC_ALL=en_US.UTF-8 R CMD check --output=baseline-standard cplxshrink_0.1.0.tar.gz`:
  Status OK, 0 errors, 0 warnings, 0 notes, including PDF manual and both test scripts.
- UTF-8 --as-cran rerun is recorded in review/logs/baseline-cran-utf8.log;
  its completed results will be included in final validation.

Full logs: review/logs/baseline-*.log; exact test counts:
review/logs/baseline-test-counts.csv; environment and remote/status snapshot:
review/logs/environment.log. R 4.4.1, aarch64-apple-darwin20, macOS 15.0.1;
Apple clang 15.0.0, MacOSX14.4 SDK. See environment log for sessionInfo and
all tested dependency versions. No dependency installation was necessary.
Core SHA-256 values: review/logs/baseline-core-sha256.txt.
Generated local object/shared-library files came only from testing and are excluded
from delivery. Source was clean before baseline tests (only review logs added).

## Existing API and model
Read README, DESCRIPTION, NAMESPACE, R/interface.R, R/canonical.R,
src/canonical.cpp, all man files and tests, inst/examples/demo.R and inst/CITATION.
cplx_shrink(y, X, prior, model, relation, hyper, iterations, warm, thin, start)
returns cplxshrink_fit. coef returns posterior mean coefficients;
coefficient_summary returns marginal real-coordinate means, SDs and intervals.
Fits hold beta, s, w, hyper, draws, sampleidx, final, iterations/warm/thin,
model settings and finite-draw diagnostics. beta has retained draws by d*p columns,
real coordinates first then imaginary coordinates for d=2; d=1 is real.
s and w are positive; hyper has zero Student columns, one lasso L column,
or p GDP lambda columns. sampleidx starts at warm+thin.

No automatic intercept, centering or scaling; any included intercept is shrunk.
Complex real-coordinate error covariance is s*[[1+Re(rho),Im(rho)],
[Im(rho),1-Re(rho)]] per observation. s is the component scale (each marginal
variance equals s only for Re(rho)=0); total complex variance is 2*s.
rho is fixed, |rho|<1. Real variance is s. Priors are Student-t, lasso, GDP;
gamma uses shape/rate, inverse gamma density x^(-a-1)*exp(-b/x).
Defaults and all six samplers remain frozen. Paper calibration differs from
package defaults. Tests validate specified reference calculations; they do not
reproduce all manuscript simulations, four-chain analyses or fMRI preprocessing.

## Development decisions
Add posterior as an Imports dependency for maintained rank-normalized R-hat,
bulk/tail ESS and MCSE; parallel (base recommended R) supplies reproducible
L'Ecuyer-CMRG streams. New functionality will wrap the existing sampler.
No sampling distributions, hyperparameter defaults or reference fixtures change.
