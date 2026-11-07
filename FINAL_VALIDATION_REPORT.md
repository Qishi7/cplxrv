# Final validation: cplxshrink 0.1.0.9000

## Provenance and status

IMPLEMENTED AND TESTED locally: all required magnitude, multiple-chain diagnostics,
prediction, summary/print and trace/density interfaces; input checks; help, examples,
README; regression and independent-oracle tests; complete source distribution.

Tested source commit: `852bff6227a8b6c45f91625eba5dc8f0236d6c15`, branch
`feature/posterior-interface-v02`. Frozen upstream main:
`d985957c4254593ae2cc5c449d6dc3d108dd1b5b` (2026-11-05 UTC), package 0.1.0.
Later handoff commits contain reports and logs only; package code remains at the
exact tested commit. No source code, fixtures or data from upstream were modified
on GitHub. No fork URL or branch URL exists from this run: fork creation/push
was unavailable (connector has no fork operation; browser signed out; no gh).
This is an explicitly disclosed deviation, using the document's local fallback.

IMPLEMENTED BUT NOT VERIFIED remotely: CI workflow for Linux/macOS/Windows R
release and Linux R 4.1. Only macOS R 4.4.1 was executed locally. Backward R 4.1
compatibility is designed and configured for CI, not independently demonstrated here.

NOT IMPLEMENTED: optional M-ReIm/paper-calibration helpers, new priors, learned rho,
spatial/temporal covariance, full manuscript replication, upstream PR/release.
No complete scientific reproducibility or convergence claim is made.

## Environment

R 4.4.1 (2024-06-14), aarch64-apple-darwin20, macOS 15.0.1. Apple clang 15.0.0
(clang-1500.3.9.4), MacOSX14.4 SDK; R itself was compiled with clang 14 and
GNU Fortran 12.2.0. Rcpp 1.1.1, RcppArmadillo 14.0.0-1, posterior 1.6.1,
testthat 3.2.3, jsonlite 1.9.1, roxygen2 7.3.2. Session information is in
review/logs/{environment,final-environment}.log. No dependencies needed installation.
All successful final commands use LC_ALL=en_US.UTF-8 to correct the inherited
unsupported C.UTF-8 locale. Full initial failure evidence is retained.

## Exact commands and actual results

Commands below are from the source root except where noted; full stdout/stderr
is in the named review/logs file. R CMD check itself also runs both test scripts.

| Check | Command | Actual result | Evidence |
|---|---|---|---|
| Unchanged baseline tests | `Rscript -e 'testthat::test_local(".", reporter="summary")'` | 151 expectations, 8 blocks, no failures/errors/warnings/skips | baseline-tests.log; baseline-tests-counted.log; baseline-test-counts.csv |
| Baseline build (parent directory) | `R CMD build baseline` | source tarball built | baseline-build.log |
| Baseline check, corrected locale (parent) | `LC_ALL=en_US.UTF-8 R CMD check --output=baseline-standard cplxshrink_0.1.0.tar.gz` | 0 errors/warnings/notes | baseline-standard-check.log |
| Phase 1/2 tests | `LC_ALL=en_US.UTF-8 Rscript -e 'testthat::test_local(".", reporter="summary")'` | 280 expectations, 13 blocks, 0 failures/errors/warnings/skips | checkpoint1-tests.log; checkpoint1-test-counts.csv |
| Final original + new tests | `LC_ALL=en_US.UTF-8 Rscript -e 'testthat::test_local(".", reporter="summary")'` | 444 expectations, 18 blocks, 0 failures/errors/warnings/skips (151 original + 293 added) | final-tests.log; final-test-counts.csv |
| Final build | `LC_ALL=en_US.UTF-8 R CMD build .` | cplxshrink_0.1.0.9000.tar.gz built | final-build.log |
| Final standard check | `LC_ALL=en_US.UTF-8 R CMD check --output=../final-standard cplxshrink_0.1.0.9000.tar.gz` | Status OK: 0 errors, 0 warnings, 0 notes; PDF manual, examples, both test scripts pass | final-check.log; final-standard/ |
| Final strict check | `LC_ALL=en_US.UTF-8 R CMD check --as-cran --output=../final-cran cplxshrink_0.1.0.9000.tar.gz` | 0 errors, 0 warnings, 4 notes | final-cran-check.log; final-cran/ |
| Install source tarball | `LC_ALL=en_US.UTF-8 R CMD INSTALL cplxshrink_0.1.0.9000.tar.gz` | succeeded | final-install.log |
| Source demo | `LC_ALL=en_US.UTF-8 Rscript examples/validate_new_features.R` | stopifnot checks pass, summary printed | demo.log |
| Installed demo | `LC_ALL=en_US.UTF-8 Rscript -e 'source(system.file("examples","validate_new_features.R",package="cplxshrink",mustWork=TRUE))'` | passed | installed-demo.log |
| Documentation sync | `LC_ALL=en_US.UTF-8 Rscript -e 'roxygen2::roxygenise()'` then `git diff --exit-code HEAD -- NAMESPACE man` | generated help identical | documentation-sync.log |
| Frozen core/fixtures | `git diff --exit-code d985957c4254593ae2cc5c449d6dc3d108dd1b5b -- src/canonical.cpp R/canonical.R src/RcppExports.cpp R/RcppExports.R tests/testthat/fixtures` | exit 0; SHA-256 lists also identical | core-diff.log; baseline-core-sha256.txt; final-core-sha256.txt |

Final strict check's four notes are precisely:
1. Incoming feasibility: maintainer/new submission; development version contains
   large components (0.1.0.9000). Development numbering is intentional.
2. Future file timestamps: "unable to verify current time".
3. Top-level files: "Files 'README.md' or 'NEWS.md' cannot be checked without
   'pandoc' being installed."
4. HTML manual validation: system validator reports R-generated HTML5 <main> as
   unrecognized, plus link/script/table/id messages. These already affect the
   pristine package; the complete detailed diagnostics are retained in the log.

These notes are not called zero. The standard check has zero notes. Initial
baseline --as-cran used invalid inherited locale and returned 1 error, 1 warning,
2 notes; the corrected pristine --as-cran run returned 0 errors, 0 warnings,
4 notes, matching these categories (without the development-version text).
An initial development check's unqualified head() note and top-level examples
note were fixed via utils::head and installed example placement. Initial test
oracle layout failure and short-chain ESS warnings are retained in separate logs.
The diagnostics wrapper now reports posterior warnings in each variable's note.
No certified-engine change was used to fix any issue.

## APIs and independent reference calculations

magnitude_summary(fit, level): original coefficient index/name, draw-wise mean
magnitude, SD, central quantiles, interval width, level and retained count.
Mod uses paired real/imaginary draws. +1/-1 draws produce mean magnitude 1,
not 0; hand reference magnitudes (5,0,13) and (0,2,4) give means (6,2),
SDs (sqrt(43),2), central 50% intervals (2.5,9) and (1,3). Tests verify
non-mutation and RNG preservation, finite inputs, custom levels, p=1, zeros,
and all three complex prior fits. Public real-fit magnitude summaries explicitly error.

cplx_shrink_chains mirrors the original fitting options with chains, seed and
starts. It calls the unchanged sampler serially with parallel::nextRNGStream
L'Ecuyer-CMRG streams. Saves individual fits, IDs, master seed, stream starts,
RNG kind, model/design/schedule and original retained iteration indices.
Tests replay a stream exactly; check differing streams, starts, iteration order,
all priors/model types with p=1 and p=3; reject inconsistent metadata and scales.
Caller RNG kind/state (or absence) is restored even on error.

as_draws_array returns retained iterations x chains x variables, including beta,
s, w, prior-specific rates and derived magnitudes. mcmc_diagnostics uses posterior
rank-normalized split/folded R-hat, bulk/tail ESS, and MCSE of means on explicit
iteration-by-chain matrices. An independently assembled full array passed to
posterior::summarise_draws matches all four statistics for nonconstant variables.
Separated-location chains yield R-hat >1.5; constant and undersized cases report
NA with explanation. Warnings, including capped ESS, are recorded as notes.
Draw pooling is used only for descriptive multi-chain summaries.

predict supports unstacked n_new x p original design, mean versus predictive,
raw draws versus marginal summaries, level and optional restoring seed. Independent
complex multiplication matches the exact mean prediction. Tests use 100000 fixed
parameter draws with s=1.7, rho=0.3-0.4i and target covariance
[[2.21,-0.68],[-0.68,1.19]]. Each empirical covariance must be within six
Gaussian standard errors sqrt((Sigma_ij^2+Sigma_ii*Sigma_jj)/(N-1)); zero means
and cross-observation independence are also checked. Real variance is checked
against 1.7 using six*s*sqrt(2/(N-1)). Varying-s standardized noise is checked.
Complex summaries compute real, imaginary and draw-wise magnitude separately.
Tests cover p=1/n_new=1, one draw, invalid/nonfinite input, named-column mismatch,
seed behavior and caller-data preservation. No original-scale transform is claimed.

summary and plot methods work with all model/prior types. Summaries show
model/prior/design/schedule and marginal coefficient/variable summaries; single
chain convergence remains unassessed. Trace/density methods limit parameters,
show chains separately and restore graphics settings. Headless PDF tests pass.
coefficient_summary retains its existing output columns while improving validation.

## Compatibility, assumptions and limits

All four protected files, the six sampling distributions/order, fixed-rho likelihood,
prior/hyperparameter defaults and original fixtures are byte-identical. The original
151 expectations, including independent quadrature and invariant-distribution checks,
still pass. The only fit-interface addition is design metadata (n, p, column names);
no raw data are stored. Existing coef behavior and real-block indexing are preserved.
Rejecting fewer than two draws in coefficient_summary intentionally replaces an
unsupported SD with an explicit error.

s is the component scale, not total complex variance; with Re(rho) nonzero the
marginal variances are s*(1+Re(rho)), s*(1-Re(rho)), and covariance s*Im(rho).
Gamma shape/rate and inverse-gamma scale conventions are unchanged. Intervals
use R type-7 quantiles and are marginal, not joint credible sets or selection rules.
No automatic intercept, scaling, whitening or transformation inversion is added.

Chain storage and conversion use memory proportional to retained draws x chains x
variables; array conversion copies draws. Prediction stores draws x new observations
and uses matrix multiplication. Chains are serial. No C++ optimization or performance
benchmark is claimed. Very large designs/draws require user memory planning.
Short examples and passing checks establish interface behavior, not scientific
convergence, coverage, causal effects or full manuscript reproducibility.

## Handoff and remaining decisions

The review ZIP contains installable source, all docs/tests/workflow/examples,
reports, logs, CSV validation summary, baseline patch and SHA-256 manifest.
A directly installable source tarball is also supplied. README_FOR_SUPERVISOR.md
contains independently runnable commands and exact fork/push instructions.
review/logs/changed-files.txt records exact changed/added paths; NEWS.md gives
release notes. No private data, credentials, .git directories or compiled binaries
are included. Archive extraction/integrity and reproduction evidence is in
review/logs/archive-reproduction.log and archive-validation.log.

Supervisor decisions: create/push the fork when authenticated; inspect CI results
on other OS/R versions; review APIs/statistical conventions and decide whether to
accept a future release. No upstream merge or official release was performed.
