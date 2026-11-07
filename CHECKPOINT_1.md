# Checkpoint 1: magnitude and multiple chains

Phases 1 and 2 complete. Version 0.1.0.9000, local development branch;
no fork/push available (see BASELINE_REPORT.md).

Command: `LC_ALL=en_US.UTF-8 Rscript -e 'testthat::test_local(".", reporter="summary")'`.
Actual counted result: 280 expectations, 13 test blocks, 0 failures/errors/warnings/skips.
Evidence: review/logs/checkpoint1-tests.log and checkpoint1-test-counts.csv.
All 151 original expectations still pass. Core files remain unchanged.

Added magnitude_summary, cplx_shrink_chains, as_draws_array, mcmc_diagnostics.
The wrapper calls the original cplx_shrink, using separate reproducible
L'Ecuyer-CMRG streams and caller RNG restoration. Different starts are supported.
Diagnostics preserve iteration/chain/variable axes and include beta, s, w,
prior rates, and magnitude. Undefined statistics return NA with explanations.
Single-chain diagnostics do not claim convergence.

Independent tests demonstrate opposing draws +1/-1 give mean magnitude 1
while modulus of their mean is zero; hand-calculated multi-coefficient quantiles;
all six real/complex and prior combinations; stream replay and ordering;
independently assembled posterior::summarise_draws oracle; separated-location
chains with R-hat > 1.5; invalid input and RNG preservation including absence.
An initial diagnostics test wrongly supplied a three-dimensional draws_array
slice to the low-level posterior matrix API. This gave a different chain layout;
the test was corrected to use posterior::summarise_draws on the independently
assembled full array. The implementation already passed explicit matrices.
Initial failure log is retained for transparency.

Changed: DESCRIPTION, NAMESPACE, R/interface.R (design metadata only),
.Rbuildignore, BASELINE_REPORT.md; added R/posterior.R, tests/testthat/test-posterior.R,
.gitignore, man/{magnitude_summary,cplx_shrink_chains,as_draws_array,mcmc_diagnostics}.Rd,
this report and logs. Roxygen 7.3.2 regenerated documentation.
posterior >=1.6.0 is a new maintained diagnostics dependency; parallel,
graphics and grDevices are standard R libraries for the planned interfaces.

Baseline UTF-8 --as-cran check completed: 0 errors, 0 warnings, 4 notes:
incoming new submission; unable to verify current time; missing pandoc for README;
old system HTML validator does not recognize R-generated <main> markup.
Baseline standard check remains 0/0/0. No baseline statistical failure occurred.

Limits: serial chains; retained draws stored in memory and copied for conversion;
short runs are software demonstrations, not convergence evidence. Real magnitudes
are available as absolute values in diagnostics; public magnitude_summary requires
a complex fit. Optional paper-calibration/M-ReIm helpers are outside required scope.
Proceeding to prediction and usability without waiting for approval, as requested.
