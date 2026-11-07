# Independent review of cplxshrink 0.1.0.9000

Development package; not an official v0.2 release. The supplied ZIP unpacks
one folder named cplxshrink. It includes all source, tests, help, CI and reports.
The owner's existing repository Qishi7/cplxrv is the publication destination.
This installs as cplxshrink; the old cplxrv interface is replaced. Prior repository
history is preserved. Historical baseline/checkpoint reports retain the initial
local-run access limitations; see PUBLISHING_REPORT.md for the upload status.

## Minimal reproduction (macOS/Linux shell, R >= 4.1 and C++ toolchain)

From the directory containing this extracted cplxshrink folder:

```sh
# en_US.UTF-8 is installed on the tested Mac; choose an installed UTF-8 locale elsewhere.
export LC_ALL=en_US.UTF-8
Rscript -e 'install.packages(c("Rcpp", "RcppArmadillo", "posterior", "testthat", "jsonlite", "roxygen2"), repos="https://cloud.r-project.org")'
R CMD build cplxshrink
R CMD check cplxshrink_0.1.0.9000.tar.gz
Rscript -e 'testthat::test_local("cplxshrink", reporter="summary")'
R CMD INSTALL cplxshrink_0.1.0.9000.tar.gz
Rscript cplxshrink/examples/validate_new_features.R
Rscript -e 'source(system.file("examples", "validate_new_features.R", package="cplxshrink", mustWork=TRUE))'
```

Dependencies were already installed during this run, so the installation command
is supplied for fresh machines rather than claimed as executed. Reproduction commands
for build/check/test/install/demo were executed against this source. On Windows,
install Rtools matching R and use R/Rscript in a shell with those tools available;
set a supported UTF-8 locale or omit the export when the default is valid.
No local Windows/Linux execution is claimed; the workflow provides those jobs.

Expect a successful source build, Status: OK on a configured standard check,
444 testthat expectations (18 test blocks) without failures, and a small summary
ending "Validation demo completed. Short chains do not establish convergence."
Exact diagnostic/random values may differ with R, posterior or toolchain versions.
The demo is intentionally too short to establish convergence.

To regenerate help from the source root:

```sh
cd cplxshrink
Rscript -e 'roxygen2::roxygenise()'
```

roxygen2 7.3.2 was used; do not regenerate Rcpp exports. R CMD check also runs
tests/stage25-smoke.R. --as-cran is an additional strict check; see the report
for the precise environmental notes. Reports, CI and source examples are in
this ZIP; R CMD build excludes them from the installed package as appropriate.
The demo is also installed at inst/examples/validate_new_features.R.

## Integrity and patch

From inside the extracted folder, `shasum -a 256 -c SHA256SUMS` verifies every
payload file except the manifest itself. The ZIP has no .git directory, compiled
objects or research data. The source tarball supplied separately is directly
installable. Frozen core hashes are in review/logs/baseline-core-sha256.txt and
review/logs/final-core-sha256.txt. development.patch applies at upstream
`d985957c4254593ae2cc5c449d6dc3d108dd1b5b` with `git apply development.patch`.

## GitHub installation and provenance

```r
install.packages("remotes")
remotes::install_github("Qishi7/cplxrv")
library(cplxshrink)
```

This existing repository replaces the need to create a new fork, as explicitly
requested by the owner. The old main commit is
7ef7f9e3decb6294e6774d08d9444768f6ad9f58. No force-push or history deletion is
needed. PUBLISHING_REPORT.md records the replacement and its checks.
The upstream manuscript repository remains unchanged.
