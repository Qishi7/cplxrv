# cplxrv

**Scale Mixtures of Complex Gaussian and Bayesian Shrinkage**

An R package implementing scale mixtures of complex Gaussian distributions for Bayesian shrinkage and variable selection in complex-valued regression models.

## Installation

You can install the development version of cplxrv from GitHub:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install cplxrv
devtools::install_github("Qishi7/cplxrv")
```

## Overview

This package extends scale mixtures of normal distributions to the complex domain, providing:

- **Complex Scale Mixtures**: Complex Student's t, complex Laplace, and complex Generalized Double Pareto (GDP) distributions
- **Bayesian Shrinkage**: MCMC algorithms for complex-valued Bayesian linear regression with shrinkage priors
- **Variable Selection**: Methods for detecting relevant predictors in complex-valued data
- **High Performance**: C++ implementations using RcppArmadillo for computational efficiency

## Key Features

### Complex Distributions
- `dcnorm()`: Complex normal density
- Complex Student's t distribution
- Complex multivariate Laplace distribution  
- Complex Generalized Double Pareto (GDP) distribution

### MCMC Samplers
- `gibbs_bl_cpp_full()`: Bayesian Lasso Gibbs sampler (real-valued)
- `gibbs_cbl_cpp_full()`: **Complex** Bayesian Lasso Gibbs sampler
- `gibbs_t_cpp_full()`: Student's t Gibbs sampler (real-valued)
- `gibbs_ct_cpp_full()`: **Complex** Student's t Gibbs sampler
- `gibbs_gdp_cpp_full()`: GDP Gibbs sampler (real-valued)
- `gibbs_cgdp_cpp_full()`: **Complex** GDP Gibbs sampler

### Utility Functions
- `cplx2real_cov()`: Convert complex covariance to real representation
- `real2cplx_var()`: Convert real variance to complex form
- `set_par()`: Set up parameter configurations for MCMC
- `mse()`, `mse_sample()`: Mean squared error calculations
- `roc()`: ROC curve analysis
- `mcc()`: Matthews correlation coefficient

## Quick Start

```r
library(cplxrv)
library(Matrix)

# Set up simulation parameters
n <- 1000         # sample size
p <- 10           # number of predictors
n_mcmc <- 3000    # MCMC iterations
sig2 <- 1         # noise variance
rho <- 0          # correlation parameter (0 for circular case)

# Define sparse coefficient vectors (first 3 non-zero)
beta_re <- c(3, 1.5, 2, rep(0, p - 3))  # real coefficients
beta_im <- c(3, 1.5, 2, rep(0, p - 3))  # imaginary coefficients

# Generate complex-valued data
set.seed(123)

# Generate design matrix
X <- matrix(rnorm(n * p), n, p)

# Generate complex noise with specified correlation structure
eps_re <- rnorm(n, 0, sqrt(sig2))
eps_im <- rnorm(n, 0, sqrt(sig2))

# If rho != 0, adjust for correlation between real and imaginary parts
# For rho = 0 (circular case), real and imaginary parts are independent

# Generate response variables
y_re <- X %*% beta_re + eps_re
y_im <- X %*% beta_im + eps_im

# Convert to vectors
y_re <- as.vector(y_re)
y_im <- as.vector(y_im)

# Example 1: Complex Bayesian Lasso
par_info <- set_par(p = p, prior = "bl", is_cplx = TRUE)

result_cbl <- gibbs_cbl_cpp_full(
  n_mcmc = n_mcmc,
  start_lst = par_info$start_lst,
  name_par = par_info$name_par,
  y = c(y_re, y_im),  # stack real and imaginary parts
  X = as.matrix(Matrix::bdiag(X, X))  # block diagonal design matrix
)

# Add parameter names
colnames(result_cbl$draws) <- par_info$name_par

# Example 2: Complex Student's t
par_info_t <- set_par(p = p, prior = "t", is_cplx = TRUE)

result_ct <- gibbs_ct_cpp_full(
  n_mcmc = n_mcmc,
  start_lst = par_info_t$start_lst,
  name_par = par_info_t$name_par,
  y = c(y_re, y_im),
  X = as.matrix(Matrix::bdiag(X, X))
)

# Example 3: Compare with separate real-valued models
par_info_real <- set_par(p = p, prior = "bl", is_cplx = FALSE)

# Real part model
result_re <- gibbs_bl_cpp_full(
  n_mcmc = n_mcmc,
  start_lst = par_info_real$start_lst,
  name_par = par_info_real$name_par,
  y = y_re,
  X = X
)

# Imaginary part model  
result_im <- gibbs_bl_cpp_full(
  n_mcmc = n_mcmc,
  start_lst = par_info_real$start_lst,
  name_par = par_info_real$name_par,
  y = y_im,
  X = X
)
```

## Model Comparison Framework

The package supports three modeling approaches:

1. **Complex Model (M-Cplx)**: Joint modeling of real and imaginary parts with complex shrinkage priors
2. **Separate Real Models (M-Re-Im)**: Independent models for real and imaginary parts
3. **Magnitude Model (M-Mag)**: Real-valued model using magnitude data only

## Applications

This package is particularly useful for:

- **Complex-valued fMRI**: Brain activation detection using both magnitude and phase information
- **Signal Processing**: Analysis of complex-valued signals preserving phase relationships
- **Communications**: Processing of complex baseband signals
- **General Complex Data**: Any application involving complex-valued observations

## Advantages of Complex Modeling

Research shows that complex-valued models provide:

- **Better Variable Selection**: Higher F1 scores and more accurate selection of relevant predictors
- **Improved Parameter Estimation**: Lower MSE, especially in high-dimensional settings
- **Enhanced Prediction**: Better out-of-sample prediction performance
- **Phase Information Preservation**: Utilizes both magnitude and phase components

## Dependencies

- R (>= 3.5.0)
- Rcpp (>= 1.0.14)
- RcppArmadillo
- Matrix
- stats
- emulator
- LaplacesDemon
- mvtnorm

## Citation

If you use this package in your research, please cite:

```
Yu, C.-H. and Zhan, Q. (2025). Scale Mixtures of Complex Gaussian and Bayesian Shrinkage.
```

## Authors

- Qishi Zhan, Department of Mathematical and Statistical Sciences, Marquette University
- Cheng-Han Yu, Department of Mathematical and Statistical Sciences, Marquette University

## License

[Specify your license]

## Issues and Contributions

Please report issues or request features at: https://github.com/Qishi7/cplxrv/issues