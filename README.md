# cplxrv

**Complex Random Variables for R**

An R package for working with complex random variables, providing efficient Gibbs sampling algorithms and MCMC methods for Bayesian inference with complex-valued data.

## Installation

You can install the development version of cplxrv from GitHub:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install cplxrv
devtools::install_github("Qishi7/cplxrv")
```

## Features

- **Complex Random Variable Functions**: Tools for handling complex-valued random variables
- **Covariance Operations**: Convert between complex and real covariance matrices
- **MCMC Sampling**: Efficient Gibbs samplers for various Bayesian models:
  - Bayesian Lasso (`gibbs_bl_cpp_full`)
  - Complex Bayesian Lasso (`gibbs_cbl_cpp_full`) 
  - Generalized Double Pareto (`gibbs_gdp_cpp_full`, `gibbs_cgdp_cpp_full`)
  - Student's t models (`gibbs_t_cpp_full`, `gibbs_ct_cpp_full`)
- **High Performance**: C++ implementations using RcppArmadillo for speed

## Quick Start

```r
library(cplxrv)

# Load sample data
data("complex_sample")
data("complex_data")
data("mcmc_params")

# Basic usage example
# Prepare your data
# data should be a list with elements: y_re, y_im, X, beta_re, beta_im

# Set up parameters
p <- 10  # number of predictors
n_mcmc <- 1000  # number of MCMC iterations

# Example 1: Bayesian Lasso for real part
par_info <- set_par(p = p, prior = "bl", is_cplx = FALSE)

# Run Gibbs sampler for real part
result_re_bl <- gibbs_bl_cpp_full(
  n_mcmc = n_mcmc, 
  start_lst = par_info$start_lst, 
  name_par = par_info$name_par, 
  y = y_re,  # your response vector (real part)
  X = X      # your design matrix
)

# Add column names
colnames(result_re_bl$draws) <- par_info$name_par

# Example 2: Complex Bayesian Lasso
par_info <- set_par(p = p, prior = "bl", is_cplx = TRUE)

result_cbl <- gibbs_cbl_cpp_full(
  n_mcmc = n_mcmc, 
  start_lst = par_info$start_lst, 
  name_par = par_info$name_par, 
  y = c(y_re, y_im),  # combined real and imaginary parts
  X = as.matrix(Matrix::bdiag(X, X))  # block diagonal design matrix
)

colnames(result_cbl$draws) <- par_info$name_par

# Example 3: Student's t model
par_info <- set_par(p = p, prior = "t", is_cplx = FALSE)

result_t <- gibbs_t_cpp_full(
  n_mcmc = n_mcmc, 
  start_lst = par_info$start_lst, 
  name_par = par_info$name_par, 
  y = y_re, 
  X = X
)

colnames(result_t$draws$draws) <- par_info$name_par

# Example 4: Generalized Double Pareto
par_info <- set_par(p = p, prior = "gdp", is_cplx = FALSE)

result_gdp <- gibbs_gdp_cpp_full(
  n_mcmc = n_mcmc, 
  start_lst = par_info$start_lst, 
  name_par = par_info$name_par, 
  y = y_re, 
  X = X
)

## Complete Workflow Example

```r
library(cplxrv)
library(Matrix)

# Set up simulation parameters
n <- 1000         # sample size
p <- 10           # number of predictors
n_mcmc <- 1000    # MCMC iterations
sig2 <- 1         # noise variance
rho <- 0          # correlation parameter (0 for circular case)

# Define sparse coefficient vectors
beta_re <- c(3, 1.5, 2, rep(0, p - 3))  # sparse real coefficients
beta_im <- c(3, 1.5, 2, rep(0, p - 3))  # sparse imaginary coefficients

# Generate simulation data
set.seed(123)

# Generate design matrix X (n x p)
X <- matrix(rnorm(n * p), n, p)

# Generate complex noise
eps_re <- rnorm(n, 0, sqrt(sig2))
eps_im <- rnorm(n, 0, sqrt(sig2))

# Generate response variables
y_re <- X %*% beta_re + eps_re
y_im <- X %*% beta_im + eps_im

# Create data structure
data_sim <- list(
  X = X,
  y_re = as.vector(y_re),
  y_im = as.vector(y_im),
  beta_re = beta_re,
  beta_im = beta_im,
  n = n,
  p = p
)

# Compare different priors and models:

# 1. Bayesian Lasso (separate real and imaginary)
par_info_bl <- set_par(p = p, prior = "bl", is_cplx = FALSE)

system.time({
  # Real part
  result_re_bl <- gibbs_bl_cpp_full(
    n_mcmc = n_mcmc, 
    start_lst = par_info_bl$start_lst, 
    name_par = par_info_bl$name_par, 
    y = data_sim$y_re, X = data_sim$X
  )
  colnames(result_re_bl$draws) <- par_info_bl$name_par
  
  # Imaginary part  
  result_im_bl <- gibbs_bl_cpp_full(
    n_mcmc = n_mcmc, 
    start_lst = par_info_bl$start_lst, 
    name_par = par_info_bl$name_par, 
    y = data_sim$y_im, X = data_sim$X
  )
  colnames(result_im_bl$draws) <- par_info_bl$name_par
})

# 2. Complex Bayesian Lasso (joint modeling)
par_info_cbl <- set_par(p = p, prior = "bl", is_cplx = TRUE)

system.time({
  result_cbl <- gibbs_cbl_cpp_full(
    n_mcmc = n_mcmc, 
    start_lst = par_info_cbl$start_lst, 
    name_par = par_info_cbl$name_par, 
    y = c(data_sim$y_re, data_sim$y_im), 
    X = as.matrix(Matrix::bdiag(data_sim$X, data_sim$X))
  )
  colnames(result_cbl$draws) <- par_info_cbl$name_par
})

# 3. Student's t models
par_info_t <- set_par(p = p, prior = "t", is_cplx = FALSE)

system.time({
  result_re_t <- gibbs_t_cpp_full(
    n_mcmc = n_mcmc, 
    start_lst = par_info_t$start_lst, 
    name_par = par_info_t$name_par, 
    y = data_sim$y_re, X = data_sim$X
  )
  colnames(result_re_t$draws$draws) <- par_info_t$name_par
})

# 4. Generalized Double Pareto
par_info_gdp <- set_par(p = p, prior = "gdp", is_cplx = FALSE)

system.time({
  result_re_gdp <- gibbs_gdp_cpp_full(
    n_mcmc = n_mcmc, 
    start_lst = par_info_gdp$start_lst, 
    name_par = par_info_gdp$name_par, 
    y = data_sim$y_re, X = data_sim$X
  )
  
# Analyze results
summary(result_re_bl$draws[, 1:p])  # posterior summaries for coefficients
plot(result_re_bl$draws[, 1])       # trace plot for first coefficient

# Compare with true values
cat("True real coefficients:", data_sim$beta_re, "\n")
cat("Posterior means (BL):", colMeans(result_re_bl$draws[, 1:p]), "\n")
```

## Data Generation Example

For generating complex random variable data with specific correlation structures:

```r
# Simple data generation (as used above)
gen_cp_data <- function(n, p, beta_re, beta_im, sig2 = 1, rho = 0, seed = 123) {
  set.seed(seed)
  
  # Generate design matrix
  X <- matrix(rnorm(n * p), n, p)
  
  # Generate complex error terms
  # For rho = 0 (circular case): independent real and imaginary parts
  eps_re <- rnorm(n, 0, sqrt(sig2))
  eps_im <- rnorm(n, 0, sqrt(sig2))
  
  # Generate response variables
  y_re <- X %*% beta_re + eps_re
  y_im <- X %*% beta_im + eps_im
  
  # Return data structure
  list(
    X = X,
    y_re = as.vector(y_re),
    y_im = as.vector(y_im),
    beta_re = beta_re,
    beta_im = beta_im,
    n = n,
    p = p,
    sig2 = sig2,
    rho = rho
  )
}

# Example usage:
data_sim <- gen_cp_data(
  n = 1000, 
  p = 10, 
  beta_re = c(3, 1.5, 2, rep(0, 7)), 
  beta_im = c(3, 1.5, 2, rep(0, 7)),
  sig2 = 1,
  rho = 0,  # circular case
  seed = 123
)
```

## Sample Data

The package includes several datasets for testing and demonstration:

- `complex_sample`: A vector of 100 complex random variables
- `complex_data`: A data frame with complex numbers and their properties
- `sigma_real`: A sample 2×2 covariance matrix
- `mcmc_params`: Default MCMC parameters

## Main Functions

### Complex Variable Operations
- `cplx2real_cov()`: Convert complex covariance to real form
- `real2cplx_var()`: Convert real variance to complex form
- `cplx2real_variance()`: Convert complex variance to real form

### MCMC Samplers
- `gibbs_bl_cpp_full()`: Bayesian Lasso Gibbs sampler
- `gibbs_cbl_cpp_full()`: Complex Bayesian Lasso Gibbs sampler
- `gibbs_gdp_cpp_full()`: Generalized Double Pareto Gibbs sampler
- `gibbs_cgdp_cpp_full()`: Complex Generalized Double Pareto Gibbs sampler
- `gibbs_t_cpp_full()`: Student's t Gibbs sampler
- `gibbs_ct_cpp_full()`: Complex Student's t Gibbs sampler

### Utility Functions
- `dcnorm()`: Complex normal density
- `mse()`, `mse_sample()`: Mean squared error calculations
- `roc()`: ROC curve analysis
- `mcc()`: Matthews correlation coefficient

## Dependencies

- R (>= 3.5.0)
- Rcpp (>= 1.0.14)
- RcppArmadillo
- stats
- emulator
- LaplacesDemon
- Matrix
- mvtnorm

## Development

This package uses:
- **Rcpp** and **RcppArmadillo** for high-performance C++ implementations
- **Roxygen2** for documentation
- **usethis** for package development workflow

## License

[Add your license information here]

## Citation

If you use this package in your research, please cite:

```
[Add citation information]
```

## Issues and Contributions

Please report issues or request features at: https://github.com/Qishi7/cplxrv/issues

## Author

[Qishi Zhan qishizhan7@gmail.com]