#' Paired Complex or Scalar Real Bayesian Shrinkage
#'
#' Fit a regression without adding an intercept. To include an intercept, add
#' a column of ones to X; this coefficient is also shrunk. Complex coefficients share
#' one local scale per real/imaginary pair. Relation is fixed and coefficient
#' priors are proper. A real fit uses a single noise scale; two calls are
#' required for the independent real and imaginary strategy.
#'
#' @param y Numeric or complex response vector.
#' @param X Numeric or complex design matrix with one row per response.
#' @param prior One of `student`, `lasso`, or `gdp`.
#' @param model `complex` for paired shrinkage or `real` for scalar shrinkage.
#' @param relation Fixed complex likelihood relation with modulus less than
#'   one. Must be zero for a real model.
#' @param hyper Named list of positive hyperparameters. All priors accept
#'   `a_s`, `b_s` (defaults 0.5). Complex Student accepts `nu` (default 1);
#'   real Student accepts `a_tau`, `b_tau` (defaults 0.5). Lasso accepts `r`,
#'   `delta` (defaults 0.5); GDP accepts `alpha`, `eta` (defaults 1).
#' @param iterations Total number of iterations.
#' @param warm Number of warmup iterations.
#' @param thin Keep every `thin`th iteration after warmup.
#' @param start Optional list of initial `beta`, `w`, `s`, and prior-specific
#'   `L` (lasso) or `lambda` (GDP).
#' @return A list with real-block coefficient draws `beta` (real columns then
#'   imaginary columns), local scales `w`, noise scale `s`, rate draws `hyper`,
#'   combined `draws`, stored `sampleidx`, final state, model metadata, and
#'   basic finite-draw diagnostics. These basic checks do not establish MCMC
#'   convergence. Assess convergence separately using multiple chains.
#' @details
#' No centering or scaling is performed. A complex design is converted
#' internally; supply the original response and design, without stacking them.
#' Stored iterations are `warm + thin`, `warm + 2 * thin`, and so on,
#' up to `iterations`. Warmup draws are discarded.
#' In real coordinates the likelihood covariance is `s * C`, where complex
#' `C` has blocks `(1+Re(relation))*I`, `Im(relation)*I`, and
#' `(1-Re(relation))*I`. Thus complex covariance is `2*s*I` and relation is
#' `2*s*relation*I`. Coefficient pair covariance is `s*w[j]*I`.
#' The noise scale `s` is the manuscript's real-component variance sigma-squared;
#' `a_s` and `b_s` are its inverse-gamma hyperparameters a_sigma and b_sigma.
#' Lasso `delta` is the manuscript's delta_L, the rate for the squared rate L.
#' Gamma parameters use shape and rate; inverse gamma has density proportional
#' to `x^(-a-1)*exp(-b/x)`. Student local scales are inverse gamma. Lasso
#' local scales have gamma shape `(d+1)/2`, rate `L/2`, with
#' `L ~ Gamma(r,delta)`. GDP uses rate `lambda[j]^2/2` and
#' `lambda[j] ~ Gamma(alpha,eta)`, with `d=2` complex and `d=1` real.
#' The GDP update draws lambda with local scales marginalized and restores
#' local scales immediately, before any coefficient/scale update or storage.
#' @export
#' @examples
#' set.seed(42)
#' X <- matrix(rnorm(60), 20, 3)
#' truth <- c(1+1i, 0, 0)
#' y <- drop(X %*% truth) + rnorm(20) + 1i*rnorm(20)
#' fit <- cplx_shrink(y, X, prior="gdp", iterations=150, warm=50)
#' coefficient_summary(fit)
cplx_shrink <- function(y, X, prior = c("student", "lasso", "gdp"),
                        model = c("complex", "real"), relation = 0,
                        hyper = list(), iterations = 2000L, warm = 500L,
                        thin = 1L, start = NULL) {
  prior <- match.arg(prior)
  model <- match.arg(model)
  allowed <- c("a_s", "b_s", switch(prior,
    student = if (model == "complex") "nu" else c("a_tau", "b_tau"),
    lasso = c("r", "delta"),
    gdp = c("alpha", "eta")
  ))
  if (!is.list(hyper) || (length(hyper) && (is.null(names(hyper)) ||
    anyNA(names(hyper)) || any(!nzchar(names(hyper))) || anyDuplicated(names(hyper)) ||
    any(!names(hyper) %in% allowed)))) {
    stop("hyper must contain only named parameters for this model and prior")
  }
  start_allowed <- c("beta", "w", "s", switch(prior,
    lasso = "L",
    gdp = "lambda",
    student = character()
  ))
  if (!is.null(start) && (!is.list(start) || (length(start) &&
    (is.null(names(start)) || anyNA(names(start)) || any(!nzchar(names(start))) ||
      anyDuplicated(names(start)) || any(!names(start) %in% start_allowed))))) {
    stop("start contains an unnamed, duplicate, or unused field")
  }
  out <- do.call(smcg_sample, c(list(
    y = y, X = X, family = prior,
    complex_model = model == "complex", rho = relation, iterations = iterations,
    warm = warm, thin = thin, start = start
  ), hyper))
  out$design <- list(n = nrow(as.matrix(X)), p = ncol(as.matrix(X)),
                     columns = colnames(as.matrix(X)))
  out$diagnostics <- list(
    draws_finite = all(is.finite(out$draws)),
    retained_draws = nrow(out$draws), convergence_assessed = FALSE
  )
  class(out) <- c("cplxshrink_fit", "list")
  out
}

#' Summarize Posterior Coefficients
#' @param fit A fit returned by [cplx_shrink()].
#' @param level Marginal central interval probability, strictly between 0 and 1.
#' @return Data frame of real-coordinate posterior means, SDs, and marginal
#'   interval endpoints. These marginal intervals are not the joint
#'   predictor-selection regions. For complex posterior mean coefficients,
#'   use `coef(fit)`. At least two finite draws are required.
#' @export
coefficient_summary <- function(fit, level = 0.95) {
  .validate_fit(fit, 2L)
  out <- .summarize_matrix(fit$beta, level)
  data.frame(parameter = if (is.null(colnames(fit$beta)))
    paste0("beta", seq_len(ncol(fit$beta))) else colnames(fit$beta),
    out[, c("mean", "sd", "lower", "upper"), drop = FALSE], row.names = NULL)
}

#' Posterior Mean Coefficients
#'
#' Return coefficients in the same order as the columns of the supplied design.
#' No intercept is added. Complex fits combine the posterior means of the
#' real and imaginary coordinates; real fits return real posterior means.
#'
#' @param object A fit returned by [cplx_shrink()].
#' @param ... Unused; included for compatibility with [stats::coef()].
#' @return A named vector of posterior mean coefficients, complex for a complex
#'   fit and numeric for a real fit.
#' @export
coef.cplxshrink_fit <- function(object, ...) {
  means <- colMeans(object$beta)
  if (object$model$d == 2L) {
    p <- length(means) %/% 2L
    means <- means[seq_len(p)] + 1i * means[p + seq_len(p)]
    names(means) <- paste0("beta", seq_len(p))
  }
  means
}
