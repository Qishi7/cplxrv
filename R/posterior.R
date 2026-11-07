# Internal validation shared by posterior interfaces; never modifies draws.
.validate_level <- function(level) {
  if (!is.numeric(level) || is.complex(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1) stop("invalid interval level")
}
.validate_fit <- function(fit, min_draws = 1L) {
  if (!inherits(fit, "cplxshrink_fit")) stop("fit must be a cplxshrink_fit")
  d <- fit$model$d
  if (length(d) != 1L || is.na(d) || !d %in% c(1L, 2L)) stop("invalid model dimension")
  b <- fit$beta
  if (!is.matrix(b) || !is.numeric(b) || is.complex(b) ||
      !ncol(b) || ncol(b) %% d || nrow(b) < min_draws || any(!is.finite(b))) {
    stop("beta must contain finite real draws, nonempty coefficient columns and at least ",
         min_draws, " retained draws")
  }
  list(n = nrow(b), p = ncol(b) %/% d, d = d)
}
.beta_complex <- function(fit) {
  z <- .validate_fit(fit)
  j <- seq_len(z$p)
  if (z$d == 1L) fit$beta else
    fit$beta[, j, drop = FALSE] + 1i * fit$beta[, z$p + j, drop = FALSE]
}
.magnitudes <- function(fit) {
  out <- Mod(.beta_complex(fit))
  if (any(!is.finite(out))) stop("nonfinite coefficient magnitudes")
  colnames(out) <- paste0("magnitude", seq_len(ncol(out)))
  out
}
.summarize_matrix <- function(x, level) {
  .validate_level(level)
  if (!is.matrix(x) || nrow(x) < 2L || !ncol(x) || any(!is.finite(x)))
    stop("summary requires finite values and at least 2 retained draws")
  q <- vapply(seq_len(ncol(x)), function(j) stats::quantile(x[, j],
    c((1-level)/2, (1+level)/2), names = FALSE), numeric(2))
  data.frame(mean = colMeans(x), sd = apply(x, 2, stats::sd),
    lower = q[1, ], upper = q[2, ], width = q[2, ] - q[1, ],
    level = level, n_draws = nrow(x), row.names = NULL)
}

#' Draw-wise Posterior Coefficient Magnitudes
#'
#' Computes the modulus of each retained complex coefficient draw, then
#' summarizes those moduli. This is the paper's draw-wise magnitude, not the
#' modulus of the posterior mean coefficient. Intervals are central marginal
#' magnitude intervals, not joint credible regions or activation decisions.
#' @param fit A complex fit from [cplx_shrink()].
#' @param level Central interval probability strictly between zero and one.
#' @return Data frame with coefficient index, name, mean, SD, lower/upper
#'   quantiles (R type 7), width, level and retained draw count. Names use
#'   original design columns where available, otherwise beta1, beta2, etc.
#' @details Requires at least two finite draws and one coefficient. Real fits,
#'   empty/malformed/nonfinite draws and invalid levels produce errors.
#'   Neither stored draws nor the RNG state is changed.
#' @export
#' @examples
#' set.seed(7)
#' X <- matrix(rnorm(12), 12, 1)
#' fit <- cplx_shrink(X[,1] + 1i*rnorm(12), X, iterations=30, warm=10)
#' magnitude_summary(fit)
magnitude_summary <- function(fit, level = 0.95) {
  z <- .validate_fit(fit, 2L)
  .validate_level(level)
  if (z$d != 2L) stop("magnitude_summary requires a complex model; real fits are unsupported")
  nm <- fit$design$columns
  if (is.null(nm)) nm <- paste0("beta", seq_len(z$p))
  cbind(data.frame(index = seq_len(z$p), coefficient = nm),
        .summarize_matrix(.magnitudes(fit), level))
}

# Restore even the absence of a caller RNG state, including on error.
.with_rng <- function(code) {
  kind <- RNGkind()
  had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had) old <- get(".Random.seed", .GlobalEnv, inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (had) assign(".Random.seed", old, .GlobalEnv) else
      if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  })
  code()
}
.validate_seed <- function(seed) {
  if (!is.numeric(seed) || is.complex(seed) || length(seed) != 1L ||
      !is.finite(seed) || seed < 0 || seed > .Machine$integer.max || seed != floor(seed))
    stop("seed must be a nonnegative integer within R's integer range")
}

#' Fit Independent Reproducible MCMC Chains
#'
#' Calls the unchanged [cplx_shrink()] serially with separate L'Ecuyer-CMRG
#' streams from [parallel::nextRNGStream()]. Caller RNG kind and state are
#' restored, including on error. No parallel workers are launched.
#' @inheritParams cplx_shrink
#' @param chains Positive integer number of chains (default four).
#' @param seed Nonnegative integer master seed. Stream starting states are saved.
#' @param starts NULL, or a list of length chains containing NULL or a start
#'   list for each chain. Different initial values are encouraged.
#' @return A cplxshrink_chains list containing fits, chain_ids, seed, streams,
#'   rng_kind, common model/design/schedule metadata and retained sampleidx.
#'   Chains use identical data and model parameters. Memory grows linearly
#'   with retained draws, chains and coefficient count.
#' @export
#' @examples
#' X <- matrix(seq(-1, 1, length.out=10), 10, 1)
#' chains <- cplx_shrink_chains(X[,1], X, model="real", seed=42,
#'                            iterations=30, warm=10)
#' mcmc_diagnostics(chains)
cplx_shrink_chains <- function(y, X, prior = c("student", "lasso", "gdp"),
    model = c("complex", "real"), relation = 0, hyper = list(),
    iterations = 2000L, warm = 500L, thin = 1L,
    chains = 4L, seed = 1L, starts = NULL) {
  prior <- match.arg(prior); model <- match.arg(model)
  if (!is.numeric(chains) || is.complex(chains) || length(chains) != 1L ||
      !is.finite(chains) || chains < 1 || chains != floor(chains) ||
      chains > .Machine$integer.max) stop("chains must be a positive integer")
  .validate_seed(seed)
  if (is.null(starts)) starts <- rep(list(NULL), chains)
  if (!is.list(starts) || length(starts) != chains ||
      any(!vapply(starts, function(x) is.null(x) || is.list(x), logical(1))))
    stop("starts must be a list with one start list or NULL per chain")
  .with_rng(function() {
    set.seed(seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
    streams <- vector("list", chains)
    streams[[1L]] <- get(".Random.seed", .GlobalEnv)
    if (chains > 1L) for (k in 2:chains)
      streams[[k]] <- parallel::nextRNGStream(streams[[k-1L]])
    fits <- lapply(seq_len(chains), function(k) {
      assign(".Random.seed", streams[[k]], .GlobalEnv)
      cplx_shrink(y, X, prior=prior, model=model, relation=relation, hyper=hyper,
        iterations=iterations, warm=warm, thin=thin, start=starts[[k]])
    })
    structure(list(fits=fits, chain_ids=seq_len(chains), seed=seed,
      streams=streams, rng_kind=RNGkind(), model=fits[[1L]]$model,
      design=fits[[1L]]$design, iterations=iterations, warm=warm, thin=thin,
      sampleidx=fits[[1L]]$sampleidx), class="cplxshrink_chains")
  })
}

#' Convert Chains to a Posterior Draws Array
#' @param x An object returned by [cplx_shrink_chains()].
#' @return A posterior draws_array with axes retained iteration, chain,
#'   variable. Variables include beta coordinates, s, w, prior-specific rates,
#'   and draw-wise magnitudes (absolute coefficients for real models).
#'   The original iteration indices remain in x$sampleidx; array indices are
#'   consecutive retained draw positions. Inconsistent chains produce errors.
#' @export
#' @examples
#' X <- matrix(1:6, 6, 1)
#' z <- cplx_shrink_chains(1:6, X, model="real", iterations=12, warm=4)
#' dim(as_draws_array(z))
as_draws_array <- function(x) {
  if (!inherits(x, "cplxshrink_chains") || !is.list(x$fits) || !length(x$fits))
    stop("x must contain cplxshrink chains")
  first <- x$fits[[1L]]; .validate_fit(first)
  mats <- lapply(x$fits, function(fit) {
    z <- .validate_fit(fit)
    if (!identical(fit$model, first$model) ||
        !identical(fit$design, first$design) ||
        !identical(fit$sampleidx, first$sampleidx) ||
        !identical(c(fit$iterations,fit$warm,fit$thin), c(first$iterations,first$warm,first$thin)) ||
        !identical(dim(fit$beta), dim(first$beta))) stop("inconsistent chain metadata or dimensions")
    if (!is.matrix(fit$w) || !identical(dim(fit$w), c(z$n,z$p)) ||
        length(fit$s) != z$n || !is.matrix(fit$hyper) || nrow(fit$hyper) != z$n ||
        ncol(fit$hyper) != switch(fit$model$family, student=0L, lasso=1L, gdp=z$p) ||
        any(!is.finite(c(fit$s,fit$w,fit$hyper))) || any(c(fit$s,fit$w,fit$hyper) <= 0))
      stop("invalid scale or hyperparameter draws")
    mat <- cbind(fit$beta, s=as.vector(fit$s), fit$w, fit$hyper, .magnitudes(fit))
    if (is.null(colnames(mat)) || anyNA(colnames(mat)) ||
        any(!nzchar(colnames(mat))) || anyDuplicated(colnames(mat))) stop("invalid variable names")
    mat
  })
  if (any(!vapply(mats, function(m) identical(colnames(m),colnames(mats[[1L]])), logical(1))))
    stop("inconsistent chain variable names")
  a <- array(NA_real_, c(nrow(mats[[1L]]), length(mats), ncol(mats[[1L]])),
    dimnames=list(iteration=seq_len(nrow(mats[[1L]])), chain=seq_along(mats),
                  variable=colnames(mats[[1L]])))
  for (k in seq_along(mats)) a[,k,] <- mats[[k]]
  posterior::as_draws_array(a)
}

#' Rank-normalized Multiple-chain Diagnostics
#' @param x A cplxshrink_chains object, or a single cplxshrink_fit.
#' @return Data frame of variable, rank-normalized split/folded rhat,
#'   ess_bulk, ess_tail, mcse_mean, and note. Uses posterior directly with
#'   iteration and chain axes preserved. With fewer than two chains or four
#'   draws per chain all statistics are NA with an explanation. Constant
#'   chains or other undefined statistics are NA and explained per variable.
#'   Warnings from posterior (such as capped ESS) are retained in the note column.
#' @details No convergence classification is made. Short examples cannot
#'   establish convergence. MCSE describes the posterior mean of each listed
#'   variable, including derived magnitudes.
#' @export
#' @examples
#' X <- matrix(1:6, 6, 1)
#' z <- cplx_shrink_chains(1:6, X, model="real", iterations=30, warm=10)
#' mcmc_diagnostics(z)
mcmc_diagnostics <- function(x) {
  if (inherits(x, "cplxshrink_fit")) x <- structure(list(fits=list(x)), class="cplxshrink_chains")
  a <- as_draws_array(x)
  vars <- dimnames(a)[[3L]]
  out <- data.frame(variable=vars, rhat=NA_real_, ess_bulk=NA_real_,
    ess_tail=NA_real_, mcse_mean=NA_real_, note="")
  if (dim(a)[2L] < 2L || dim(a)[1L] < 4L) {
    out$note <- "Requires at least two chains and four retained draws per chain"
    return(out)
  }
  funs <- list(posterior::rhat, posterior::ess_bulk, posterior::ess_tail, posterior::mcse_mean)
  for (j in seq_along(vars)) {
    m <- matrix(a[,,j], nrow=dim(a)[1L], ncol=dim(a)[2L])
    warnings <- character()
    vals <- withCallingHandlers(
      vapply(funs, function(f) f(m), numeric(1)),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    out$note[j] <- paste(unique(warnings), collapse="; ")
    if (any(!is.finite(vals))) {
      out$note[j] <- paste(c(warnings, "Undefined statistic (e.g. constant within-chain draws)"), collapse="; ")
      vals[!is.finite(vals)] <- NA_real_
    }
    out[j, 2:5] <- as.list(vals)
  }
  out
}
