#' Predict Conditional Means or New Responses
#'
#' Uses unstacked design rows in the original coefficient order. No intercept,
#' centering, scaling or whitening is added or inverted.
#' @param object A fit from [cplx_shrink()].
#' @param newdata Finite numeric or complex matrix, n_new by p, in original
#'   design column order. Real fits require a real matrix. When both training
#'   and new matrices have column names they must match exactly and in order.
#' @param type "mean" returns conditional linear predictor draws; "predictive"
#'   adds new independent observation errors conditional on each posterior draw.
#' @param summary One logical value; FALSE returns draws, TRUE summaries.
#' @param level Central marginal interval probability.
#' @param seed NULL uses and advances the current RNG for predictive draws.
#'   An integer seed makes predictions reproducible and restores the caller's
#'   RNG state/kind. Mean predictions do not use the RNG.
#' @param ... Unused, for S3 compatibility.
#' @return With summary=FALSE, a retained-draws by n_new matrix, complex for
#'   complex fits. With summary=TRUE, one row per observation for real fits,
#'   or three rows per observation (real, imaginary, magnitude) for complex
#'   fits, with mean, SD, central lower/upper quantiles, width, level and count.
#'   Complex magnitude is computed draw-wise, not from the mean prediction.
#' @details Complex errors have real-coordinate covariance
#'   s * matrix(c(1+Re(rho), Im(rho), Im(rho), 1-Re(rho)), 2, 2).
#'   Total complex error variance is 2*s; real fits have variance s.
#'   Intervals are marginal, not joint complex credible sets. Summaries require
#'   two draws; raw predictions support a single draw. Errors must be finite
#'   and scales strictly positive. There is no temporal/spatial covariance model.
#' @export
#' @examples
#' set.seed(11)
#' X <- matrix(rnorm(12), 12, 1)
#' f <- cplx_shrink(X[,1] + 1i*rnorm(12), X, iterations=30, warm=10)
#' predict(f, X[1:2,,drop=FALSE], type="mean", summary=TRUE)
#' predict(f, X[1:2,,drop=FALSE], type="predictive", seed=12, summary=TRUE)
predict.cplxshrink_fit <- function(object, newdata, type = c("mean", "predictive"),
    summary = FALSE, level = 0.95, seed = NULL, ...) {
  z <- .validate_fit(object)
  type <- match.arg(type); .validate_level(level)
  if (!is.logical(summary) || length(summary) != 1L || is.na(summary)) stop("summary must be one logical value")
  if (!is.null(seed)) .validate_seed(seed)
  if (!is.matrix(newdata) || !(is.numeric(newdata) || is.complex(newdata)) ||
      !nrow(newdata) || ncol(newdata) != z$p || any(!is.finite(newdata)))
    stop("newdata must be a finite unstacked matrix with ", z$p, " columns and at least one row")
  if (z$d == 1L && is.complex(newdata)) stop("real fits require real newdata")
  if (!is.null(object$design$columns) && !is.null(colnames(newdata)) &&
      !identical(object$design$columns, colnames(newdata))) stop("newdata column names/order must match the original design")
  mu <- .beta_complex(object) %*% t(newdata)
  if (type == "predictive") {
    s <- as.vector(object$s); rho <- object$model$rho
    if (!is.numeric(s) || is.complex(s) || length(s) != z$n || any(!is.finite(s) | s <= 0)) stop("invalid positive finite noise scale draws")
    if (!(is.numeric(rho) || is.complex(rho)) || length(rho) != 1L ||
        !is.finite(rho) || Mod(rho) >= 1 || (z$d == 1L && rho != 0)) stop("invalid fixed relation rho")
    generate <- function() {
      z1 <- matrix(stats::rnorm(length(mu)), nrow=z$n)
      if (z$d == 1L) return(mu + sqrt(s)*z1)
      z2 <- matrix(stats::rnorm(length(mu)), nrow=z$n)
      a <- Re(rho); u <- Im(rho)
      mu + sqrt(s*(1+a))*z1 + 1i*sqrt(s)*
        (u/sqrt(1+a)*z1 + sqrt((1-a*a-u*u)/(1+a))*z2)
    }
    if (is.null(seed)) mu <- generate() else mu <- .with_rng(function() {
      set.seed(seed, kind="Mersenne-Twister", normal.kind="Inversion", sample.kind="Rejection")
      generate()
    })
  }
  if (any(!is.finite(mu))) stop("nonfinite predictions; check design and coefficient scales")
  colnames(mu) <- if (is.null(rownames(newdata))) paste0("observation",seq_len(nrow(newdata))) else rownames(newdata)
  if (!summary) return(mu)
  parts <- if (z$d == 2L) list(real=Re(mu), imaginary=Im(mu), magnitude=Mod(mu)) else list(real=mu)
  out <- lapply(names(parts), function(nm) cbind(
    data.frame(observation=seq_len(nrow(newdata)), component=nm),
    .summarize_matrix(parts[[nm]],level)))
  do.call(rbind,out)
}
