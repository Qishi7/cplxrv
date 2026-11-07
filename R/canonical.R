# Canonical API accepts an unstacked n x p design; complex responses/designs
# are converted to the [Re; Im] representation. No intercept is added.
smcg_sample <- function(y, X, family = c("student", "lasso", "gdp"),
                        complex_model = is.complex(y) || is.complex(X),
                        rho = 0, rho_b = 0, learn_rho = FALSE, griddy = FALSE,
                        iterations = 2000L, warm = 500L, thin = 1L,
                        nu = 1, a_tau = 0.5, b_tau = 0.5,
                        r = 0.5, delta = 0.5, alpha = 1, eta = 1,
                        a_s = 0.5, b_s = 0.5, start = NULL) {
  # 1. Validate model settings, data, and the storage schedule.
  family <- match.arg(family)
  scalar <- function(x) length(x) == 1L && is.finite(x)
  if (!is.logical(complex_model) || length(complex_model) != 1L || is.na(complex_model)) {
    stop("complex_model must be one logical value")
  }
  if (!identical(learn_rho, FALSE)) stop("unknown-rho learning is disabled")
  if (!identical(griddy, FALSE)) stop("GDP griddy learning is disabled")
  if (!scalar(rho_b) || rho_b != 0) stop("coefficient priors require rho_b=0")
  if (!scalar(rho) || Mod(rho) >= 1 || (!complex_model && rho != 0)) stop("invalid fixed rho")
  ints <- c(iterations, warm, thin)
  if (length(ints) != 3L || any(!is.finite(ints)) || any(ints != floor(ints)) ||
    any(ints > .Machine$integer.max) || iterations < 1 || warm < 0 || warm >= iterations || thin < 1 ||
    floor((iterations - warm) / thin) < 1) {
    stop("invalid iterations/warm/thin")
  }
  pos <- list(nu, a_tau, b_tau, r, delta, alpha, eta, a_s, b_s)
  if (any(!vapply(pos, function(x) scalar(x) && is.numeric(x) && !is.complex(x) && x > 0, logical(1)))) {
    stop("hyperparameters must be positive finite real scalars")
  }
  X <- as.matrix(X)
  y <- as.vector(y)
  if (!(is.numeric(X) || is.complex(X)) || !(is.numeric(y) || is.complex(y)) || !ncol(X) || !nrow(X) ||
    nrow(X) != length(y) || any(!is.finite(X)) || any(!is.finite(y))) {
    stop("invalid data")
  }
  if (!complex_model && (is.complex(y) || is.complex(X))) stop("real sampler requires real data")

  # 2. Construct the structured real representation internally.
  p <- ncol(X)
  d <- if (complex_model) 2L else 1L
  H <- if (complex_model) rbind(cbind(Re(X), -Im(X)), cbind(Im(X), Re(X))) else X
  yy <- if (complex_model) c(Re(y), Im(y)) else y

  # 3. Initialize coefficients, local scales, noise, and rates.
  if (is.null(start)) start <- list()
  if (!is.list(start)) stop("start must be a list")
  unknown <- setdiff(names(start), c("beta", "w", "s", "L", "lambda"))
  if (length(unknown)) stop("unknown start field: ", paste(unknown, collapse = ", "))
  b <- if (is.null(start$beta)) rep(0, d * p) else start$beta
  if (complex_model && is.complex(b) && length(b) == p) b <- c(Re(b), Im(b))
  if (!is.numeric(b) || is.complex(b) || length(b) != d * p || any(!is.finite(b))) stop("invalid beta start")
  w <- if (is.null(start$w)) rep(1, p) else start$w
  lam <- if (is.null(start$lambda)) rep(1, p) else start$lambda
  s <- if (is.null(start$s)) 1 else start$s
  L <- if (is.null(start$L)) 1 else start$L
  valid_vec <- function(x, k) is.numeric(x) && !is.complex(x) && length(x) == k && all(is.finite(x) & x > 0)
  if (!valid_vec(w, p) || !valid_vec(lam, p) || !valid_vec(s, 1) || !valid_vec(L, 1)) stop("invalid positive start")
  at <- if (complex_model) nu else a_tau
  bt <- if (complex_model) nu else b_tau

  # 4. Call the certified MCMC engine without changing its update order.
  out <- canonical_cpp(
    as.numeric(yy), H, d, family, Re(rho), Im(rho), as.integer(iterations),
    as.integer(warm), as.integer(thin), as.numeric(b), w, s, L, lam,
    at, bt, r, delta, alpha, eta, a_s, b_s
  )

  # 5. Label and organize retained posterior draws.
  colnames(out$beta) <- if (complex_model) {
    c(paste0("beta_re", seq_len(p)), paste0("beta_im", seq_len(p)))
  } else {
    paste0("beta", seq_len(p))
  }
  colnames(out$w) <- paste0("w", seq_len(p))
  colnames(out$hyper) <- if (family == "gdp") {
    paste0("lambda", seq_len(p))
  } else if (family == "lasso") {
    "L"
  } else {
    character(0)
  }
  out$draws <- cbind(out$beta, s = as.vector(out$s), out$w, out$hyper)
  out$iterations <- as.integer(iterations)
  out$warm <- as.integer(warm)
  out$thin <- as.integer(thin)
  out$model <- list(
    family = family, d = d, rho = rho, rho_b = 0, nu = nu, a_tau = at, b_tau = bt,
    r = r, delta = delta, alpha = alpha, eta = eta, a_s = a_s, b_s = b_s, convention = "c=2"
  )

  # 6. Return draws, the final state, and model conventions.
  class(out) <- c("smcgcert_draws", "list")
  out
}
sample_complex_student <- function(y, X, ...) smcg_sample(y, X, family = "student", complex_model = TRUE, ...)
sample_real_student <- function(y, X, ...) smcg_sample(y, X, family = "student", complex_model = FALSE, ...)
sample_complex_lasso <- function(y, X, ...) smcg_sample(y, X, family = "lasso", complex_model = TRUE, ...)
sample_real_lasso <- function(y, X, ...) smcg_sample(y, X, family = "lasso", complex_model = FALSE, ...)
sample_complex_gdp <- function(y, X, ...) smcg_sample(y, X, family = "gdp", complex_model = TRUE, ...)
sample_real_gdp <- function(y, X, ...) smcg_sample(y, X, family = "gdp", complex_model = FALSE, ...)
