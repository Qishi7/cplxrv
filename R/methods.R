#' Summarize a Regression Fit
#' @param object A cplxshrink_fit object.
#' @param level Central marginal interval probability.
#' @param ... Unused, for S3 compatibility.
#' @return A summary.cplxshrink_fit list with model, design, schedule,
#'   coefficient and magnitude summaries. Single-chain convergence is not assessed.
#' @export
#' @examples
#' X <- matrix(1:6, 6, 1)
#' f <- cplx_shrink(1:6, X, model="real", iterations=20, warm=5)
#' summary(f)
summary.cplxshrink_fit <- function(object, level = 0.95, ...) {
  z <- .validate_fit(object, 2L); .validate_level(level)
  b <- .summarize_matrix(object$beta,level)
  b <- cbind(parameter=colnames(object$beta),b)
  structure(list(model=object$model, design=object$design,
    schedule=c(iterations=object$iterations,warm=object$warm,thin=object$thin,retained=z$n),
    coefficients=b, magnitude=if(z$d==2L) magnitude_summary(object,level) else NULL,
    convergence="Not assessed from a single chain"),class="summary.cplxshrink_fit")
}

#' Print a Fit Summary
#' @param x A summary returned by summary.cplxshrink_fit.
#' @param ... Passed to print.data.frame.
#' @return x invisibly.
#' @export
print.summary.cplxshrink_fit <- function(x, ...) {
  cat("cplxshrink",if(x$model$d==2L) "complex" else "real",x$model$family,"fit\n")
  if (!is.null(x$design)) cat("Design:", x$design$n,"observations,",x$design$p,"coefficients\n")
  print(x$schedule)
  print(x$coefficients,...)
  if (!is.null(x$magnitude)) { cat("Draw-wise magnitudes:\n"); print(x$magnitude,...) }
  cat(x$convergence,"\n")
  invisible(x)
}

#' Summarize Multiple Chains
#' @param object A cplxshrink_chains object.
#' @param level Central marginal interval probability.
#' @param ... Unused, for S3 compatibility.
#' @return A summary.cplxshrink_chains list with model/design/schedule,
#'   pooled marginal summaries of all variables and separate-chain diagnostics.
#'   Pooling is used only for descriptive summaries, never for diagnostics.
#' @export
#' @examples
#' X <- matrix(1:6, 6, 1)
#' z <- cplx_shrink_chains(1:6, X, model="real", iterations=20, warm=5)
#' summary(z)
summary.cplxshrink_chains <- function(object, level = 0.95, ...) {
  .validate_level(level); a <- as_draws_array(object)
  mat <- matrix(a, nrow=dim(a)[1]*dim(a)[2], ncol=dim(a)[3],
    dimnames=list(NULL,dimnames(a)[[3]]))
  structure(list(model=object$fits[[1]]$model, design=object$fits[[1]]$design,
    schedule=c(chains=dim(a)[2],iterations=object$fits[[1]]$iterations,
      warm=object$fits[[1]]$warm,thin=object$fits[[1]]$thin,retained_per_chain=dim(a)[1]),
    parameters=cbind(parameter=colnames(mat),.summarize_matrix(mat,level)),
    diagnostics=mcmc_diagnostics(object)),class="summary.cplxshrink_chains")
}

#' Print Multiple-chain Summary
#' @param x A summary returned by summary.cplxshrink_chains.
#' @param ... Passed to print.data.frame.
#' @return x invisibly.
#' @export
print.summary.cplxshrink_chains <- function(x, ...) {
  cat("cplxshrink multiple-chain", x$model$family, "summary\n")
  print(x$schedule); print(x$parameters,...); print(x$diagnostics,...)
  cat("Assess diagnostics in context; no convergence classification is made.\n")
  invisible(x)
}

.plot_draws <- function(a, sampleidx, parameters, type, dots) {
  vars <- dimnames(a)[[3L]]
  if (is.null(parameters)) parameters <- utils::head(vars,4L)
  if (!is.character(parameters) || !length(parameters) || anyNA(parameters) ||
      any(!parameters %in% vars) || anyDuplicated(parameters) || length(parameters)>16L)
    stop("parameters must select 1 to 16 distinct variable names")
  if (dim(a)[1L] < 2L) stop("plot requires at least 2 retained draws per chain")
  old <- graphics::par(mfrow=c(ceiling(length(parameters)/2),min(2,length(parameters))))
  on.exit(graphics::par(old))
  colors <- grDevices::hcl.colors(dim(a)[2L],"Dark 3")
  for (nm in parameters) {
    j <- match(nm,vars)
    m <- matrix(a[,,j],nrow=dim(a)[1L],ncol=dim(a)[2L])
    if (type=="trace") {
      do.call(graphics::matplot,c(list(x=sampleidx,y=m,type="l",lty=1,col=colors,
        xlab="Original iteration",ylab=nm,main=nm),dots))
    } else {
      dens <- lapply(seq_len(ncol(m)),function(k) stats::density(m[,k]))
      xr <- range(vapply(dens,function(d) range(d$x),numeric(2)))
      yr <- c(0,max(vapply(dens,function(d) max(d$y),numeric(1))))
      do.call(graphics::plot,c(list(x=dens[[1]],xlim=xr,ylim=yr,col=colors[1],
        xlab=nm,main=nm),dots))
      if(ncol(m)>1L) for(k in 2:ncol(m)) graphics::lines(dens[[k]],col=colors[k])
    }
    if(ncol(m)>1L) graphics::legend("topright",legend=paste("chain",seq_len(ncol(m))),
      col=colors,lty=1,bty="n",cex=.65)
  }
  invisible(a)
}

#' Plot Posterior Traces or Densities
#' @param x A cplxshrink_fit object.
#' @param parameters Character vector of 1 to 16 variable names. Default:
#'   first four variables. Includes beta coordinates, s, w, rates, magnitudes.
#' @param type "trace" or "density".
#' @param ... Additional plot arguments excluding labels, layout and colors
#'   supplied by the method.
#' @return Draws array invisibly. Graphics parameters are restored.
#' @details A trace is not a convergence guarantee. Base graphics works with
#'   PDF devices in headless environments. At most 16 parameters per call.
#' @export
#' @examples
#' X <- matrix(1:6, 6, 1)
#' f <- cplx_shrink(1:6, X, model="real", iterations=20, warm=5)
#' plot(f, parameters="beta1")
plot.cplxshrink_fit <- function(x, parameters = NULL, type = c("trace","density"), ...) {
  type <- match.arg(type)
  a <- as_draws_array(structure(list(fits=list(x)),class="cplxshrink_chains"))
  .plot_draws(a,x$sampleidx,parameters,type,list(...))
}

#' Plot Multiple Chains Separately
#' @param x A cplxshrink_chains object.
#' @inheritParams plot.cplxshrink_fit
#' @return Draws array invisibly. Each chain is shown in its own color;
#'   density curves are computed separately for each chain.
#' @export
#' @examples
#' X <- matrix(1:6, 6, 1)
#' z <- cplx_shrink_chains(1:6, X, model="real", iterations=20, warm=5)
#' plot(z, parameters="beta1", type="density")
plot.cplxshrink_chains <- function(x, parameters = NULL, type = c("trace","density"), ...) {
  type <- match.arg(type); a <- as_draws_array(x)
  .plot_draws(a,x$fits[[1]]$sampleidx,parameters,type,list(...))
}
