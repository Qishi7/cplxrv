test_that("coef pairs coordinates in design-column order", {
  fit <- structure(list(
    beta = rbind(c(1, -2, 10, 3), c(3, 0, 30, 5)),
    model = list(d = 2L)
  ), class = c("cplxshrink_fit", "list"))
  expect_identical(coef(fit), c(beta1 = 2 + 20i, beta2 = -1 + 4i))
  fit$beta <- fit$beta[, 1:2]
  colnames(fit$beta) <- c("beta1", "beta2")
  fit$model$d <- 1L
  expect_identical(coef(fit), c(beta1 = 2, beta2 = -1))
  fit$beta <- matrix(c(1, 10), nrow = 1)
  fit$model$d <- 2L
  expect_identical(coef(fit), c(beta1 = 1 + 10i))
})

test_that("coef is available for each fitted model and leaves draws and RNG intact", {
  X <- matrix(seq(-1, 1, length.out = 24), 8, 3)
  y <- seq(-0.4, 0.5, length.out = 8)
  for (model in c("real", "complex")) {
    yy <- if (model == "complex") y + 1i * rev(y) else y
    XX <- if (model == "complex") X + 1i * X / 3 else X
    for (prior in c("student", "lasso", "gdp")) {
      set.seed(420)
      fit <- cplx_shrink(yy, XX, prior = prior, model = model,
                         iterations = 40, warm = 3, thin = 2)
      saved <- fit
      rng <- .Random.seed
      actual <- coef(fit)
      expect_identical(is.complex(actual), model == "complex")
      expect_length(actual, 3)
      expect_identical(fit, saved)
      expect_identical(.Random.seed, rng)
    }
  }
})
