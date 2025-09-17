mse <- function(pt_est, test) {
  mean((pt_est - test) ^ 2)
}
