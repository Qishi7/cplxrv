mse_sample <- function(sample_mat, test) {
  # apply(pred_mat, 1, function(z) {
  #     mean((z - test) ^ 2)
  # })
  sapply(1:ncol(sample_mat), function(k) {
    mean((sample_mat[, k] - test[k]) ^ 2)
  })
}