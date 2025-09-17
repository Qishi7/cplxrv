short_dist <- function(cm) {
  results <- list()
  
  for (i in seq_along(cm)) {
    conf_mat <- cm[[i]]
    num_thresholds <- ncol(conf_mat)
    distances <- numeric(num_thresholds)
    
    for (j in 1:num_thresholds) {
      tn <- conf_mat[1, j]
      fp <- conf_mat[2, j]
      fn <- conf_mat[3, j]
      tp <- conf_mat[4, j]
      total <- tn + fp + fn + tp
      
      if (total == 0) {
        distances[j] <- Inf
        next
      }
      
      correct_rate <- tp / (tp + fn)
      false_rate <- fp / (fp + tn)
      distances[j] <- sqrt((1 - correct_rate) ^ 2 + false_rate ^ 2)
    }
    
    min_dist <- min(distances, na.rm = TRUE)
    min_cis <- which(distances == min_dist)
    
    results[[i]] <- list(min_dist = min_dist, min_cis = min_cis)
  }
  
  return(results)
}
