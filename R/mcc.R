mcc <- function(tnfpfntp) {
  TN <- tnfpfntp[1]
  FP <- tnfpfntp[2]
  FN <- tnfpfntp[3]
  TP <- tnfpfntp[4]
  
  A <- TP * TN - FP * FN
  S1 <- TP + FP
  S2 <- TP + FN
  S3 <- TN + FP
  S4 <- TN + FN
  
  if (S1 == 0 | S2 == 0 | S3 == 0 | S4 == 0) {
    B = 1
  } else {
    B <- S1 * S2 * S3 * S4
  }
  A / sqrt(B)
}