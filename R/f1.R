f1 <- function(tnfpfntp) {
  TN <- tnfpfntp[1]
  FP <- tnfpfntp[2]
  FN <- tnfpfntp[3]
  TP <- tnfpfntp[4]
  if (2 * TP + FP + FN == 0) {
    return (0)
  } else {
    return(2 * TP / (2 * TP + FP + FN))
  }
}
