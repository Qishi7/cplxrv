roc <- function(tnfpfntp, ...) {
  TN <- tnfpfntp[1]
  FP <- tnfpfntp[2]
  FN <- tnfpfntp[3]
  TP <- tnfpfntp[4]
  
  TPR <- TP / (FN + TP)
  FPR <- FP / (TN + FP)
  
  plot(FPR, TPR, xlim = c(0, 1), ylim = c(0, 1), type = "b",
       main = "ROC", ...)
}