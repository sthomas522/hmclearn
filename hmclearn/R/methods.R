
#' @export
summary.hmclearn <- function(x, burnin=100, probs=c(0.05, 0.25, 0.5, 0.75, 0.95), ...) {
  cat("Summary of HMC simulation\n\n")

  thetaDF <- x$thetaDF[-c(1:burnin), ]
  t(apply(thetaDF, 2, quantile, probs=probs, ...))
}

#' @export
plot.hmclearn <- function(x, ...) {
  diagplots(x, ...)
}

#' @export
coef.hmclearn <- function(x, burnin=100, prob=0.5, ...) {
  thetaDF <- x$thetaDF[-c(1:burnin), ]
  apply(thetaDF, 2, quantile, probs=prob, ...)
}


