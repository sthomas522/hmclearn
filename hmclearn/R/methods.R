
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


#' @export
predict.hmclearn <- function(object, X, fam = "linear", burnin=100, nsamp=NULL, ...) {

  thetaDF <- object$thetaDF[-c(1:burnin), ]
  k <- ncol(thetaDF)

  if (is.null(nsamp)) {
    nsamp <- nrow(thetaDF)
  }
  index <- 1:nsamp
  sampvals <- sample(index, size=nsamp, replace = FALSE)

  if (!is.matrix(X)) {
    X <- matrix(X, nrow=1)
  }

  preds <- NULL
  if (fam == "linear") {
    preds <- sapply(sampvals, function(xx, dat=thetaDF, bparam=1:(k-1), sigparam=k) {
      betaval <- as.numeric(dat[xx, bparam])
      sigval <- sqrt(exp(dat[xx, sigparam]))
      xval <- X[sample(1:nrow(X), size=1), ]
      z <- rnorm(1, mean=xval%*%betaval, sd=sigval)
      z
    })
  } else if (fam == "binomial") {
    preds <- sapply(sampvals, function(xx, dat=thetaDF, bparam=1:k) {
      betaval <- as.numeric(dat[xx, bparam])
      xval <- X[sample(1:nrow(X), size=1), ]
      bx <- exp(sum(xval%*%betaval))
      z <- rbinom(1, 1, bx / (1 + bx))
      z
    })
  } else if (fam == "poisson") {
    preds <- sapply(sampvals, function(xx, dat=thetaDF, bparam=1:k) {
      betaval <- as.numeric(dat[xx, bparam])
      xval <- X[sample(1:nrow(X), size=1), ]
      z <- exp(sum(xval%*%betaval))
      z
    })
  }
  return(preds)
}
