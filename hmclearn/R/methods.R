

# combArray <- function(x, burnin, marg=3) {
#   xnew <- x[-c(1:burnin), 1:dim(x)[2], 1:dim(x)[3]]
#   if (length(dim(xnew)) > 2) {
#     xnew <- aperm(xnew, c(3, 1, 2))
#     xnew <- matrix(xnew, prod(dim(xnew)[1:2]), dim(xnew)[3])
#   } else {
#     xnew <- apply(xnew, 2, as.numeric)
#   }
#   print(class(xnew))
#   print(apply(xnew, 2, class))
#   xnewDF <- as.data.frame(xnew)
#   return(xnewDF)
# }

combMatrix <- function(x, burnin) {
  if (is.null(burnin)) {
    return(x)
  }
  xnew <- lapply(x, function(zz) {
    zz <- zz[-c(1:burnin), ]
    zz
  })
  xnew
}

#' @export
summary.hmclearn <- function(x, burnin=NULL, probs=c(0.05, 0.25, 0.5, 0.75, 0.95), ...) {
  cat("Summary of HMC simulation\n\n")

  # remove burnin
  thetaCombined <- combMatrix(x$thetaCombined, burnin=burnin)
  thetaCombined <- do.call(rbind, thetaCombined)

  # quantiles
  res <- t(apply(thetaCombined, 2, quantile, probs=probs, ...))
  if (!is.null(x$varnames)) {
    row.names(res) <- x$varnames
  }

  # rhat calc
  rhat <- psrf(x, burnin=burnin)

  res <- cbind(res, rhat)
  colnames(res)[ncol(res)] <- "rhat"
  res
}

#' @export
plot.hmclearn <- function(x, burnin=1, ...) {
  # diagplots(x, ...)
  thetaCombined <- combMatrix(x$thetaCombined, burnin=burnin)
  bayesplot::mcmc_hist(thetaCombined, )
}

#' @export
coef.hmclearn <- function(x, burnin=1, prob=0.5, ...) {
  thetaCombined <- combMatrix(x$thetaCombined, burnin=burnin)
  thetaCombined <- do.call(rbind, thetaCombined)
  apply(thetaCombined, 2, quantile, probs=prob, ...)
}


#' @export
predict.hmclearn <- function(object, X, fam = "linear", burnin=1, nsamp=NULL, ...) {

  thetaCombined <- combMatrix(object$thetaCombined, burnin=burnin)
  thetaCombined <- do.call(rbind, thetaCombined)
  k <- ncol(thetaCombined)

  if (is.null(nsamp)) {
    nsamp <- nrow(thetaCombined)
  }
  index <- 1:nsamp
  sampvals <- sample(index, size=nsamp, replace = FALSE)

  if (!is.matrix(X)) {
    X <- matrix(X, nrow=1)
  }

  preds <- NULL
  if (fam == "linear") {
    preds <- sapply(sampvals, function(xx, dat=thetaCombined, bparam=1:(k-1), sigparam=k) {
      betaval <- as.numeric(dat[xx, bparam])
      sigval <- sqrt(exp(dat[xx, sigparam]))
      xval <- X[sample(1:nrow(X), size=1), ]
      z <- rnorm(1, mean=xval%*%betaval, sd=sigval)
      z
    })
  } else if (fam == "binomial") {
    preds <- sapply(sampvals, function(xx, dat=thetaCombined, bparam=1:k) {
      betaval <- as.numeric(dat[xx, bparam])
      xval <- X[sample(1:nrow(X), size=1), ]
      bx <- exp(sum(xval%*%betaval))
      z <- rbinom(1, 1, bx / (1 + bx))
      z
    })
  } else if (fam == "poisson") {
    preds <- sapply(sampvals, function(xx, dat=thetaCombined, bparam=1:k) {
      betaval <- as.numeric(dat[xx, bparam])
      xval <- X[sample(1:nrow(X), size=1), ]
      z <- exp(sum(xval%*%betaval))
      z
    })
  }
  return(preds)
}


#' @export
psrf <- function(object, ...) {
  UseMethod("psrf")
}

#' @export
psrf.hmclearn <- function(object, burnin=NULL, ...) {

  data <- combMatrix(object$thetaCombined, burnin=burnin)

  if (object$chains == 1) {
    # stop("Multiple chains needed to calculate Potential Scale Reduction Factor")
    rhat <- rep(NA, ncol(data[[1]]))
    return(rhat)
  }

  data <- combMatrix(object$thetaCombined, burnin=burnin)

  # number of samples after burnin
  Nraw <- object$N
  if (!is.null(burnin)) {
    N <- Nraw - burnin
  } else {
    N <- Nraw
  }

  # number of chains
  M <- object$chains

  # means for each chain
  thetaBarm <- lapply(data, colMeans)
  thetaBarAll <- colMeans(do.call(rbind, thetaBarm))

  # between-chain var
  diff2 <- lapply(thetaBarm, function(x1, x2=thetaBarAll) {
    (x1 - x2)^2
  })
  B <- N / (M-1) * base::colSums(do.call(rbind, diff2))

  # within-chain var
  sampvarByChain <- lapply(data, function(xx) {
    apply(xx, 2, var)
  })
  W <- colMeans(do.call(rbind, sampvarByChain))

  # variance estimator
  varest <- (N-1)/N * W + 1/N*B

  rhat <- sqrt(varest/W)
  return(rhat)

}

