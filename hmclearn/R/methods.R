
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


#' Summarizing HMC Model Fits
#'
#' summary method for class \code{hmclearn}
#'
#' @param x an object of class \code{hmclearn}, usually a result of a call to \code{mh} or \code{hmc}
#' @param burnin optional numeric parameter for the number of initial MCMC samples to omit from the summary
#' @param probs quantiles to summarize the posterior distribution
#' @param ... additional arguments to pass to \code{quantile}
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
plot.hmclearn <- function(x, burnin=NULL, ...) {
  # diagplots(x, ...)
  thetaCombined <- combMatrix(x$thetaCombined, burnin=burnin)
  bayesplot::mcmc_hist(thetaCombined, )
}

#' @export
coef.hmclearn <- function(x, burnin=NULL, prob=0.5, ...) {
  thetaCombined <- combMatrix(x$thetaCombined, burnin=burnin)
  thetaCombined <- do.call(rbind, thetaCombined)
  apply(thetaCombined, 2, quantile, probs=prob, ...)
}


#' @export
predict.hmclearn <- function(object, y, X, fam = "linear", burnin=NULL, draws=NULL, ...) {

  thetaCombined <- combMatrix(object$thetaCombined, burnin=burnin)
  thetaCombined <- do.call(rbind, thetaCombined)
  k <- ncol(thetaCombined)

  if (is.null(draws)) {
    draws <- min(500, nrow(thetaCombined))
  }
  index <- 1:nrow(thetaCombined)
  sampvals <- sample(index, size=draws, replace = FALSE)

  # for each value of theta sample all N values of y

  if (!is.matrix(X)) {
    X <- matrix(X, nrow=1)
  }

  preds <- NULL
  if (fam == "linear") {
    preds <- t(sapply(sampvals, function(xx, dat=thetaCombined, bparam=1:(k-1), sigparam=k) {
      betaval <- as.numeric(dat[xx, bparam])
      sigval <- sqrt(exp(dat[xx, sigparam]))
      bx <- X %*% betaval
      v <- sapply(bx, function(zz, sd=sigval, n=1) {
        rnorm(n, mean=zz, sd=sd)
      })
      v
    }))
  } else if (fam == "binomial") {
    preds <- t(sapply(sampvals, function(xx, dat=thetaCombined, bparam=1:k) {
      betaval <- as.numeric(dat[xx, bparam])
      bx <- exp(X%*%betaval)
      v <- sapply(bx, function(zz, n=1) {
        rbinom(n, 1, zz / (1 + zz))
      })
      v
    }))
  } else if (fam == "poisson") {
    preds <- t(sapply(sampvals, function(xx, dat=thetaCombined, bparam=1:k) {
      betaval <- as.numeric(dat[xx, bparam])
      bx <- exp(X%*%betaval)
      v <- sapply(bx, function(zz, n=1) {
        rpois(n, bx)
      })
      v
    }))
  }

  retval <- list(y = y,
                 yrep = preds,
                 X = X)
  class(retval) <- c("hmclearnpred", "list")

  return(retval)
}

#' #' @export
pp_check <- function(object, ...) {
  UseMethod("pp_check")
}

#' @export
pp_check.hmclearnpred <- function(object,
                                  type=c("overlaid",
                                         "multiple",
                                         "bars",
                                         "stat"), ...) {

  y <- object[["y"]]
  yrep <- object[["yrep"]]
  switch(match.arg(type),
         multiple = bayesplot::ppc_hist(y, yrep[1:min(8, nrow(yrep)),, drop = FALSE]),
         overlaid = bayesplot::ppc_dens_overlay(y, yrep),
         bars = bayesplot::ppc_bars(y, yrep, ...),
         stat = bayesplot::ppc_stat(y, yrep, ...)
  )

}
