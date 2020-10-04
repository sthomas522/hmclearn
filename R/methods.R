
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
#' @param object an object of class \code{hmclearn}, usually a result of a call to \code{mh} or \code{hmc}
#' @param burnin optional numeric parameter for the number of initial MCMC samples to omit from the summary
#' @param probs quantiles to summarize the posterior distribution
#' @param ... additional arguments to pass to \code{quantile}
#' @returns Returns a matrix with posterior quantiles and the posterior scale reduction factor statistic for each parameter.
#' @references Gelman, A., et. al. (2013) \emph{Bayesian Data Analysis}.  Chapman and Hall/CRC.
#' @references Gelman, A. and Rubin, D. (1992) \emph{Inference from Iterative Simulation Using Multiple Sequences}.  Statistical Science 7(4) 457-472.
#' @export
#' @examples
#' # Linear regression example
#' set.seed(521)
#' X <- cbind(1, matrix(rnorm(300), ncol=3))
#' betavals <- c(0.5, -1, 2, -3)
#' y <- X%*%betavals + rnorm(100, sd=.2)
#'
#' f1 <- hmc(N = 500,
#'           theta.init = c(rep(0, 4), 1),
#'           epsilon = 0.01,
#'           L = 10,
#'           logPOSTERIOR = linear_posterior,
#'           glogPOSTERIOR = g_linear_posterior,
#'           varnames = c(paste0("beta", 0:3), "log_sigma_sq"),
#'           param=list(y=y, X=X), parallel=FALSE, chains=1)
#'
#' summary(f1)
summary.hmclearn <- function(object, burnin=NULL, probs=c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975), ...) {
  cat("Summary of MCMC simulation\n\n")

  # remove burnin
  thetaCombined <- combMatrix(object$thetaCombined, burnin=burnin)
  thetaCombined <- do.call(rbind, thetaCombined)

  # quantiles
  res <- t(apply(thetaCombined, 2, quantile, probs=probs, ...))
  if (!is.null(object$varnames)) {
    row.names(res) <- object$varnames
  }

  # rhat calc
  if (object$chains > 1) {
    rhat <- psrf(object, burnin=burnin)

    res <- cbind(res, rhat)
    colnames(res)[ncol(res)] <- "rhat"
  }

  res
}


print.hmclearn <- function(obj) {
  if (obj$algorithm == "MH") {
    alg <- "Metropolis Hastings"
  } else if (obj$algorithm == "HMC") {
    alg <- "Hamiltonian Monte Carlo"
  }

  cat(paste(alg, "\n\n"))
  cat(paste(obj$N, "simulations from", obj$chains, "chains.\n\n"))
  cat(paste("Average acceptance rate of", round(mean(obj$accept/obj$N), 3)), "\n\n")
  cat("Quantiles with no burnin period\n\n")

  # remove burnin
  thetaCombined <- combMatrix(obj$thetaCombined, burnin=NULL)
  thetaCombined <- do.call(rbind, thetaCombined)

  # quantiles
  res <- t(apply(thetaCombined, 2, quantile, probs=c(0.25, 0.5, 0.75)))
  if (!is.null(obj$varnames)) {
    row.names(res) <- obj$varnames
  }
  print(res)

}

#' Extract Model Coefficients
#'
#' Method for \code{hmclearn} objects created by \code{mh} and \code{hmc} functions.  Extracts the specified quantile of the posterior.
#'
#' @param object an object of class \code{hmclearn}, usually a result of a call to \code{mh} or \code{hmc}
#' @param burnin optional numeric parameter for the number of initial MCMC samples to omit from the summary
#' @param prob quantile to extract coefficients
#' @param ... additional arguments to pass to \code{quantile}
#' @return Numeric vector of parameter point estimates based on the given \code{prob}, with a default of the median estimate.
#' @export
#' @examples
#' # Linear regression example
#' set.seed(521)
#' X <- cbind(1, matrix(rnorm(300), ncol=3))
#' betavals <- c(0.5, -1, 2, -3)
#' y <- X%*%betavals + rnorm(100, sd=.2)
#'
#' f1 <- hmc(N = 500,
#'           theta.init = c(rep(0, 4), 1),
#'           epsilon = 0.01,
#'           L = 10,
#'           logPOSTERIOR = linear_posterior,
#'           glogPOSTERIOR = g_linear_posterior,
#'           varnames = c(paste0("beta", 0:3), "log_sigma_sq"),
#'           param=list(y=y, X=X), parallel=FALSE, chains=1)
#'
#' summary(f1)
#' coef(f1)
coef.hmclearn <- function(object, burnin=NULL, prob=0.5, ...) {
  thetaCombined <- combMatrix(object$thetaCombined, burnin=burnin)
  thetaCombined <- do.call(rbind, thetaCombined)
  apply(thetaCombined, 2, quantile, probs=prob, ...)
}

#' Model Predictions for HMC or MH
#'
#' \code{predict} generates simulated data from the posterior predictive distribution.
#' This simulated data can be used for posterior predictive check diagnostics from the \code{bayesplot} package
#'
#' @param object an object of class \code{hmclearn}, usually a result of a call to \code{mh} or \code{hmc}
#' @param X design matrix, either from fitting the model or new data
#' @param fam generalized linear model family.  Currently "linear", "binomial", and "poisson" are supported
#' @param burnin optional numeric parameter for the number of initial MCMC samples to omit from the summary
#' @param draws Number of simulated values from the posterior conditioned on \code{X}
#' @param ... additional parameters, currently unsupported
#' @return An object of class \code{hmclearnpred}.
#' @section Elements of \code{hmclearnpred} objects:
#' \describe{
#'   \item{\code{y}}{
#'   Median simulated values for each observation in \code{X}
#'   }
#'   \item{\code{yrep}}{
#'   Matrix of simulated values where each row is a draw from the posterior predictive distribution
#'   }
#'   \item{\code{X}}{
#'   Numeric design matrix
#'   }
#' }
#' @references Gabry, Jonah and Mahr, Tristan (2019).  \emph{bayesplot:  Plotting for Bayesian Models}.  \url{https://mc-stan.org/bayesplot/}
#' @references Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., and Gelman, A (2019).  \emph{Visualization in Bayesian Workflow}.  Journal of the Royal Statistical Society: Series A. Vol 182.  Issue 2.  p.389-402.
#' @export
#'
#' @examples
#' # Linear regression example
#' set.seed(521)
#' X <- cbind(1, matrix(rnorm(300), ncol=3))
#' betavals <- c(0.5, -1, 2, -3)
#' y <- X%*%betavals + rnorm(100, sd=.2)
#'
#' f1 <- hmc(N = 500,
#'           theta.init = c(rep(0, 4), 1),
#'           epsilon = 0.01,
#'           L = 10,
#'           logPOSTERIOR = linear_posterior,
#'           glogPOSTERIOR = g_linear_posterior,
#'           varnames = c(paste0("beta", 0:3), "log_sigma_sq"),
#'           param=list(y=y, X=X), parallel=FALSE, chains=1)
#'
#' summary(f1)
#'
#' p <- predict(f1, X)
#' predvals <- p$y
#' plot(predvals, y, xlab="predicted", ylab="actual")
#'
#' X2 <- cbind(1, matrix(rnorm(30), ncol=3))
#' p2 <- predict(f1, X2)
#' p2$y
#'
predict.hmclearn <- function(object, X, fam = "linear", burnin=NULL, draws=NULL, ...) {

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

  retval <- list(y = apply(preds, 2, stats::median),
                 yreps = preds,
                 X = X)
  class(retval) <- c("hmclearnpred", "list")

  return(retval)
}

#'
#' #' Posterior (or prior) predictive checks based on \code{bayesplot} package
#' #'
#' #' Facilitates the display of several predictive checks from \code{bayesplot}
#' #'
#' #' @param object an object of class \code{hmclearn}, usually a result of a call to \code{mh} or \code{hmc}
#' #' @param type character identification of \code{pp_check} plot.  Available types are "overlaid", "multiple", "bars", and "stat".
#' #' @param ... optional additional argument passed to the \code{bayesplot} functions
#' #'
#' #' @section Plot/Data Descriptions from \code{bayesplot}:
#' #' \describe{
#' #'   \item{overlaid:  `ppc_dens_overlay()`}{
#' #'    Kernel density or empirical CDF estimates of each dataset (row) in
#' #'    `yrep` are overlaid, with the distribution of `y` itself on top
#' #'    (and in a darker shade).
#' #'   }
#' #'   \item{multiple: `ppc_hist()`}{
#' #'    A separate histogram, shaded frequency polygon, smoothed kernel density
#' #'    estimate, or box and whiskers plot is displayed for `y` and each
#' #'    dataset (row) in `yrep`. For these plots `yrep` should therefore
#' #'    contain only a small number of rows. See the **Examples** section.
#' #'   }
#' #'   \item{`ppc_bars()`}{
#' #'    For discrete distributions.  Bar plot of `y` with `yrep` medians and uncertainty intervals
#' #'    superimposed on the bars.
#' #'   }
#' #'   \item{`ppc_stat()`}{
#' #'    A histogram of the distribution of a test statistic computed by applying
#' #'    `stat` to each dataset (row) in `yrep`. The value of the statistic in the
#' #'    observed data, `stat(y)`, is overlaid as a vertical line.
#' #'   }
#' #' }
#' #' @references Gabry, Jonah and Mahr, Tristan (2019).  \emph{bayesplot:  Plotting for Bayesian Models}.  \url{https://mc-stan.org/bayesplot/}
#' #' @references Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., and Gelman, A (2019).  \emph{Visualization in Bayesian Workflow}.  Journal of the Royal Statistical Society: Series A. Vol 182.  Issue 2.  p.389-402.
#' #' @export
#' pp_check <- function(object, type, ...) {
#'   UseMethod("pp_check")
#' }
#'
#' #' Posterior (or prior) predictive checks based on \code{bayesplot} package
#' #'
#' #' Facilitates the display of several predictive checks from \code{bayesplot}
#' #'
#' #' @param object an object of class \code{hmclearn}, usually a result of a call to \code{mh} or \code{hmc}
#' #' @param type character identification of \code{pp_check} plot.  Available types are "overlaid", "multiple", "bars", and "stat".
#' #' @param ... optional additional argument passed to the \code{bayesplot} functions
#' #'
#' #' @section Plot/Data Descriptions from \code{bayesplot}:
#' #' \describe{
#' #'   \item{overlaid:  `ppc_dens_overlay()`}{
#' #'    Kernel density or empirical CDF estimates of each dataset (row) in
#' #'    `yrep` are overlaid, with the distribution of `y` itself on top
#' #'    (and in a darker shade).
#' #'   }
#' #'   \item{multiple: `ppc_hist()`}{
#' #'    A separate histogram, shaded frequency polygon, smoothed kernel density
#' #'    estimate, or box and whiskers plot is displayed for `y` and each
#' #'    dataset (row) in `yrep`. For these plots `yrep` should therefore
#' #'    contain only a small number of rows. See the **Examples** section.
#' #'   }
#' #'   \item{`ppc_bars()`}{
#' #'    For discrete distributions.  Bar plot of `y` with `yrep` medians and uncertainty intervals
#' #'    superimposed on the bars.
#' #'   }
#' #'   \item{`ppc_stat()`}{
#' #'    A histogram of the distribution of a test statistic computed by applying
#' #'    `stat` to each dataset (row) in `yrep`. The value of the statistic in the
#' #'    observed data, `stat(y)`, is overlaid as a vertical line.
#' #'   }
#' #' }
#' #' @references Gabry, Jonah and Mahr, Tristan (2019).  \emph{bayesplot:  Plotting for Bayesian Models}.  \url{https://mc-stan.org/bayesplot/}
#' #' @references Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., and Gelman, A (2019).  \emph{Visualization in Bayesian Workflow}.  Journal of the Royal Statistical Society: Series A. Vol 182.  Issue 2.  p.389-402.
#' #' @export
#' pp_check.hmclearnpred <- function(object,
#'                                   type=c("overlaid",
#'                                          "multiple",
#'                                          "bars",
#'                                          "stat"), ...) {
#'
#'   y <- object[["y"]]
#'   yrep <- object[["yrep"]]
#'   switch(match.arg(type),
#'          multiple = bayesplot::ppc_hist(y, yrep[1:min(8, nrow(yrep)),, drop = FALSE]),
#'          overlaid = bayesplot::ppc_dens_overlay(y, yrep),
#'          bars = bayesplot::ppc_bars(y, yrep, ...),
#'          stat = bayesplot::ppc_stat(y, yrep, ...)
#'   )
#'
#' }
