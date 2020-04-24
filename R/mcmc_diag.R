

#' Calculates Potential Scale Reduction Factor (psrf), also called the Rhat statistic,
#' from models fit via \code{mh} or \code{hmc}
#'
#' Gelman and Rubin's diagnostic assesses the mix of multiple MCMC chain with different initial parameter values
#' Values close to 1 indicate that the posterior simulation has sufficiently converged, while
#' values above 1 indicate that additional samples may be necessary to ensure convergence.  A general
#' guideline suggests that values less than 1.05 are good, between 1.05 and 1.10 are ok, and above 1.10
#' have not converged well.
#'
#' @param object an object of class \code{hmclearn}, usually a result of a call to \code{mh} or \code{hmc}
#' @param burnin optional numeric parameter for the number of initial MCMC samples to omit from the summary
#' @param ... currently unused
#'
#' @references Gelman, A. and Rubin, D. (1992) \emph{Inference from Iterative Simulation Using Multiple Sequences}.  Statistical Science 7(4) 457-472.
#' @references Gelman, A., et. al. (2013) \emph{Bayesian Data Analysis}.  Chapman and Hall/CRC.
#' @references Gabry, Jonah and Mahr, Tristan (2019).  \emph{bayesplot:  Plotting for Bayesian Models}.  \url{https://mc-stan.org/bayesplot}
#' @return Numeric vector of Rhat statistics for each parameter
#' @export
#'
#' @examples
#' # poisson regression example
#' set.seed(7363)
#' X <- cbind(1, matrix(rnorm(40), ncol=2))
#' betavals <- c(0.8, -0.5, 1.1)
#' lmu <- X %*% betavals
#' y <- sapply(exp(lmu), FUN = rpois, n=1)
#'
#' f <- hmc(N = 1000,
#'           theta.init = rep(0, 3),
#'           epsilon = 0.01,
#'           L = 10,
#'           logPOSTERIOR = poisson_posterior,
#'           glogPOSTERIOR = g_poisson_posterior,
#'           varnames = paste0("beta", 0:2),
#'           param = list(y=y, X=X),
#'           parallel=FALSE, chains=2)
#'
#' psrf(f, burnin=100)

psrf <- function(object, burnin, ...) {
  UseMethod("psrf")
}

#' Calculates Potential Scale Reduction Factor (psrf), also called the Rhat statistic,
#' from models fit via \code{mh} or \code{hmc}
#'
#' Gelman and Rubin's diagnostic assesses the mix of multiple MCMC chain with different initial parameter values
#' Values close to 1 indicate that the posterior simulation has sufficiently converged, while
#' values above 1 indicate that additional samples may be necessary to ensure convergence.  A general
#' guideline suggests that values less than 1.05 are good, between 1.05 and 1.10 are ok, and above 1.10
#' have not converged well.
#'
#' @param object an object of class \code{hmclearn}, usually a result of a call to \code{mh} or \code{hmc}
#' @param burnin optional numeric parameter for the number of initial MCMC samples to omit from the summary
#' @param ... currently unused
#'
#' @references Gelman, A. and Rubin, D. (1992) \emph{Inference from Iterative Simulation Using Multiple Sequences}.  Statistical Science 7(4) 457-472.
#' @references Gelman, A., et. al. (2013) \emph{Bayesian Data Analysis}.  Chapman and Hall/CRC.
#' @references Gabry, Jonah and Mahr, Tristan (2019).  \emph{bayesplot:  Plotting for Bayesian Models}.  \url{https://mc-stan.org/bayesplot}
#' @return Numeric vector of Rhat statistics for each parameter
#' @export
#'
#' @examples
#' # poisson regression example
#' set.seed(7363)
#' X <- cbind(1, matrix(rnorm(40), ncol=2))
#' betavals <- c(0.8, -0.5, 1.1)
#' lmu <- X %*% betavals
#' y <- sapply(exp(lmu), FUN = rpois, n=1)
#'
#' f <- hmc(N = 1000,
#'           theta.init = rep(0, 3),
#'           epsilon = 0.01,
#'           L = 10,
#'           logPOSTERIOR = poisson_posterior,
#'           glogPOSTERIOR = g_poisson_posterior,
#'           varnames = paste0("beta", 0:2),
#'           param = list(y=y, X=X),
#'           parallel=FALSE, chains=2)
#'
#' psrf(f, burnin=100)
#'
psrf.hmclearn <- function(object, burnin=NULL, ...) {

  data <- combMatrix(object$thetaCombined, burnin=burnin)

  if (object$chains == 1) {
    # stop("Multiple chains needed to calculate Potential Scale Reduction Factor")
    rhat <- rep(NA, ncol(data[[1]]))
    return(rhat)
  }

  # number of samples after burnin
  Nraw <- object$N
  if (!is.null(burnin)) {
    N <- Nraw - burnin
  } else {
    N <- Nraw
  }

  # variance estimation
  allvar <- varest(data, N)

  rhat <- sqrt(allvar$varest/allvar$W)
  return(rhat)

}



# variance estimator
# data: list of matrices
varest <- function(data, N) {
  M <- length(data)

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
  if (M > 1) {
    varest <- (N-1)/N * W + 1/N*B
  } else {
    varest <- W
  }


  # return variance estimator between and within sequence variance
  retval <- list(B=B,
                 W=W,
                 varest=varest)
  return(retval)

}

#' Effective sample size calculation
#'
#' Calculates an estimate of the adjusted MCMC sample size per parameter
#' adjusted for autocorrelation.
#'
#' @param object an object of class \code{hmclearn}, usually a result of a call to \code{mh} or \code{hmc}
#' @param burnin optional numeric parameter for the number of initial MCMC samples to omit from the summary
#' @param lagmax maximum lag to extract for determining effective sample sizes
#' @param ... currently unused
#'
#' @references Gelman, A., et. al. (2013) \emph{Bayesian Data Analysis}.  Chapman and Hall/CRC.  Section 11.5
#' @return Numeric vector with effective sample sizes for each parameter in the model
#' @export
#' @examples
#' # poisson regression example
#' set.seed(7363)
#' X <- cbind(1, matrix(rnorm(40), ncol=2))
#' betavals <- c(0.8, -0.5, 1.1)
#' lmu <- X %*% betavals
#' y <- sapply(exp(lmu), FUN = rpois, n=1)
#'
#' f <- hmc(N = 1000,
#'           theta.init = rep(0, 3),
#'           epsilon = c(0.03, 0.02, 0.015),
#'           L = 10,
#'           logPOSTERIOR = poisson_posterior,
#'           glogPOSTERIOR = g_poisson_posterior,
#'           varnames = paste0("beta", 0:2),
#'           param = list(y=y, X=X),
#'           parallel=FALSE, chains=2)
#'
#' neff(f, burnin=100)
neff <- function(object, burnin=NULL, lagmax=NULL, ...) {
  UseMethod("neff")
}

#' Effective sample size calculation
#'
#' Calculates an estimate of the adjusted MCMC sample size per parameter
#' adjusted for autocorrelation.
#'
#' @param object an object of class \code{hmclearn}, usually a result of a call to \code{mh} or \code{hmc}
#' @param burnin optional numeric parameter for the number of initial MCMC samples to omit from the summary
#' @param lagmax maximum lag to extract for determining effective sample sizes
#' @param ... currently unused
#'
#' @references Gelman, A., et. al. (2013) \emph{Bayesian Data Analysis}.  Chapman and Hall/CRC.  Section 11.5
#' @return Numeric vector with effective sample sizes for each parameter in the model
#' @export
#' @examples
#' # poisson regression example
#' set.seed(7363)
#' X <- cbind(1, matrix(rnorm(40), ncol=2))
#' betavals <- c(0.8, -0.5, 1.1)
#' lmu <- X %*% betavals
#' y <- sapply(exp(lmu), FUN = rpois, n=1)
#'
#' f <- hmc(N = 1000,
#'           theta.init = rep(0, 3),
#'           epsilon = c(0.03, 0.02, 0.015),
#'           L = 10,
#'           logPOSTERIOR = poisson_posterior,
#'           glogPOSTERIOR = g_poisson_posterior,
#'           varnames = paste0("beta", 0:2),
#'           param = list(y=y, X=X),
#'           parallel=FALSE, chains=2)
#'
#' neff(f, burnin=100)
neff.hmclearn <- function(object, burnin=NULL, lagmax=NULL, ...) {
  data <- combMatrix(object$thetaCombined, burnin=burnin)
  M <- length(data)

  # number of samples after burnin
  Nraw <- object$N
  if (!is.null(burnin)) {
    N <- Nraw - burnin
  } else {
    N <- Nraw
  }

  # variance estimator
  vest <- varest(data, N)$varest

  if (is.null(lagmax)) {
    lagmax <- min(500, N-2)
  }

  # Vt for each chain
  Vtchains <- lapply(data, function(xx) {
    apply(xx, 2, variog, lagmax=lagmax)
  })

  # elementwise mean over chains
  Vt <- Reduce("+", Vtchains) / M

  # estimated autocorrelations
  rhat <- 1 -Vt %*% diag(1 / vest / 2)

  # partial sums
  psum <- apply(rhat, 2, partialsum)

  neff <- M*N / (1 + 2*psum)

  floor(neff)

}
#
#   # estimate lagmax
#   # lagmax <- ceiling(10 * log(N))
#   # lagmax <- 2000
#   lagmax <- min(500, N-2)
#
#   x <- data[[1]][, 1]
#   r = vector()
#   for (k in 0:lagmax){
#     r[k+1] = cor(x[1:(N-k)],x[(1+k):N])
#   }
#
#   G = r[2:(lagmax+1)] + r[1:lagmax]
#
#   # estimate autocorrelation
#   tauf = -r[1] + 2*G[1]
#   for (kk in 1:(lagmax-1)){
#     if (G[kk+1]< G[kk] & G[kk+1]>0){
#       tauf = tauf + 2*G[kk+1]
#     } else {
#       break
#     }
#   }
#
#   # Vt for each chain
#   Vtchains <- lapply(data, function(xx) {
#     t(sapply(1:lagmax, function(t) {
#       apply(xx, 2, variog, lagval=t)
#     }))
#   })
#
#   # elementwise mean over chains
#   Vt <- Reduce("+", Vtchains) / M
#
#   # estimated autocorrelations
#   rhat <- 1 -Vt %*% diag(1 / vest / 2)
#
#   M*N / (1 + 2* colSums(rhat))

partialsum <- function(rhat, lagmax=NULL) {
  if (is.null(lagmax)) {
    lagmax <- length(rhat)
  }

  chk <- rhat + c(0, rhat[2:length(rhat)])
  maxt <- min(lagmax, which(chk < 0))
  if (maxt %% 2 == 0) {
    maxt <- maxt + 1
  }

  sum(rhat[1:maxt], na.rm=T)
}


# Variogram 11.6 partial in BDA3
variog <- function(x, lagmax=NULL) {
  n <- length(x)

  if (is.null(lagmax)) {
    lagvals <- seq_len(min(500, n-2))
  } else {
    lagvals <- seq_len(lagmax)
  }

  res <- sapply(lagvals, function(zz) {
    sum(diff(x, lag=zz)^2) / (n - zz)
  })
  res

}


# ess <- function(x) {
#   N <- length(x)
#   V <- purrr::map_dbl(seq_len(N - 1),
#                function(t) {
#                  mean(diff(x, lag = t) ^ 2, na.rm = TRUE)
#                })
#   rho <- purrr::head_while(1 - V / var(x), ~ . > 0)
#   N / (1 + sum(rho))
# }



