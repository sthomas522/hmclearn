
# get pre-stored logistic regression MCMCpack data


#' @importFrom stats aggregate
#' @importFrom stats quantile
#' @importFrom stats rbinom
#' @importFrom stats reshape
#' @importFrom stats rnorm
#' @importFrom stats rpois
#' @importFrom stats runif
#' @importFrom stats var
#' @importFrom utils tail
#' @importFrom mvtnorm dmvnorm
#' @importFrom MASS mvrnorm

# random effects covariance
create_Uj <- function(uj, neg=TRUE) {
  q <- length(uj)
  if (q == 1) return(0)
  nr <- q
  nc <- q*(q-1)/2
  Uj <- matrix(0, nrow = nr, ncol = nc)
  for (kk in 2:nr)
  { Uj[kk, ((kk-1)*(kk-2)/2 + 1):(kk*(kk-1)/2)] <- uj[1:(kk-1)] }
  if (neg) {
    return(-Uj)
  } else {
    return(Uj)
  }
}


#' Diagnostic plots for \code{hmclearn}
#'
#' Plots histograms of the posterior estimates.  Optionally, displays the 'actual'
#' values given a simulated dataset.
#'
#' @param object an object of class \code{hmclearn}, usually a result of a call to \code{mh} or \code{hmc}
#' @param burnin optional numeric parameter for the number of initial MCMC samples to omit from the summary
#' @param plotfun integer 1 or 2 indicating which plots to display.  1 shows trace plots.  2 shows a histogram
#' @param comparison.theta optional numeric vector of true parameter values
#' @param cols optional integer index indicating which parameters to display
#' @param ... currently unused
#' @return Returns a customized \code{ggplot} object
#' @export
#' @examples
#' # Linear regression example
#' set.seed(522)
#' X <- cbind(1, matrix(rnorm(300), ncol=3))
#' betavals <- c(0.5, -1, 2, -3)
#' y <- X%*%betavals + rnorm(100, sd=.2)
#'
#' f <- hmc(N = 1000,
#'           theta.init = c(rep(0, 4), 1),
#'           epsilon = 0.01,
#'           L = 10,
#'           logPOSTERIOR = linear_posterior,
#'           glogPOSTERIOR = g_linear_posterior,
#'           varnames = c(paste0("beta", 0:3), "log_sigma_sq"),
#'           param=list(y=y, X=X), parallel=FALSE, chains=1)
#'
#' diagplots(f, burnin=300, comparison.theta=c(betavals, 2*log(.2)))
#'
diagplots <- function(object, burnin=NULL, plotfun=2, comparison.theta=NULL, cols=NULL, ...) {
  UseMethod("diagplots")
}

#' Diagnostic plots for \code{hmclearn}
#'
#' Plots histograms of the posterior estimates.  Optionally, displays the 'actual'
#' values given a simulated dataset.
#'
#' @param object an object of class \code{hmclearn}, usually a result of a call to \code{mh} or \code{hmc}
#' @param burnin optional numeric parameter for the number of initial MCMC samples to omit from the summary
#' @param plotfun integer 1 or 2 indicating which plots to display.  1 shows trace plots.  2 shows a histogram
#' @param comparison.theta optional numeric vector of parameter values to compare to the Bayesian estimates
#' @param cols optional integer index indicating which parameters to display
#' @param ... currently unused
#' @return Returns a customized \code{ggplot} object
#' @export
#' @examples
#' # Linear regression example
#' set.seed(522)
#' X <- cbind(1, matrix(rnorm(300), ncol=3))
#' betavals <- c(0.5, -1, 2, -3)
#' y <- X%*%betavals + rnorm(100, sd=.2)
#'
#' f <- hmc(N = 1000,
#'           theta.init = c(rep(0, 4), 1),
#'           epsilon = 0.01,
#'           L = 10,
#'           logPOSTERIOR = linear_posterior,
#'           glogPOSTERIOR = g_linear_posterior,
#'           varnames = c(paste0("beta", 0:3), "log_sigma_sq"),
#'           param=list(y=y, X=X), parallel=FALSE, chains=1)
#'
#' diagplots(f, burnin=300, comparison.theta=c(betavals, 2*log(.2)))
#'
diagplots.hmclearn <- function(object, burnin=NULL, plotfun=2, comparison.theta=NULL, cols=NULL, ...) {

  data <- combMatrix(object$thetaCombined, burnin=burnin)
  data <- do.call(rbind, data)

  if (is.null(cols)) {
    cols <- 1:ncol(data)
  }

  thetaCombinedsubs <- data[, cols]

  # if (!is.null(burnin)) {
  #   thetaCombinedsubs <- data[-c(1:burnin), cols]
  # } else {
  #   thetaCombinedsubs <- data[, cols]
  # }

  pdata <- as.data.frame(thetaCombinedsubs)
  pdata$t <- 1:nrow(pdata)
  pdata <- reshape(pdata,
                   varying = list(1:(ncol(pdata)-1)),
                   v.names = "value",
                   idvar = "t",
                   timevar = "coefficient",
                   times = colnames(pdata)[-ncol(pdata)],
                   direction = "long")
  pdata$true.mu <- rep(comparison.theta, each=nrow(thetaCombinedsubs))
  pdata$coefficient <- as.factor(pdata$coefficient)

  # k <- ncol(data)
  k <- length(cols)

  # return list

  # line plots of simulation
  p1 <- NULL
  if (1 %in% plotfun) {
    p1 <- ggplot2::ggplot(pdata, ggplot2::aes_string(x="t", y="value", colour="coefficient")) + ggplot2::geom_line()
    p1 <- p1 + ggplot2::facet_wrap(~ coefficient, ncol=trunc(sqrt(k)), scales="free_y")
    p1 <- p1 + ggplot2::theme_bw()
    p1
  }

  # histograms
  p2 <- NULL
  if (2 %in% plotfun) {
    p2 <- ggplot2::ggplot(pdata, ggplot2::aes_string(x="value", y="..density..", fill="coefficient",
                                              colour="coefficient")) +
      ggplot2::geom_histogram(bins=40)

    if (!is.null(comparison.theta)) {
      p2 <- p2 + ggplot2::geom_vline(data=aggregate(pdata[4], pdata[2], mean),
                                     mapping=ggplot2::aes_string(xintercept = "true.mu"), colour="red")
    }

    p2 <- p2 + ggplot2::facet_wrap(~ coefficient, ncol=trunc(sqrt(k)), scales="free")

    p2 <- p2 + ggplot2::theme_bw()
    p2
  }

  allplots <- NULL
  if (1 %in% plotfun) {
    allplots <- c(allplots, list(trace=p1))
  }
  if (2 %in% plotfun) {
    allplots <- c(allplots, list(histogram=p2))
  }

  allplots



}


