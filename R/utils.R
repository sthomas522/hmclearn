
#' @importFrom stats aggregate
#' @importFrom stats quantile
#' @importFrom stats rbinom
#' @importFrom stats reshape
#' @importFrom stats rnorm
#' @importFrom stats rpois
#' @importFrom stats runif
#' @importFrom stats var
#' @importFrom utils tail

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
#' @param actual.mu optional numeric vector of true parameter values
#' @param cols optional integer index indicating which parameters to display
#' @param ... currently unused
#' @export
diagplots <- function(object, burnin=NULL, plotfun=2, actual.mu=NULL, cols=NULL, ...) {
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
#' @param actual.mu optional numeric vector of true parameter values
#' @param cols optional integer index indicating which parameters to display
#' @param ... currently unused
#' @export
diagplots.hmclearn <- function(object, burnin=NULL, plotfun=2, actual.mu=NULL, cols=NULL, ...) {

  data <- combMatrix(object$thetaCombined, burnin=burnin)
  data <- do.call(rbind, data)

  if (is.null(cols)) {
    cols <- 1:ncol(data)
  }

  if (!is.null(burnin)) {
    thetaCombinedsubs <- data[-c(1:burnin), cols]
  } else {
    thetaCombinedsubs <- data[, cols]
  }

  pdata <- as.data.frame(thetaCombinedsubs)
  pdata$t <- 1:nrow(pdata)
  pdata <- reshape(pdata,
                   varying = list(1:(ncol(pdata)-1)),
                   v.names = "value",
                   idvar = "t",
                   timevar = "coefficient",
                   times = colnames(pdata)[-ncol(pdata)],
                   direction = "long")
  pdata$true.mu <- rep(actual.mu, each=nrow(thetaCombinedsubs))
  pdata$coefficient <- as.factor(pdata$coefficient)

  k <- ncol(data)

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

    if (!is.null(actual.mu)) {
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




