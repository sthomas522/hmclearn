

#' @export
diagplots <- function(result, actual.mu=NULL, burnin=100, cols=NULL) {

  if (is.null(cols)) {
    cols <- 1:ncol(result$thetaCombined)
  }

  thetaCombinedsubs <- result$thetaCombined[-c(1:burnin), cols]
  pdata <- thetaCombinedsubs
  pdata$t <- 1:nrow(pdata)
  pdata <- reshape(pdata,
                   varying = list(1:(ncol(pdata)-1)),
                   v.names = "value",
                   idvar = "t",
                   timevar = "coefficient",
                   times = colnames(pdata)[-ncol(pdata)],
                   direction = "long")
  pdata$true.mu <- rep(actual.mu, each=nrow(thetaCombinedsubs))

  k <- ncol(result$thetaCombined)

  # return list

  # line plots of simulation
  p1 <- ggplot2::ggplot(pdata, ggplot2::aes(t, value, colour=factor(coefficient))) + ggplot2::geom_line()
  p1 <- p1 + ggplot2::facet_wrap(~ coefficient, ncol=trunc(sqrt(k)), scales="free_y")
  p1 <- p1 + ggplot2::theme_bw()
  p1

  # histograms
  p2 <- ggplot2::ggplot(pdata, ggplot2::aes(x=value, y=..density.., fill=factor(coefficient),
                                            colour=factor(coefficient))) +
    ggplot2::geom_histogram(bins=40)

  if (!is.null(actual.mu)) {
    p2 <- p2 + ggplot2::geom_vline(data=aggregate(pdata[4], pdata[2], mean),
                                   mapping=ggplot2::aes(xintercept = true.mu), colour="red")
  }

  p2 <- p2 + ggplot2::facet_wrap(~ coefficient, ncol=trunc(sqrt(k)), scales="free")

  p2 <- p2 + ggplot2::theme_bw()
  p2


  list(p1, p2)
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
  varest <- (N-1)/N * W + 1/N*B

  # return variance estimator between and within sequence variance
  retval <- list(B=B,
                 W=W,
                 varest=varest)
  return(retval)

}

#' @export
neff <- function(object, ...) {
  UseMethod("neff")
}

#' @export
neff.hmclearn <- function(object, burnin=NULL, ...) {
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

  # estimate lagmax
  lagmax <- ceiling(10 * log(N))
  # lagmax <- 2000

  # Vt for each chain
  Vtchains <- lapply(data, function(xx) {
    t(sapply(1:lagmax, function(t) {
      apply(xx, 2, variog, lagval=t)
    }))
  })

  # elementwise mean over chains
  Vt <- Reduce("+", Vtchains) / M

  # estimated autocorrelations
  rhat <- 1 -Vt %*% diag(1 / vest / 2)

}

# Variogram 11.6 partial in BDA3
variog <- function(x, lagval=1) {
  n <- length(x)
  # startval <- lagval+1
  # vt <- 0
  # for (i in startval:n) {
  #   # print(paste("i", i, (x[i] - x[i-lagval])^2))
  #   vt <- vt + (x[i] - x[i-lagval])^2
  # }
  # vt / (n-lagval)

  sum(diff(x, lag=lagval)^2) / (n - lagval)

}






