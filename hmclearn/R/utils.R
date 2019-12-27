

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

