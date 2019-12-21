
#######################################################################################
# metropolis hastings functions
#

# metropolis algorithm
# N:  number of iterations
# param.init:  initial value in chain
# qPROP:  function to generate proposal
# qFUN:  probability for proposal function.
#    first argument is where to evaluate, and second argument is the conditional parameter
# pFUN:  function to calculate probabilities of proposals (e.g. log likelihood)
# ... additional parameters for pFUN
#' @export
mh <- function(N, paramInit, qPROP, qFUN, pdFUN, nu, ...) {
  paramSim <- list()
  paramSim[[1]] <- paramInit
  accept <- 0
  for (j in 2:N) {
    u <- runif(1)
    paramProposal <- qPROP(paramSim[[j-1]], nu)
    lnum <- pdFUN(paramProposal, ...) + qFUN(paramSim[[j-1]], paramProposal, nu)
    lden <- pdFUN(paramSim[[j-1]], ...) + qFUN(paramProposal, paramSim[[j-1]], nu)
    l.alpha <- pmin(0, lnum - lden)
    if (l.alpha > log(u)) {
      paramSim[[j]] <- paramProposal
      accept <- accept + 1
    } else {
      paramSim[[j]] <- paramSim[[j-1]]
    }
  }

  # create dataframe from simulation
  thetaDF <- as.data.frame(do.call(rbind, paramSim))

  obj <- list(paramSim = paramSim,
       thetaDF = thetaDF,
       accept = accept)

  class(obj) <- c("hmclearn", "list")
  return(obj)
}


# q(theta1 | theta2, nu) where nu is a kxk positive definite matrix
#' @export
qfun <- function(theta1, theta2, nu) {
  k <- length(theta1)
  nu <- diag(nu, k, k)
  mvtnorm::dmvnorm(theta1, theta2, nu, log=TRUE)
}

# sample from proposal density
#' @export
qprop <- function(theta1, nu) {
  k <- length(theta1)
  require(MASS)
  nu <- diag(nu, k, k)
  MASS::mvrnorm(1, theta1, nu)
}

# nu is a vector
#' @export
qprop_all <- function(theta1, nu) {
  nu <- diag(nu)
  MASS::mvrnorm(1, theta1, nu)
}

#' @export
qfun_all <- function(theta1, theta2, nu) {
  nu <- diag(nu)
  mvtnorm::dmvnorm(theta1, theta2, nu, log=TRUE)
}





###############################################
# functions to support
# Hamiltonian Monte Carlo
#
#

# leapfrog integrator from Hamiltonian dynamics
# theta:  parameter of interest
# r:  momentum variable
# epsilon:  step size parameter
# logDENS:  log of joint density of parameter of interest
#   (log likelihood)
# ... additional parameters to pass to logDENS
#' @export
leapfrog <- function(theta_lf, r, epsilon, logPOSTERIOR, glogPOSTERIOR, Minv, constrain,
                     lastSTEP=FALSE, ...) {

  # gradient of log posterior for old theta
  g.ld <- glogPOSTERIOR(theta_lf,  ...)

  # first momentum update
  r.new <- as.numeric(r + epsilon/2*g.ld)

  # theta update
  theta.new <- theta_lf + as.numeric(epsilon* Minv %*% r.new)

  # check positive
  switch_sign <- constrain & theta.new < 0
  r.new[switch_sign] <- -r.new[switch_sign]
  theta.new[switch_sign] <- -theta.new[switch_sign]

  # gradient of log posterior for new theta
  g.ld.new <- glogPOSTERIOR(theta.new, ...)

  # if not on last step, second momentum update
  if (!lastSTEP) {
    r.new <- r.new + epsilon/2*g.ld.new
  }
  list(theta.new=theta.new,
       r.new=as.numeric(r.new))
}



# theta.init:  initial values of theta
# Nstep:  number of steps (called L in paper)
# N:  number of times to repeat
# epsilon:  step size
# logPOSTERIOR:  log of joint density of parameter of interest
# ...:  additional parameters to pass to logPOSTERIOR


#' Fit a generic model using Hamiltonian Monte Carlo (HMC)
#'
#' This function runs the HMC algorithm on a generic model provided
#' the \code{logPOSTERIOR} and gradient \code{glogPOSTERIOR} functions.
#' All parameters specified in ... are passed to these two functions.
#' The tuning parameters \code{epsilon} and \code{L} are passed to the
#' Leapfrog algorithm.
#'
#' @param N Number of MCMC samples
#' @param theta.init Vector of initial values for the parameters
#' @param epsilon Step-size parameter for \code{leapfrog}
#' @param L Number of \code{leapfrog} steps parameter
#' @param logPOSTERIOR Function to calculate and return the log posterior given a vector of values of \code{theta}
#' @param glogPOSTERIOR Function to calculate and return the gradient of the log posterior given a vector of values of  \code{theta}
#' @param randlength Logical to determine whether to apply some randomness to the number of leapfrog steps tuning parameter \code{L}
#' @param Mdiag Optional vector of the diagonal of the mass matrix \code{M}.  Defaults to unit diagonal.
#' @param constrain Optional vector of which parameters in \code{theta} accept positive values only.  Default is that all parameters accept all real numbers
#' @param verbose Logical to determine whether to display the progress of the HMC algorithm
#' @return Object of class \code{hmclearn}
#'
#' @section Elements for \code{hmclearn} objects:
#' \describe{
#'   \item{\code{N}}{
#'   Number of MCMC samples
#'   }
#'   \item{\code{theta}}{
#'   List of length \code{N} of the sampled values of \code{theta}
#'   }
#'   \item{\code{thetaDF}}{
#'   Sampled values in dataframe form
#'   }
#'   \item{\code{r}}{
#'   List of length \code{N} of the sampled momenta
#'   }
#'   \item{\code{theta.all}}{
#'   List of all parameter values of \code{theta} sampled prior to accept/reject step
#'   }
#'   \item{\code{r.all}}{
#'   List of all values of the momenta \code{r} sampled prior to accept/reject
#'   }
#'   \item{\code{accept}}{
#'   Number of accepted proposals.  The ratio \code{accept} / \code{N} is the acceptance rate
#'   }
#'   \item{\code{accept_v}}{
#'   Vector of length \code{N} indicating which samples were accepted
#'   }
#'   \item{\code{M_mx}}{
#'   Mass matrix used in the HMC algorithm
#'   }
#' }
#'
#' @section Available \code{logPOSTERIOR} and \code{glogPOSTERIOR} functions:
#' \describe{
#'   \item{\code{linear_posterior}}{
#'   Linear regression:  log posterior
#'   }
#'   \item{\code{g_linear_posterior}}{
#'   Linear regression:  gradient of the log posterior
#'   }
#'   \item{\code{logistic_posterior}}{
#'   Logistic regression:  log posterior
#'   }
#'   \item{\code{g_logistic_posterior}}{
#'   Logistic regression:  gradient of the log posterior
#'   }
#'   \item{\code{poisson_posterior}}{
#'   Poisson (count) regression:  log posterior
#'   }
#'   \item{\code{g_poisson_posterior}}{
#'   Poisson (count) regression: gradient of the log posterior
#'   }
#' }
#'
#' @author Samuel Thomas \email{samthoma@@iu.edu}, Wanzhu Tu \email{wtu@iu.edu}
#' @references \emph{HMC in R} paper
#' @references Thomas, S., Li, X., and Tu, W.  2019.  \emph{Hamiltonian Monte Carlo}.  Wiley
#' @references Neal, Radford. 2011. \emph{MCMC Using Hamiltonian Dynamics.} In Handbook of Markov Chain Monte Carlo, edited by Steve Brooks, Andrew Gelman, Galin L. Jones, and Xiao-Li Meng, 116–62. Chapman; Hall/CRC.
#' @keywords hamiltonian monte carlo
hmc <- function(N=10000, theta.init, epsilon=1e-2, L=10, logPOSTERIOR, glogPOSTERIOR, varnames=NULL,
                randlength=FALSE, Mdiag=NULL, constrain=NULL, verbose=FALSE, ...) {

  p <- length(theta.init) # number of parameters

  # mass matrix
  mu.p <- rep(0, p)

  # epsilon values
  eps_orig <- epsilon
  if (length(epsilon) == 1) {
    eps_orig <- rep(epsilon, p)
  }

  # epsilon matrix
  eps_vals <- matrix(rep(eps_orig, N), ncol=N, byrow=F)

  # number of steps
  L_vals <- rep(L, N)

  # randomize epsilon and L
  if (randlength) {
    randvals <- replicate(N, runif(p, -0.1*eps_orig, 0.1*eps_orig), simplify=T)
    eps_vals <- eps_vals + randvals
    L_vals <- round(runif(N, 0.5*L, 2.0*L))
  }

  # invert covariance M for leapfrog
  if (is.null(Mdiag)) {
    M_mx <- diag(p)
    Minv <- M_mx
  } else {
    M_mx <- diag(Mdiag)
    Minv <- diag(1 / Mdiag)
  }

  # print(diag(Minv))

  # store all momentum and theta values
  iter.all <- 1
  theta.all <- list()
  r.all <- list()

  # store theta and momentum (usually not of interest)
  theta <- list()
  theta[[1]] <- theta.init
  r <- list()
  r[[1]] <- NA
  accept <- 0
  accept_v <- vector()
  accept_v[1] <- 1
  for (jj in 2:N) {
    theta[[jj]] <- theta.new <- theta[[jj-1]]
    r0 <- MASS::mvrnorm(1, mu.p, M_mx)
    r.new <- r[[jj]] <- r0

    for (i in 1:L_vals[jj]) {
      theta.all[[iter.all]] <- theta.new
      r.all[[iter.all]] <- r.new
      iter.all <- iter.all + 1

      lstp <- i == L_vals[jj]
      lf <- leapfrog(theta_lf = theta.new, r = r.new, epsilon = eps_vals[, jj], logPOSTERIOR = logPOSTERIOR,
                     glogPOSTERIOR = glogPOSTERIOR,
                     Minv=Minv, constrain=constrain, lastSTEP=lstp, ...)

      theta.new <- lf$theta.new
      r.new <- lf$r.new

    }
    if (verbose) print(jj)

    # standard metropolis-hastings update
    u <- runif(1)

    # use log transform for ratio due to low numbers
    num <- logPOSTERIOR(theta.new,  ...) - 0.5*(r.new %*% r.new)
    den <- logPOSTERIOR(theta[[jj-1]], ...) - 0.5*(r0 %*% r0)

    log.alpha <- pmin(0, num - den)

    if (log(u) < log.alpha) {
      theta[[jj]] <- theta.new
      r[[jj]] <- -r.new
      accept <- accept + 1
      accept_v <- c(accept_v, 1)
    } else {
      theta[[jj]] <- theta[[jj-1]]
      r[[jj]] <- r[[jj-1]]
      accept_v <- c(accept_v, 0)
    }

  }

  # create dataframe from simulation
  thetaDF <- as.data.frame(do.call(rbind, theta))

  if (!is.null(varnames)) {
    colnames(thetaDF) <- varnames
  }

  obj <- list(N=N,
       theta=theta,
       thetaDF = thetaDF,
       r=r,
       theta.all = theta.all,
       r.all = r.all,
       accept=accept,
       accept_v = accept_v,
       M=M_mx)

  class(obj) <- c("hmclearn", "list")
  return(obj)
}

#' @export
diagplots <- function(result, actual.mu=NULL, burnin=100, cols=NULL) {

  if (is.null(cols)) {
    cols <- 1:ncol(result$thetaDF)
  }

  thetaDFsubs <- result$thetaDF[-c(1:burnin), cols]
  pdata <- thetaDFsubs
  pdata$t <- 1:nrow(pdata)
  pdata <- reshape(pdata,
                   varying = list(1:(ncol(pdata)-1)),
                   v.names = "value",
                   idvar = "t",
                   timevar = "coefficient",
                   times = colnames(pdata)[-ncol(pdata)],
                   direction = "long")
  pdata$true.mu <- rep(actual.mu, each=nrow(thetaDFsubs))

  k <- ncol(result$thetaDF)

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
summary.hmclearn <- function(x, burnin=100, probs=c(0.05, 0.25, 0.5, 0.75, 0.95), ...) {
  cat("Summary of HMC simulation\n\n")

  thetaDF <- x$thetaDF[-c(1:burnin), ]
  t(apply(thetaDF, 2, quantile, probs=probs, ...))
}

#' @export
plot.hmclearn <- function(x, ...) {
  diagplots(x, ...)
}

