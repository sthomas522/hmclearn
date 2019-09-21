
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
  list(paramSim = paramSim,
       accept = accept)
}


# q(theta1 | theta2, nu) where nu is a kxk positive definite matrix
qfun <- function(theta1, theta2, nu) {
  k <- length(theta1)
  nu <- diag(nu, k, k)
  mvtnorm::dmvnorm(theta1, theta2, nu, log=TRUE)
}

# sample from proposal density
qprop <- function(theta1, nu) {
  k <- length(theta1)
  require(MASS)
  nu <- diag(nu, k, k)
  MASS::mvrnorm(1, theta1, nu)
}

# nu is a vector
qprop_all <- function(theta1, nu) {
  nu <- diag(nu)
  MASS::mvrnorm(1, theta1, nu)
}

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
leapfrog <- function(theta_lf, r, epsilon, logPOSTERIOR, glogPOSTERIOR, y, X, Z, Minv, constrain,
                     lastSTEP=FALSE, ...) {

  # gradient of log posterior for old theta
  g.ld <- glogPOSTERIOR(theta_lf, y=y, X=X, Z=Z, ...)

  # first momentum update
  r.new <- r + epsilon/2*g.ld

  # theta update
  # note diagonal matrix update
  # theta.new <- theta_lf + as.numeric(epsilon*r.new/ diag(M_mx))
  theta.new <- theta_lf + as.numeric(epsilon* Minv %*% r.new)

  # check positive
  switch_sign <- constrain & theta.new < 0
  r.new[switch_sign] <- -r.new[switch_sign]
  theta.new[switch_sign] <- -theta.new[switch_sign]

  # gradient of log posterior for new theta
  g.ld.new <- glogPOSTERIOR(theta.new, y=y, X=X, Z=Z, ...)

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
hmc <- function(N, theta.init, epsilon, L, logPOSTERIOR, glogPOSTERIOR, y, X, Mdiag=NULL,
                verbose=FALSE, ...) {
  p <- length(theta.init)
  mu.p <- rep(0, p)

  if (is.null(Mdiag)) {
    Mdiag <- rep(1, length(theta.init))
  }

  # store theta and momentum (usually not of interest)
  theta <- list()
  theta[[1]] <- theta.init
  r <- list()
  r[[1]] <- NA
  accept <- 0
  for (m in 2:N) {
    theta[[m]] <- theta.new <- theta[[m-1]]
    # r0 <- MASS::mvrnorm(1, mu.p, diag(p))
    r0 <- mapply(rnorm, mean=mu.p, sd=Mdiag, MoreArgs=list(n=1))
    r.new <- r[[m]] <- r0
    for (i in 1:L) {
      lf <- leapfrog(theta.new, r.new, epsilon, logPOSTERIOR, glogPOSTERIOR, y, X, ...)
      theta.new <- lf$theta.new
      r.new <- lf$r.new
    }

    if (verbose) print(m)

    # standard metropolis-hastings update
    u <- runif(1)

    # use log transform for ratio due to low numbers
    log.alpha <- pmin(0, logPOSTERIOR(theta.new, y=y, X=X, ...) - 0.5*(r.new %*% r.new) -
                        logPOSTERIOR(theta[[m-1]], y=y, X=X, ...) - 0.5*(r0 %*% r0))

    if (log.alpha < log(u)) {
      theta[[m]] <- theta.new
      r[[m]] <- -r.new
      accept <- accept + 1
    } else {
      theta[[m]] <- theta[[m-1]]
      r[[m]] <- r[[m-1]]
    }

  }
  list(theta=theta,
       r=r,
       accept = accept)
}

