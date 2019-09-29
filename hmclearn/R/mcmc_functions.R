
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
leapfrog <- function(theta_lf, r, epsilon, logPOSTERIOR, glogPOSTERIOR, Minv, constrain,
                     lastSTEP=FALSE, ...) {

  # gradient of log posterior for old theta
  g.ld <- glogPOSTERIOR(theta_lf,  ...)

  # first momentum update
  r.new <- r + epsilon/2*g.ld

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
hmc <- function(N, theta.init, epsilon, L, logPOSTERIOR, glogPOSTERIOR,
                randlength=FALSE, Mdiag=NULL, constrain=FALSE, verbose=FALSE, ...) {

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
  }
  Minv <- qr.solve(M_mx)
  # print(diag(Minv))

  # store theta and momentum (usually not of interest)
  theta <- list()
  theta[[1]] <- theta.init
  r <- list()
  r[[1]] <- NA
  accept <- 0
  for (jj in 2:N) {
    theta[[jj]] <- theta.new <- theta[[jj-1]]
    r0 <- MASS::mvrnorm(1, mu.p, M_mx)
    r.new <- r[[jj]] <- r0

    for (i in 1:L_vals[jj]) {
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
    } else {
      theta[[jj]] <- theta[[jj-1]]
      r[[jj]] <- r[[jj-1]]
    }

  }
  list(theta=theta,
       r=r,
       accept=accept,
       M=M_mx)
}

