#' Fitter function for Metropolis-Hastings (MH)
#'
#' This is the basic computing function for MH and should not be called directly except by experienced users.
#'
#' @param N Number of MCMC samples
#' @param theta.init Vector of initial values for the parameters
#' @param qPROP Function to generate proposal
#' @param qFUN Probability for proposal function.  First argument is where to evaluate, and second argument is the conditional parameter
#' @param logPOSTERIOR Function to calculate and return the log posterior given a vector of values of \code{theta}
#' @param nu Single value or vector parameter passed to \code{qPROP} or \code{qFUN} for the proposal density
#' @param varnames Optional vector of theta parameter names
#' @param param(...) List of additional parameters for \code{logPOSTERIOR}
#' @return List for \code{mh}
#'
#' @export
mh.fit <- function(N, theta.init, qPROP, qFUN, logPOSTERIOR, nu=1e-3,
               varnames=NULL, ...) {
  paramSim <- list()
  paramSim[[1]] <- theta.init
  accept <- 0
  accept_v <- vector()
  accept_v[1] <- 1
  for (j in 2:N) {
    u <- runif(1)
    paramProposal <- qPROP(paramSim[[j-1]], nu)
    lnum <- logPOSTERIOR(paramProposal, ...) + qFUN(paramSim[[j-1]], paramProposal, nu)
    lden <- logPOSTERIOR(paramSim[[j-1]], ...) + qFUN(paramProposal, paramSim[[j-1]], nu)
    l.alpha <- pmin(0, lnum - lden)
    if (l.alpha > log(u)) {
      paramSim[[j]] <- paramProposal
      accept <- accept + 1
      accept_v <- c(accept_v, 1)
    } else {
      paramSim[[j]] <- paramSim[[j-1]]
      accept_v <- c(accept_v, 0)
    }
  }

  # create dataframe from simulation
  thetaCombined <- as.data.frame(do.call(rbind, paramSim))

  if (!is.null(varnames)) {
    colnames(thetaCombined) <- varnames
  }

  # obj <- list(paramSim = paramSim,
  #      thetaCombined = thetaCombined,
  #      accept = accept)

  obj <- list(N=N,
              theta=paramSim,
              thetaCombined = thetaCombined,
              r=NULL,
              theta.all = paramSim,
              r.all = NULL,
              accept=accept,
              accept_v = accept_v,
              M=NULL,
              algorithm="MH")

  # class(obj) <- c("hmclearn", "list")
  return(obj)
}

withGlobals <- function(FUN, lst){
  environment(FUN) <- list2env(lst)
  FUN
}


mhpar <- function(paramlst, ...) {
  do.call(mh.fit, paramlst)
}

#' Fit a generic model using Metropolis-Hastings (MH)
#'
#' This function runs the MH algorithm on a generic model provided
#' the \code{logPOSTERIOR} function.
#' All parameters specified within the list \code{param} are passed to these the posterior function.
#'
#' @param N Number of MCMC samples
#' @param theta.init Vector of initial values for the parameters
#' @param qPROP Function to generate proposal
#' @param qFUN Probability for proposal function.  First argument is where to evaluate, and second argument is the conditional parameter
#' @param logPOSTERIOR Function to calculate and return the log posterior given a vector of values of \code{theta}
#' @param nu Single value or vector parameter passed to \code{qPROP} or \code{qFUN} for the proposal density
#' @param varnames Optional vector of theta parameter names
#' @param ... Additional parameters for \code{logPOSTERIOR}
#' @return Object of class \code{hmclearn}
#'
#' @section Elements for \code{hmclearn} objects:
#' \describe{
#'   \item{\code{N}}{
#'   Number of MCMC samples
#'   }
#'   \item{\code{theta}}{
#'   Nested list of length \code{N} of the sampled values of \code{theta} for each chain
#'   }
#'   \item{\code{thetaCombined}}{
#'   List of dataframes containing sampled values, one for each chain
#'   }
#'   \item{\code{r}}{
#'   NULL for Metropolis-Hastings
#'   }
#'   \item{\code{theta.all}}{
#'   Nested list of all parameter values of \code{theta} sampled prior to accept/reject step for each
#'   }
#'   \item{\code{r.all}}{
#'   NULL for Metropolis-Hastings
#'   }
#'   \item{\code{accept}}{
#'   Number of accepted proposals.  The ratio \code{accept} / \code{N} is the acceptance rate
#'   }
#'   \item{\code{accept_v}}{
#'   Vector of length \code{N} indicating which samples were accepted
#'   }
#'   \item{\code{M}}{
#'   NULL for Metropolis-Hastings
#'   }
#'   \item{\code{algorithm}}{
#'   \code{MH} for Metropolis-Hastings
#'   }
#'   \item{\code{varnames}}{
#'   Optional vector of parameter names
#'   }
#'   \item{\code{chains}}{
#'   Number of MCMC chains
#'   }
#' }
#'
#' @section Available \code{logPOSTERIOR} functions:
#' \describe{
#'   \item{\code{linear_posterior}}{
#'   Linear regression:  log posterior
#'   }
#'   \item{\code{logistic_posterior}}{
#'   Logistic regression:  log posterior
#'   }
#'   \item{\code{poisson_posterior}}{
#'   Poisson (count) regression:  log posterior
#'   }
#'   \item{\code{lmm_posterior}}{
#'   Linear mixed effects model:  log posterior
#'   }
#'   \item{\code{glmm_bin_posterior}}{
#'   Logistic mixed effects model:  log posterior
#'   }
#'   \item{\code{glmm_poisson_posterior}}{
#'   Poisson mixed effects model:  log posterior
#'   }
#' }
#'
#' @author Samuel Thomas \email{samthoma@@iu.edu}, Wanzhu Tu \email{wtu@iu.edu}
#' @export
mh <- function(N, theta.init, qPROP, qFUN, logPOSTERIOR, nu=1e-3,
                   varnames=NULL, param = list(...),
               chains=1, parallel=FALSE, ...) {

  allparam <- c(list(N=N,
                     theta.init=theta.init,
                     qPROP=qPROP,
                     qFUN=qFUN,
                     logPOSTERIOR=logPOSTERIOR,
                     nu=nu,
                     varnames=varnames),
                param)

  if (parallel) {
    no_cores <- pmin(parallel::detectCores(), chains)
    cl <- parallel::makeCluster(no_cores)

    allparamParallel <- replicate(no_cores, allparam, FALSE)
    parallel::clusterExport(cl, varlist=c("mhpar", "mh.fit", "leapfrog"), envir=environment())

    res <- parallel::parLapply(cl=cl, X=allparamParallel, fun="mhpar")
    parallel::stopCluster(cl)

    # store array
    thetaCombined <- lapply(res, function(xx) as.matrix(xx$thetaCombined))


    obj <- list(N=N,
                theta = lapply(res, function(xx) xx$theta),
                # thetaCombined = sapply(thetaCombined, as.matrix, simplify="array"),
                thetaCombined = thetaCombined,
                r = NULL,
                theta.all = lapply(res, function(xx) xx$theta),
                r.all = NULL,
                accept = sapply(res, function(xx) xx$accept),
                accept_v = lapply(res, function(xx) xx$accept_v),
                M=NULL,
                algorithm = "MH",
                varnames = varnames,
                chains = no_cores)
    class(obj) <- c("hmclearn", "list")
    return(obj)

  } else {

    allparamParallel <- replicate(chains, allparam, FALSE)

    res <- lapply(allparamParallel, mhpar)

    # store array
    thetaCombined <- lapply(res, function(xx) as.matrix(xx$thetaCombined))


    obj <- list(N=N,
                theta = lapply(res, function(xx) xx$theta),
                # thetaCombined = sapply(thetaCombined, as.matrix, simplify="array"),
                thetaCombined = thetaCombined,
                r = NULL,
                theta.all = lapply(res, function(xx) xx$theta),
                r.all = NULL,
                accept = sapply(res, function(xx) xx$accept),
                accept_v = lapply(res, function(xx) xx$accept_v),
                M=NULL,
                algorithm = "MH",
                varnames = varnames,
                chains = chains)
    class(obj) <- c("hmclearn", "list")
    return(obj)

    # # res <- do.call(mh.fit, allparam )
    # res <- mhpar(allparam)
    # # res$thetaCombined <- array(res$thetaCombined,
    # #                            dim=c(dim(res$thetaCombined), 1))
    # res$thetaCombined <- list(as.matrix(res$thetaCombined))
    # res$varnames <- varnames
    # res$chains <- 1
    # return(res)
  }
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
  # info <- paste("r", length(r), "epsilon", length(epsilon), "g.ld", length(g.ld))
  # print(info)

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

#' Fitter function for Hamiltonian Monte Carlo (HMC)
#'
#' This is the basic computing function for HMC and should not be called directly except by experienced users.
#'
#' @param N Number of MCMC samples
#' @param theta.init Vector of initial values for the parameters
#' @param epsilon Step-size parameter for \code{leapfrog}
#' @param L Number of \code{leapfrog} steps parameter
#' @param logPOSTERIOR Function to calculate and return the log posterior given a vector of values of \code{theta}
#' @param glogPOSTERIOR Function to calculate and return the gradient of the log posterior given a vector of values of  \code{theta}
#' @param varnames Optional vector of theta parameter names
#' @param randlength Logical to determine whether to apply some randomness to the number of leapfrog steps tuning parameter \code{L}
#' @param Mdiag Optional vector of the diagonal of the mass matrix \code{M}.  Defaults to unit diagonal.
#' @param constrain Optional vector of which parameters in \code{theta} accept positive values only.  Default is that all parameters accept all real numbers
#' @param verbose Logical to determine whether to display the progress of the HMC algorithm
#' @param ... additional parameters for \code{logPOSTERIOR} and \code{glogPOSTERIOR}
#' @return List for \code{hmc}
#'
#' @export
hmc.fit <- function(N, theta.init, epsilon, L, logPOSTERIOR, glogPOSTERIOR, varnames=NULL,
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
  thetaCombined <- as.data.frame(do.call(rbind, theta))

  if (!is.null(varnames)) {
    colnames(thetaCombined) <- varnames
  }

  obj <- list(N=N,
       theta=theta,
       thetaCombined = thetaCombined,
       r=r,
       theta.all = theta.all,
       r.all = r.all,
       accept=accept,
       accept_v = accept_v,
       M=M_mx,
       algorithm="HMC")

  # class(obj) <- c("hmclearn", "list")
  return(obj)
}


hmcpar <- function(paramlst, ...) {
  do.call(hmc.fit, paramlst)
}

#' Fit a generic model using Hamiltonian Monte Carlo (HMC)
#'
#' This function runs the HMC algorithm on a generic model provided
#' the \code{logPOSTERIOR} and gradient \code{glogPOSTERIOR} functions.
#' All parameters specified within the list \code{param}are passed to these two functions.
#' The tuning parameters \code{epsilon} and \code{L} are passed to the
#' Leapfrog algorithm.
#'
#' @param N Number of MCMC samples
#' @param theta.init Vector of initial values for the parameters
#' @param epsilon Step-size parameter for \code{leapfrog}
#' @param L Number of \code{leapfrog} steps parameter
#' @param logPOSTERIOR Function to calculate and return the log posterior given a vector of values of \code{theta}
#' @param glogPOSTERIOR Function to calculate and return the gradient of the log posterior given a vector of values of  \code{theta}
#' @param varnames Optional vector of theta parameter names
#' @param randlength Logical to determine whether to apply some randomness to the number of leapfrog steps tuning parameter \code{L}
#' @param Mdiag Optional vector of the diagonal of the mass matrix \code{M}.  Defaults to unit diagonal.
#' @param constrain Optional vector of which parameters in \code{theta} accept positive values only.  Default is that all parameters accept all real numbers
#' @param verbose Logical to determine whether to display the progress of the HMC algorithm
#' @param param(...) List of additional parameters for \code{logPOSTERIOR} and \code{glogPOSTERIOR}
#' @return Object of class \code{hmclearn}
#'
#' @section Elements for \code{hmclearn} objects:
#' \describe{
#'   \item{\code{N}}{
#'   Number of MCMC samples
#'   }
#'   \item{\code{theta}}{
#'   Nested list of length \code{N} of the sampled values of \code{theta} for each chain
#'   }
#'   \item{\code{thetaCombined}}{
#'   List of dataframes containing sampled values, one for each chain
#'   }
#'   \item{\code{r}}{
#'   List of length \code{N} of the sampled momenta
#'   }
#'   \item{\code{theta.all}}{
#'   Nested list of all parameter values of \code{theta} sampled prior to accept/reject step for each
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
#'   \item{\code{M}}{
#'   Mass matrix used in the HMC algorithm
#'   }
#'   \item{\code{varnames}}{
#'   Optional vector of parameter names
#'   }
#'   \item{\code{chains}}{
#'   Number of MCMC chains
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
#'   \item{\code{lmm_posterior}}{
#'   Linear mixed effects model:  log posterior
#'   }
#'   \item{\code{g_lmm_posterior}}{
#'   Linear mixed effects model:  gradient of the log posterior
#'   }
#'   \item{\code{glmm_bin_posterior}}{
#'   Logistic mixed effects model:  log posterior
#'   }
#'   \item{\code{g_glmm_bin_posterior}}{
#'   Logistic mixed effects model:  gradient of the log posterior
#'   }
#'   \item{\code{glmm_poisson_posterior}}{
#'   Poisson mixed effects model:  log posterior
#'   }
#'   \item{\code{g_glmm_poisson_posterior}}{
#'   Poisson mixed effects model:  gradient of the log posterior
#'   }
#' }
#'
#' @author Samuel Thomas \email{samthoma@@iu.edu}, Wanzhu Tu \email{wtu@iu.edu}
#' @references \emph{HMC in R} paper
#' @references Thomas, S., Li, X., and Tu, W.  2019.  \emph{Hamiltonian Monte Carlo}.  Wiley
#' @references Neal, Radford. 2011. \emph{MCMC Using Hamiltonian Dynamics.} In Handbook of Markov Chain Monte Carlo, edited by Steve Brooks, Andrew Gelman, Galin L. Jones, and Xiao-Li Meng, 116–62. Chapman; Hall/CRC.
#' @keywords hamiltonian monte carlo
#' @export
hmc <- function(N=10000, theta.init, epsilon=1e-2, L=10, logPOSTERIOR, glogPOSTERIOR,
                randlength=FALSE, Mdiag=NULL, constrain=NULL, verbose=FALSE, varnames=NULL,
                param = list(...),
               chains=1, parallel=FALSE, ...) {

  allparam <- c(list(N=N,
                     theta.init=theta.init,
                     epsilon=epsilon,
                     L=L,
                     logPOSTERIOR=logPOSTERIOR,
                     glogPOSTERIOR=glogPOSTERIOR,
                     randlength=randlength,
                     Mdiag=Mdiag,
                     constrain=constrain,
                     verbose=verbose,
                     varnames=varnames),
                param)

  if (parallel) {
    no_cores <- pmin(parallel::detectCores(), chains)
    cl <- parallel::makeCluster(no_cores)

    allparamParallel <- replicate(no_cores, allparam, FALSE)
    parallel::clusterExport(cl, varlist=c("hmcpar", "hmc.fit", "leapfrog"), envir=environment())

    res <- parallel::parLapply(cl=cl, X=allparamParallel, fun="hmcpar")
    parallel::stopCluster(cl)

    thetaCombined <- lapply(res, function(xx) as.matrix(xx$thetaCombined))


    # store array
    obj <- list(N=N,
                theta = lapply(res, function(xx) xx$theta),
                # thetaCombined = sapply(thetaCombined, as.matrix, simplify="array"),
                thetaCombined = thetaCombined,
                r = lapply(res, function(xx) xx$r),
                theta.all = lapply(res, function(xx) xx$theta),
                r.all =lapply(res, function(xx) xx$r.all),
                accept = sapply(res, function(xx) xx$accept),
                accept_v = lapply(res, function(xx) xx$accept_v),
                M=lapply(res, function(xx) xx$M_Mx),
                algorithm = "HMC",
                varnames = varnames,
                chains = no_cores)
    class(obj) <- c("hmclearn", "list")
    return(obj)

  } else {

    allparamParallel <- replicate(chains, allparam, FALSE)

    res <- lapply(allparamParallel, hmcpar)

    # store array
    thetaCombined <- lapply(res, function(xx) as.matrix(xx$thetaCombined))


    obj <- list(N=N,
                theta = lapply(res, function(xx) xx$theta),
                # thetaCombined = sapply(thetaCombined, as.matrix, simplify="array"),
                thetaCombined = thetaCombined,
                r = NULL,
                theta.all = lapply(res, function(xx) xx$theta),
                r.all = NULL,
                accept = sapply(res, function(xx) xx$accept),
                accept_v = lapply(res, function(xx) xx$accept_v),
                M=NULL,
                algorithm = "HMC",
                varnames = varnames,
                chains = chains)
    class(obj) <- c("hmclearn", "list")
    return(obj)

    # # res <- do.call(mh.fit, allparam )
    # res <- hmcpar(allparam)
    # # res$thetaCombined <- array(res$thetaCombined,
    # #                            dim=c(dim(res$thetaCombined), 1))
    # res$thetaCombined <- list(as.matrix(res$thetaCombined))
    # res$varnames <- varnames
    # res$chains <- 1
    # return(res)
  }
}
