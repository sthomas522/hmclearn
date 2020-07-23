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
#' @param param List of additional parameters for \code{logPOSTERIOR}
#' @param ... Additional parameters for \code{logPOSTERIOR}
#' @return List for \code{mh}
#'
#' @section Elements in \code{mh} list:
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
#' }
#'
#' @examples
#' # Logistic regression example
#' X <- cbind(1, seq(-100, 100, by=0.25))
#' betavals <- c(-0.9, 0.2)
#' lodds <- X %*% betavals
#' prob1 <- as.numeric(1 / (1 + exp(-lodds)))
#'
#' set.seed(9874)
#' y <- sapply(prob1, function(xx) {
#'   sample(c(0, 1), 1, prob=c(1-xx, xx))
#' })
#'
#' f1 <- mh.fit(N = 2000,
#'          theta.init = rep(0, 2),
#'          nu = c(0.03, 0.001),
#'          qPROP = qprop,
#'          qFUN = qfun,
#'          logPOSTERIOR = logistic_posterior,
#'          varnames = paste0("beta", 0:1),
#'          y=y, X=X)
#'
#' f1$accept / f1$N
#' @export
mh.fit <- function(N, theta.init, qPROP, qFUN, logPOSTERIOR, nu=1e-3,
               varnames=NULL, param=list(), ...) {
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
#' @param param List of additional parameters for \code{logPOSTERIOR} and \code{glogPOSTERIOR}
#' @param chains Number of MCMC chains to run
#' @param parallel Logical to set whether multiple MCMC chains should be run in parallel
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
#' @examples
#' # Linear regression example
#' set.seed(521)
#' X <- cbind(1, matrix(rnorm(300), ncol=3))
#' betavals <- c(0.5, -1, 2, -3)
#' y <- X%*%betavals + rnorm(100, sd=.2)
#'
#' f1_mh <- mh(N = 3e3,
#'          theta.init = c(rep(0, 4), 1),
#'          nu <- c(rep(0.001, 4), 0.1),
#'          qPROP = qprop,
#'          qFUN = qfun,
#'          logPOSTERIOR = linear_posterior,
#'          varnames = c(paste0("beta", 0:3), "log_sigma_sq"),
#'          param=list(y=y, X=X), parallel=FALSE, chains=1)
#'
#' summary(f1_mh, burnin=1000)
#'
#'
#' @author Samuel Thomas \email{samthoma@@iu.edu}, Wanzhu Tu \email{wtu@iu.edu}
#' @export
mh <- function(N, theta.init, qPROP, qFUN, logPOSTERIOR, nu=1e-3,
                   varnames=NULL, param = list(),
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

  }
}


#' Multivariate Normal Density of Theta1 | Theta2
#'
#' Provided for Random Walk Metropolis algorithm
#'
#' @param theta1 Vector of current quantiles
#' @param theta2 Vector for mean parameter
#' @param nu Either a single numeric value for the covariance matrix, or a vector for the diagonal
#'
#' @return Multivariate normal density vector log-transformed
#' @references Alan Genz, Frank Bretz, Tetsuhisa Miwa, Xuefei
#' Mi, Friedrich Leisch, Fabian Scheipl and Torsten Hothorn (2019).  \emph{mvtnorm: Multivariate Normal and t Distributions}
#' @export
#'
#' @examples
#' qfun(0, 0, 1)
#' log(1/sqrt(2*pi))
#'
qfun <- function(theta1, theta2, nu) {
  k <- length(theta1)
  nu <- diag(nu, k, k)
  mvtnorm::dmvnorm(theta1, theta2, nu, log=TRUE)
}

#' Simulate from Multivariate Normal Density for Metropolis Algorithm
#'
#' Provided for Random Walk Metropolis algorithm
#'
#' @param theta1 Vector of current quantiles
#' @param nu Either a single numeric value for the covariance matrix, or a vector for the diagonal
#'
#' @references B. D. Ripley (1987) \emph{Stochastic Simulation}. Wiley.  Page 98
#' @references Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S.} Fourth edition. Springer.
#' @return Returns a single numeric simulated value from a Normal distribution or vector of length \code{theta1}.
#' \code{length(mu)} matrix with one sample in each row.
#'
#' @export
#'
#' @examples
#' s <- replicate(1000, qprop(0, 1))
#' summary(s)
#' hist(s, col='light blue')
qprop <- function(theta1, nu) {
  k <- length(theta1)
  nu <- diag(nu, k, k)
  MASS::mvrnorm(1, theta1, nu)
}

#' Leapfrog Algorithm for Hamiltonian Monte Carlo
#'
#' Runs a single iteration of the leapfrog algorithm.  Typically called directly from \code{hmc}
#'
#' @param theta_lf starting parameter vector
#' @param r starting momentum vector
#' @param epsilon Step-size parameter for \code{leapfrog}
# #' @param logPOSTERIOR Function to calculate and return the log posterior given a vector of values of \code{theta}
#' @param glogPOSTERIOR Function to calculate and return the gradient of the log posterior given a vector of values of \code{theta}
#' @param Minv Inverse Mass matrix
#' @param constrain Optional vector of which parameters in \code{theta} accept positive values only.  Default is that all parameters accept all real numbers
#' @param lastSTEP Boolean indicating whether to calculate the last half-step of the momentum update
#' @param ... Additional parameters passed to glogPOSTERIOR
#' @references Neal, Radford. 2011. \emph{MCMC Using Hamiltonian Dynamics.} In Handbook of Markov Chain Monte Carlo, edited by Steve Brooks, Andrew Gelman, Galin L. Jones, and Xiao-Li Meng, 116–62. Chapman; Hall/CRC.
#' @return List containing two elements:  \code{theta.new} the ending value of theta and \code{r.new} the ending value of the momentum
#' @export
#'
#' @examples
#' set.seed(321)
#' X <- cbind(1, rnorm(10))
#' y <- rnorm(10)
#' p <- runif(3) - 0.5
#' leapfrog(rep(0,3), p, 0.01, g_linear_posterior,
#'          diag(3), FALSE, X=X, y=y)
leapfrog <- function(theta_lf, r, epsilon, glogPOSTERIOR, Minv, constrain,
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
#' @param ... Additional parameters for \code{logPOSTERIOR} and \code{glogPOSTERIOR}
#' @return List for \code{hmc}
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
#'   \item{\code{algorithm}}{
#'   \code{HMC} for Hamiltonian Monte Carlo
#'   }
#' }
#' @references Neal, Radford. 2011. \emph{MCMC Using Hamiltonian Dynamics.} In Handbook of Markov Chain Monte Carlo, edited by Steve Brooks, Andrew Gelman, Galin L. Jones, and Xiao-Li Meng, 116–62. Chapman; Hall/CRC.
#' @references Betancourt, Michael. 2017.  \emph{A Conceptual Introduction to Hamiltonian Monte Carlo}.
#' @references Thomas, S., Tu, W. 2020. \emph{Learning Hamiltonian Monte Carlo in R}.
# #' @references Thomas, S., Tu, W. 2020. \emph{Hamiltonian Monte Carlo}. In Wiley StatsRef: Statistics Reference Online (eds N. Balakrishnan, T. Colton, B. Everitt, W. Piegorsch, F. Ruggeri and J.L. Teugels).
#' @examples
#' # Logistic regression example
#' X <- cbind(1, seq(-100, 100, by=0.25))
#' betavals <- c(-0.9, 0.2)
#' lodds <- X %*% betavals
#' prob1 <- as.numeric(1 / (1 + exp(-lodds)))
#'
#' set.seed(9874)
#' y <- sapply(prob1, function(xx) {
#'   sample(c(0, 1), 1, prob=c(1-xx, xx))
#' })
#'
#' f1 <- hmc.fit(N = 500,
#'           theta.init = rep(0, 2),
#'           epsilon = c(0.1, 0.002),
#'           L = 10,
#'           logPOSTERIOR = logistic_posterior,
#'           glogPOSTERIOR = g_logistic_posterior,
#'           y=y, X=X)
#'
#' f1$accept / f1$N
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
      lf <- leapfrog(theta_lf = theta.new, r = r.new, epsilon = eps_vals[, jj], # logPOSTERIOR = logPOSTERIOR,
                     glogPOSTERIOR = glogPOSTERIOR,
                     Minv=Minv, constrain=constrain, lastSTEP=lstp, ...)

      theta.new <- lf$theta.new
      r.new <- lf$r.new

    }
    if (verbose) print(jj)

    # standard metropolis-hastings update
    u <- runif(1)

    # use log transform for ratio due to low numbers
    num <- logPOSTERIOR(theta.new,  ...) - 0.5*(r.new %*% Minv %*% r.new)
    den <- logPOSTERIOR(theta[[jj-1]], ...) - 0.5*(r0 %*% Minv %*% r0)

    # alpha = min[1, exp(num - den) ]
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
#' @param param List of additional parameters for \code{logPOSTERIOR} and \code{glogPOSTERIOR}
#' @param chains Number of MCMC chains to run
#' @param parallel Logical to set whether multiple MCMC chains should be run in parallel
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
#'   \item{\code{algorithm}}{
#'   \code{HMC} for Hamiltonian Monte Carlo
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
#' @examples
#' # Linear regression example
#' set.seed(521)
#' X <- cbind(1, matrix(rnorm(300), ncol=3))
#' betavals <- c(0.5, -1, 2, -3)
#' y <- X%*%betavals + rnorm(100, sd=.2)
#'
#' fm1_hmc <- hmc(N = 500,
#'           theta.init = c(rep(0, 4), 1),
#'           epsilon = 0.01,
#'           L = 10,
#'           logPOSTERIOR = linear_posterior,
#'           glogPOSTERIOR = g_linear_posterior,
#'           varnames = c(paste0("beta", 0:3), "log_sigma_sq"),
#'           param=list(y=y, X=X), parallel=FALSE, chains=1)
#'
#' summary(fm1_hmc, burnin=100)
#'
#'
#' # poisson regression example
#' set.seed(7363)
#' X <- cbind(1, matrix(rnorm(40), ncol=2))
#' betavals <- c(0.8, -0.5, 1.1)
#' lmu <- X %*% betavals
#' y <- sapply(exp(lmu), FUN = rpois, n=1)
#'
#' fm2_hmc <- hmc(N = 500,
#'           theta.init = rep(0, 3),
#'           epsilon = 0.01,
#'           L = 10,
#'           logPOSTERIOR = poisson_posterior,
#'           glogPOSTERIOR = g_poisson_posterior,
#'           varnames = paste0("beta", 0:2),
#'           param = list(y=y, X=X),
#'           parallel=FALSE, chains=1)
#'
#' summary(fm2_hmc, burnin=100)
#'
#'
#' @author Samuel Thomas \email{samthoma@@iu.edu}, Wanzhu Tu \email{wtu@iu.edu}
#' @references Neal, Radford. 2011. \emph{MCMC Using Hamiltonian Dynamics.} In Handbook of Markov Chain Monte Carlo, edited by Steve Brooks, Andrew Gelman, Galin L. Jones, and Xiao-Li Meng, 116–62. Chapman; Hall/CRC.
#' @references Betancourt, Michael. 2017.  \emph{A Conceptual Introduction to Hamiltonian Monte Carlo}.
#' @references Thomas, S., Tu, W. 2020. \emph{Learning Hamiltonian Monte Carlo in R}.
# #' @references Thomas, S., Tu, W. 2020. \emph{Hamiltonian Monte Carlo}. In Wiley StatsRef: Statistics Reference Online (eds N. Balakrishnan, T. Colton, B. Everitt, W. Piegorsch, F. Ruggeri and J.L. Teugels).
#' @keywords hamiltonian monte carlo
#' @export
hmc <- function(N=10000, theta.init, epsilon=1e-2, L=10, logPOSTERIOR, glogPOSTERIOR,
                randlength=FALSE, Mdiag=NULL, constrain=NULL, verbose=FALSE, varnames=NULL,
                param = list(),
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
                     varnames=varnames,
                     ... = ...),
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

  }
}

