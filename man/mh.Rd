% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_functions.R
\name{mh}
\alias{mh}
\title{Fit a generic model using Metropolis-Hastings (MH)}
\usage{
mh(
  N,
  theta.init,
  qPROP,
  qFUN,
  logPOSTERIOR,
  nu = 0.001,
  varnames = NULL,
  param = list(),
  chains = 1,
  parallel = FALSE,
  ...
)
}
\arguments{
\item{N}{Number of MCMC samples}

\item{theta.init}{Vector of initial values for the parameters}

\item{qPROP}{Function to generate proposal}

\item{qFUN}{Probability for proposal function.  First argument is where to evaluate, and second argument is the conditional parameter}

\item{logPOSTERIOR}{Function to calculate and return the log posterior given a vector of values of \code{theta}}

\item{nu}{Single value or vector parameter passed to \code{qPROP} or \code{qFUN} for the proposal density}

\item{varnames}{Optional vector of theta parameter names}

\item{param}{List of additional parameters for \code{logPOSTERIOR} and \code{glogPOSTERIOR}}

\item{chains}{Number of MCMC chains to run}

\item{parallel}{Logical to set whether multiple MCMC chains should be run in parallel}

\item{...}{Additional parameters for \code{logPOSTERIOR}}
}
\value{
Object of class \code{hmclearn}
}
\description{
This function runs the MH algorithm on a generic model provided
the \code{logPOSTERIOR} function.
All parameters specified within the list \code{param} are passed to these the posterior function.
}
\section{Elements for \code{hmclearn} objects}{

\describe{
  \item{\code{N}}{
  Number of MCMC samples
  }
  \item{\code{theta}}{
  Nested list of length \code{N} of the sampled values of \code{theta} for each chain
  }
  \item{\code{thetaCombined}}{
  List of dataframes containing sampled values, one for each chain
  }
  \item{\code{r}}{
  NULL for Metropolis-Hastings
  }
  \item{\code{theta.all}}{
  Nested list of all parameter values of \code{theta} sampled prior to accept/reject step for each
  }
  \item{\code{r.all}}{
  NULL for Metropolis-Hastings
  }
  \item{\code{accept}}{
  Number of accepted proposals.  The ratio \code{accept} / \code{N} is the acceptance rate
  }
  \item{\code{accept_v}}{
  Vector of length \code{N} indicating which samples were accepted
  }
  \item{\code{M}}{
  NULL for Metropolis-Hastings
  }
  \item{\code{algorithm}}{
  \code{MH} for Metropolis-Hastings
  }
  \item{\code{varnames}}{
  Optional vector of parameter names
  }
  \item{\code{chains}}{
  Number of MCMC chains
  }
}
}

\section{Available \code{logPOSTERIOR} functions}{

\describe{
  \item{\code{linear_posterior}}{
  Linear regression:  log posterior
  }
  \item{\code{logistic_posterior}}{
  Logistic regression:  log posterior
  }
  \item{\code{poisson_posterior}}{
  Poisson (count) regression:  log posterior
  }
  \item{\code{lmm_posterior}}{
  Linear mixed effects model:  log posterior
  }
  \item{\code{glmm_bin_posterior}}{
  Logistic mixed effects model:  log posterior
  }
  \item{\code{glmm_poisson_posterior}}{
  Poisson mixed effects model:  log posterior
  }
}
}

\examples{
# Linear regression example
set.seed(521)
X <- cbind(1, matrix(rnorm(300), ncol=3))
betavals <- c(0.5, -1, 2, -3)
y <- X\%*\%betavals + rnorm(100, sd=.2)

f1_mh <- mh(N = 3e3,
         theta.init = c(rep(0, 4), 1),
         nu <- c(rep(0.001, 4), 0.1),
         qPROP = qprop,
         qFUN = qfun,
         logPOSTERIOR = linear_posterior,
         varnames = c(paste0("beta", 0:3), "log_sigma_sq"),
         param=list(y=y, X=X), parallel=FALSE, chains=1)

summary(f1_mh, burnin=1000)


}
\author{
Samuel Thomas \email{samthoma@iu.edu}, Wanzhu Tu \email{wtu@iu.edu}
}
