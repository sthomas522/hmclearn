
# proposal function in form expected by mh
pfun_logistic <- function(PARAM, ...) {
  d <- length(PARAM)
  theta <- PARAM
  logistic_posterior(theta=theta, ...)
}

# proposal function in form expected by mh
pfun_poisson <- function(PARAM, ...) {
  d <- length(PARAM)
  theta <- PARAM
  poisson_posterior(theta=theta, ...)
}

# proposal function in form expected by mh
pfun_lmm <- function(PARAM, ...) {
  d <- length(PARAM)
  theta <- PARAM
  lmm_posterior(theta=theta, ...)
}

# proposal function in form expected by mh
pfun_glmm_bin <- function(PARAM, ...) {
  d <- length(PARAM)
  theta <- PARAM
  glmm_bin_posterior(theta=theta, ...)
}

# proposal function in form expected by mh
pfun_glmm_poisson <- function(PARAM, ...) {
  d <- length(PARAM)
  theta <- PARAM
  glmm_poisson_posterior(theta=theta, ...)
}

#' Sample log posterior and gradient functions for select generalized linear models
#' and mixed effect models
#'
#' These functions can be used to fit common generalized linear models and mixed effect models.
#' See the accompanying vignettes for details on the derivations of the log posterior and gradient.
#' In addition, these functions can be used as templates to build custom models to fit using HMC.
#'
#' @name hmclearn-glm-posterior
#'
#' @param theta vector of parameters.  See details below for the order of parameters for each model
#' @param y numeric vector for the dependent variable for all models
#' @param X numeric design matrix of fixed effect parameters for all models
#' @param Z numeric design matrix of random effect parameters for all mixed effects models
#' @param a hyperparameter for the Inverse Gamma shape parameter for \eqn{\sigma_\epsilon} in linear regression models
#' @param b hyperparameter for the Inverse Gamma scale parameter for \eqn{\sigma_\epsilon} in linear regression models
#' @param sig2beta diagonal covariance of prior for linear predictors is multivariate normal with mean 0 for linear regression and linear mixed effect models.
#' @param n number of observations for standard glm models, or number of subjects for all mixed effect models
#' @param nrandom number of random effects covariance parameters for all mixed effects models
#' @param d number of observations per subject for mixed effects models, but an input for linear mixed effect models only.
# #' @param nueps hyperparameter \eqn{\nu} for the half-t prior of the error parameter for linear mixed effects model \eqn{\sigma_\epsilon}
#' @param nuxi hyperparameter \eqn{\nu} for the half-t prior of the random effects diagonal for all mixed effects models \eqn{\xi}
#' @param nugamma hyperparameter \eqn{\nu} for the half-t prior of the log transformed error for linear mixed effects model \eqn{\gamma}
# #' @param Aeps hyperparameter \eqn{A} for the half-t prior of the error parameter for linear mixed effects model \eqn{\sigma_\epsilon}
#' @param Axi hyperparameter \eqn{A} for the half-t prior of the random effects diagonal for all mixed effects models\eqn{\xi}
#' @param Agamma hyperparameter \eqn{A} for the half-t prior of the log transformed error for linear mixed effects model \eqn{\gamma}
#' @section Generalized Linear Models with available posterior and gradient functions:
#' \describe{
#'   \item{`linear_posterior(theta, y, X, a=1e-4, b=1e-4, sig2beta = 1e3)`}{
#'    The log posterior function for linear regression
#'    \deqn{f(y | X, \beta, \sigma) = \frac{1}{(2\pi\sigma^2)^{n/2}}\exp{\left(-\frac{1}{2\sigma^2} (y - X\beta)^T(y-X\beta) \right)}}
#'    with priors \eqn{p(\sigma^2) \sim IG(a, b)} and \eqn{\beta \sim N(0, \sigma_\beta^2 I)}.  The variance term is log transformed \eqn{\gamma = \log\sigma}
#'    The input parameter vector \eqn{theta} is of length \eqn{k}.  The first \eqn{k-1} parameters are for \eqn{\beta}, and the last parameter is \eqn{\gamma}
#'    Note that the Inverse Gamma prior can be problematic for certain applications with low variance, such as hierarchical models.  See Gelman (2006)
#'   }
#'   \item{`g_linear_posterior(theta, y, X, a = 1e-04, b = 1e-04, sig2beta=1e3)`}{
#'    Gradient of the log posterior for a linear regression model with Normal prior for the linear parameters and Inverse Gamma for the error term.
#'    \deqn{f(y | X, \beta, \sigma) = \frac{1}{(2\pi\sigma^2)^{n/2}}\exp{\left(-\frac{1}{2\sigma^2} (y - X\beta)^T(y-X\beta) \right)}}
#'    with priors \eqn{p(\sigma^2) \sim IG(a, b)} and \eqn{\beta \sim N(0, \sigma_\beta^2 I)}.  The variance term is log transformed \eqn{\gamma = \log\sigma}
#'    The input parameter vector \eqn{theta} is of length \eqn{k}.  The first \eqn{k-1} parameters are for \eqn{\beta}, and the last parameter is \eqn{\gamma}
#'    Note that the Inverse Gamma prior can be problematic for certain applications with low variance, such as hierarchical models.  See Gelman (2006)
#'   }
#'   \item{`logistic_posterior(theta, y, X, sig2beta=1e3) `}{
#'    Log posterior for a logistic regression model with Normal prior for the linear parameters.
#'    The likelihood function for logistic regression
#'    \deqn{f(\beta| X, y) = \prod_{i=1}^{n} \left(\frac{1}{1+e^{-X_i\beta}}\right)^{y_i} \left(\frac{e^{-X_i\beta}}{1+e^{-X_i\beta}}\right)^{1-y_i}}
#'    with priors \eqn{\beta \sim N(0, \sigma_\beta^2 I)}.
#'    The input parameter vector \eqn{theta} is of length \eqn{k}, containing parameter values for \eqn{\beta}
#'   }
#'   \item{`g_logistic_posterior(theta, y, X, sig2beta=1e3) `}{
#'    Gradient of the log posterior for a logistic regression model with Normal prior for the linear parameters.
#'    The likelihood function for logistic regression
#'    \deqn{f(\beta| X, y) = \prod_{i=1}^{n} \left(\frac{1}{1+e^{-X_i\beta}}\right)^{y_i} \left(\frac{e^{-X_i\beta}}{1+e^{-X_i\beta}}\right)^{1-y_i}}
#'    with priors \eqn{\beta \sim N(0, \sigma_\beta^2 I)}.
#'    The input parameter vector \eqn{theta} is of length \eqn{k}, containing parameter values for \eqn{\beta}
#'   }
#'   \item{`poisson_posterior(theta, y, X, sig2beta=1e3) `}{
#'    Log posterior for a Poisson regression model with Normal prior for the linear parameters.
#'    The likelihood function for poisson regression
#'    \deqn{f(\beta| y, X) = \prod_{i=1}^n \frac{e^{-e^{X_i\beta}}e^{y_iX_i\beta}}{y_i!}}
#'    with priors \eqn{\beta \sim N(0, \sigma_\beta^2 I)}.
#'    The input parameter vector \eqn{theta} is of length \eqn{k}, containing parameter values for \eqn{\beta}
#'   }
#'   \item{`g_poisson_posterior(theta, y, X, sig2beta=1e3) `}{
#'    Gradient of the log posterior for a Poisson regression model with Normal prior for the linear parameters.
#'    The likelihood function for poisson regression
#'    \deqn{f(\beta| y, X) = \prod_{i=1}^n \frac{e^{-e^{X_i\beta}}e^{y_iX_i\beta}}{y_i!}}
#'    with priors \eqn{\beta \sim N(0, \sigma_\beta^2 I)}.
#'    The input parameter vector \eqn{theta} is of length \eqn{k}, containing parameter values for \eqn{\beta}
#'   }
#'  }
#' @section Generalized Linear Mixed Effect with available posterior and gradient functions:
#' \describe{
#'   \item{`lmm_posterior(theta, y, X, Z, n, d, nrandom = 1, nueps = 1, nuxi = 1, Aeps = 25, Axi = 25, sig2beta = 1e3) `}{
#'    The log posterior function for linear mixed effects regression
#'    \deqn{f(y | \beta, u, \sigma_\epsilon) \propto (\sigma_\epsilon^2)^{-nd/2} e^{-\frac{1}{2\sigma_\epsilon^2}(y - X\beta - Zu)^T (y - X\beta - Zu)}}
#'    with priors \eqn{\beta \sim N(0, \sigma_\beta^2 I)}, \eqn{\sigma_\epsilon \sim half-t(A_\epsilon, nu_\epsilon)}, \eqn{\lambda \sim half-t}.
#'    The vector \eqn{\xi} is the diagonal of the covariance \eqn{G} log transformed hyperprior where \eqn{u \sim N(0, G}, \eqn{\xi = \log\lambda} and \eqn{A_\xi, \nu_\xi} are parameters for the transformed distribution
#'    The standard deviation of the error is log transformed, where \eqn{\gamma = \log \sigma_\epsilon} and \eqn{\sigma_\epsilon \sim half-t}. The parameters for \eqn{\gamma} are \eqn{A_\gamma, \nu_\gamma}
#'    The input parameter vector \eqn{theta} is of length \eqn{k}.  The order of parameters for the vector is \eqn{\beta, \tau, \gamma, \xi}.
#'   }
#'   \item{`g_lmm_posterior(theta, y, X, Z, n, d, nrandom = 1, nueps = 1, nuxi = 1, Aeps = 25, Axi = 25, sig2beta = 1e3)`}{
#'    Gradient of the log posterior for a linear mixed effects regression model
#'    \deqn{f(y | \beta, u, \sigma_\epsilon) \propto (\sigma_\epsilon^2)^{-n/2} e^{-\frac{1}{2\sigma_\epsilon^2}(y - X\beta - Zu)^T (y - X\beta - Zu)}}
#'    with priors \eqn{\beta \sim N(0, \sigma_\beta^2 I)}, \eqn{\sigma_\epsilon \sim half-t(A_\epsilon, nu_\epsilon)}, \eqn{\lambda \sim half-t}.
#'    The vector \eqn{\xi} is the diagonal of the covariance \eqn{G} log transformed hyperprior where \eqn{u \sim N(0, G}, \eqn{\xi = \log\lambda} and \eqn{A_\xi, \nu_\xi} are parameters for the transformed distribution
#'    The standard deviation of the error is log transformed, where \eqn{\gamma = \log \sigma_\epsilon} and \eqn{\sigma_\epsilon \sim half-t}. The parameters for \eqn{\gamma} are \eqn{A_\gamma, \nu_\gamma}
#'    The input parameter vector \eqn{theta} is of length \eqn{k}.  The order of parameters for the vector is \eqn{\beta, \tau, \gamma, \xi}
#'   }
#'   \item{`glmm_bin_posterior(theta, y, X, Z, n, nrandom = 1, nuxi = 1, Axi = 25, sig2beta=1e3)`}{
#'    The log posterior function for logistic mixed effects regression
#'    \deqn{f(y | X, Z, \beta, u) = \prod_{i=1}^n\prod_{j=1}^d \left(\frac{1}{1 + e^{-X_{i}\beta - Z_{ij}u_i}}\right)^{y_{ij}} \left(\frac{e^{-X_i\beta - Z_{ij}u_i}}{1 + e^{-X_{i}\beta - Z_{ij}u_i}}\right)^{1-y_{ij}} }
#'    with priors \eqn{\beta \sim N(0, \sigma_\beta^2 I)}, \eqn{\sigma_\epsilon \sim half-t(A_\epsilon, nu_\epsilon)}, \eqn{\lambda \sim half-t(A_\lambda, nu_\lambda )}.
#'    The vector \eqn{\lambda} is the diagonal of the covariance \eqn{G} hyperprior where \eqn{u \sim N(0, G}, \eqn{\xi = \log\lambda} and \eqn{A_\xi, \nu_\xi} are parameters for the transformed distribution
#'    The input parameter vector \eqn{theta} is of length \eqn{k}.  The order of parameters for the vector is \eqn{\beta, \tau, \xi}
#'   }
#'   \item{`g_glmm_bin_posterior(theta, y, X, Z, n, nrandom = 1, nuxi = 1, Axi = 25, sig2beta = 1e3) `}{
#'    Gradient of the log posterior function for logistic mixed effects regression
#'    \deqn{f(y | X, Z, \beta, u) = \prod_{i=1}^n\prod_{j=1}^m \left(\frac{1}{1 + e^{-X_{i}\beta - Z_{ij}u_i}}\right)^{y_{ij}} \left(\frac{e^{-X_i\beta - Z_{ij}u_i}}{1 + e^{-X_{i}\beta - Z_{ij}u_i}}\right)^{1-y_{ij}} }
#'    with priors \eqn{\beta \sim N(0, \sigma_\beta^2 I)}, \eqn{\sigma_\epsilon \sim half-t(A_\epsilon, nu_\epsilon)}, \eqn{\lambda \sim half-t(A_\lambda, nu_\lambda )}.
#'    The vector \eqn{\lambda} is the diagonal of the covariance \eqn{G} hyperprior where \eqn{u \sim N(0, G}, \eqn{\xi = \log\lambda} and \eqn{A_\xi, \nu_\xi} are parameters for the transformed distribution
#'    The input parameter vector \eqn{theta} is of length \eqn{k}.  The order of parameters for the vector is \eqn{\beta, \tau, \xi}
#'   }
#'   \item{`glmm_poisson_posterior(theta, y, X, Z, n, nrandom = 1, nuxi = 1, Axi = 25, sig2beta = 1e3) `}{
#'    Log posterior for a Poisson mixed effect regression
#'    \deqn{f(y | X, Z, \beta, u) = \prod_{i=1}^n \prod_{j=1}^m \frac{e^{-e^{X_i\beta + Z_{ij}u_{ij}}}e^{y_i(X_i\beta + Z_{ij}u_{ij})}}{y_i!} }
#'    with priors \eqn{\beta \sim N(0, \sigma_\beta^2 I)}, \eqn{\sigma_\epsilon \sim half-t(A_\epsilon, nu_\epsilon)}, \eqn{\lambda \sim half-t(A_\lambda, nu_\lambda )}.
#'    The vector \eqn{\lambda} is the diagonal of the covariance \eqn{G} hyperprior where \eqn{u \sim N(0, G}, \eqn{\xi = \log\lambda} and \eqn{A_\xi, \nu_\xi} are parameters for the transformed distribution
#'    The input parameter vector \eqn{theta} is of length \eqn{k}.  The order of parameters for the vector is \eqn{\beta, \tau, \xi}
#'   }
#'   \item{`g_glmm_poisson_posterior(theta, y, X, Z, n, nrandom = 1, nuxi = 1, Axi = 25, sig2beta = 1e3) `}{
#'    Gradient of the log posterior for a Poisson mixed effect regression
#'    \deqn{f(y | X, Z, \beta, u) = \prod_{i=1}^n \prod_{j=1}^m \frac{e^{-e^{X_i\beta + Z_{ij}u_{ij}}}e^{y_i(X_i\beta + Z_{ij}u_{ij})}}{y_i!} }
#'    with priors \eqn{\beta \sim N(0, \sigma_\beta^2 I)}, \eqn{\sigma_\epsilon \sim half-t(A_\epsilon, nu_\epsilon)}, \eqn{\lambda \sim half-t(A_\lambda, nu_\lambda )}.
#'    The vector \eqn{\lambda} is the diagonal of the covariance \eqn{G} hyperprior where \eqn{u \sim N(0, G}, \eqn{\xi = \log\lambda} and \eqn{A_\xi, \nu_\xi} are parameters for the transformed distribution
#'    The input parameter vector \eqn{theta} is of length \eqn{k}.  The order of parameters for the vector is \eqn{\beta, \tau, \xi}
#'   }
#'  }
#' @return Numeric vector for the log posterior or gradient of the log posterior
#' @references Gelman, A. (2006). \emph{Prior distributions for variance parameters in hierarchical models (comment on article by Browne and Draper)}. Bayesian analysis, 1(3), 515-534.
#' @references Chan, J. C. C., & Jeliazkov, I. (2009). \emph{MCMC estimation of restricted covariance matrices}. Journal of Computational and Graphical Statistics, 18(2), 457-480.
#' @references Betancourt, M., & Girolami, M. (2015). \emph{Hamiltonian Monte Carlo for hierarchical models}. Current trends in Bayesian methodology with applications, 79, 30.
#' @examples
#' # Linear regression example
#' set.seed(521)
#' X <- cbind(1, matrix(rnorm(300), ncol=3))
#' betavals <- c(0.5, -1, 2, -3)
#' y <- X%*%betavals + rnorm(100, sd=.2)
#'
#' f1_hmc <- hmc(N = 500,
#'           theta.init = c(rep(0, 4), 1),
#'           epsilon = 0.01,
#'           L = 10,
#'           logPOSTERIOR = linear_posterior,
#'           glogPOSTERIOR = g_linear_posterior,
#'           varnames = c(paste0("beta", 0:3), "log_sigma_sq"),
#'           param=list(y=y, X=X), parallel=FALSE, chains=1)
#'
#' summary(f1_hmc, burnin=100)
#'
#'
#' # poisson regression example
#' set.seed(7363)
#' X <- cbind(1, matrix(rnorm(40), ncol=2))
#' betavals <- c(0.8, -0.5, 1.1)
#' lmu <- X %*% betavals
#' y <- sapply(exp(lmu), FUN = rpois, n=1)
#'
#' f2_hmc <- hmc(N = 500,
#'           theta.init = rep(0, 3),
#'           epsilon = 0.01,
#'           L = 10,
#'           logPOSTERIOR = poisson_posterior,
#'           glogPOSTERIOR = g_poisson_posterior,
#'           varnames = paste0("beta", 0:2),
#'           param = list(y=y, X=X),
#'           parallel=FALSE, chains=1)
#'
NULL

#' @rdname hmclearn-glm-posterior
#' @export
linear_posterior <- function(theta, y, X, a=1e-4, b=1e-4, sig2beta=1e3) {
  k <- length(theta)
  beta_param <- as.numeric(theta[1:(k-1)])
  gamma_param <- theta[k]

  n <- nrow(X)
  result <- -(n/2+a)*gamma_param - exp(-gamma_param)/2 * t(y - X%*%beta_param) %*%
    (y - X%*%beta_param) - b*exp(-gamma_param) - 1/2* t(beta_param) %*% beta_param / sig2beta
  return(result)
}

#' @rdname hmclearn-glm-posterior
#' @export
g_linear_posterior <- function(theta, y, X, a=1e-4, b=1e-4, sig2beta=1e3) {
  k <- length(theta)
  beta_param <- as.numeric(theta[1:(k-1)])
  gamma_param <- theta[k]
  n <- nrow(X)

  grad_beta <- exp(-gamma_param)  * t(X) %*% (y - X%*%beta_param)  - beta_param / sig2beta
  grad_gamma <- -(n/2 + a) + exp(-gamma_param)/2 * t(y - X%*%beta_param) %*%
    (y - X%*%beta_param) + b*exp(-gamma_param)
  c(as.numeric(grad_beta), as.numeric(grad_gamma))
}

#' @rdname hmclearn-glm-posterior
#' @export
logistic_posterior <- function(theta, y, X, sig2beta=1e3) {
  k <- length(theta)
  beta_param <- as.numeric(theta)
  onev <- rep(1, length(y))

  ll_bin <- t(beta_param) %*% t(X) %*% (y - 1) -
    t(onev) %*% log(1 + exp(-X %*% beta_param))

  result <- ll_bin - 1/2* t(beta_param) %*% beta_param / sig2beta

  return(result)
}

#' @rdname hmclearn-glm-posterior
#' @export
g_logistic_posterior <- function(theta, y, X, sig2beta=1e3) {
  n <- length(y)
  k <- length(theta)
  beta_param <- as.numeric(theta)

  result <- t(X) %*% ( y - 1  + exp(-X %*% beta_param) /
                         (1 + exp(-X %*% beta_param))) -beta_param/sig2beta

  return(result)
}

#' @rdname hmclearn-glm-posterior
#' @export
poisson_posterior <- function(theta, y, X, sig2beta=1e3) {
  k <- length(theta)
  beta_param <- theta
  onev <- rep(1, length(y))

  ll_pois <- t(y) %*% X %*% beta_param - onev %*% exp( X %*% beta_param)

  result <- ll_pois - 1/2* t(beta_param) %*% beta_param / sig2beta

  return(result)
}

#' @rdname hmclearn-glm-posterior
#' @export
g_poisson_posterior <- function(theta, y, X, sig2beta=1e3) {
  n <- length(y)
  k <- length(theta)

  t(X) %*% (y - exp(X %*% theta)) - theta / sig2beta
}

#' @rdname hmclearn-glm-posterior
#' @export
lmm_posterior <- function(theta, y, X, Z, n, d, nrandom=1,
                          nugamma=1, nuxi=1, Agamma=25, Axi=25, sig2beta=1e3) {
  Z <- as.matrix(Z)
  p <- ncol(X)

  # extract parameters from theta vector
  beta_param <- as.numeric(theta[1:p])
  tau_param <- theta[(p+1):(p+n*nrandom)]

  # log of sig2eps
  gamma_param <- theta[p+n*nrandom+1]

  # diagonal of G matrix
  xi_param <- theta[(p+n*nrandom+2):(p+n*nrandom+nrandom+1)]

  # reconstruct G LDLT decomposition
  LDhalf <- diag(exp(xi_param), nrandom, nrandom)
  LDhalf_block <- kronecker(diag(n), LDhalf)

  # u is deterministic function of xi and tau
  u_param <- LDhalf_block %*% tau_param

  log_likelihood <- -n*d*gamma_param - exp(-2*gamma_param)/2 *
    t(y - X%*%beta_param - Z%*%u_param) %*% (y - X%*%beta_param - Z%*%u_param)
  log_beta_prior <- - 1/2* t(beta_param) %*% beta_param / sig2beta
  log_tau_prior <- -1/2 * t(tau_param) %*% tau_param
  log_gamma_prior <- -(nugamma + 1)/2 * log(1 + 1/nugamma*exp(2*gamma_param)/Agamma^2) + gamma_param
  log_xi_prior <- -(nuxi+1)/2 * log(1 + 1/nuxi*exp(2*xi_param)/Axi^2) + xi_param

  result <- log_likelihood + log_beta_prior + log_tau_prior +
    log_gamma_prior + sum(log_xi_prior)
  return(as.numeric(result))
}


#' @rdname hmclearn-glm-posterior
#' @export
g_lmm_posterior <- function(theta, y, X, Z, n, d, nrandom=1,
                            nugamma=1, nuxi=1, Agamma=25, Axi=25, sig2beta=1e3) {
  Z <- as.matrix(Z)
  p <- ncol(X)

  # extract parameters from theta vector
  beta_param <- as.numeric(theta[1:p])
  tau_param <- theta[(p+1):(p+n*nrandom)]

  # log of sig2eps
  gamma_param <- theta[p+n*nrandom+1]

  # diagonal of G matrix
  xi_param <- theta[(p+n*nrandom+2):(p+n*nrandom+nrandom+1)]

  # reconstruct G LDLT decomposition
  Dhalf <- diag(exp(xi_param), nrandom, nrandom)
  Dhalf_block <- kronecker(diag(n), Dhalf)

  L <- diag(nrandom)

  L_block <- kronecker(diag(n), L)

  LDhalf <- L %*% Dhalf
  LDhalf_block <- kronecker(diag(n), LDhalf)

  # u is deterministic function of xi, a, and tau
  u_param <- LDhalf_block %*% tau_param

  # gradient
  g_beta <- exp(-2*gamma_param) * t(X) %*% (y - X%*%beta_param - Z%*%u_param) -
    beta_param / sig2beta
  g_tau <- exp(-2*gamma_param) * Dhalf_block %*% t(L_block) %*% t(Z) %*%
    (y - X%*%beta_param - Z%*%u_param) - tau_param
  g_gamma <- -n*d + exp(-2*gamma_param) * t(y - X%*%beta_param - Z%*%u_param) %*%
    (y - X%*%beta_param - Z%*%u_param) -
    (nugamma + 1) / (1 + nugamma*Agamma^2 * exp(-2*gamma_param)) + 1

  # gradient for xi using matrix algebra
  zero_v <- rep(0, nrandom)
  g_xi <- sapply(seq_along(1:nrandom), function(jj) {
    zv <- zero_v
    zv[jj] <- 1
    bd <- kronecker(diag(n), diag(zv))
    tot <- (y - X%*%beta_param - Z%*%u_param) %*%
      t(Z %*% L_block %*% bd %*% Dhalf_block %*% tau_param)
    exp(-2*gamma_param) * sum(diag(tot))
  })
  g_xi <- g_xi - (nuxi + 1) / (1 + nuxi*Axi^2 * exp(-2*xi_param)) + 1

  g_all <- c(as.numeric(g_beta),
             as.numeric(g_tau),
             as.numeric(g_gamma),
             g_xi)

  return(g_all)
}

#' @rdname hmclearn-glm-posterior
#' @export
glmm_bin_posterior <- function(theta, y, X, Z, n, nrandom=1,
                               nuxi=1, Axi=25, sig2beta=1e3) {
  Z <- as.matrix(Z)
  p <- ncol(X)

  # extract parameters from theta vector
  beta_param <- as.numeric(theta[1:p])
  tau_param <- theta[(p+1):(p+n*nrandom)]

  # diagonal of G matrix
  xi_param <- theta[(p+n*nrandom+1):(p+n*nrandom+nrandom)]

  # reconstruct G LDLT decomposition
  Dhalf <- diag(exp(xi_param), nrandom, nrandom)

  L <- diag(nrandom)

  LDhalf <- L %*% Dhalf
  LDhalf_block <- kronecker(diag(n), LDhalf)

  # u is deterministic function of xi and tau
  u_param <- LDhalf_block %*% tau_param

  XZbetau <- X %*% beta_param + Z %*% u_param

  log_likelihood <- t(y-1) %*% XZbetau - t(rep(1, length(y))) %*% log(1 + exp(-XZbetau))
  log_beta_prior <- - 1/2*t(beta_param) %*% beta_param/sig2beta
  log_tau_prior <- -1/2 * t(tau_param) %*% tau_param
  log_xi_prior <- -(nuxi + 1)/2 * log(1 + 1/nuxi * exp(2*xi_param) / Axi^2)

  result <- log_likelihood + log_beta_prior + log_tau_prior + sum(log_xi_prior)
  return(as.numeric(result))
}


#' @rdname hmclearn-glm-posterior
#' @export
g_glmm_bin_posterior <- function(theta, y, X, Z, n, nrandom=1,
                                 nuxi=1, Axi=25, sig2beta=1e3) {
  Z <- as.matrix(Z)
  p <- ncol(X)

  # extract parameters from theta vector
  beta_param <- as.numeric(theta[1:p])
  tau_param <- theta[(p+1):(p+n*nrandom)]

  # diagonal of G matrix
  xi_param <- theta[(p+n*nrandom+1):(p+n*nrandom+nrandom)]

  # reconstruct G LDLT decomposition
  Dhalf <- diag(exp(xi_param), nrandom, nrandom)

  L <- diag(nrandom)

  LDhalf <- L %*% Dhalf
  LDhalf_block <- kronecker(diag(n), LDhalf)

  # u is deterministic function of xi and tau
  u_param <- LDhalf_block %*% tau_param

  XZbetau <- X %*% beta_param + Z %*% u_param

  # block L and xi for xi gradient
  L_block <- kronecker(diag(n), L)
  Dhalf_block <- kronecker(diag(n), Dhalf)

  # gradient
  gradprop <- 1 / (1 + exp(X %*% beta_param + Z %*% u_param))

  # g_beta <- -t(1 - y) %*% X + t(gradprop) %*% X - t(beta_param/sig2beta)
  g_beta <- t(X) %*% (y - 1 + gradprop) - beta_param/sig2beta

  # tau gradient
  # g_tau <- (-t(1-y) + t(gradprop)) %*% Z %*% LDhalf_block - tau_param
  g_tau <- LDhalf_block %*% t(Z) %*% (y - 1 + gradprop) - tau_param

  # gradient for xi using matrix algebra
  zero_v <- rep(0, nrandom)
  g_xi <- sapply(seq_along(1:nrandom), function(jj) {
    zv <- zero_v
    zv[jj] <- 1
    bd <- kronecker(diag(n), diag(zv, nrandom, nrandom))
    # (-t(1-y) + t(gradprop)) %*% Z %*% L_block %*% bd %*% Dhalf_block %*% tau_param
    t(tau_param) %*% Dhalf_block %*% bd %*% L_block %*% t(Z) %*% (y - 1 + gradprop)
  })
  g_xi <- g_xi - (nuxi + 1) / (1 + nuxi*Axi^2 * exp(-2*xi_param))

  g_all <- c(as.numeric(g_beta),
             as.numeric(g_tau),
             g_xi)

}

#' @rdname hmclearn-glm-posterior
#' @export
glmm_poisson_posterior <- function(theta, y, X, Z, n, nrandom=1,
                                   nuxi=1, Axi=25, sig2beta=1e3) {
  Z <- as.matrix(Z)
  p <- ncol(X)

  # extract parameters from theta vector
  beta_param <- theta[1:p]
  tau_param <- theta[(p+1):(p+n*nrandom)]

  # diagonal of G matrix
  xi_param <- theta[(p+n*nrandom+1):(p+n*nrandom+nrandom)]

  # reconstruct G LDLT decomposition
  Dhalf <- diag(exp(xi_param), nrandom, nrandom)

  L <- diag(nrandom)

  LDhalf <- L %*% Dhalf
  LDhalf_block <- kronecker(diag(n), LDhalf)

  # u is deterministic function of xi and tau
  u_param <- LDhalf_block %*% tau_param

  XZbetau <- X %*% beta_param + Z %*% u_param

  # log_likelihood <- -sum(exp(XZbetau)) + y %*% XZbetau
  onev <- rep(1, length(y))
  log_likelihood <- -t(onev) %*% exp(XZbetau) + y %*% XZbetau
  log_beta_prior <- - 1/2*t(beta_param)%*% beta_param/sig2beta
  log_tau_prior <- -1/2 * t(tau_param) %*% tau_param
  log_xi_prior <- -(nuxi + 1)/2 * log(1 + 1/nuxi * exp(2*xi_param) / Axi^2)

  result <- log_likelihood + log_beta_prior + log_tau_prior + sum(log_xi_prior)
  return(as.numeric(result))
}

#' @rdname hmclearn-glm-posterior
#' @export
g_glmm_poisson_posterior <- function(theta, y, X, Z, n, nrandom=1,
                                     nuxi=1, Axi=25, sig2beta=1e3) {
  Z <- as.matrix(Z)
  p <- ncol(X)

  # extract parameters from theta vector
  beta_param <- theta[1:p]
  tau_param <- theta[(p+1):(p+n*nrandom)]

  # diagonal of G matrix
  xi_param <- theta[(p+n*nrandom+1):(p+n*nrandom+nrandom)]

  # reconstruct G LDLT decomposition
  Dhalf <- diag(exp(xi_param), nrandom, nrandom)

  L <- diag(nrandom)

  LDhalf <- L %*% Dhalf
  LDhalf_block <- kronecker(diag(n), LDhalf)

  # u is deterministic function of xi and tau
  u_param <- LDhalf_block %*% tau_param

  XZbetau <- X %*% beta_param + Z %*% u_param

  # block L and xi for xi gradient
  L_block <- kronecker(diag(n), L)
  Dhalf_block <- kronecker(diag(n), Dhalf)

  # gradient
  g_beta <- t(X) %*% (-exp(XZbetau) + y)- (beta_param)/sig2beta

  # tau gradient
  g_tau <- t(LDhalf_block) %*% t(Z) %*% (-exp(XZbetau) + y) - tau_param

  # gradient for xi using matrix algebra
  zero_v <- rep(0, nrandom)
  g_xi <- sapply(seq_along(1:nrandom), function(jj) {
    zv <- zero_v
    zv[jj] <- 1
    bd <- kronecker(diag(n), diag(zv, nrandom, nrandom))
      t(L_block %*% bd %*% Dhalf_block %*% tau_param) %*% t(Z) %*% (-exp(XZbetau) + y)
  })
  g_xi <- g_xi - (nuxi + 1) / (1 + nuxi*Axi^2 * exp(-2*xi_param)) + 1

  g_all <- c(as.numeric(g_beta),
             as.numeric(g_tau),
             g_xi)

  return(g_all)
}


