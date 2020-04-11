
# proposal function in form expected by mh
pfun_logistic <- function(PARAM, ...) {
  d <- length(PARAM)
  theta <- PARAM
  logistic_posterior(theta=theta, ...)
}

log_lik_bin <- function(beta_param, y, X) {
  result <- X %*% beta_param * (y - 1) - log(1 + exp(-X %*% beta_param))
  sum(result)
}

# proposal function in form expected by mh
pfun_poisson <- function(PARAM, ...) {
  d <- length(PARAM)
  theta <- PARAM
  poisson_posterior(theta=theta, ...)
}

# loglikelihood poisson regression
# note that factorial omitted
log_lik_poisson <- function(beta_param, y, X) {
  y %*% X %*% beta_param - sum(exp( X %*% beta_param)) # - sum(lfactorial(y))
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
#' @name hmclearn-glm-posterior
#'
#' @param theta vector of parameters.  Stored as a single vector in order fixed effect, random effect, log-transformed diagonal \eqn{\lambda}, and off-diagonal of \code{G} vector \code{a}
#' @param y numeric vector for the dependent variable
#' @param X numeric design matrix of fixed effect parameters
#' @param Z numeric design matrix of random effect parameters
#' @param a hyperprior for the Inverse Gamma shape parameter
#' @param b hyperprior for the Inverse Gamma scale parameter
#' @param B prior for linear predictors is multivariate Normal with mean 0 with diagonal covariance B^-1
#' @param m number of random effect linear parameters
#' @param q number of random effects covariance parameters
#' @param A hyperprior numeric vector for the random effects off-diagonal \code{a}
#' @param nueps hyperprior for the half-t prior of the error parameter \eqn{\nu}
#' @param nulambda hyperprior for the half-t prior of the random effects diagonal \eqn{\lambda}
#' @param Aeps hyperprior for the half-t prior of the error parameter \eqn{A_\epsilon}
#' @param Alambda hyperprior for the half-t prior of the random effects diagonal \eqn{A_\lambda}
#' @section Generalized Linear Models with available posterior and gradient functions:
#' \describe{
#'   \item{`linear_posterior(theta, y, X, a=1e-4, b=1e-4, B=0.001)`}{
#'    The log posterior function for linear regression
#'    \deqn{p(y | X, \beta; \sigma^2) = \frac{1}{(2\pi\sigma^2)^{n/2}}\exp{\left(-\frac{1}{2\sigma^2} (y - X\beta)^T(y-X\beta) \right)}}
#'    with priors \eqn{p(\sigma^2) \sim IG(a, b)} and \eqn{\beta \sim N(0, BI)}.  The variance term is log transformed \eqn{\gamma = \log\sigma^2}
#'    The input parameter vector \code{theta} is of length \code{k}.  The first \code{k-1} parameters are for \eqn{\beta}, and the last parameter is \eqn{\gamma}
#'    Note that the Inverse Gamma prior can be problematic for certain applications with low variance.  See Gelman (2006)
#'   }
#'   \item{`g_linear_posterior(theta, y, X, a = 1e-04, b = 1e-04, B = 0.001)`}{
#'    Gradient of the log posterior for a linear regression model with Normal prior for the linear parameters and Inverse Gamma for the error term.
#'    \deqn{p(y | X, \beta; \sigma^2) = \frac{1}{(2\pi\sigma^2)^{n/2}}\exp{\left(-\frac{1}{2\sigma^2} (y - X\beta)^T(y-X\beta) \right)}}
#'    with priors \eqn{p(\sigma^2) \sim IG(a, b)} and \eqn{\beta \sim N(0, BI)}.  The variance term is log transformed \eqn{\gamma = \log\sigma^2}
#'    The input parameter vector \code{theta} is of length \code{k}.  The first \code{k-1} parameters are for \eqn{\beta}, and the last parameter is \eqn{\gamma}
#'    Note that the Inverse Gamma prior can be problematic for certain applications with low variance.  See Gelman (2006)
#'   }
#'   \item{`logistic_posterior(theta, y, X, B = 0.01) `}{
#'    Log posterior for a logistic regression model with Normal prior for the linear parameters.
#'    The likelihood function for logistic regression
#'    \deqn{L(\beta; X, y) = \prod_{i=1}^{n} \left(\frac{1}{1+e^{-X_i\beta}}\right)^{y_i} \left(\frac{e^{-X_i\beta}}{1+e^{-X_i\beta}}\right)^{1-y_i}}
#'    with priors \eqn{\beta \sim N(0, BI)}.
#'    The input parameter vector \code{theta} is of length \code{k}, containing parameter values for \eqn{\beta}
#'   }
#'   \item{`g_logistic_posterior(theta, y, X, B = 0.01) `}{
#'    Gradient of the log posterior for a logistic regression model with Normal prior for the linear parameters.
#'    The likelihood function for logistic regression
#'    \deqn{L(\beta; X, y) = \prod_{i=1}^{n} \left(\frac{1}{1+e^{-X_i\beta}}\right)^{y_i} \left(\frac{e^{-X_i\beta}}{1+e^{-X_i\beta}}\right)^{1-y_i}}
#'    with priors \eqn{\beta \sim N(0, BI)}.
#'    The input parameter vector \code{theta} is of length \code{k}, containing parameter values for \eqn{\beta}
#'   }
#'   \item{`glmm_poisson_posterior(theta, y, X, Z, m, q = 1, A = 10000, nulambda = 1, Alambda = 25, B = 10000) `}{
#'    Log posterior for a Poisson mixed effect regression
#'    \deqn{L(\beta; y, X) = \prod_{i=1}^n \prod_{j=1}^m \frac{e^{-e^{X_i\beta + Z_{ij}u_{ij}}}e^{y_i(X_i\beta + Z_{ij}u_{ij})}}{y_i!} }
#'    with priors \eqn{\beta \sim N(0, BI)}, \eqn{\sigma_\epsilon \sim half-t(A_\epsilon, nu_\epsilon)}, \eqn{\lambda \sim half-t(A_\lambda, nu_\lambda )}.
#'    The vector \eqn{\lambda} is the diagonal of the covariance \code{G} hyperprior where \eqn{u \sim N(0, G}.  The off-diagonal hyperpriors are stored in a vector \eqn{a \sim N(0, A}.  See Chan, Jeliazkov (2009) for details.
#'    The input parameter vector \code{theta} is of length \code{k}.  The first \code{k-1} parameters are for \eqn{\beta}, and the last parameter is \eqn{\gamma}
#'   }
#'   \item{`g_glmm_poisson_posterior(theta, y, X, Z, m, q = 1, A = 10000, nulambda = 1, Alambda = 25, B = 10000) `}{
#'    Gradient of the log posterior for a Poisson mixed effect regression
#'    \deqn{L(\beta; y, X) = \prod_{i=1}^n \prod_{j=1}^m \frac{e^{-e^{X_i\beta + Z_{ij}u_{ij}}}e^{y_i(X_i\beta + Z_{ij}u_{ij})}}{y_i!} }
#'    with priors \eqn{\beta \sim N(0, BI)}, \eqn{\sigma_\epsilon \sim half-t(A_\epsilon, nu_\epsilon)}, \eqn{\lambda \sim half-t(A_\lambda, nu_\lambda )}.
#'    The vector \eqn{\lambda} is the diagonal of the covariance \code{G} hyperprior where \eqn{u \sim N(0, G}.  The off-diagonal hyperpriors are stored in a vector \eqn{a \sim N(0, A}.  See Chan, Jeliazkov (2009) for details.
#'    The input parameter vector \code{theta} is of length \code{k}.  The first \code{k-1} parameters are for \eqn{\beta}, and the last parameter is \eqn{\gamma}
#'   }
#'  }
#' @section Generalized Linear Mixed Effect with available posterior and gradient functions:
#' \describe{
#'   \item{`lmm_posterior(theta, y, X, Z, m, q = 1, A = 10000, nueps = 1, nulambda = 1, Aeps = 25, Alambda = 25, B = 0.001) `}{
#'    The log posterior function for linear mixed effects regression
#'    \deqn{p(y | \beta, u, \sigma_\epsilon^2) \propto (\sigma_\epsilon^2)^{-n/2} e^{-\frac{1}{2\sigma_\epsilon^2}(y - X\beta - Zu)^T (y - X\beta - Zu)}}
#'    with priors \eqn{\beta \sim N(0, BI)}, \eqn{\sigma_\epsilon \sim half-t(A_\epsilon, nu_\epsilon)}, \eqn{\lambda \sim half-t(A_\lambda, nu_\lambda )}.
#'    The vector \eqn{\lambda} is the diagonal of the covariance \code{G} hyperprior where \eqn{u \sim N(0, G}.  The off-diagonal hyperpriors are stored in a vector \eqn{a \sim N(0, A}.  See Chan, Jeliazkov (2009) for details.
#'    The input parameter vector \code{theta} is of length \code{k}.  The first \code{k-1} parameters are for \eqn{\beta}, and the last parameter is \eqn{\gamma}
#'   }
#'   \item{`g_lmm_posterior(theta, y, X, Z, m, q = 1, A = 10000, nueps = 1, nulambda = 1, Aeps = 25, Alambda = 25, B = 0.001)`}{
#'    Gradient of the log posterior for a linear mixed effects regression model
#'    \deqn{p(y | \beta, u, \sigma_\epsilon^2) \propto (\sigma_\epsilon^2)^{-n/2} e^{-\frac{1}{2\sigma_\epsilon^2}(y - X\beta - Zu)^T (y - X\beta - Zu)}}
#'    with priors \eqn{\beta \sim N(0, BI)}, \eqn{\sigma_\epsilon \sim half-t(A_\epsilon, nu_\epsilon)}, \eqn{\lambda \sim half-t(A_\lambda, nu_\lambda )}.
#'    The vector \eqn{\lambda} is the diagonal of the covariance \code{G} hyperprior where \eqn{u \sim N(0, G}.  The off-diagonal hyperpriors are stored in a vector \eqn{a \sim N(0, A}.  See Chan, Jeliazkov (2009) for details.
#'    The input parameter vector \code{theta} is of length \code{k}.  The first \code{k-1} parameters are for \eqn{\beta}, and the last parameter is \eqn{\gamma}
#'   }
#'   \item{`glmm_bin_posterior(theta, y, X, Z, m, q = 1, A = 10000, nulambda = 1, Alambda = 25, B = 10000)`}{
#'    The log posterior function for logistic mixed effects regression
#'    \deqn{p(y | X, Z, \beta, u) = \prod_{i=1}^n\prod_{j=1}^m \left(\frac{1}{1 + e^{-X_{i}\beta - Z_{ij}u_i}}\right)^{y_{ij}} \left(\frac{e^{-X_i\beta - Z_{ij}u_i}}{1 + e^{-X_{i}\beta - Z_{ij}u_i}}\right)^{1-y_{ij}} }
#'    with priors \eqn{\beta \sim N(0, BI)}, \eqn{\sigma_\epsilon \sim half-t(A_\epsilon, nu_\epsilon)}, \eqn{\lambda \sim half-t(A_\lambda, nu_\lambda )}.
#'    The vector \eqn{\lambda} is the diagonal of the covariance \code{G} hyperprior where \eqn{u \sim N(0, G}.  The off-diagonal hyperpriors are stored in a vector \eqn{a \sim N(0, A}.  See Chan, Jeliazkov (2009) for details.
#'    The input parameter vector \code{theta} is of length \code{k}.  The first \code{k-1} parameters are for \eqn{\beta}, and the last parameter is \eqn{\gamma}
#'   }
#'   \item{`g_glmm_bin_posterior(theta, y, X, Z, m, q = 1, A = 10000, nulambda = 1, Alambda = 25, B = 10000) `}{
#'    Gradient of the log posterior function for logistic mixed effects regression
#'    \deqn{p(y | X, Z, \beta, u) = \prod_{i=1}^n\prod_{j=1}^m \left(\frac{1}{1 + e^{-X_{i}\beta - Z_{ij}u_i}}\right)^{y_{ij}} \left(\frac{e^{-X_i\beta - Z_{ij}u_i}}{1 + e^{-X_{i}\beta - Z_{ij}u_i}}\right)^{1-y_{ij}} }
#'    with priors \eqn{\beta \sim N(0, BI)}, \eqn{\sigma_\epsilon \sim half-t(A_\epsilon, nu_\epsilon)}, \eqn{\lambda \sim half-t(A_\lambda, nu_\lambda )}.
#'    The vector \eqn{\lambda} is the diagonal of the covariance \code{G} hyperprior where \eqn{u \sim N(0, G}.  The off-diagonal hyperpriors are stored in a vector \eqn{a \sim N(0, A}.  See Chan, Jeliazkov (2009) for details.
#'    The input parameter vector \code{theta} is of length \code{k}.  The first \code{k-1} parameters are for \eqn{\beta}, and the last parameter is \eqn{\gamma}
#'   }
#'   \item{`poisson_posterior(theta, y, X, B = 0.01) `}{
#'    Log posterior for a Poisson regression model with Normal prior for the linear parameters.
#'    The likelihood function for poisson regression
#'    \deqn{L(\beta; y, X) = \prod_{i=1}^n \frac{e^{-e^{X_i\beta}}e^{y_iX_i\beta}}{y_i!}}
#'    with priors \eqn{\beta \sim N(0, BI)}.
#'    The input parameter vector \code{theta} is of length \code{k}, containing parameter values for \eqn{\beta}
#'   }
#'   \item{`g_poisson_posterior(theta, y, X, B = 0.01) `}{
#'    Gradient of the log posterior for a Poisson regression model with Normal prior for the linear parameters.
#'    The likelihood function for poisson regression
#'    \deqn{L(\beta; y, X) = \prod_{i=1}^n \frac{e^{-e^{X_i\beta}}e^{y_iX_i\beta}}{y_i!}}
#'    with priors \eqn{\beta \sim N(0, BI)}.
#'    The input parameter vector \code{theta} is of length \code{k}, containing parameter values for \eqn{\beta}
#'   }
#'  }
#' @return numeric value for the log posterior or gradient of the log posterior
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
linear_posterior <- function(theta, y, X, a=1e-4, b=1e-4, B=0.001) {
  k <- length(theta)
  beta_param <- as.numeric(theta[1:(k-1)])
  gamma_param <- theta[k]

  nu <- diag(1/B, k-1, k-1)
  inv.nu <- diag(B, k-1, k-1)

  n <- nrow(X)
  result <- -(n/2 + a)*gamma_param - exp(-gamma_param)/2 * t(y - X%*%beta_param) %*%
    (y - X%*%beta_param) - b*exp(-gamma_param) - 1/2* t(beta_param) %*% inv.nu %*% beta_param
  return(result)
}

#' @rdname hmclearn-glm-posterior
#' @export
g_linear_posterior <- function(theta, y, X, a=1e-4, b=1e-4, B=0.001) {
  k <- length(theta)
  beta_param <- as.numeric(theta[1:(k-1)])
  gamma_param <- theta[k]
  n <- nrow(X)

  nu <- diag(1/B, k-1, k-1)
  inv.nu <- diag(B, k-1, k-1)

  grad_beta <- exp(-gamma_param)  * t(X) %*% (y - X%*%beta_param)  - (inv.nu %*% beta_param)
  grad_gamma <- -(n/2 + a) + exp(-gamma_param)/2 * t(y - X%*%beta_param) %*%
    (y - X%*%beta_param) + b*exp(-gamma_param)
  c(as.numeric(grad_beta), as.numeric(grad_gamma))
}

#' @rdname hmclearn-glm-posterior
#' @export
logistic_posterior <- function(theta, y, X, B=0.01) {
  k <- length(theta)
  beta_param <- as.numeric(theta)

  nu <- diag(1/B, k, k)
  inv.nu <- diag(B, k, k)

  result <- log_lik_bin(beta_param=beta_param, y=y, X=X) +
    -1/2*log(det(2*pi*nu)) - 1/2* t(beta_param) %*% inv.nu %*% beta_param

  return(result)
}

#' @rdname hmclearn-glm-posterior
#' @export
g_logistic_posterior <- function(theta, y, X, B=0.01) {
  n <- length(y)
  k <- length(theta)
  beta <- as.numeric(theta)

  nu <- diag(1/B, k, k)
  inv.nu <- diag(B, k, k)

  result <- t(y-1) %*% X %*% diag(1, k, k) +
    t(exp(-X %*% theta) / (1 + exp(-X %*% theta))) %*% X - t(inv.nu %*% theta)
  return(result)
}

#' @rdname hmclearn-glm-posterior
#' @export
poisson_posterior <- function(theta, y, X, B=0.01) {
  k <- length(theta)
  beta_param <- theta

  nu <- diag(1/B, k, k)
  inv.nu <- diag(B, k, k)

  result <- log_lik_poisson(beta_param=beta_param, y=y, X=X) -
    - 1/2* t(beta_param) %*% inv.nu %*% beta_param

  return(result)
}

#' @rdname hmclearn-glm-posterior
#' @export
g_poisson_posterior <- function(theta, y, X, B=0.01) {
  n <- length(y)
  k <- length(theta)

  nu <- diag(1/B, k, k)
  inv.nu <- diag(B, k, k)

  y %*% X - crossprod(exp(X %*% theta), X) - crossprod(theta, inv.nu)
}

#' @rdname hmclearn-glm-posterior
#' @export
lmm_posterior <- function(theta, y, X, Z, m, q=1, A = 1e4, nueps=1, nulambda=1, Aeps=25, Alambda=25, B=0.001) {
  Z <- as.matrix(Z)
  p <- ncol(X)
  n <- nrow(X)

  # extract parameters from theta vector
  beta_param <- as.numeric(theta[1:p])
  tau_param <- theta[(p+1):(p+m*q)]

  # log of sig2eps
  gamma_param <- theta[p+m*q+1]

  # diagonal of G matrix
  xi_param <- theta[(p+m*q+2):(p+m*q+q+1)]

  # off diagonal of G matrix
  a_param <- tail(theta, q*(q-1)/2)
  # a_param <- 0

  nu <- diag(1/B, p, p)
  inv.nu <- diag(B, p, p)

  # reconstruct L D^1/2
  Dhalf <- diag(exp(xi_param))

  L <- diag(q)
  L[lower.tri(L, diag=FALSE)] <- a_param

  LDhalf <- L %*% Dhalf

  Ainv <- diag(1/A, length(a_param), length(a_param))

  LDhalf_block <- kronecker(diag(m), LDhalf)

  # u is deterministic function of xi, a, and tau
  u_param <- LDhalf_block %*% tau_param

  log_likelihood <- -n*gamma_param - exp(-2*gamma_param)/2 *
    t(y - X%*%beta_param - Z%*%u_param) %*% (y - X%*%beta_param - Z%*%u_param)
  log_beta_prior <- - 1/2* t(beta_param) %*% inv.nu %*% beta_param
  log_tau_prior <- -1/2 * t(tau_param) %*% tau_param
  log_gamma_prior <- -(nueps + 1)/2 * log(1 + 1/nueps*exp(2*gamma_param)/Aeps^2) + gamma_param
  log_xi_prior <- -(nulambda+1)/2 * log(1 + 1/nulambda*exp(2*xi_param)/Alambda^2) + xi_param
  log_a_prior <- -1/2*(t(a_param) %*% Ainv %*% a_param)

  result <- log_likelihood + log_beta_prior + log_tau_prior + log_gamma_prior + sum(log_xi_prior) + log_a_prior
  return(as.numeric(result))
}

#' @rdname hmclearn-glm-posterior
#' @export
g_lmm_posterior <- function(theta, y, X, Z, m, q=1, A = 1e4, nueps=1, nulambda=1, Aeps=25, Alambda=25, B=0.001) {
  Z <- as.matrix(Z)
  p <- ncol(X)
  n <- nrow(X)

  # extract parameters from theta vector
  beta_param <- as.numeric(theta[1:p])
  tau_param <- theta[(p+1):(p+m*q)]

  # log of sig2eps
  gamma_param <- theta[p+m*q+1]

  # diagonal of G matrix
  xi_param <- theta[(p+m*q+2):(p+m*q+q+1)]

  # off diagonal of G matrix
  a_param <- tail(theta, q*(q-1)/2)
  # a_param <- 0

  nu <- diag(1/B, p, p)
  inv.nu <- diag(B, p, p)

  # reconstruct L D^1/2
  Dhalf <- diag(exp(xi_param))
  Dhalf_block <- kronecker(diag(m), Dhalf)

  L <- diag(q)
  L[lower.tri(L, diag=FALSE)] <- a_param
  L_block <- kronecker(diag(m), L)

  LDhalf <- L %*% Dhalf

  Ainv <- diag(1/A, length(a_param), length(a_param))

  LDhalf_block <- kronecker(diag(m), LDhalf)

  # u is deterministic function of xi, a, and tau
  u_param <- LDhalf_block %*% tau_param

  # gradient
  g_beta <- exp(-2*gamma_param) * t(X) %*% (y - X%*%beta_param - Z%*%u_param) - (inv.nu %*% beta_param)
  g_tau <- exp(-2*gamma_param) * Dhalf_block %*% t(L_block) %*% t(Z) %*%
    (y - X%*%beta_param - Z%*%u_param) - tau_param
  g_gamma <- -(n-1) + exp(-2*gamma_param) * t(y - X%*%beta_param - Z%*%u_param) %*%
    (y - X%*%beta_param - Z%*%u_param) -
    (nueps + 1) / (1 + nueps*Aeps^2 * exp(-2*gamma_param))

  # gradient for xi using matrix algebra
  zero_v <- rep(0, q)
  g_xi <- sapply(seq_along(1:q), function(jj) {
    zv <- zero_v
    zv[jj] <- 1
    bd <- kronecker(diag(m), diag(zv))
    tot <- (y - X%*%beta_param - Z%*%u_param) %*% t(Z %*% L_block %*% bd %*% Dhalf_block %*% tau_param)
    exp(-2*gamma_param) * sum(diag(tot))
  })
  g_xi <- g_xi - (nulambda + 1) / (1 + nulambda*Alambda^2 * exp(-2*xi_param)) + 1

  # gradient for a
  tau_tilde <- exp(xi_param) * a_param
  T_tilde <- create_Uj(tau_tilde, neg=FALSE)

  # gradient for a
  tau_lst <- split(tau_param, ceiling(seq_along(tau_param)/q))
  Tj <- lapply(tau_lst, function(tauvals) {
    Tj <- create_Uj(Dhalf %*% tauvals, neg=FALSE)
  })
  # Tja <- as.matrix(bdiag(Tja))
  Tj <- do.call(rbind, Tj)

  g_a1 <- exp(-2*gamma_param) * t(Tj) %*% t(Z) %*%
    (y - X%*%beta_param - Z%*%Dhalf_block%*%tau_param - Z%*%Tj%*%a_param)
  g_a2 <- -Ainv %*% a_param
  g_a <- g_a1 + g_a2

  # g_a1 <- sapply(tau_lst, function(tauvals) {
  #   Tj <- create_Uj(tauvals)
  #   exp(-2*gamma_param) * t(Tj) %*% Dinv %*% (uvals - Uj * a_param)
  # })
  # g_a <- - Ainv %*% a_param + sum(g_a1)
  # g_a <- 0


  g_all <- c(as.numeric(g_beta),
             as.numeric(g_tau),
             as.numeric(g_gamma),
             as.numeric(g_xi),
             as.numeric(g_a))
  return(g_all)
}

#' @rdname hmclearn-glm-posterior
#' @export
glmm_bin_posterior <- function(theta, y, X, Z, m, q=1, A = 1e4,
                               nulambda=1, Alambda=25, B=1e4) {
  Z <- as.matrix(Z)
  p <- ncol(X)
  n <- nrow(X)

  # prior covariance for beta
  Sig_beta <- diag(B, p, p)
  Sig_inv_beta <- diag(1/B, p, p)

  # extract parameters from theta vector
  beta_param <- theta[1:p]
  tau_param <- theta[(p+1):(p+m*q)]

  # diagonal of G matrix
  xi_param <- theta[(p+m*q+1):(p+m*q+q)]

  # off diagonal of G matrix
  if (q > 1) {
    a_param <- tail(theta, q*(q-1)/2)
    Ainv <- diag(1/A, length(a_param), length(a_param))
    log_a_prior <- -1/2*(t(a_param) %*% Ainv %*% a_param)
  } else {
    log_a_prior <- 0
    a_param <- 0
  }

  # reconstruct G LDLT decomposition
  Dhalf <- diag(exp(xi_param), q, q)

  L <- diag(q)
  L[lower.tri(L, diag=FALSE)] <- a_param

  LDhalf <- L %*% Dhalf
  LDhalf_block <- kronecker(diag(m), LDhalf)

  # u is deterministic function of xi, a, and tau
  u_param <- LDhalf_block %*% tau_param

  XZbetau <- X %*% beta_param + Z %*% u_param

  log_likelihood <- sum( -(1-y) * XZbetau - log(1 + exp(-XZbetau)))
  log_beta_prior <- -1/2 * p*log(B) - 1/2*t(beta_param) %*% Sig_inv_beta%*% beta_param
  log_tau_prior <- -1/2 * t(tau_param) %*% tau_param
  log_xi_prior <- -(nulambda + 1)/2 * log(1 + 1/nulambda * exp(2*xi_param) / Alambda^2)

  result <- log_likelihood + log_beta_prior + log_tau_prior + sum(log_xi_prior) + log_a_prior
  return(as.numeric(result))
}

#' @rdname hmclearn-glm-posterior
#' @export
g_glmm_bin_posterior <- function(theta, y, X, Z, m, q=1, A = 1e4,
                                 nulambda=1, Alambda=25, B=1e4) {
  Z <- as.matrix(Z)
  p <- ncol(X)
  n <- nrow(X)

  # prior covariance for beta
  Sig_beta <- diag(B, p, p)
  Sig_inv_beta <- diag(1/B, p, p)

  # extract parameters from theta vector
  beta_param <- theta[1:p]
  tau_param <- theta[(p+1):(p+m*q)]

  # diagonal of G matrix
  xi_param <- theta[(p+m*q+1):(p+m*q+q)]

  # off diagonal of G matrix
  if (q > 1) {
    a_param <- tail(theta, q*(q-1)/2)
    Ainv <- diag(1/A, length(a_param), length(a_param))
  } else {
    log_a_prior <- 0
    a_param <- 0
  }

  # reconstruct G LDLT decomposition
  Dhalf <- diag(exp(xi_param), q, q)

  L <- diag(q)
  L[lower.tri(L, diag=FALSE)] <- a_param

  LDhalf <- L %*% Dhalf
  LDhalf_block <- kronecker(diag(m), LDhalf)

  # u is deterministic function of xi, a, and tau
  u_param <- LDhalf_block %*% tau_param

  XZbetau <- X %*% beta_param + Z %*% u_param

  Ainv <- diag(1/A, length(a_param), length(a_param))

  # block L and xi for xi gradient
  L_block <- kronecker(diag(m), L)
  Dhalf_block <- kronecker(diag(m), Dhalf)

  # gradient
  gradprop <- 1 / (1 + exp(X %*% beta_param + Z %*% u_param))

  g_beta <- -t(1 - y) %*% X + t(gradprop) %*% X - t(Sig_inv_beta %*% beta_param)

  # tau gradient
  g_tau <- (-t(1-y) + t(gradprop)) %*% Z %*% LDhalf_block - tau_param

  # gradient for xi using matrix algebra
  zero_v <- rep(0, q)
  g_xi <- sapply(seq_along(1:q), function(jj) {
    zv <- zero_v
    zv[jj] <- 1
    bd <- kronecker(diag(m), diag(zv, q, q))
    (-t(1-y) + t(gradprop)) %*% Z %*% L_block %*% bd %*% Dhalf_block %*% tau_param
  })
  g_xi <- g_xi - (nulambda + 1) / (1 + nulambda*Alambda^2 * exp(-2*xi_param))

  # gradient for a
  if (q > 1) {
    tau_lst <- split(tau_param, ceiling(seq_along(tau_param)/q))
    Tj <- lapply(tau_lst, function(tauvals) {
      Tj <- create_Uj(Dhalf %*% tauvals, neg=FALSE)
    })
    Tj <- do.call(rbind, Tj)

    g_a <- -(t(1-y) + t(gradprop)) %*% Z %*% Tj -Ainv %*% a_param

    g_all <- c(as.numeric(g_beta),
               as.numeric(g_tau),
               g_xi,
               as.numeric(g_a))
  } else {
    g_all <- c(as.numeric(g_beta),
               as.numeric(g_tau),
               g_xi)
  }

  return(g_all)
}


#' @rdname hmclearn-glm-posterior
#' @export
glmm_poisson_posterior <- function(theta, y, X, Z, m, q=1, A = 1e4,
                                   nulambda=1, Alambda=25, B=1e4) {
  Z <- as.matrix(Z)
  p <- ncol(X)
  n <- nrow(X)
  q <- 1

  # prior covariance for beta
  Sig_beta <- diag(B, p, p)
  Sig_inv_beta <- diag(1/B, p, p)

  # extract parameters from theta vector
  beta_param <- theta[1:p]
  tau_param <- theta[(p+1):(p+m*q)]

  # diagonal of G matrix
  xi_param <- theta[(p+m*q+1):(p+m*q+q)]

  # off diagonal of G matrix
  if (q > 1) {
    a_param <- tail(theta, q*(q-1)/2)
    Ainv <- diag(1/A, length(a_param), length(a_param))
    log_a_prior <- -1/2*(t(a_param) %*% Ainv %*% a_param)
  } else {
    log_a_prior <- 0
    a_param <- 0
  }

  # reconstruct G LDLT decomposition
  Dhalf <- diag(exp(xi_param), q, q)

  L <- diag(q)
  L[lower.tri(L, diag=FALSE)] <- a_param

  LDhalf <- L %*% Dhalf
  LDhalf_block <- kronecker(diag(m), LDhalf)

  # u is deterministic function of xi, a, and tau
  u_param <- LDhalf_block %*% tau_param

  XZbetau <- X %*% beta_param + Z %*% u_param

  log_likelihood <- -sum(exp(XZbetau)) + y %*% XZbetau
  log_beta_prior <- -1/2 * p*log(B) - 1/2*t(beta_param) %*% Sig_inv_beta%*% beta_param
  log_tau_prior <- -1/2 * t(tau_param) %*% tau_param
  log_xi_prior <- -(nulambda + 1)/2 * log(1 + 1/nulambda * exp(2*xi_param) / Alambda^2)

  result <- log_likelihood + log_beta_prior + log_tau_prior + sum(log_xi_prior) + log_a_prior
  return(as.numeric(result))
}

#' @rdname hmclearn-glm-posterior
#' @export
g_glmm_poisson_posterior <- function(theta, y, X, Z, m, q=1, A = 1e4,
                                     nulambda=1, Alambda=25, B=1e4) {
  Z <- as.matrix(Z)
  p <- ncol(X)
  n <- nrow(X)
  q <- 1

  # prior covariance for beta
  Sig_beta <- diag(B, p, p)
  Sig_inv_beta <- diag(1/B, p, p)

  # extract parameters from theta vector
  beta_param <- theta[1:p]
  tau_param <- theta[(p+1):(p+m*q)]

  # diagonal of G matrix
  xi_param <- theta[(p+m*q+1):(p+m*q+q)]

  # off diagonal of G matrix
  if (q > 1) {
    a_param <- tail(theta, q*(q-1)/2)
    Ainv <- diag(1/A, length(a_param), length(a_param))
    log_a_prior <- -1/2*(t(a_param) %*% Ainv %*% a_param)
  } else {
    log_a_prior <- 0
    a_param <- 0
  }

  # reconstruct G LDLT decomposition
  Dhalf <- diag(exp(xi_param), q, q)

  L <- diag(q)
  L[lower.tri(L, diag=FALSE)] <- a_param

  LDhalf <- L %*% Dhalf
  LDhalf_block <- kronecker(diag(m), LDhalf)

  # u is deterministic function of xi, a, and tau
  u_param <- LDhalf_block %*% tau_param

  XZbetau <- X %*% beta_param + Z %*% u_param

  Ainv <- diag(1/A, length(a_param), length(a_param))

  # block L and xi for xi gradient
  L_block <- kronecker(diag(m), L)
  Dhalf_block <- kronecker(diag(m), Dhalf)

  # gradient
  g_beta <- -t(exp(XZbetau) - y) %*% X - t(Sig_inv_beta %*% beta_param)

  # tau gradient
  g_tau <- -t(exp(XZbetau) - y) %*% Z %*% LDhalf_block - tau_param

  # gradient for xi using matrix algebra
  zero_v <- rep(0, q)
  g_xi <- sapply(seq_along(1:q), function(jj) {
    zv <- zero_v
    zv[jj] <- 1
    bd <- kronecker(diag(m), diag(zv, q, q))
    -t(exp(XZbetau) - y) %*% Z %*% L_block %*% bd %*% Dhalf_block %*% tau_param
  })
  g_xi <- g_xi - (nulambda + 1) / (1 + nulambda*Alambda^2 * exp(-2*xi_param))

  if (q > 1) {
    tau_lst <- split(tau_param, ceiling(seq_along(tau_param)/q))
    Tj <- lapply(tau_lst, function(tauvals) {
      Tj <- create_Uj(Dhalf %*% tauvals, neg=FALSE)
    })
    # Tja <- as.matrix(bdiag(Tja))
    Tj <- do.call(rbind, Tj)

    g_a <- -t(exp(XZbetau)) %*% Z %*% Tj + t(y)%*%Z%*%Tj -Ainv %*% a_param

    g_all <- c(as.numeric(g_beta),
               as.numeric(g_tau),
               g_xi,
               as.numeric(g_a))
  } else {
    g_all <- c(as.numeric(g_beta),
               as.numeric(g_tau),
               g_xi)
  }
  return(g_all)
}

