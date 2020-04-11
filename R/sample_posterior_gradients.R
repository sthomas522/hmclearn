
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

#' Sample log posterior and gradient functions for select generalized linear models
#'
#' @name hmclearn-glm-posterior
#'
#' @param theta vector of parameters.  Stored as a single vector in order fixed effect, random effect, log-transformed diagonal \eqn{\lambda}, and off-diagonal of \code{G} vector \code{a}
#' @param y numeric vector for the dependent variable
#' @param X numeric design matrix of fixed effect parameters
#' @param a hyperprior for the Inverse Gamma shape parameter
#' @param b hyperprior for the Inverse Gamma scale parameter
#' @param B prior for linear predictors is multivariate Normal with mean 0 with diagonal covariance B^-1
#'
#' @section Models with available posterior and gradient functions:
#' \describe{
#'   \item{`linear_posterior(theta, y, X, a=1e-4, b=1e-4, B=0.001)`}{
#'    The likelihood function for linear regression
#'    \deqn{p(y | X, \beta; \sigma^2) = \frac{1}{(2\pi\sigma^2)^{n/2}}\exp{\left(-\frac{1}{2\sigma^2} (y - X\beta)^T(y-X\beta) \right)}}
#'    with priors \eqn{p(\sigma^2) \sim IG(a, b)} and \eqn{\beta \sim N(0, BI)}.  The variance term is log transformed \eqn{\gamma = \log\sigma^2}
#'    The input parameter vector \code{theta} is of length \code{k}.  The first \code{k-1} parameters are for \eqn{\beta}, and the last parameter is \eqn{\gamma}
#'    Note that the Inverse Gamma prior can be problematic for certain applications with low variance.  See Gelman (2006)
#'    @return numeric value for the log posterior
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
