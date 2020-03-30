
# proposal function in form expected by mh
pfun_poisson <- function(PARAM, ...) {
  d <- length(PARAM)
  theta <- PARAM
  poisson_posterior(theta=theta, ...)
}

#' Poisson regression log posterior
#'
#' Compute the log posterior of a poisson regression model.
#' Priors are multivariate Normal for the linear predictors
#'
#' @param theta vector of linear predictor parameters
#' @param y numeric vector for the dependent variable
#' @param X numeric design matrix of independent variables
#' @param B prior for linear predictors is multivariate Normal with mean 0 with diagonal covariance B^-1
#' @details The likelihood function for poisson regression
#' @details \deqn{L(\beta; y, X) = \prod_{i=1}^n \frac{e^{-e^{X_i\beta}}e^{y_iX_i\beta}}{y_i!}}
#' @details with priors \eqn{\beta \sim N(0, BI)}.
#' @details The input parameter vector \code{theta} is of length \code{k}, containing parameter values for \eqn{\beta}
#' @return numeric value for the log posterior
#' @export
poisson_posterior <- function(theta, y, X, B=0.01) {
  k <- length(theta)
  beta_param <- theta

  nu <- diag(1/B, k, k)
  inv.nu <- diag(B, k, k)
  # result <- -(y + alpha1 - 1)* log(1 + exp(-X %*% beta)) + (1 - y + beta1 - 1) *
  #   log(exp(-X %*% beta) / (1 + exp(-X %*% beta)))

  # loglikelihood + log beta + log sig2
  result <- log_lik_poisson(beta_param=beta_param, y=y, X=X) -
    - 1/2* t(beta_param) %*% inv.nu %*% beta_param

  return(result)
}

# loglikelihood poisson regression
# note that factorial omitted
log_lik_poisson <- function(beta_param, y, X) {
  y %*% X %*% beta_param - sum(exp( X %*% beta_param)) # - sum(lfactorial(y))
}

#' Gradient of Poisson regression log posterior
#'
#' Compute the gradient of the log posterior of a poisson regression model.
#' Priors are multivariate Normal for the linear predictors
#'
#' @param theta vector of linear predictor parameters
#' @param y numeric vector for the dependent variable
#' @param X numeric design matrix of independent variables
#' @param B prior for linear predictors is multivariate Normal with mean 0 with diagonal covariance B^-1
#' @details The likelihood function for poisson regression
#' @details \deqn{L(\beta; y, X) = \prod_{i=1}^n \frac{e^{-e^{X_i\beta}}e^{y_iX_i\beta}}{y_i!}}
#' @details with priors \eqn{\beta \sim N(0, BI)}.
#' @details The input parameter vector \code{theta} is of length \code{k}, containing parameter values for \eqn{\beta}
#' @return numeric vector for the gradient of the log posterior
#' @export
g_poisson_posterior <- function(theta, y, X, B=0.01) {
  n <- length(y)
  k <- length(theta)

  nu <- diag(1/B, k, k)
  inv.nu <- diag(B, k, k)

  # gradient
  y %*% X - crossprod(exp(X %*% theta), X) - crossprod(theta, inv.nu)
}

