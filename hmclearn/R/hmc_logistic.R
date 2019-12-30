
# proposal function in form expected by mh
pfun_logistic <- function(PARAM, ...) {
  d <- length(PARAM)
  theta <- PARAM
  logistic_posterior(theta=theta, ...)
}

#' Logistic regression log posterior
#'
#' Compute the log posterior of a logistic regression model.
#' Priors are multivariate Normal for the linear predictors
#'
#' @param theta vector of linear predictor parameters
#' @param y numeric vector for the dependent variable
#' @param X numeric design matrix of independent variables
#' @param B prior for linear predictors is multivariate Normal with mean 0 with diagonal covariance B^-1
#' @details The likelihood function for logistic regression
#' @details \deqn{L(\beta; X, y) = \prod_{i=1}^{n} \left(\frac{1}{1+e^{-X_i\beta}}\right)^{y_i} \left(\frac{e^{-X_i\beta}}{1+e^{-X_i\beta}}\right)^{1-y_i}}
#' @details with priors \eqn{\beta \sim N(0, BI)}.
#' @details The input parameter vector \code{theta} is of length \code{k}, containing parameter values for \eqn{\beta}
#' @return numeric value for the log posterior
#' @export
logistic_posterior <- function(theta, y, X, B=0.01) {
  k <- length(theta)
  beta_param <- as.numeric(theta)

  nu <- diag(1/B, k, k)
  inv.nu <- diag(B, k, k)
  # result <- -(y + alpha1 - 1)* log(1 + exp(-X %*% beta)) + (1 - y + beta1 - 1) *
  #   log(exp(-X %*% beta) / (1 + exp(-X %*% beta)))

  # loglikelihood + log beta + log sig2
  result <- log_lik_bin(beta_param=beta_param, y=y, X=X) +
    -1/2*log(det(2*pi*nu)) - 1/2* t(beta_param) %*% inv.nu %*% beta_param

  return(result)
}

log_lik_bin <- function(beta_param, y, X) {
  # n <- length(y)
  # result <- -y* log(1 + exp(-X %*% beta)) + (1 - y) * log(exp(-X %*% beta) / (1 + exp(-X %*% beta))  )
  # result <- y* X %*% beta - X %*% beta - log(1 + exp(-X %*% beta))
  result <- X %*% beta_param * (y - 1) - log(1 + exp(-X %*% beta_param))
  sum(result)
}


#' Gradient of logistic regression log posterior
#'
#' Compute the gradient of the log posterior of a logistic regression model.
#' Priors are multivariate Normal for the linear predictors
#'
#' @param theta vector of linear predictor parameters
#' @param y numeric vector for the dependent variable
#' @param X numeric design matrix of independent variables
#' @param B prior for linear predictors is multivariate Normal with mean 0 with diagonal covariance B^-1
#' @details The likelihood function for logistic regression
#' @details \deqn{L(\beta; X, y) = \prod_{i=1}^{n} \left(\frac{1}{1+e^{-X_i\beta}}\right)^{y_i} \left(\frac{e^{-X_i\beta}}{1+e^{-X_i\beta}}\right)^{1-y_i}}
#' @details with priors \eqn{\beta \sim N(0, BI)}.
#' @details The input parameter vector \code{theta} is of length \code{k}, containing parameter values for \eqn{\beta}
#' @return numeric vector for the gradient of the log posterior
#' @export
g_logistic_posterior <- function(beta, y, X, B=0.01) {
  n <- length(y)
  k <- length(beta)
  beta <- as.numeric(beta)

  nu <- diag(1/B, k, k)
  inv.nu <- diag(B, k, k)

  result <- t(y-1) %*% X %*% diag(1, k, k) +
    t(exp(-X %*% beta) / (1 + exp(-X %*% beta))) %*% X - t(inv.nu %*% beta)
  return(result)
}
