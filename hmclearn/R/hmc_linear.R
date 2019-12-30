

#' Linear regression log posterior
#'
#' Compute the log posterior of a linear regression model.
#' Priors are multivariate Normal for the predictor variables
#' and Inverse Gamma for the error term squared.
#'
#' @param theta vector of parameters.  All but the last parameter should be linear predictor parameter values.  The last parameter should be the log of the squared error parameter
#' @param y numeric vector for the dependent variable
#' @param X numeric design matrix of independent variables
#' @param a hyperprior for the Inverse Gamma shape parameter
#' @param b hyperprior for the Inverse Gamma scale parameter
#' @param B inverse of the diagonal of the Covariance hyperprior for the linear predictors
#' @details The likelihood function for linear regression
#' @details \deqn{p(y | X, \beta; \sigma^2) = \frac{1}{(2\pi\sigma^2)^{n/2}}\exp{\left(-\frac{1}{2\sigma^2} (y - X\beta)^T(y-X\beta) \right)}}
#' @details with priors \eqn{p(\sigma^2) \sim IG(a, b)} and \eqn{\beta \sim N(0, BI)}.  The variance term is log transformed \eqn{\gamma = \log\sigma^2}
#' @details The input parameter vector \code{theta} is of length \code{k}.  The first \code{k-1} parameters are for \eqn{\beta}, and the last parameter is \eqn{\gamma}
#' @details Note that the Inverse Gamma prior can be problematic for certain applications with low variance.  See Gelman (2006)
#' @return numeric value for the log posterior
#' @references Gelman, A. (2006). \emph{Prior distributions for variance parameters in hierarchical models (comment on article by Browne and Draper)}. Bayesian analysis, 1(3), 515-534.
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

#' Gradient of linear regression log posterior
#'
#' Compute the gradient of the log posterior of a linear regression model.
#' Priors are multivariate Normal for the predictor variables
#' and Inverse Gamma for the error term squared.
#'
#' @param theta vector of parameters.  All but the last parameter should be linear predictor parameter values.  The last parameter should be the log of the squared error parameter
#' @param y numeric vector for the dependent variable
#' @param X numeric design matrix of independent variables
#' @param a hyperprior for the Inverse Gamma shape parameter
#' @param b hyperprior for the Inverse Gamma scale parameter
#' @param B inverse of the diagonal of the Covariance hyperprior for the linear predictors
#' @details The likelihood function for linear regression
#' @details \deqn{p(y | X, \beta; \sigma^2) = \frac{1}{(2\pi\sigma^2)^{n/2}}\exp{\left(-\frac{1}{2\sigma^2} (y - X\beta)^T(y-X\beta) \right)}}
#' @details with priors \eqn{p(\sigma^2) \sim IG(a, b)} and \eqn{\beta \sim N(0, BI)}.  The variance term is log transformed \eqn{\gamma = \log\sigma^2}
#' @details The input parameter vector \code{theta} is of length \code{k}.  The first \code{k-1} parameters are for \eqn{\beta}, and the last parameter is \eqn{\gamma}
#' @details Note that the Inverse Gamma prior can be problematic for certain applications with low variance.  See Gelman (2006)
#' @return numeric vector for the gradient of the log posterior
#' @references Gelman, A. (2006). \emph{Prior distributions for variance parameters in hierarchical models (comment on article by Browne and Draper)}. Bayesian analysis, 1(3), 515-534.
#' @export
g_linear_posterior <- function(theta, y, X, a=1e-4, b=1e-4, B=0.001) {
  # parameters
  k <- length(theta)
  beta_param <- as.numeric(theta[1:(k-1)])
  gamma_param <- theta[k]
  n <- nrow(X)

  nu <- diag(1/B, k-1, k-1)
  inv.nu <- diag(B, k-1, k-1)

  # gradient of beta
  grad_beta <- exp(-gamma_param)  * t(X) %*% (y - X%*%beta_param)  - (inv.nu %*% beta_param)
  grad_gamma <- -(n/2 + a) + exp(-gamma_param)/2 * t(y - X%*%beta_param) %*%
    (y - X%*%beta_param) + b*exp(-gamma_param)
  c(as.numeric(grad_beta), as.numeric(grad_gamma))
}


