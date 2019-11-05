#' @export
linear_posterior <- function(theta, y, X, a=1e-4, b=1e-4) {
  k <- length(theta)
  beta_param <- as.numeric(theta[1:(k-1)])
  gamma_param <- theta[k]

  n <- nrow(X)
  -(n/2 + a)*gamma_param - exp(-gamma_param)/2 * t(y - X%*%beta_param) %*%
    (y - X%*%beta_param) - b*exp(-gamma_param)
}

#' @export
g_linear_posterior <- function(theta, y, X, a=1e-4, b=1e-4) {
  # parameters
  k <- length(theta)
  beta_param <- as.numeric(theta[1:(k-1)])
  gamma_param <- theta[k]
  n <- nrow(X)

  # gradient of beta
  grad_beta <- exp(-gamma_param)  * t(X) %*% (y - X%*%beta_param)
  grad_gamma <- -(n/2 + a) + exp(-gamma_param)/2 * t(y - X%*%beta_param) %*%
    (y - X%*%beta_param) + b*exp(-gamma_param)
  c(as.numeric(grad_beta), as.numeric(grad_gamma))
}


