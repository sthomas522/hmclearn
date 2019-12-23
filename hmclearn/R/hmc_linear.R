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


