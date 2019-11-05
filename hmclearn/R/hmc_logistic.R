
# proposal function in form expected by mh
pfun_logistic <- function(PARAM, ...) {
  d <- length(PARAM)
  theta <- PARAM
  logistic_posterior(theta=theta, ...)
}

# prior for beta is mean 0 with diagonal covariance B^-1
logistic_posterior <- function(theta, y, X, B=0.01) {
  k <- length(theta)
  beta_param <- theta

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


# gradient of the log posterior for hmc
g_logistic_posterior <- function(beta, y, X, B=0.01) {
  n <- length(y)
  k <- length(beta)

  nu <- diag(1/B, k, k)
  inv.nu <- diag(B, k, k)

  result <- t(y-1) %*% X %*% diag(1, k, k) +
    t(exp(-X %*% beta) / (1 + exp(-X %*% beta))) %*% X - t(inv.nu %*% beta)
  return(result)
}
