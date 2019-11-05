
# proposal function in form expected by mh
pfun_poisson <- function(PARAM, ...) {
  d <- length(PARAM)
  theta <- PARAM
  poisson_posterior(theta=theta, ...)
}

# prior for beta is mean 0 with diagonal covariance B^-1
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

# gradient of the log posterior for hmc
g_poisson_posterior <- function(beta_param, y, X, B=0.01) {
  n <- length(y)
  k <- length(beta_param)

  nu <- diag(1/B, k, k)
  inv.nu <- diag(B, k, k)

  # gradient
  y %*% X - crossprod(exp(X %*% beta_param), X) - crossprod(beta_param, inv.nu)
}

