# proposal function in form expected by mh
pfun_glmm_poisson <- function(PARAM, ...) {
  d <- length(PARAM)
  theta <- PARAM
  glmm_poisson_posterior(theta=theta, ...)
}

# prior for beta is mean 0 with diagonal covariance B^-1
#' @export
glmm_poisson_posterior <- function(theta, y, X, Z, m=10, q=1, A = 1e4, B=1e4,
                                   nuxi=1, Axi=25) {
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

  log_likelihood <- -sum(exp(XZbetau)) + y %*% XZbetau
  log_beta_prior <- -1/2 * p*log(B) - 1/2*t(beta_param) %*% Sig_inv_beta%*% beta_param
  log_tau_prior <- -1/2 * t(tau_param) %*% tau_param
  log_xi_prior <- -(nuxi + 1)/2 * log(1 + 1/nuxi * exp(2*xi_param) / Axi^2)

  result <- log_likelihood + log_beta_prior + log_tau_prior + sum(log_xi_prior) + log_a_prior
  return(as.numeric(result))
}

#' @export
g_glmm_poisson_posterior <- function(theta, y, X, Z, m=10, q=1, A = 1e4, B=1e4,
                                     nuxi=1, Axi=25) {
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
  g_xi <- g_xi - (nuxi + 1) / (1 + nuxi*Axi^2 * exp(-2*xi_param))

  # TODO:  gradient for a
  if (q > 1) {
    u_lst <- split(u_param, ceiling(seq_along(u_param)/q))
    g_a1 <- sapply(u_lst, function(uvals) {
      Uj <- create_Uj(uvals)
      as.numeric(t(Uj) %*% Dinv %*% (uvals - Uj %*% a_param))
    })
    g_a <- - Ainv %*% a_param + rowSums(g_a1)

    g_all <- c(as.numeric(g_beta),
               as.numeric(g_u),
               g_xi,
               as.numeric(g_a))
  } else {
    g_all <- c(as.numeric(g_beta),
               as.numeric(g_tau),
               g_xi)
  }

  return(g_all)

}
