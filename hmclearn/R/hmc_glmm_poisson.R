# proposal function in form expected by mh
pfun_glmm_poisson <- function(PARAM, ...) {
  d <- length(PARAM)
  theta <- PARAM
  glmm_poisson_posterior(theta=theta, ...)
}

#' Poisson Mixed Effects model log posterior
#'
#' Compute the log posterior of a Poisson mixed effects regression model.
#' Priors are multivariate Normal for the fixed effects
#'
#' @param theta vector of parameters.  Stored as a single vector in order fixed effect, random effect, log-transformed diagonal \eqn{\lambda}, and off-diagonal of \code{G} vector \code{a}
#' @param y numeric vector for the dependent variable
#' @param X numeric design matrix of fixed effect parameters
#' @param Z numeric design matrix of random effect parameters
#' @param m number of random effect parameters
#' @param A hyperprior numeric vector for the random effects off-diagonal \code{a}
#' @param nulambda hyperprior for the half-t prior of the random effects diagonal \eqn{\lambda}
#' @param Alambda hyperprior for the half-t prior of the random effects diagonal \eqn{A_\lambda}
#' @param B prior for linear predictors is multivariate Normal with mean 0 with diagonal covariance B^-1
#' @details The likelihood function for Poisson mixed effect regression
#' @details \deqn{L(\beta; y, X) = \prod_{i=1}^n \prod_{j=1}^m \frac{e^{-e^{X_i\beta + Z_{ij}u_{ij}}}e^{y_i(X_i\beta + Z_{ij}u_{ij})}}{y_i!} }
#' @details with priors \eqn{\beta \sim N(0, BI)}, \eqn{\sigma_\epsilon \sim half-t(A_\epsilon, nu_\epsilon)}, \eqn{\lambda \sim half-t(A_\lambda, nu_\lambda )}.
#' @details The vector \eqn{\lambda} is the diagonal of the covariance \code{G} hyperprior where \eqn{u \sim N(0, G}.  The off-diagonal hyperpriors are stored in a vector \eqn{a \sim N(0, A}.  See Chan, Jeliazkov (2009) for details.
#' @details The input parameter vector \code{theta} is of length \code{k}.  The first \code{k-1} parameters are for \eqn{\beta}, and the last parameter is \eqn{\gamma}
#' @return numeric value for the log posterior
#' @references Gelman, A. (2006). \emph{Prior distributions for variance parameters in hierarchical models (comment on article by Browne and Draper)}. Bayesian analysis, 1(3), 515-534.
#' @references Chan, J. C. C., & Jeliazkov, I. (2009). \emph{MCMC estimation of restricted covariance matrices}. Journal of Computational and Graphical Statistics, 18(2), 457-480.
#' @references Betancourt, M., & Girolami, M. (2015). \emph{Hamiltonian Monte Carlo for hierarchical models}. Current trends in Bayesian methodology with applications, 79, 30.
#' @export
glmm_poisson_posterior <- function(theta, y, X, Z, m=10, A = 1e4,
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

#' Gradient of a Poisson Mixed Effects model log posterior
#'
#' Compute the gradient of the log posterior of a Poisson mixed effects regression model.
#' Priors are multivariate Normal for the fixed effects
#'
#' @param theta vector of parameters.  Stored as a single vector in order fixed effect, random effect, log-transformed diagonal \eqn{\lambda}, and off-diagonal of \code{G} vector \code{a}
#' @param y numeric vector for the dependent variable
#' @param X numeric design matrix of fixed effect parameters
#' @param Z numeric design matrix of random effect parameters
#' @param m number of random effect parameters
#' @param A hyperprior numeric vector for the random effects off-diagonal \code{a}
#' @param nulambda hyperprior for the half-t prior of the random effects diagonal \eqn{\lambda}
#' @param Alambda hyperprior for the half-t prior of the random effects diagonal \eqn{A_\lambda}
#' @param B prior for linear predictors is multivariate Normal with mean 0 with diagonal covariance B^-1
#' @details The likelihood function for Poisson mixed effect regression
#' @details \deqn{L(\beta; y, X) = \prod_{i=1}^n \prod_{j=1}^m \frac{e^{-e^{X_i\beta + Z_{ij}u_{ij}}}e^{y_i(X_i\beta + Z_{ij}u_{ij})}}{y_i!} }
#' @details with priors \eqn{\beta \sim N(0, BI)}, \eqn{\sigma_\epsilon \sim half-t(A_\epsilon, nu_\epsilon)}, \eqn{\lambda \sim half-t(A_\lambda, nu_\lambda )}.
#' @details The vector \eqn{\lambda} is the diagonal of the covariance \code{G} hyperprior where \eqn{u \sim N(0, G}.  The off-diagonal hyperpriors are stored in a vector \eqn{a \sim N(0, A}.  See Chan, Jeliazkov (2009) for details.
#' @details The input parameter vector \code{theta} is of length \code{k}.  The first \code{k-1} parameters are for \eqn{\beta}, and the last parameter is \eqn{\gamma}
#' @return numeric vector for the gradient of the log posterior
#' @references Gelman, A. (2006). \emph{Prior distributions for variance parameters in hierarchical models (comment on article by Browne and Draper)}. Bayesian analysis, 1(3), 515-534.
#' @references Chan, J. C. C., & Jeliazkov, I. (2009). \emph{MCMC estimation of restricted covariance matrices}. Journal of Computational and Graphical Statistics, 18(2), 457-480.
#' @references Betancourt, M., & Girolami, M. (2015). \emph{Hamiltonian Monte Carlo for hierarchical models}. Current trends in Bayesian methodology with applications, 79, 30.
#' @export
g_glmm_poisson_posterior <- function(theta, y, X, Z, m=10, A = 1e4,
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
