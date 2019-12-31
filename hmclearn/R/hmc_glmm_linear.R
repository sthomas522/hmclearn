
# proposal function in form expected by mh
pfun_lmm <- function(PARAM, ...) {
  d <- length(PARAM)
  theta <- PARAM
  lmm_posterior(theta=theta, ...)
}

#' Linear Mixed Effects model log posterior
#'
#' Compute the log posterior of a linear mixed effects regression model.
#' Priors are multivariate Normal for the fixed effects
#' and half-t for the error
#'
#' @param theta vector of parameters.  Stored as a single vector in order fixed effect, random effect, log-transformed error, log-transformed diagonal \eqn{\lambda}, and off-diagonal of \code{G} vector \code{a}
#' @param y numeric vector for the dependent variable
#' @param X numeric design matrix of fixed effect parameters
#' @param Z numeric design matrix of random effect parameters
#' @param m number of random effect parameters
#' @param A hyperprior numeric vector for the random effects off-diagonal \code{a}
#' @param nueps hyperprior for the half-t prior of the error parameter \eqn{\nu}
#' @param nulambda hyperprior for the half-t prior of the random effects diagonal \eqn{\lambda}
#' @param Aeps hyperprior for the half-t prior of the error parameter \eqn{A_\epsilon}
#' @param Alambda hyperprior for the half-t prior of the random effects diagonal \eqn{A_\lambda}
#' @param B prior for linear predictors is multivariate Normal with mean 0 with diagonal covariance B^-1
#' @details The likelihood function for linear mixed effect regression
#' @details \deqn{p(y | \beta, u, \sigma_\epsilon^2) &\propto (\sigma_\epsilon^2)^{-n/2} e^{-\frac{1}{2\sigma_\epsilon^2}(y - X\beta - Zu)^T (y - X\beta - Zu)}}
#' @details with priors \eqn{\beta \sim N(0, BI)}, \eqn{\sigma_\epsilon \sim half-t(A_\epsilon, nu_\epsilon)}, \eqn{\lambda \sim half-t(A_\lambda, nu_\lambda )}.
#' @details The vector \eqn{\lambda} is the diagonal of the covariance \code{G} hyperprior where \eqn{u \sim N(0, G}.  The off-diagonal hyperpriors are stored in a vector \eqn{a \sim N(0, A}.  See Chan, Jeliazkov (2009) for details.
#' @details The input parameter vector \code{theta} is of length \code{k}.  The first \code{k-1} parameters are for \eqn{\beta}, and the last parameter is \eqn{\gamma}
#' @return numeric value for the log posterior
#' @references Gelman, A. (2006). \emph{Prior distributions for variance parameters in hierarchical models (comment on article by Browne and Draper)}. Bayesian analysis, 1(3), 515-534.
#' @references Chan, J. C. C., & Jeliazkov, I. (2009). \emph{MCMC estimation of restricted covariance matrices}. Journal of Computational and Graphical Statistics, 18(2), 457-480.
#' @references Betancourt, M., & Girolami, M. (2015). \emph{Hamiltonian Monte Carlo for hierarchical models}. Current trends in Bayesian methodology with applications, 79, 30.
#' @export
lmm_posterior <- function(theta, y, X, Z, m, A = 1e4, nueps=1, nulambda=1, Aeps=25, Alambda=25, B=0.001) {
  Z <- as.matrix(Z)
  p <- ncol(X)
  n <- nrow(X)
  q <- 1

  # extract parameters from theta vector
  beta_param <- as.numeric(theta[1:p])
  tau_param <- theta[(p+1):(p+m*q)]

  # log of sig2eps
  gamma_param <- theta[p+m*q+1]

  # diagonal of G matrix
  xi_param <- theta[(p+m*q+2):(p+m*q+q+1)]

  # off diagonal of G matrix
  a_param <- tail(theta, q*(q-1)/2)
  # a_param <- 0

  nu <- diag(1/B, p, p)
  inv.nu <- diag(B, p, p)

  # reconstruct L D^1/2
  Dhalf <- diag(exp(xi_param))

  L <- diag(q)
  L[lower.tri(L, diag=FALSE)] <- a_param

  LDhalf <- L %*% Dhalf

  Ainv <- diag(1/A, length(a_param), length(a_param))

  LDhalf_block <- kronecker(diag(m), LDhalf)

  # u is deterministic function of xi, a, and tau
  u_param <- LDhalf_block %*% tau_param

  log_likelihood <- -n*gamma_param - exp(-2*gamma_param)/2 *
    t(y - X%*%beta_param - Z%*%u_param) %*% (y - X%*%beta_param - Z%*%u_param)
  log_beta_prior <- - 1/2* t(beta_param) %*% inv.nu %*% beta_param
  log_tau_prior <- -1/2 * t(tau_param) %*% tau_param
  log_gamma_prior <- -(nueps + 1)/2 * log(1 + 1/nueps*exp(2*gamma_param)/Aeps^2) + gamma_param
  log_xi_prior <- -(nulambda+1)/2 * log(1 + 1/nulambda*exp(2*xi_param)/Alambda^2) + xi_param
  log_a_prior <- -1/2*(t(a_param) %*% Ainv %*% a_param)

  result <- log_likelihood + log_beta_prior + log_tau_prior + log_gamma_prior + sum(log_xi_prior) + log_a_prior
  return(as.numeric(result))
}


#' Gradient of linear Mixed Effects model log posterior
#'
#' Compute the gradient of the log posterior of a linear mixed effects regression model.
#' Priors are multivariate Normal for the fixed effects
#' and half-t for the error
#'
#' @param theta vector of parameters.  Stored as a single vector in order fixed effect, random effect, log-transformed error, log-transformed diagonal \eqn{\lambda}, and off-diagonal of \code{G} vector \code{a}
#' @param y numeric vector for the dependent variable
#' @param X numeric design matrix of fixed effect parameters
#' @param Z numeric design matrix of random effect parameters
#' @param m number of random effect parameters
#' @param A hyperprior numeric vector for the random effects off-diagonal \code{a}
#' @param nueps hyperprior for the half-t prior of the error parameter \eqn{\nu}
#' @param nulambda hyperprior for the half-t prior of the random effects diagonal \eqn{\lambda}
#' @param Aeps hyperprior for the half-t prior of the error parameter \eqn{A_\epsilon}
#' @param Alambda hyperprior for the half-t prior of the random effects diagonal \eqn{A_\lambda}
#' @param B prior for linear predictors is multivariate Normal with mean 0 with diagonal covariance B^-1
#' @details The likelihood function for linear mixed effect regression
#' @details \deqn{p(y | \beta, u, \sigma_\epsilon^2) &\propto (\sigma_\epsilon^2)^{-n/2} e^{-\frac{1}{2\sigma_\epsilon^2}(y - X\beta - Zu)^T (y - X\beta - Zu)}}
#' @details with priors \eqn{\beta \sim N(0, BI)}, \eqn{\sigma_\epsilon \sim half-t(A_\epsilon, nu_\epsilon)}, \eqn{\lambda \sim half-t(A_\lambda, nu_\lambda )}.
#' @details The vector \eqn{\lambda} is the diagonal of the covariance \code{G} hyperprior where \eqn{u \sim N(0, G}.  The off-diagonal hyperpriors are stored in a vector \eqn{a \sim N(0, A}.  See Chan, Jeliazkov (2009) for details.
#' @details The input parameter vector \code{theta} is of length \code{k}.  The first \code{k-1} parameters are for \eqn{\beta}, and the last parameter is \eqn{\gamma}
#' @return numeric vector of the gradient of the log posterior
#' @references Gelman, A. (2006). \emph{Prior distributions for variance parameters in hierarchical models (comment on article by Browne and Draper)}. Bayesian analysis, 1(3), 515-534.
#' @references Chan, J. C. C., & Jeliazkov, I. (2009). \emph{MCMC estimation of restricted covariance matrices}. Journal of Computational and Graphical Statistics, 18(2), 457-480.
#' @references Betancourt, M., & Girolami, M. (2015). \emph{Hamiltonian Monte Carlo for hierarchical models}. Current trends in Bayesian methodology with applications, 79, 30.
#' @export
g_lmm_posterior <- function(theta, y, X, Z, m, A = 1e4, nueps=1, nulambda=1, Aeps=25, Alambda=25, B=0.001) {
  Z <- as.matrix(Z)
  p <- ncol(X)
  n <- nrow(X)
  q <- 1

  # extract parameters from theta vector
  beta_param <- as.numeric(theta[1:p])
  tau_param <- theta[(p+1):(p+m*q)]

  # log of sig2eps
  gamma_param <- theta[p+m*q+1]

  # diagonal of G matrix
  xi_param <- theta[(p+m*q+2):(p+m*q+q+1)]

  # off diagonal of G matrix
  a_param <- tail(theta, q*(q-1)/2)
  # a_param <- 0

  nu <- diag(1/B, p, p)
  inv.nu <- diag(B, p, p)

  # reconstruct L D^1/2
  Dhalf <- diag(exp(xi_param))
  Dhalf_block <- kronecker(diag(m), Dhalf)

  L <- diag(q)
  L[lower.tri(L, diag=FALSE)] <- a_param
  L_block <- kronecker(diag(m), L)

  LDhalf <- L %*% Dhalf

  Ainv <- diag(1/A, length(a_param), length(a_param))

  LDhalf_block <- kronecker(diag(m), LDhalf)

  # u is deterministic function of xi, a, and tau
  u_param <- LDhalf_block %*% tau_param

  # gradient
  g_beta <- exp(-2*gamma_param) * t(X) %*% (y - X%*%beta_param - Z%*%u_param) - (inv.nu %*% beta_param)
  g_tau <- exp(-2*gamma_param) * Dhalf_block %*% t(L_block) %*% t(Z) %*%
    (y - X%*%beta_param - Z%*%u_param) - tau_param
  g_gamma <- -(n-1) + exp(-2*gamma_param) * t(y - X%*%beta_param - Z%*%u_param) %*%
    (y - X%*%beta_param - Z%*%u_param) -
    (nueps + 1) / (1 + nueps*Aeps^2 * exp(-2*gamma_param))

  # gradient for xi using matrix algebra
  zero_v <- rep(0, q)
  g_xi <- sapply(seq_along(1:q), function(jj) {
    zv <- zero_v
    zv[jj] <- 1
    bd <- kronecker(diag(m), diag(zv))
    tot <- (y - X%*%beta_param - Z%*%u_param) %*% t(Z %*% L_block %*% bd %*% Dhalf_block %*% tau_param)
    exp(-2*gamma_param) * sum(diag(tot))
  })
  g_xi <- g_xi - (nulambda + 1) / (1 + nulambda*Alambda^2 * exp(-2*xi_param)) + 1

  # gradient for a
  tau_tilde <- exp(xi_param) * a_param
  T_tilde <- create_Uj(tau_tilde, neg=FALSE)

  # gradient for a
  tau_lst <- split(tau_param, ceiling(seq_along(tau_param)/q))
  Tj <- lapply(tau_lst, function(tauvals) {
    Tj <- create_Uj(Dhalf %*% tauvals, neg=FALSE)
  })
  # Tja <- as.matrix(bdiag(Tja))
  Tj <- do.call(rbind, Tj)

  g_a1 <- exp(-2*gamma_param) * t(Tj) %*% t(Z) %*%
    (y - X%*%beta_param - Z%*%Dhalf_block%*%tau_param - Z%*%Tj%*%a_param)
  g_a2 <- -Ainv %*% a_param
  g_a <- g_a1 + g_a2

  # g_a1 <- sapply(tau_lst, function(tauvals) {
  #   Tj <- create_Uj(tauvals)
  #   exp(-2*gamma_param) * t(Tj) %*% Dinv %*% (uvals - Uj * a_param)
  # })
  # g_a <- - Ainv %*% a_param + sum(g_a1)
  # g_a <- 0


  g_all <- c(as.numeric(g_beta),
             as.numeric(g_tau),
             as.numeric(g_gamma),
             as.numeric(g_xi),
             as.numeric(g_a))
  return(g_all)

}

