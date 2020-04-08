
# Linear regression example
set.seed(521)
X <- cbind(1, matrix(rnorm(300), ncol=3))
y <- 0.5*X[, 1] -1*X[, 2] + 2*X[, 3] -3*X[, 4] +
  rnorm(100, sd=.2)

f1 <- hmc(N = 500,
          theta.init = c(rep(0, 4), 1),
          epsilon = 0.01,
          L = 10,
          logPOSTERIOR = linear_posterior,
          glogPOSTERIOR = g_linear_posterior,
          varnames = c(paste0("beta", 0:3), "log_sigma_sq"),
          param=list(y=y, X=X), parallel=F, chains=1)

medparam1 <- as.vector(summary(f1)[, 3])

expect_equal(round(medparam1, 8),
             c(0.5325537,
               -1.0118904,
               2.0164956,
               -2.9807522,
               -3.0392322))

# Logistic regression example
X <- cbind(1, seq(-100, 100, by=0.25))
betavals <- c(-0.9, 0.2)
lodds <- X %*% betavals
prob1 <- as.numeric(1 / (1 + exp(-lodds)))

set.seed(9874)
y <- sapply(prob1, function(xx) {
  sample(c(0, 1), 1, prob=c(1-xx, xx))
})

f2 <- hmc(N = 500,
          theta.init = rep(0, 2),
          epsilon = c(0.1, 0.002),
          L = 10,
          logPOSTERIOR = logistic_posterior,
          glogPOSTERIOR = g_logistic_posterior,
          randlength=TRUE,
          varnames = paste0("beta", 0:1),
          param = list(y=y, X=X),
          parallel=F, chains=1)

medparam2 <- as.vector(summary(f2, burnin=100)[, 3])
medparam2

expect_equal(round(medparam2, 7),
             c(-0.8497155, 0.1892506))

# poisson regression example





