test_that("hmc testing", {
  # Linear regression example
  set.seed(521)
  X <- cbind(1, matrix(rnorm(300), ncol=3))
  betavals <- c(0.5, -1, 2, -3)
  y <- X%*%betavals + rnorm(100, sd=.2)

  f1 <- hmc(N = 500,
            theta.init = c(rep(0, 4), 1),
            epsilon = 0.01,
            L = 10,
            logPOSTERIOR = linear_posterior,
            glogPOSTERIOR = g_linear_posterior,
            varnames = c(paste0("beta", 0:3), "log_sigma_sq"),
            param=list(y=y, X=X), parallel=FALSE, chains=1)

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
            parallel=FALSE, chains=1)

  medparam2 <- as.vector(summary(f2, burnin=100)[, 3])
  medparam2

  expect_equal(round(medparam2, 7),
               c(-0.8497155, 0.1892506))

  # poisson regression example
  set.seed(7363)
  X <- cbind(1, matrix(rnorm(40), ncol=2))
  betavals <- c(0.8, -0.5, 1.1)
  lmu <- X %*% betavals
  y <- sapply(exp(lmu), FUN = rpois, n=1)

  f3 <- hmc(N = 500,
            theta.init = rep(0, 3),
            epsilon = 0.01,
            L = 10,
            logPOSTERIOR = poisson_posterior,
            glogPOSTERIOR = g_poisson_posterior,
            varnames = paste0("beta", 0:2),
            param = list(y=y, X=X),
            parallel=FALSE, chains=1)

  medparam3 <- as.vector(summary(f3, burnin=100)[, 3])
  medparam3

  expect_equal(round(medparam3, 7),
               c(0.7770322, -0.5192680, 1.2087953))

})

