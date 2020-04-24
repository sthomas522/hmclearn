test_that("mh testing", {
  # Linear regression example
  set.seed(521)
  X <- cbind(1, matrix(rnorm(300), ncol=3))
  betavals <- c(0.5, -1, 2, -3)
  y <- X%*%betavals + rnorm(100, sd=.2)

  f1 <- mh(N = 5e3,
           theta.init = c(rep(0, 4), 1),
           nu <- c(rep(0.001, 4), 0.1),
           qPROP = qprop,
           qFUN = qfun,
           logPOSTERIOR = linear_posterior,
           varnames = c(paste0("beta", 0:3), "log_sigma_sq"),
           param=list(y=y, X=X), parallel=FALSE, chains=2)

  medparam1 <- as.vector(summary(f1, burnin=1000)[, 3])

  expect_equal(round(medparam1, 6),
               c(0.535630,
                 -1.008264,
                 2.016371,
                 -2.981911,
                 -3.057343))


  # Logistic regression example
  X <- cbind(1, seq(-100, 100, by=0.25))
  betavals <- c(-0.9, 0.2)
  lodds <- X %*% betavals
  prob1 <- as.numeric(1 / (1 + exp(-lodds)))

  set.seed(9874)
  y <- sapply(prob1, function(xx) {
    sample(c(0, 1), 1, prob=c(1-xx, xx))
  })

  f2 <- mh(N = 2000,
            theta.init = rep(0, 2),
            nu = c(0.03, 0.001),
            qPROP = qprop,
            qFUN = qfun,
            logPOSTERIOR = logistic_posterior,
            varnames = paste0("beta", 0:1),
            param = list(y=y, X=X),
            parallel=FALSE, chains=2)

  medparam2 <- as.vector(summary(f2, burnin=200)[, 3])

  expect_equal(round(medparam2, 6),
               c(-0.872098, 0.188425))

  # poisson regression example
  set.seed(7363)
  X <- cbind(1, matrix(rnorm(40), ncol=2))
  betavals <- c(0.8, -0.5, 1.1)
  lmu <- X %*% betavals
  y <- sapply(exp(lmu), FUN = rpois, n=1)

  f3 <- mh(N = 5e3,
            theta.init = rep(0, 3),
            nu = rep(0.01, 3),
            qPROP = qprop,
            qFUN = qfun,
            logPOSTERIOR = poisson_posterior,
            varnames = paste0("beta", 0:2),
            param = list(y=y, X=X),
            parallel=FALSE, chains=2)

  medparam3 <- as.vector(summary(f3, burnin=1000)[, 3])

  expect_equal(round(medparam3, 6),
               c(0.824951, -0.515598, 1.183672))

})

