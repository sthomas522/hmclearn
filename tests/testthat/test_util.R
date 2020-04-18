test_that("utility and diagnostic functions", {
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
            param=list(y=y, X=X), parallel=FALSE, chains=2)

  # psrf
  p1 <- as.vector(round(psrf(f1), 6))

  test1 <- c(1.008833, 1.003896, 1.001133, 1.000041, 1.000087)

  expect_equal(p1, test1)

  #neff

  p2 <- neff(f1)

  expect_equal(c(7, rep(6, 4)), p2)

})
