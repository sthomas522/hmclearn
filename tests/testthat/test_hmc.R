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

  medparam1 <- as.vector(summary(f1)[, 4])

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

  medparam2 <- as.vector(summary(f2, burnin=100)[, 4])

  expect_equal(round(medparam2, 7),
               c(-0.8504805, 0.1892643))

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

  medparam3 <- as.vector(summary(f3, burnin=100)[, 4])

  expect_equal(round(medparam3, 7),
               c(0.7756224, -0.5167261, 1.2088527))

  # linear regression
  data(warpbreaks)

  X <- model.matrix(breaks ~ wool*tension, data=warpbreaks)
  y <- warpbreaks$breaks

  N <- 1e3
  eps_vals <- c(rep(2e-1, 6), 2e-2)

  set.seed(143)
  f4_hmc <- hmc(N, theta.init = c(rep(0, 6), 1),
                 epsilon = eps_vals, L = 20,
                 logPOSTERIOR = linear_posterior,
                 glogPOSTERIOR = g_linear_posterior,
                 varnames = c(colnames(X), "log_sigma_sq"),
                 param=list(y=y, X=X), parallel = FALSE, chains = 2)

  # test values
  c4_hmc <- as.vector(round(coef(f4_hmc), 6))

  test4 <- c(42.579371, -12.955274, -17.135850, -17.007076, 17.180934,
             7.381530, 4.814274)

  expect_equal(c4_hmc, test4)

  p4_hmc <- as.vector(round(psrf(f4_hmc), 6))

  test4b <- c(1.000371, 1.000456, 0.999587,
              1.002521, 1.007281, 1.000957, 0.999622)

  expect_equal(p4_hmc, test4b)

  # logistic regression
  data("PimaIndiansDiabetes", package = "mlbench")

  y <- ifelse(PimaIndiansDiabetes$diabetes == 'pos', 1, 0)
  X <- cbind(1, as.matrix(PimaIndiansDiabetes[, -which(colnames(PimaIndiansDiabetes) == "diabetes")]))
  colnames(X)[1] <- "int"

  N <- 1e3
  eps_vals <- c(5e-2, 2e-3, 2e-4, 1e-3, 1e-3,
                1e-4, 1e-3, 3e-2, 4e-4)

  set.seed(412)
  f5_hmc <- hmc(N = N, theta.init = rep(0, 9),
                 epsilon = eps_vals, L = 10,
                 logPOSTERIOR = logistic_posterior,
                 glogPOSTERIOR = g_logistic_posterior,
                 param=list(y = y, X=X),
                 parallel=FALSE, chains=2)

  # test values
  c5_hmc <- as.vector(round(coef(f5_hmc), 6))

  test5 <- c(-8.486601, 0.127246, 0.035835, -0.013833,
             0.000857, -0.001223, 0.089886, 0.950046, 0.013887)

  expect_equal(c5_hmc, test5)

  p5_hmc <- as.vector(round(psrf(f5_hmc), 6))

  test5b <- c(1.011002, 1.002309, 1.000229, 0.999595,
              1.002694, 1.007651, 1.003839, 0.999561, 1.017046)

  expect_equal(p5_hmc, test5b)

  # poisson regression
  data("AMSsurvey", package="carData")

  # design matrix
  X <- model.matrix(count ~ type + sex + citizen, data=AMSsurvey)

  # independent variable is count data
  y <- AMSsurvey$count
  p <- ncol(X)

  N <- 2e3

  eps_vals <- c(2.2e-3, 2e-3, 2e-3, 2e-3, 2e-3, 3e-3,
                2e-3, 2e-3)

  f6_hmc <- hmc(N, theta.init = rep(0, p),
                 epsilon = eps_vals, L = 20,
                 logPOSTERIOR = poisson_posterior,
                 glogPOSTERIOR=g_poisson_posterior,
                 varnames = colnames(X),
                 param=list(y = y, X=X), parallel=FALSE, chains=2)

  # test values
  c6_hmc <- as.vector(round(coef(f6_hmc), 6))

  test6 <- c(3.564045, 0.458186, 0.311374, -0.188827,
             0.535655, -0.850916, 0.743093, -0.128503)

  expect_equal(c6_hmc, test6)

  p6_hmc <- as.vector(round(psrf(f6_hmc), 6))

  test6b <- c(1.008315, 1.011858, 1.021309, 1.005917,
              1.017213, 1.003739, 1.013696, 1.000277)

  expect_equal(p6_hmc, test6b)

})



