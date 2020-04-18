test_that("paper examples", {
  # Linear regression example
  y <- warpbreaks$breaks
  X <- model.matrix(breaks ~ wool*tension, data=warpbreaks)
  N <- 1e4
  set.seed(143)

  eps_vals <- c(rep(2e-1, 6), 2e-2)

  t1 <- Sys.time()
  fm1_hmc <- hmc(N, theta.init = c(rep(0, 6), 1),
                 epsilon = eps_vals, L = 20,
                 logPOSTERIOR = linear_posterior,
                 glogPOSTERIOR = g_linear_posterior,
                 varnames = c(colnames(X), "log_sigma_sq"),
                 param=list(y=y, X=X))
  t2 <- Sys.time()

  c1_hmc <- as.vector(round(coef(fm1_hmc), 6))

  test1 <- c(43.013150, -14.012657, -18.327812, -18.018012,
             18.005013, 7.972607, 4.798094)

  expect_equal(c1_hmc, test1)

})
