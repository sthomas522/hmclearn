test_that("paper examples", {
  # Linear regression example
  y <- warpbreaks$breaks
  X <- model.matrix(breaks ~ wool*tension, data=warpbreaks)
  N <- 2e3
  set.seed(143)

  eps_vals <- c(rep(2e-1, 6), 2e-2)

  t1 <- Sys.time()
  fm1_hmc <- hmc(N, theta.init = c(rep(0, 6), 1),
                 epsilon = eps_vals, L = 20,
                 logPOSTERIOR = linear_posterior,
                 glogPOSTERIOR = g_linear_posterior,
                 varnames = c(colnames(X), "log_sigma_sq"),
                 param=list(y=y, X=X), chains=2,
                 parallel=FALSE)
  t2 <- Sys.time()

  # test values
  c1_hmc <- as.vector(round(coef(fm1_hmc), 6))

  test1 <- c(42.950208, -13.744052, -18.029510, -17.597738,
             17.869968, 7.675399, 4.800039)

  expect_equal(c1_hmc, test1)

  p1_hmc <- as.vector(round(psrf(fm1_hmc), 6))

  test1b <- c(0.999787, 1.000431, 1.008846, 1.005493, 1.008580, 0.999912,
              0.999788)

  expect_equal(p1_hmc, test1b)

  #######################################
  # logistic regression
  birthwt2 <- MASS::birthwt

  # label race variable
  birthwt2$race2 <- factor(birthwt2$race, labels = c("white", "black", "other"))

  # reduce to indicator variable for positive number of premature labors
  birthwt2$ptd <- ifelse(birthwt2$ptl > 0, 1, 0)

  # reduce to three levels
  birthwt2$ftv2 <- factor(ifelse(birthwt2$ftv > 2, 2, birthwt2$ftv),
                          labels = c("0", "1", "2+"))

  # create design matrix
  X <- model.matrix(low ~ age + lwt + race2 + smoke + ptd + ht + ui + ftv2,
                    data = birthwt2)
  y <- birthwt2$low

  N <- 2e3

  # use different epsilon values for continuous and dichotomous variables
  continuous_ind <- c(FALSE, TRUE, TRUE, rep(FALSE, 8))
  eps_vals <- ifelse(continuous_ind, 1e-3, 5e-2)

  t1 <- Sys.time()
  set.seed(143)
  fm2_hmc <- hmc(N, theta.init = rep(0, 11),
                 epsilon = eps_vals, L = 10,
                 logPOSTERIOR = logistic_posterior,
                 glogPOSTERIOR = g_logistic_posterior,
                 param=list(y=y, X=X),
                 varnames = colnames(X),
                 chains=2, parallel=FALSE)
  t2 <- Sys.time()

  # test values
  c2_hmc <- as.vector(round(coef(fm2_hmc), 6))

  test2 <- c(0.972062, -0.040255, -0.016622, 1.194040, 0.730470,
             0.748092, 1.435274, 2.021169, 0.674046, -0.483113,
             0.140353)

  expect_equal(c2_hmc, test2)

  p2_hmc <- as.vector(round(psrf(fm2_hmc), 6))

  test2b <- c(1.004376, 1.003049, 1.005980, 1.002901, 1.002973,
              1.000847, 1.000158, 1.010873, 0.999796, 1.004207,
              1.001912)

  expect_equal(p2_hmc, test2b)

  ############################################
  # Poisson mixed effects model
  #

  library(lme4)

  # Gopher data:  Gdat dataframe
  data(Gdat)

  ##########
  # block diagonal
  Zi.lst <- split(rep(1, nrow(Gdat)), Gdat$Site)
  Zi.lst <- lapply(Zi.lst, as.matrix)
  Z <- Matrix::bdiag(Zi.lst)
  Z <- as.matrix(Z)
  X <- model.matrix(~ factor(year), data=Gdat)
  X <- cbind(X, Gdat$prev)
  colnames(X)[ncol(X)] <- "prev"
  colnames(X) <- make.names(colnames(X))
  colnames(X)[1] <- "intercept"
  y <- Gdat$shells
  p <- ncol(X)

  N <- 2e3

  set.seed(412)
  initvals <- c(rep(0, 4),
                rep(0, 10),
                0)

  # eps_vals <- c(3e-2, 3e-2, 3e-2, 1e-3, rep(1e-1, 10), 5e-2)
  eps_vals <- c(3e-2, 3e-2, 3e-2, 1e-3, rep(1e-1, 10), 3e-2)

  t1.hmc <- Sys.time()
  fm3_hmc <- hmc(N = N, theta.init = initvals, epsilon = eps_vals, L = 10,
                 logPOSTERIOR  = glmm_poisson_posterior,
                 glogPOSTERIOR  = g_glmm_poisson_posterior,
                 varnames=c(colnames(X), paste0("u", 1:ncol(Z)), "lambda"),
                 param=list(y = y, X=X, Z=Z, m=10, nulambda=1, Alambda=25),
                 chains=2, parallel=FALSE)
  t2.hmc <- Sys.time()


  # test values
  c3_hmc <- as.vector(round(coef(fm3_hmc), 6))

  test3 <- c(-0.130375, -0.656025, -0.358779, 0.023281, -0.792656, -0.146328, -0.558491,
             0.747858, -0.089237, 1.151311, 0.239391, -0.186775, 0.975497, -1.040002,
             -0.355951)
  expect_equal(c3_hmc, test3)

  p3_hmc <- as.vector(round(psrf(fm3_hmc), 6))

  test3b <- c(0.999790, 1.004770, 1.009279, 1.002720,
              0.999764, 0.999925, 1.000083, 1.007014,
              1.001045, 1.000027, 1.000493, 1.002373,
              0.999822, 0.999778, 1.002191)

  expect_equal(p3_hmc, test3b)

})
