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

  # test values
  c1_hmc <- as.vector(round(coef(fm1_hmc), 6))

  test1 <- c(43.013150, -14.012657, -18.327812, -18.018012,
             18.005013, 7.972607, 4.798094)

  expect_equal(c1_hmc, test1)

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

  N <- 1e4

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
                 varnames = colnames(X))
  t2 <- Sys.time()

  # test values
  c2_hmc <- as.vector(round(coef(fm2_hmc), 6))

  test2 <- c(0.759534, -0.034473, -0.016576, 1.238999, 0.785729,
             0.777453, 1.421063, 2.027768, 0.709706, -0.480946,
             0.151321)

  expect_equal(c2_hmc, test2)

  ############################################
  # Poisson mixed effects model
  #

  library(lme4)

  # Gopher data:  Gdat dataframe
  data(gopherdat2)

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

  # shorten simulation for testing
  # N <- 1e5
  N <- 1e3

  set.seed(412)
  initvals <- c(rep(0, 4),
                rnorm(10, mean=0, sd=1e-3),
                0)

  M_vals <- c(1e-3, 1e-3, 1e-3, 1,
              rep(1e-3, 10),
              1e-3)


  eps_vals <- c(rep(1e-3, 3), 1e-4, rep(2e-4, 11))

  t1.hmc <- Sys.time()
  res <- hmc(N = N, theta.init = initvals, epsilon = eps_vals, L = 10,
             logPOSTERIOR  = glmm_poisson_posterior,
             Mdiag = M_vals,
             varnames=c(colnames(X), paste0("u", 1:ncol(Z)), "lambda"),
             glogPOSTERIOR  = g_glmm_poisson_posterior,
             param=list(y = y, X=X, Z=Z, m=10, nulambda=1, Alambda=25))
  t2.hmc <- Sys.time()
  # t2.hmc - t1.hmc
  # res$accept/N

  # test values
  c3_hmc <- as.vector(round(coef(res), 6))

  test3 <- c(0.195788, -0.660140, -0.319137, 0.015836, -1.135854,
    -0.195184, -0.551731, 0.574802, -0.163860, 0.548948,
    0.173786, -0.202746, 0.444318, -0.699641, 0.267427)

  expect_equal(c3_hmc, test3)


})
