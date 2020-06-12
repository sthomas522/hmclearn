test_that("paper examples", {
  # Linear regression example
  y <- warpbreaks$breaks
  X <- model.matrix(breaks ~ wool*tension, data=warpbreaks)
  # N <- 2e3
  N <- 2e2
  set.seed(143)

  eps_vals <- c(rep(2e-1, 6), 2e-2)

  fm1_hmc <- hmc(N, theta.init = c(rep(0, 6), 1),
                 epsilon = eps_vals, L = 20,
                 logPOSTERIOR = linear_posterior,
                 glogPOSTERIOR = g_linear_posterior,
                 varnames = c(colnames(X), "log_sigma_sq"),
                 param=list(y=y, X=X), chains=2,
                 parallel=FALSE)

  # test values
  c1_hmc <- as.vector(round(coef(fm1_hmc), 6))

  test1 <- c(43.466057, -7.679181, -12.104186, -12.363732,
             15.987310, 4.629290, 4.967581)

  expect_equal(c1_hmc, test1)

  p1_hmc <- as.vector(round(psrf(fm1_hmc), 6))

  test1b <- c(0.999640, 0.997576, 1.004052, 1.008374, 1.033657, 1.003823,
              0.998089)

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

  # N <- 2e3
  N <- 3e2

  # use different epsilon values for continuous and dichotomous variables
  continuous_ind <- c(FALSE, TRUE, TRUE, rep(FALSE, 8))
  eps_vals <- ifelse(continuous_ind, 1e-3, 5e-2)

  set.seed(143)
  fm2_hmc <- hmc(N, theta.init = rep(0, 11),
                 epsilon = eps_vals, L = 10,
                 logPOSTERIOR = logistic_posterior,
                 glogPOSTERIOR = g_logistic_posterior,
                 param=list(y=y, X=X),
                 varnames = colnames(X),
                 chains=2, parallel=FALSE)

  # test values
  c2_hmc <- as.vector(round(coef(fm2_hmc), 6))

  test2 <- c(-0.864417, -0.003891, -0.007844, 0.913143, 0.161931,
             0.270894, 0.720779, 0.486956, 0.200137, -0.575171,
             0.066911)

  expect_equal(c2_hmc, test2)

  p2_hmc <- as.vector(round(psrf(fm2_hmc), 6))

  test2b <- c(1.174999, 1.003712, 1.198086, 0.999307, 1.056735,
              1.276597, 1.119958, 1.062370, 1.302114, 1.026340,
              1.000074)

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

  N <- 1e2


  initvals <- c(rep(0, 4),
                rep(0, 10),
                0)

  # eps_vals <- c(3e-2, 3e-2, 3e-2, 1e-3, rep(1e-1, 10), 5e-2)
  eps_vals <- c(3e-2, 3e-2, 3e-2, 1e-3, rep(1e-1, 10), 3e-2)

  set.seed(412)
  fm3_hmc <- hmc(N = N, theta.init = initvals, epsilon = eps_vals, L = 10,
                 logPOSTERIOR  = glmm_poisson_posterior,
                 glogPOSTERIOR  = g_glmm_poisson_posterior,
                 varnames=c(colnames(X), paste0("u", 1:ncol(Z)), "xi"),
                 param=list(y = y, X=X, Z=Z, n=10, nuxi=1, Axi=25),
                 chains=2, parallel=FALSE)


  # test values
  c3_hmc <- as.vector(round(coef(fm3_hmc), 6))

  test3 <- c(-0.144120, -0.567178, -0.216619, 0.020996, -0.911150,
             -0.305638, -0.659170, 0.693371, -0.096723, 1.135807,
             0.350674, -0.131506, 0.892042, -1.035583, -0.386174)
  expect_equal(c3_hmc, test3)

  p3_hmc <- as.vector(round(psrf(fm3_hmc), 6))

  test3b <- c(1.019176, 0.994988, 1.006717, 1.002968, 0.995020,
              0.996206, 1.014475, 0.999683, 1.002627, 1.001639,
              1.050004, 1.032394, 1.042095, 1.028268, 1.002090)

  expect_equal(p3_hmc, test3b)

})
