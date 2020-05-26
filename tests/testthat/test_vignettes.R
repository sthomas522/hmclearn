test_that("vignette examples", {

  ################################################
  # Linear regression
  ################################################

  ScotsRaces <- MASS::hills

  X <- cbind(1, as.matrix(ScotsRaces[, -which(colnames(ScotsRaces) == "time")]))

  # add interaction
  X <- cbind(X, X[, "dist"] * X[, "climb"])
  colnames(X)[ncol(X)] <- "climb_distance"
  colnames(X)[1] <- "intercept"

  y <- ScotsRaces$time

  # reduce from 1e4 for testing
  N <- 1e3

  eps_vals <- c(1e-1, 1e-2, 1e-4, 5e-6, 3e-3)

  set.seed(412)
  f_hmc <- hmc(N=N, theta.init = rep(0, 5),
               epsilon = eps_vals, L = 20,
               logPOSTERIOR = linear_posterior,
               glogPOSTERIOR = g_linear_posterior,
               param=list(y = y, X=X),
               varnames = c(colnames(X), "log_sigma2"),
               parallel=FALSE, chains=1)

  # test values
  c1_hmc <- as.vector(round(coef(f_hmc), 6))

  test1 <- c(10.882923, 4.636202, -0.105222, -0.000477, 12.031474)

  expect_equal(c1_hmc, test1)

  ################################################
  # Logistic regression
  ################################################

  data(Endometrial)

  # data prep
  Endometrial$PI2 <- with(Endometrial, (PI - mean(PI)) / sd(PI))
  Endometrial$EH2 <- with(Endometrial, (EH - mean(EH)) / sd(EH))
  Endometrial$NV2 <- Endometrial$NV - 0.5

  X <- cbind(1, as.matrix(Endometrial[, which(colnames(Endometrial) %in% c("PI2", "EH2", "NV2"))]))
  y <- Endometrial$HG

  colnames(X) <- c("(Intercept)", "PI2", "EH2", "NV2")

  # reduce from 1e4 for testing
  N <- 1e3

  set.seed(412)
  f_hmc <- hmc(N = N, theta.init = rep(1, 4),
               epsilon = 0.1, L = 20,
               logPOSTERIOR = logistic_posterior,
               glogPOSTERIOR = g_logistic_posterior,
               varnames = colnames(X),
               param=list(y = y, X=X, sig2beta=10),
               parallel=FALSE, chains=1)

  # test values
  c2_hmc <- as.vector(round(coef(f_hmc), 6))

  test2 <- c(0.630469, -0.409110, -2.104589, 3.817009)

  expect_equal(c2_hmc, test2)

  ################################################
  # Poisson regression
  ################################################

  data(Drugs)

  # design matrix
  X <- model.matrix(count ~ A + C + M + A:C + A:M + C:M , data=Drugs)
  X <- X[, 1:ncol(X)]

  # independent variable is count data
  y <- Drugs$count

  # reduce from 1e4 for testing
  N <- 1e3

  eps_vals <- c(rep(5e-4, 2), 1e-3, 2e-3, 1e-3, 2e-3, 5e-4)

  set.seed(412)
  f_hmc <- hmc(N = N, theta.init = rep(0, 7),
               epsilon = eps_vals, L = 50,
               logPOSTERIOR = poisson_posterior,
               glogPOSTERIOR = g_poisson_posterior,
               varnames = colnames(X),
               parallel=FALSE, chains=1,
               param=list(y=y, X=X))

  # test values
  c3_hmc <- as.vector(round(coef(f_hmc), 6))

  test3 <- c(5.626569, 0.506499, -1.836742, -5.566284,
             2.011625, 3.194251, 2.745701)

  expect_equal(c3_hmc, test3)

  ################################################
  # Logistic mixed effects regression
  ################################################

  Contraception <- mlmRev::Contraception

  Contraception$liv2 <- ifelse(Contraception$livch == "0", 0, 1)

  ##########
  # block diagonal
  Zi.lst <- split(rep(1, nrow(Contraception)), Contraception$district)
  Zi.lst <- lapply(Zi.lst, as.matrix)
  Z <- Matrix::bdiag(Zi.lst)
  Z <- as.matrix(Z)

  urban <- ifelse(Contraception$urban == "Y", 1, 0)

  X <- cbind(1, Contraception$age, Contraception$age^2, urban, Contraception$liv2)
  colnames(X) <- c("int", "age", "age_sq", "urban", "liv2")
  y <- ifelse(Contraception$use == "Y", 1, 0)

  xqr <- qr(X)
  Q <- qr.Q(xqr)
  R <- qr.R(xqr)

  n <- nrow(X)
  X2 <- Q * sqrt(n-1)
  Rstar <- R / sqrt(n-1)
  Rstar_inv <- solve(Rstar)
  colnames(X2) <- c("int", "age", "age_sq", "urban", "liv2")

  # new intercept from QR decomposition
  diagval <- X2[1,1]

  ##########
  # block diagonal
  Zi.lst <- split(rep(diagval, nrow(Contraception)), Contraception$district)
  Zi.lst <- lapply(Zi.lst, as.matrix)
  Z2 <- Matrix::bdiag(Zi.lst)
  Z2 <- as.matrix(Z2)

  # reduce from 2e3 for testing
  N <- 5e2

  set.seed(412)
  initvals<- c(rep(0, 5), # fixed effects
               rnorm(60, mean=0, sd=0.1), # random intercepts
               0) # variance of random intercepts


  vnames <- c(colnames(X),
              paste0("tau_int", 1:60),
              "xi1")

  epsvals <- c(5e-2, rep(1e-2, 4), rep(5e-2, 61))

  set.seed(412)
  f_hmc <- hmc(N = N, theta.init = initvals,
               epsilon = epsvals, L = 10,
               logPOSTERIOR = glmm_bin_posterior,
               glogPOSTERIOR = g_glmm_bin_posterior,
               varnames = vnames,
               parallel=FALSE, chains=1,
               param=list(y = y, X=X2, Z=Z2, n=60, nrandom=1, sig2beta=5,
                          nuxi=1, Axi=25)  )

  # test values
  c4_hmc <- as.vector(round(coef(f_hmc), 6))

  test4 <- c(0.505508, 0.085187, -0.517873, -0.298942, 0.305182,
             1.574743, 0.214264, -0.083402, -0.405973, -0.063414,
             0.573493, 0.517492, -0.118635, 0.200436, 0.780959,
             1.776997, 0.317879, -0.080616, -1.221131, 0.36502,
             -1.135727, 0.385947, 0.179772, -0.05121, -0.277666,
             -0.060179, 0.976638, 0.651959, 1.225312, -0.174169,
             0.217565, 1.12392, 0.725287, 0.437492, -0.777459,
             -0.358952, 1.021771, 0.443896, -1.354307, -0.492309,
             0.479458, -0.313183, 0.716553, -0.889986, -0.194919,
             -0.338574, -0.330028, -0.848203, 0.748297, 0.283666,
             -1.162581, -0.271286, -0.710266, 0.551858, -0.425514,
             -0.058836, -0.374909, 0.397397, 1.003512, -1.079143,
             1.226214, -0.362492, 0.554207, 0.848197, 1.058152,
             -0.690198)

  expect_equal(c4_hmc, test4)

  ################################################
  # Linear mixed effects regression
  ################################################

  library(lme4)
  library(Matrix)
  data(sleepstudy)

  # dependent variable
  y <- sleepstudy$Reaction

  yi.lst <- split(sleepstudy$Reaction, sleepstudy$Subject)

  # fixed effects
  ss2 <- sleepstudy
  ss2$int <- 1
  ss2 <- ss2[, c(4, 1:3)] # rearrange columns to store in list
  Xi.lst <- split(ss2[, which(colnames(ss2) %in% c("Days", "int"))],
                  ss2$Subject)
  Xi.lst <- lapply(Xi.lst, as.matrix)

  X <- as.matrix(do.call(rbind, Xi.lst))

  # random effects
  n <- length(unique(sleepstudy$Subject))
  d <- length(y)/n
  nrandom <- 1

  ##########
  # intercept
  Z <- kronecker(diag(n), matrix(rep(1, d), ncol=1))
  # Zi.lst <- replicate(m, matrix(rep(1, n/m), ncol=1), simplify=FALSE)
  # Z <- as.matrix(bdiag(Zi.lst))

  # reduce from 2e3 for testing
  N <- 5e2

  theta.init <- c(0, 1, # beta
                  rep(0, 18), # tau
                  3, # gamma (log sig2eps)
                  1) # xi

  vnames <- c(paste0("beta", 0:1),
              paste0("tau_int", 1:18),
              "sigeps", "xi")

  eps_vals <- c(5e-1, 5e-2,
                rep(3e-2, 18),
                6e-3,
                5e-2)

  set.seed(41132)
  f_hmc <- hmc(N = N, theta.init = theta.init,
               epsilon = eps_vals, L = 10,
               logPOSTERIOR = lmm_posterior,
               glogPOSTERIOR = g_lmm_posterior,
               varnames = vnames,
               param=list(y = y, X=X, Z=Z, n=n, d=d, nrandom=1,
                          nugamma=4, Agamma=1,
                          nuxi=1, Axi=1, sig2beta=1e5),
               parallel=FALSE, chains=1)

  # test values
  c5_hmc <- as.vector(round(coef(f_hmc), 6))

  test5 <- c(251.357622, 10.338225, 1.126278, -2.165279, -1.726435,
             0.168601, 0.312531, 0.224049, 0.465446, -0.194339,
             -1.210437, 1.986091, -0.527670, 0.404418, -0.217797,
             0.999496, 0.168943, -0.178436, -0.081471, 0.612903,
             3.443459, 3.489608)

  expect_equal(c5_hmc, test5)
})
