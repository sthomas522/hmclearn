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
  N <- 1e2

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

  test1 <- c(5.390459, 2.032531, -0.212157, -0.003633, 17.819437)

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
  N <- 1e2

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

  test2 <- c(0.786272, -0.338544, -2.153504, 3.991616)

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
  N <- 1e2

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

  test3 <- c(4.602609, 1.264835, -0.357982, -3.021825,
             0.960862, 1.867834, 1.493367)

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
  N <- 20

  set.seed(412)
  initvals<- c(rep(0, 5), # fixed effects
               rep(0, 60), # random intercepts
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

  test4 <- c(0.54491, 0.129053, -0.495055, -0.341037, 0.30362,
            1.385695, -0.169534, -1.015453, -0.841186, -0.006765,
            0.48774, -0.050163, -0.173102, 0.621903, 1.629583,
            1.170337, -0.027646, 0.267821, -1.028982, 0.70248,
            -0.903477, -0.026878, 0.052465, -0.336058, -0.246761,
            -0.693597, 0.994767, -0.014793, 1.063293, -0.407184,
            -0.088851, 0.707301, 0.619891, 0.072133, -0.243582,
            -0.248633, 0.777454, 0.132293, -0.908939, -0.272922,
            0.706662, -0.301936, 0.249533, -0.822124, 0.157421,
            -0.372467, -0.515733, -1.057183, 1.093029, 0.003732,
            -1.220364, 0.023378, -1.236613, -0.082753, -1.171446,
            0.563121, 0.086442, 0.655766, 0.193117, -1.133754,
            1.18816, -0.464122, 0.344491, 1.034123, 1.199272,
            -0.623314)

  # test4 <- c(0.505508, 0.085253, -0.517873, -0.298942, 0.305182,
  #            1.574587, 0.214333, -0.083195, -0.405973, -0.063414,
  #            0.573493, 0.517492, -0.118635, 0.200475, 0.780959,
  #            1.778907, 0.317892, -0.08589, -1.221131, 0.36502,
  #            -1.135727, 0.385947, 0.179772, -0.050349, -0.277637,
  #            -0.060179, 0.994767, 0.651959, 1.225312, -0.177366,
  #            0.218792, 1.12392, 0.725356, 0.437492, -0.777459,
  #            -0.358952, 1.021771, 0.448855, -1.354307, -0.492309,
  #            0.47946, -0.313183, 0.716556, -0.889983, -0.194919,
  #            -0.338574, -0.330158, -0.848203, 0.748297, 0.283672,
  #            -1.162581, -0.271286, -0.710266, 0.551858, -0.425433,
  #            -0.058832, -0.374909, 0.397397, 1.005636, -1.078609,
  #            1.226558, -0.362492, 0.554207, 0.848197, 1.058174,
  #            -0.690198)

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
  N <- 20

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

  test5 <- c(216.014829, 12.53231, 0.174431, -0.042903,
             0.512965, 0.051781, -0.641552, 0.854632,
             -0.131122, -0.582717, 1.928844, 0.210722,
             1.409785, -0.177323, -0.423615, 0.198745,
             0.309729, -0.569999, 0.340722, 1.320643,
             9.570696, -0.118472)

  # test5 <- c(251.357622, 10.338225, 1.126278, -2.165279, -1.726435,
  #            0.168601, 0.312531, 0.224049, 0.465446, -0.194339,
  #            -1.210437, 1.986091, -0.527670, 0.404418, -0.217797,
  #            0.999496, 0.168943, -0.178436, -0.081471, 0.612903,
  #            3.443459, 3.489608)

  expect_equal(c5_hmc, test5)
})
