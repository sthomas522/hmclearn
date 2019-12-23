# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}

# test all glm families
if (1 == 0) {
  ###################################################################
  # Linear regression
  ###################################################################
  X <- model.matrix(breaks ~ wool*tension, data=warpbreaks)
  y <- warpbreaks$breaks

  N <- 10000
  set.seed(143)

  # eps_vals <- c(rep(2e-1, 6), 2e-2)
  eps_vals <- c(rep(2e-2, 6), 2e-3)
  # eps_vals <- c(rep(2.8e-1, 6), 2e-2)

  t1 <- Sys.time()
  set.seed(321)
  fm1_hmc <- hmc(N, theta.init = c(rep(0, 6), 1), epsilon = eps_vals, L = 20,
                 logPOSTERIOR = linear_posterior,
                 glogPOSTERIOR = g_linear_posterior,
                 varnames = c(colnames(X), "log_sigma_sq"), y=y, X=X)
  t2 <- Sys.time()
  t2 - t1

  fm1_hmc$accept / N
  summary(fm1_hmc)

  fm1_pred <- predict(fm1_hmc, X=X)

  ###################################################################
  # logistic regression
  ###################################################################

  library(mlbench)
  data(BreastCancer)

  bc <- BreastCancer[complete.cases(BreastCancer), ]

  X <- model.matrix(Class ~ Cl.thickness + Cell.size + Cell.shape,
                 #     Marg.adhesion + Epith.c.size + Bare.nuclei +
                #      Bl.cromatin + Normal.nucleoli + Mitoses,
                    data = bc)
  y <- ifelse(bc$Class == "benign", 0, 1)

  p <- ncol(X)

  N <- 10000

  t1 <- Sys.time()
  set.seed(321)
  fm2_hmc <- hmc(N, theta.init = rep(0, p),
                 epsilon = 1e-1, L=20,
                 logPOSTERIOR = logistic_posterior,
                 glogPOSTERIOR = g_logistic_posterior,
                 randlength = TRUE,
                 varnames = colnames(X), y=y, X=X)
  t2 <- Sys.time()

  fm2_hmc$accept / N

  fm2_pred <- predict(fm2_hmc, X=X, fam="binomial")

  ###################################################################
  # poisson regression
  ###################################################################

  library(carData)
  data(AMSsurvey)

  # design matrix
  X <- model.matrix(count ~ type + sex + citizen, data=AMSsurvey)

  # independent variable is count data
  y <- AMSsurvey$count
  p <- ncol(X)

  N <- 10000

  fm3_hmc <- hmc(N, theta.init = rep(0, p), epsilon = 2e-3, L = 20,
                     logPOSTERIOR = poisson_posterior,
                 glogPOSTERIOR=g_poisson_posterior,
                     y = y, X=X)
  fm3_hmc$accept / N

  fm3_pred <- predict(fm3_hmc, X=X, fam="poisson")

  ###################################################################
  # Linear Mixed effects model
  ###################################################################

  N <- 1e5

  library(lme4)
  # dependent variable
  y <- sleepstudy$Reaction

  # fixed effects
  ss2 <- sleepstudy
  ss2$int <- 1
  ss2 <- ss2[, c(4, 1:3)] # rearrange columns to store in list
  Xi.lst <- split(ss2[, which(colnames(ss2) %in% c("Days", "int"))],
                  ss2$Subject)
  Xi.lst <- lapply(Xi.lst, as.matrix)

  X <- as.matrix(do.call(rbind, Xi.lst))

  # random effects
  m <- length(unique(sleepstudy$Subject))
  d <- length(unique(sleepstudy$Days))

  ##########
  # intercept and slope
  Zi <- cbind(1, sort(unique(sleepstudy$Days)))
  q <- ncol(Zi)
  Zi.lst <- replicate(m, cbind(1, sort(unique(sleepstudy$Days))), simplify=FALSE)
  Z <- bdiag(Zi.lst)
  Z <- as.matrix(Z)


  thetaInit <-c(0,0, # beta
             rnorm(36), # tau
             # rnorm(36, mean=0, sd=sqrt(10)),
             5, # gamma (log sig2eps)
             c(2, 2, 0)) # xi and a (log G diagonal and a off-diagonal)


  # 2min
  t1 <- Sys.time()
  fm4_hmc <- hmc(N = N, theta.init = thetaInit,
                 epsilon = 3e-3, L = 20,
             logPOSTERIOR = lmm_posterior,
             glogPOSTERIOR = g_lmm_posterior,
            varnames = c("beta1", "beta2",
                         paste0("u", 1:36), "gamma", "xi", "a1", "a2"),
             y = y, X=X, Z=Z, m=m, q=q,
             A = 1e4, nueps=1, nulambda=1, Aeps=25, Alambda=25)
  t2 <- Sys.time()
  t2 - t1

  fm4_hmc$accept / N


  ###################################################################
  # Logistic Mixed effects model
  ###################################################################
  library(mlmRev)

  data(Contraception)
  Contraception$liv2 <- ifelse(Contraception$livch == "0", 0, 1)

  ##########
  # block diagonal
  Zi.lst <- split(rep(1, nrow(Contraception)), Contraception$district)
  Zi.lst <- lapply(Zi.lst, as.matrix)
  Z <- bdiag(Zi.lst)
  Z <- as.matrix(Z)

  urban <- ifelse(Contraception$urban == "Y", 1, 0)

  X <- cbind(1, Contraception$age, Contraception$age^2, urban, Contraception$liv2)
  colnames(X) <- c("int", "age", "age_sq", "urban", "liv2")
  y <- ifelse(Contraception$use == "Y", 1, 0)

  # qr decomposition
  xqr <- qr(X)
  Q <- qr.Q(xqr)
  R <- qr.R(xqr)

  n <- nrow(X)
  X2 <- Q * sqrt(n-1)
  Rstar <- R / sqrt(n-1)
  Rstar_inv <- solve(Rstar)
  colnames(X2) <- c("int", "age", "age_sq", "urban", "liv2")


  N <- 1e4

  set.seed(412)
  thetaInit <- c(rep(0, 5), # fixed effects
                 rnorm(60, mean=0, sd=0.1), # random intercepts
                 0) # variance of random intercepts

  fm5_hmc <- hmc(N = N, theta.init = thetaInit,
             epsilon = 5e-3, L = 10,
             logPOSTERIOR = glmm_bin_posterior,
             glogPOSTERIOR = g_glmm_bin_posterior,
             y = y, X=X2, Z=Z, m=60, q=1, B=5, nuxi=1, Axi=25)

  ###################################################################
  # Poisson Mixed effects model
  ###################################################################

  # https://ms.mcmaster.ca/~bolker/R/misc/foxchapter/bolker_chap.html

  data(gopherdat2)

  ##########
  # block diagonal
  Zi.lst <- split(rep(1, nrow(Gdat)), Gdat$Site)
  Zi.lst <- lapply(Zi.lst, as.matrix)
  Z <- bdiag(Zi.lst)
  Z <- as.matrix(Z)

  X <- model.matrix(~ factor(year), data=Gdat)
  X <- cbind(X, Gdat$prev)
  colnames(X)[ncol(X)] <- "prev"
  colnames(X) <- make.names(colnames(X))
  colnames(X)[1] <- "intercept"

  y <- Gdat$shells
  p <- ncol(X)

  Mruns <- 5e5

  set.seed(412)
  thetaInit <- c(rep(0, 4),
                             rnorm(10, mean=0, sd=1e-3),
                             0)


  M_vals <- c(1e-3, 1e-3, 1e-3, 1,
              rep(1e-3, 10),
              1e-3)

  t1.hmc <- Sys.time()

  fm6_hmc <- hmc(N = Mruns, theta.init = thetaInit,
                 epsilon = 1e-4, L = 10,
             logPOSTERIOR = glmm_poisson_posterior,
              Mdiag = M_vals,
             glPOSTERIOR = g_glmm_poisson_posterior,
             y = y, X=X, Z=Z, m=10, q=1, nuxi=1, Axi=25)

  t2.hmc <- Sys.time()
  t2.hmc - t1.hmc
  fm6_hmc$accept/Mruns

  ###################################################################
  # Linear Mixed effects model
  ###################################################################
#
#   library(MCMCglmm)
#   data("BTdata")
#
#   # sort BTdata by random effects level
#   BTdata <- BTdata[order(BTdata$dam, BTdata$sex), ]
#
#   # dependent variable
#   y <- BTdata$tarsus
#   n <- length(y)
#
#   yi.lst <- split(BTdata$tarsus, BTdata$dam)
#
#   levels(BTdata$sex) <- c("UNK", "Fem", "Male")
#
#   X <- model.matrix(~ sex, data=BTdata)
#
#   # random effects
#   m <- length(unique(BTdata$dam))
#
#   ##########
#   # block diagonal
#   Zi.lst <- split(data.frame(X), BTdata$dam)
#   Zi.lst <- lapply(Zi.lst, as.matrix)
#   Z <- bdiag(Zi.lst)
#   Z <- as.matrix(Z)
#
#   Mruns <- 4e4
#   # Mruns <- 100
#   set.seed(41121)
#   thetaInit <- initvals <- c(0,0, 0, # beta
#                              rnorm(3*106, sd=0.1), # tau
#                              # rep(0, 3*106),
#                              # rnorm(36, mean=0, sd=sqrt(10)),
#                              0, # gamma (log sigeps)
#                              c(0, 0, 0),
#                              c(0, 0, 0)) # xi and a (log G diagonal and a off-diagonal)
#
#
#   M_vals <- c(1, 1, 0.1,
#               rep(1, 3*106),
#               1,
#               rep(0.1, 3),
#               rep(0.1, 3))
#
#   set.seed(41132)
#   # tuning: get mass matrix
#   t1.hmc <- Sys.time()
#   # profvis({
#   res <- hmc(N = Mruns, theta.init = thetaInit, epsilon = 3e-3, L = 10,
#              logPOSTERIOR = lmm_posterior,
#              glogPOSTERIOR = g_lmm_posterior,
#              Mdiag = M_vals,
#              y = y, X=X, Z=Z, m=106, q=3, nulambda=4, Alambda=1, A=0.1)
#   # })
#   t2.hmc <- Sys.time()
#   t2.hmc - t1.hmc
#
#   res$accept/Mruns


}




# hmc example
if (1 == 0) {
  library(ggplot2)
  library(gganimate)
  library(hmclearn)
  library(magick)

  fm1 <- lm(breaks ~ wool*tension, data=warpbreaks)
  summary(fm1)
  X <- model.matrix(breaks ~ wool*tension, data=warpbreaks)
  y <- warpbreaks$breaks

  N <- 10000
  set.seed(143)

  # eps_vals <- c(rep(2e-1, 6), 2e-2)
  eps_vals <- c(rep(2e-2, 6), 2e-3)
  # eps_vals <- c(rep(2.8e-1, 6), 2e-2)

  t1 <- Sys.time()
  set.seed(321)
  fm1_hmc <- hmc(N, theta.init = c(rep(0, 6), 1), epsilon = eps_vals, L = 20,
                 logPOSTERIOR = linear_posterior,
                 glogPOSTERIOR = g_linear_posterior,
                 varnames = c(colnames(X), "log_sigma_sq"), y=y, X=X)
  t2 <- Sys.time()
  t2 - t1

  fm1_hmc$accept / N

  theta.all <- as.data.frame(do.call(rbind, fm1_hmc$theta.all))
  r.all <- as.data.frame(do.call(rbind, fm1_hmc$r.all))

  # create a contour plot for V2
  # get median theta vals
  theta.median <- apply(fm1_hmc$thetaDF, 2, median)

  cdata <- expand.grid(V1 = theta.median[1],
                       theta = seq(-30, -5, by=0.1),
                       V3 = theta.median[3],
                       V4 = theta.median[4],
                       V5 = theta.median[5],
                       V6 = theta.median[6],
                       V7 = theta.median[7],
                       p = seq(-4.2, 4.2, by=0.1))

  linpost <- function(theta_p, y, x, ...) {
    k <- length(theta_p)
    theta <- theta_p[1:(k-1)]
    p <- theta_p[k]
    res <- linear_posterior(theta=theta, y=y, X=x, ...) - p^2/2
    res
  }

  cdata$zval <- apply(X=cdata, MARGIN=1, FUN=linpost, x=X, y=y)

  v <- ggplot(cdata, aes(x=p, y=theta, z=zval))
  v <- v + geom_contour()
  v

  Lval <- 20

  # dataframe of accepted proposals
  thetaDF <- fm1_hmc$thetaDF
  accept_v <- fm1_hmc$accept_v

  # vector of woolB (V2) proposals
  V2accept <- thetaDF$woolB
  V2accept_all <- rep(V2accept, each=Lval)
  V2accept_all <- V2accept_all[1:(length(V2accept_all) - Lval)]

  # V2all[Lval*i + 1]
  V2all <- theta.all$V2
  pV2all <- r.all$V2

  burnin_start <- 1000
  basenum <- burnin_start * Lval

  pdata <- NULL
  tempDF2 <- NULL
  kk <- 0

  for (jj in 1:200) {

    if (jj > 1) {
      tempDF2 <- pdata[pdata$keepval, ]
      tempDF2$timeval <- max(tempDF2$timeval)
      tempDF2 <- tempDF2[!duplicated(tempDF2), ]
      tempDF2$accept <- 1
    }

    tempDF <- data.frame(p = pV2all[(basenum + (jj-1)*Lval+1):(basenum + jj*Lval+1)],
                         theta = V2all[(basenum + (jj-1)*Lval+1):(basenum + jj*Lval+1)],
                         theta_accept = V2accept_all[(basenum + (jj-1)*Lval+1):(basenum + jj*Lval+1)],
                         timeval = 1:(Lval+1),
                         keepval = 1:(Lval+1) %% 21 == 0)

    # adjust momentum for visual
    tempDF$p[Lval+1] <- 2*tempDF$p[Lval] - tempDF$p[Lval-1]

    tempDF$accept <- as.integer(tempDF$keepval)

    # inner loop expand pdata
    temp3 <- NULL
    for (kk in 1:nrow(tempDF)) {
      xx1 <- tempDF[1:kk, ]
      xx1$timeval <- kk + (Lval+1)*(jj-1)

      xx2 <- tempDF2
      if (!is.null(tempDF2)) {
        xx2$timeval <- kk + (Lval+1)*(jj-1)
      }
      temp3 <- rbind(temp3, xx2, xx1)
    }

    # tempDF$keepval <- tempDF$keepval + (jj-1)*Lval

    pdata <- rbind(pdata, temp3)
    print(jj)

  }

  row.names(pdata) <- 1:nrow(pdata)
  pdata$accept[pdata$keepval] <- 1
  pdata$accept <- factor(pdata$accept)

  color.codes<-as.character(c("#b2c3a3", "#6f0022"))
  color.names<-c("blue", "red")
  #
  #   p <- ggplot(pdata, aes(x=p, y=theta, colour=accept, shape=accept))
  #   p <- p + geom_point(size=2) + theme_bw()
  #   p <- p + scale_colour_manual(values=setNames(color.codes, c("0", "1")))
  #   p <- p + stat_contour(data=cdata, aes(x=p, y=theta, z=zval))
  #   p <- p + transition_time(timeval)
  #   animate(p, renderer = av_renderer('~/webmfiles/test.webm'), width = 1280,
  #                                     height = 720, res = 104, duration = 120)

  # plot with contour
  v <- ggplot(cdata, aes(x=p, y=theta, z=zval))
  v <- v + geom_contour()
  v <- v + geom_point(data = pdata, size=2, inherit.aes=FALSE,
                      aes(x=p, y=theta, colour=accept, shape=accept))
  v <- v + theme_bw()
  v <- v + scale_colour_manual(values=setNames(color.codes, c("0", "1")))
  v <- v + transition_time(timeval)
  animate(v, renderer = av_renderer('~/webmfiles/test2.webm'), width = 1280,
          height = 720, res = 104, duration = 120)


  # # animation plot
  # basenum <- 1000*Lval
  # pdata <- NULL
  #
  # library(gganimate)
  # for (jj in 1:50) {
  #   tempDF <- data.frame(p = r.all$V2[basenum:(basenum+jj*Lval)],
  #                        theta = theta.all$V2[basenum:(basenum+jj*Lval)],
  #                        timevar =  jj)
  #   pdata <- rbind(pdata, tempDF)
  # }
  #
  #
  # for (k in 1:50) {
  #
  #   x1 <- theta.all$V2[(basenum+1):(basenum+Lval)]
  #   y1 <- r.all$V2[(basenum+1):(basenum+Lval)]
  #   if (k == 1) {
  #     plot(x1, y1, type='o', ylim = c(-2.5, 2.5),
  #          xlim = c(-40, 0))
  #   } else {
  #     lines(x1, y1, type='o')
  #   }
  #
  # }


  # pdata <- data.frame(p = r.all$V2[basenum:(basenum+1000)],
  #                     theta = theta.all$V2[basenum:(basenum+1000)])
  # pdata$timevar <- 1:nrow(pdata)
  #
  # pdata$selected <- pmin(pdata$timevar %% 10, 1)

  p <- ggplot(pdata, aes(x=p, y=theta))
  p <- p + geom_point()
  p <- p + transition_time(timevar)
  # animate(p, renderer = file_renderer(dir='~/animation/',
  #                                     overwrite = T))
  animate(p, renderer = av_renderer())

  library(gapminder)
  p <- ggplot(
    subset(gapminder, year==2007),
    aes(x = gdpPercap, y=lifeExp, size = pop, colour = country)
  ) +
    geom_point(show.legend = FALSE, alpha = 0.7) +
    scale_color_viridis_d() +
    scale_size(range = c(2, 12)) +
    scale_x_log10() +
    labs(x = "GDP per capita", y = "Life expectancy")

  p

  p <- p + transition_time(year) +
    labs(title = "Year: {frame_time}")

  animate(p, renderer = av_renderer())

  ##############################################
  # Trace plot
  ##############################################

  eps_vals <- c(rep(2.8e-1, 6), 2e-2)
  # eps_vals <- c(rep(2e-2, 6), 2e-3)

  t1 <- Sys.time()
  set.seed(321)
  fm2_hmc <- hmc(N, theta.init = c(rep(0, 6), 1), epsilon = eps_vals, L = 20,
                 logPOSTERIOR = linear_posterior,
                 glogPOSTERIOR = g_linear_posterior, y=y, X=X,
                 varnames = c(colnames(X), "log_sigma_sq"))
  t2 <- Sys.time()
  t2 - t1


  savepath <- '/Users/samuelthomas/Documents/anim'
  savename <- 'testanim.webm'

  pdata2 <- fm2_hmc$thetaDF
  pdata2$index <- 1:nrow(pdata2)
  n1 <- nrow(pdata2)

  # use thetaDF
  p2 <- ggplot(pdata2[1:4000, ], aes(x=index, y=woolB))
  p2 <- p2 + geom_line(color='steelblue') + theme_bw() + theme(text=element_text(size=18))
  p2 <- p2 + ggtitle('Trace Plot - HMC') + xlab('Simulation')
  p2 <- p2 + ylim(-40, 40) + geom_hline(yintercept=-16.3, color='red')
  p2 <- p2 + transition_reveal(index)

  savename_gif <- 'testanim.gif'

  a_anim <- animate(p2, duration=60, width=600, height=400,
                    renderer = gifski_renderer(paste(savepath, savename_gif, sep='/')))
  # renderer = av_renderer(paste(savepath, savename, sep='/')))

  ##############################################
  # Histogram
  ##############################################

  t1b <- Sys.time()
  pdata3 <- NULL
  iter <- c(1:1000, seq(1200, 10000, by=200))
  for (jj in 1:4000) {
    tempDF <- pdata2[1:jj, ]
    tempDF$timeval <- jj
    pdata3 <- rbind(pdata3, tempDF)
    print(jj)
  }
  t2b <- Sys.time()
  t2b - t1b

  savename2 <- 'testanim3.webm'

  savename2_gif <- 'testanim3.gif'
  p3 <- ggplot(pdata3, aes(x=woolB))
  p3 <- p3 + geom_histogram(colour="#B47846", fill="#B47846", binwidth=2) + theme_bw()
  p3 <- p3 + ggtitle('Histogram - HMC') + theme(text=element_text(size=18)) + xlim(-40, 40)
  p3 <- p3 + coord_flip()
  p3 <- p3 + transition_time(timeval)
  b_anim <- animate(p3, duration=60, width=600, height=400,
                    # renderer = av_renderer(paste(savepath, savename2, sep='/')))
                    renderer = gifski_renderer(paste(savepath, savename2_gif, sep='/')))



  ##############################################
  # Combine gifs
  ##############################################

  a_mgif <- image_read(a_anim)
  b_mgif <- image_read(b_anim)


  new_gif <- image_append(c(a_mgif[1], b_mgif[1]))
  for(i in 2:100){
    combined <- image_append(c(a_mgif[i], b_mgif[i]))
    new_gif <- c(new_gif, combined)
  }

  # save video
  savename3 <- 'hmc_animation2.webm'
  savename3_gif <- 'hmc_animation2.gif'
  image_write_video(new_gif, path=paste(savepath, savename3, sep='/'), framerate=1)
  image_write_gif(new_gif, path=paste(savepath, savename3_gif, sep='/'))




  # MH version
  N2 <- 1e5
  t1 <- Sys.time()
  set.seed(321)
  fm2_mh <- mh(N2, paramInit=rep(0, 7), qPROP = qprop_all,
               qFUN = qfun, pdFUN = linear_posterior,
               nu=c(rep(5.0, 6), 0.1),
               y =y, X = X)
  t2 <- Sys.time()
  t2 - t1

  fm2_mh$accept/N2

  plot(fm2_mh, burnin=2e4)

  # use thetaDF
  pdata4 <- fm2_mh$thetaDF
  pdata4$index <- 1:nrow(pdata4)
  colnames(pdata4)[2] <- "woolB"
  n1 <- nrow(pdata4)

  p4 <- ggplot(pdata4[1:4000, ], aes(x=index, y=woolB))
  p4 <- p4 + geom_line(color='steelblue') + theme_bw() + theme(text=element_text(size=18))
  p4 <- p4 + ggtitle('Trace Plot - Metropolis Hastings') + xlab('Simulation')
  p4 <- p4 + ylim(-40, 40)
  p4 <- p4 + geom_hline(yintercept=-16.3, color='red')
  p4 <- p4 + transition_reveal(index)

  c_anim <- animate(p4, duration=60, width=600, height=400,
                    renderer = gifski_renderer(paste(savepath, "mh1.gif", sep='/')))
  # renderer = av_renderer(paste(savepath, savename, sep='/')))

  ##############################################
  # Histogram
  ##############################################

  t1b <- Sys.time()
  pdata5 <- NULL
  for (jj in 1:4000) {
    tempDF <- pdata4[1:jj, ]
    tempDF$timeval <- jj
    pdata5 <- rbind(pdata5, tempDF)
    print(jj)
  }
  t2b <- Sys.time()
  t2b - t1b


  p5 <- ggplot(pdata5, aes(x=woolB))
  p5 <- p5 + geom_histogram(colour="#B47846", fill="#B47846", binwidth=2) + theme_bw()
  p5 <- p5 + ggtitle('Histogram - Metropolis Hastings') + theme(text=element_text(size=18))
  p5 <- p5 + xlim(-40, 40)
  p5 <- p5 + coord_flip()
  p5 <- p5 + transition_time(timeval)

  d_anim <- animate(p5, duration=60, width=600, height=400,
                    # renderer = av_renderer(paste(savepath, savename2, sep='/')))
                    renderer = gifski_renderer(paste(savepath, "mh2.gif", sep='/')))

  ##############################################
  # Combine gifs
  ##############################################

  c_mgif <- image_read(c_anim)
  d_mgif <- image_read(d_anim)

  mh_gif <- image_append(c(c_mgif[1], d_mgif[1]))
  for(i in 2:100){
    combined <- image_append(c(c_mgif[i], d_mgif[i]))
    mh_gif <- c(mh_gif, combined)
  }

  save(pdata, pdata2, pdata3, pdata3, pdata4, pdata5,
       a_anim, a_mgif, b_anim, b_mgif, c_anim, c_mgif, d_anim, d_mgif, new_gif, mh_gif,
       file = "animfiles.RData")

  # stack vertically
  hmc_mh_gif <- image_append(c(new_gif[1], mh_gif[1]), stack = T)
  for(i in 2:100){
    combined <- image_append(c(new_gif[i], mh_gif[i]), stack=T)
    hmc_mh_gif <- c(hmc_mh_gif, combined)
  }

  hmc_mh_gif


  }

