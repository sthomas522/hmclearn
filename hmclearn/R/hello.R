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

# hmc example
if (1 == 0) {
  library(ggplot2)
  library(gganimate)
  fm1 <- lm(breaks ~ wool*tension, data=warpbreaks)
  summary(fm1)
  X <- model.matrix(breaks ~ wool*tension, data=warpbreaks)
  y <- warpbreaks$breaks

  N <- 10000
  set.seed(143)

  eps_vals <- c(rep(2e-1, 6), 2e-2)
  # eps_vals <- c(rep(2e-2, 6), 2e-3)

  t1 <- Sys.time()
  set.seed(321)
  fm1_hmc <- hmc(N, theta.init = c(rep(0, 6), 1), epsilon = eps_vals, L = 20,
                 logPOSTERIOR = linear_posterior,
                 glogPOSTERIOR = g_linear_posterior, y=y, X=X,
                 varnames = c(colnames(X), "log_sigma_sq"))
  t2 <- Sys.time()
  t2 - t1

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

  eps_vals <- c(rep(2e-1, 6), 2e-2)
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
  p2 <- ggplot(pdata2[1:2000, ], aes(x=index, y=woolB))
  p2 <- p2 + geom_line(color='steelblue') + theme_bw() + theme(text=element_text(size=18))
  p2 <- p2 + ggtitle('Trace Plot') + xlab('Simulation')
  p2 <- p2 + transition_reveal(index)

  savename_gif <- 'testanim.gif'

  a_anim <- animate(p2, duration=60, width=600, height=400,
            renderer = gifski_renderer(paste(savepath, savename_gif, sep='/')))
          # renderer = av_renderer(paste(savepath, savename, sep='/')))

  ##############################################
  # Histogram
  ##############################################

  # t1b <- Sys.time()
  # pdata3 <- NULL
  # for (jj in 1:2000) {
  #   tempDF <- pdata2[1:jj, ]
  #   tempDF$timeval <- jj
  #   pdata3 <- rbind(pdata3, tempDF)
  # }
  # t2b <- Sys.time()
  # t2b - t1b
  savename2 <- 'testanim2.webm'

  savename2_gif <- 'testanim2.gif'
  p3 <- ggplot(pdata3, aes(x=woolB))
  p3 <- p3 + geom_histogram(colour="#B47846", fill="#B47846", binwidth=2) + theme_bw()
  p3 <- p3 + ggtitle('Histogram') + theme(text=element_text(size=18))
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
  savename3 <- 'hmc_animation.webm'
  savename3_gif <- 'hmc_animation.gif'
  image_write_video(new_gif, path=paste(savepath, savename3, sep='/'), framerate=1)
  image_write_gif(new_gif, path=paste(savepath, savename3_gif, sep='/'))

}

