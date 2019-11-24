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
                       V2 = seq(-30, -5, by=0.1),
                       V3 = theta.median[3],
                       V4 = theta.median[4],
                       V5 = theta.median[5],
                       V6 = theta.median[6],
                       V7 = theta.median[7])


  foo <- apply(cdata, MARGIN=1, FUN=linear_posterior, X=X, y=y)

  Lval <- 20

  for (k in 1:100) {

    basenum <- 1000*Lval + Lval*k
    y1 <- theta.all$V2[(basenum+1):(basenum+Lval)]
    x1 <- r.all$V2[(basenum+1):(basenum+Lval)]
    if (k == 1) {
      plot(x1, y1, type='o', xlim = c(-2.5, 2.5),
           ylim = c(-40, 0))
    } else {
      lines(x1, y1, type='o')
    }

  }

  # animation plot
  basenum <- 1000*Lval
  pdata <- NULL

  library(gganimate)
  for (jj in 1:50) {
    tempDF <- data.frame(p = r.all$V2[basenum:(basenum+jj*Lval)],
                         theta = theta.all$V2[basenum:(basenum+jj*Lval)],
                         timevar =  jj)
    pdata <- rbind(pdata, tempDF)
  }


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


}


