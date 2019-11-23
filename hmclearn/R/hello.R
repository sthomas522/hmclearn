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

  # eps_vals <- c(rep(2e-1, 6), 2e-2)
  eps_vlas <- c(rep(2e-2, 6), 2e-3)

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

  # animation plot
  Lval <- 20
  plot()
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



}


