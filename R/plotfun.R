
#' Plotting for MCMC visualization and diagnostics provided by \code{bayesplot} package
#'
#' Plots of Rhat statistics, ratios of effective sample size to total sample
#' size, and autocorrelation of MCMC draws.
#'
#' @name hmclearn-plots
#'
#' @param object an object of class \code{hmclearn}, usually a result of a call to \code{mh} or \code{hmc}
#' @param burnin optional numeric parameter for the number of initial MCMC samples to omit from the summary
#' @param lagmax maximum lag to extract for determining effective sample sizes
#' @param ... optional additional arguments to pass to the \code{bayesplot} functions
#'
#' @section Plot Descriptions from the \code{bayesplot} package documentation:
#' \describe{
#'   \item{`mcmc_hist(object, burnin=NULL, ...)`}{
#'    Default plot called by `plot` function.  Histograms of posterior draws with all chains merged.
#'   }
#'   \item{`mcmc_dens(object, burnin=NULL, ...)`}{
#'    Kernel density plots of posterior draws with all chains merged.
#'   }
#'   \item{`mcmc_hist_by_chain(object, burnin=NULL, ...)`}{
#'    Histograms of posterior draws with chains separated via faceting.
#'   }
#'   \item{`mcmc_dens_overlay(object, burnin=NULL, ...)`}{
#'    Kernel density plots of posterior draws with chains separated but
#'    overlaid on a single plot.
#'   }
#'   \item{`mcmc_violin(object, burnin=NULL, ...)`}{
#'    The density estimate of each chain is plotted as a violin with
#'    horizontal lines at notable quantiles.
#'   }
#'   \item{`mcmc_dens_chains(object, burnin=NULL, ...)`}{
#'    Ridgeline kernel density plots of posterior draws with chains separated
#'    but overlaid on a single plot. In `mcmc_dens_overlay()` parameters
#'    appear in separate facets; in `mcmc_dens_chains()` they appear in the
#'    same panel and can overlap vertically.
#'   }
#'   \item{`mcmc_intervals(object, burnin=NULL, ...)`}{
#'    Plots of uncertainty intervals computed from posterior draws with all
#'    chains merged.
#'   }
#'   \item{`mcmc_areas(object, burnin=NULL, ...)`}{
#'    Density plots computed from posterior draws with all chains merged,
#'    with uncertainty intervals shown as shaded areas under the curves.
#'   }
#'   \item{`mcmc_scatter(object, burnin=NULL, ...)`}{
#'    Bivariate scatterplot of posterior draws. If using a very large number of
#'    posterior draws then `mcmc_hex()` may be preferable to avoid
#'    overplotting.
#'   }
#'   \item{`mcmc_hex(object, burnin=NULL, ...)`}{
#'    Hexagonal heatmap of 2-D bin counts. This plot is useful in cases where
#'    the posterior sample size is large enough that `mcmc_scatter()` suffers
#'    from overplotting.
#'   }
#'   \item{`mcmc_pairs(object, burnin=NULL, ...)`}{
#'    A square plot matrix with univariate marginal distributions along the
#'    diagonal (as histograms or kernel density plots) and bivariate
#'    distributions off the diagonal (as scatterplots or hex heatmaps).
#'
#'    For the off-diagonal plots, the default is to split the chains so that
#'    (roughly) half are displayed above the diagonal and half are below (all
#'    chains are always merged together for the plots along the diagonal). Other
#'    possibilities are available by setting the `condition` argument.
#'   }
#' \item{`mcmc_rhat(object, burnin=NULL, ...)`, `mcmc_rhat_hist(object, burnin=NULL, ...)`}{
#'   Rhat values as either points or a histogram. Values are colored using
#'   different shades (lighter is better). The chosen thresholds are somewhat
#'   arbitrary, but can be useful guidelines in practice.
#'   * _light_: below 1.05 (good)
#'   * _mid_: between 1.05 and 1.1 (ok)
#'   * _dark_: above 1.1 (too high)
#'  }
#'  \item{`mcmc_neff(object, burnin=NULL, ...)`, `mcmc_neff_hist(object, burnin=NULL, ...)`}{
#'   Ratios of effective sample size to total sample size as either points or a
#'   histogram. Values are colored using different shades (lighter is better).
#'   The chosen thresholds are somewhat arbitrary, but can be useful guidelines
#'   in practice.
#'   * _light_: between 0.5 and 1 (high)
#'   * _mid_: between 0.1 and 0.5 (good)
#'   * _dark_: below 0.1 (low)
#'  }
#'  \item{`mcmc_acf(object, burnin=NULL, ...)`, `mcmc_acf_bar(object, burnin=NULL, ...)`}{
#'   Grid of autocorrelation plots by chain and parameter. The `lags` argument
#'   gives the maximum number of lags at which to calculate the autocorrelation
#'   function. `mcmc_acf()` is a line plot whereas `mcmc_acf_bar()` is a
#'   barplot.
#'  }
#' }
#' @return These functions call various plotting functions from the \code{bayesplot} package, which returns a list including \code{ggplot2} objects.
#' @references Gabry, Jonah and Mahr, Tristan (2019).  \emph{bayesplot:  Plotting for Bayesian Models}.  \url{https://mc-stan.org/bayesplot/}
#' @references Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., and Gelman, A (2019).  \emph{Visualization in Bayesian Workflow}.  Journal of the Royal Statistical Society: Series A. Vol 182.  Issue 2.  p.389-402.
#' @references Gelman, A. and Rubin, D. (1992) \emph{Inference from Iterative Simulation Using Multiple Sequences}.  Statistical Science 7(4) 457-472.
#' @references Gelman, A., et. al. (2013) \emph{Bayesian Data Analysis}.  Chapman and Hall/CRC.
NULL

#' @rdname hmclearn-plots
#' @export
#' @examples
#' # poisson regression example
#' set.seed(7363)
#' X <- cbind(1, matrix(rnorm(40), ncol=2))
#' betavals <- c(0.8, -0.5, 1.1)
#' lmu <- X %*% betavals
#' y <- sapply(exp(lmu), FUN = rpois, n=1)
#'
#' f <- hmc(N = 1000,
#'           theta.init = rep(0, 3),
#'           epsilon = c(0.03, 0.02, 0.015),
#'           L = 10,
#'           logPOSTERIOR = poisson_posterior,
#'           glogPOSTERIOR = g_poisson_posterior,
#'           varnames = paste0("beta", 0:2),
#'           param = list(y=y, X=X),
#'           parallel=FALSE, chains=2)
#'
#' mcmc_trace(f, burnin=100)
#' mcmc_hist(f, burnin=100)
#' mcmc_intervals(f, burnin=100)
#' mcmc_rhat(f, burnin=100)
#' mcmc_violin(f, burnin=100)
mcmc_intervals <- function(object, ...) {
  UseMethod("mcmc_intervals")
}

#' @rdname hmclearn-plots
#' @export
mcmc_intervals.hmclearn <- function(object, burnin=NULL, ...) {
  data <- combMatrix(object$thetaCombined, burnin=burnin)
  bayesplot::mcmc_intervals(data, ...)
}

#' @rdname hmclearn-plots
#' @export
mcmc_areas <- function(object, ...) {
  UseMethod("mcmc_areas")
}

#' @rdname hmclearn-plots
#' @export
mcmc_areas.hmclearn <- function(object, burnin=NULL, ...) {
  data <- combMatrix(object$thetaCombined, burnin=burnin)
  bayesplot::mcmc_areas(data, ...)
}

#' @rdname hmclearn-plots
#' @export
mcmc_hist <- function(object, ...) {
  UseMethod("mcmc_hist")
}

#' @rdname hmclearn-plots
#' @export
mcmc_hist.hmclearn <- function(object, burnin=NULL, ...) {
  data <- combMatrix(object$thetaCombined, burnin=burnin)
  bayesplot::mcmc_hist(data, ...)
}

#' @rdname hmclearn-plots
#' @export
mcmc_hist_by_chain <- function(object, ...) {
  UseMethod("mcmc_hist_by_chain")
}

#' @rdname hmclearn-plots
#' @export
mcmc_hist_by_chain.hmclearn <- function(object, burnin=NULL, ...) {
  data <- combMatrix(object$thetaCombined, burnin=burnin)
  bayesplot::mcmc_hist_by_chain(data, ...)
}

#' @rdname hmclearn-plots
#' @export
mcmc_dens <- function(object, ...) {
  UseMethod("mcmc_dens")
}

#' @rdname hmclearn-plots
#' @export
mcmc_dens.hmclearn <- function(object, burnin=NULL, ...) {
  data <- combMatrix(object$thetaCombined, burnin=burnin)
  bayesplot::mcmc_dens(data, ...)
}

#' @rdname hmclearn-plots
#' @export
mcmc_scatter <- function(object, ...) {
  UseMethod("mcmc_scatter")
}

#' @rdname hmclearn-plots
#' @export
mcmc_scatter.hmclearn <- function(object, burnin=NULL, ...) {
  data <- combMatrix(object$thetaCombined, burnin=burnin)
  bayesplot::mcmc_scatter(data, ...)
}

#' @rdname hmclearn-plots
#' @export
mcmc_hex <- function(object, ...) {
  UseMethod("mcmc_hex")
}

#' @rdname hmclearn-plots
#' @export
mcmc_hex.hmclearn <- function(object, burnin=NULL, ...) {
  data <- combMatrix(object$thetaCombined, burnin=burnin)
  bayesplot::mcmc_hex(data, ...)
}

#' @rdname hmclearn-plots
#' @export
mcmc_pairs <- function(object, ...) {
  UseMethod("mcmc_pairs")
}

#' @rdname hmclearn-plots
#' @export
mcmc_pairs.hmclearn <- function(object, burnin=NULL, ...) {
  data <- combMatrix(object$thetaCombined, burnin=burnin)
  bayesplot::mcmc_pairs(data, ...)
}

#' @rdname hmclearn-plots
#' @export
mcmc_acf <- function(object, ...) {
  UseMethod("mcmc_acf")
}

#' @rdname hmclearn-plots
#' @export
mcmc_acf.hmclearn <- function(object, burnin=NULL, ...) {
  data <- combMatrix(object$thetaCombined, burnin=burnin)
  bayesplot::mcmc_acf(data, ...)
}

#' @rdname hmclearn-plots
#' @export
mcmc_acf_bar <- function(object, ...) {
  UseMethod("mcmc_acf_bar")
}

#' @rdname hmclearn-plots
#' @export
mcmc_acf_bar.hmclearn <- function(object, burnin=NULL, ...) {
  data <- combMatrix(object$thetaCombined, burnin=burnin)
  bayesplot::mcmc_acf_bar(data, ...)
}

#' @rdname hmclearn-plots
#' @export
mcmc_trace <- function(object, ...) {
  UseMethod("mcmc_trace")
}

#' @rdname hmclearn-plots
#' @export
mcmc_trace.hmclearn <- function(object, burnin=NULL, ...) {
  data <- combMatrix(object$thetaCombined, burnin=burnin)
  bayesplot::mcmc_trace(data, ...)
}

#' @rdname hmclearn-plots
#' @export
mcmc_rhat <- function(object, ...) {
  UseMethod("mcmc_rhat")
}

#' @rdname hmclearn-plots
#' @export
mcmc_rhat.hmclearn <- function(object, burnin=NULL, ...) {
  rhatvals <- psrf(object, burnin=burnin)
  names(rhatvals) <- object$varnames
  bayesplot::mcmc_rhat(rhatvals, ...) #+ bayesplot::yaxis_text(hjust=1)
}

#' @rdname hmclearn-plots
#' @export
mcmc_rhat_hist <- function(object, ...) {
  UseMethod("mcmc_rhat_hist")
}

#' @rdname hmclearn-plots
#' @export
mcmc_rhat_hist.hmclearn <- function(object, burnin=NULL, ...) {
  rhatvals <- psrf(object, burnin=burnin)
  names(rhatvals) <- object$varnames
  bayesplot::mcmc_rhat_hist(rhatvals, ...)
}

#' @rdname hmclearn-plots
#' @export
mcmc_neff <- function(object, ...) {
  UseMethod("mcmc_neff")
}

#' @rdname hmclearn-plots
#' @export
mcmc_neff.hmclearn <- function(object, burnin=NULL, lagmax=NULL, ...) {
  neffvals <- neff(object, burnin, lagmax)

  if (!is.null(burnin)) {
    N <- object$N - burnin
  } else {
    N <- object$N
  }

  ratio <- neffvals / N

  bayesplot::mcmc_neff(ratio, ...)
}

#' @rdname hmclearn-plots
#' @export
mcmc_neff_hist <- function(object, ...) {
  UseMethod("mcmc_neff_hist")
}

#' @rdname hmclearn-plots
#' @export
mcmc_neff_hist.hmclearn <- function(object, burnin=NULL, lagmax=NULL, ...) {
  neffvals <- neff(object, burnin, lagmax)

  if (!is.null(burnin)) {
    N <- object$N - burnin
  } else {
    N <- object$N
  }

  ratio <- neffvals / N

  bayesplot::mcmc_neff_hist(ratio, ...)
}

#' @rdname hmclearn-plots
#' @export
mcmc_neff_data <- function(object, ...) {
  UseMethod("mcmc_neff_data")
}

#' @rdname hmclearn-plots
#' @export
mcmc_neff_data.hmclearn <- function(object, burnin=NULL, lagmax=NULL, ...) {
  neffvals <- neff(object, burnin, lagmax)

  if (!is.null(burnin)) {
    N <- object$N - burnin
  } else {
    N <- object$N
  }

  ratio <- neffvals / N

  bayesplot::mcmc_neff_data(ratio, ...)
}

# requires multiple chains
#' @rdname hmclearn-plots
#' @export
mcmc_violin <- function(object, ...) {
  UseMethod("mcmc_violin")
}

#' @rdname hmclearn-plots
#' @export
mcmc_violin.hmclearn <- function(object, burnin=NULL, ...) {
  data <- combMatrix(object$thetaCombined, burnin=burnin)
  bayesplot::mcmc_violin(data, ...)
}


#' Plot Histograms of the Posterior Distribution
#'
#' Calls \code{mcmc_hist} from the \code{bayesplot} package to display histograms of the posterior
#'
#' @param x an object of class \code{hmclearn}, usually a result of a call to \code{mh} or \code{hmc}
#' @param burnin optional numeric parameter for the number of initial MCMC samples to omit from the summary
#' @param ... optional additional arguments to pass to the \code{bayesplot} functions
#'
#' @references Gabry, Jonah and Mahr, Tristan (2019).  \emph{bayesplot:  Plotting for Bayesian Models}.  \url{https://mc-stan.org/bayesplot/}
#' @return Calls \code{mcmc_hist} from the \code{bayesplot} package, which returns a list including a \code{ggplot2} object.
#' @export
#' @examples
#' # poisson regression example
#' set.seed(7363)
#' X <- cbind(1, matrix(rnorm(40), ncol=2))
#' betavals <- c(0.8, -0.5, 1.1)
#' lmu <- X %*% betavals
#' y <- sapply(exp(lmu), FUN = rpois, n=1)
#'
#' f <- hmc(N = 1000,
#'           theta.init = rep(0, 3),
#'           epsilon = c(0.03, 0.02, 0.015),
#'           L = 10,
#'           logPOSTERIOR = poisson_posterior,
#'           glogPOSTERIOR = g_poisson_posterior,
#'           varnames = paste0("beta", 0:2),
#'           param = list(y=y, X=X),
#'           parallel=FALSE, chains=2)
#'
#' plot(f, burnin=100)
plot.hmclearn <- function(x, burnin=NULL, ...) {
  thetaCombined <- combMatrix(x$thetaCombined, burnin=burnin)
  bayesplot::mcmc_hist(thetaCombined, ...)
}


