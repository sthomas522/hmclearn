
#' @export
mcmc_intervals <- function(object, ...) {
  UseMethod("mcmc_intervals")
}

#' @export
mcmc_intervals.hmclearn <- function(object, burnin=NULL, ...) {
  data <- object$thetaDF
  if (!is.null(burnin)) {
    data <- data[-c(1:burnin), ]
  }
  bayesplot::mcmc_intervals(data, ...)
}

#' @export
mcmc_areas <- function(object, ...) {
  UseMethod("mcmc_areas")
}

#' @export
mcmc_areas.hmclearn <- function(object, burnin=NULL, ...) {
  data <- object$thetaDF
  if (!is.null(burnin)) {
    data <- data[-c(1:burnin), ]
  }
  bayesplot::mcmc_areas(data, ...)
}


#' @export
mcmc_hist <- function(object, ...) {
  UseMethod("mcmc_hist")
}

#' @export
mcmc_hist.hmclearn <- function(object, burnin=NULL, ...) {
  data <- object$thetaDF
  if (!is.null(burnin)) {
    data <- data[-c(1:burnin), ]
  }
  bayesplot::mcmc_hist(data, ...)
}



#' @export
mcmc_dens <- function(object, ...) {
  UseMethod("mcmc_dens")
}

#' @export
mcmc_dens.hmclearn <- function(object, burnin=NULL, ...) {
  data <- object$thetaDF
  if (!is.null(burnin)) {
    data <- data[-c(1:burnin), ]
  }
  bayesplot::mcmc_dens(data, ...)
}


#' @export
mcmc_scatter <- function(object, ...) {
  UseMethod("mcmc_scatter")
}

#' @export
mcmc_scatter.hmclearn <- function(object, burnin=NULL, ...) {
  data <- object$thetaDF
  if (!is.null(burnin)) {
    data <- data[-c(1:burnin), ]
  }
  bayesplot::mcmc_scatter(data, ...)
}


#' @export
mcmc_hex <- function(object, ...) {
  UseMethod("mcmc_hex")
}

#' @export
mcmc_hex.hmclearn <- function(object, burnin=NULL, ...) {
  data <- object$thetaDF
  if (!is.null(burnin)) {
    data <- data[-c(1:burnin), ]
  }
  bayesplot::mcmc_hex(data, ...)
}


# requires multiple chains
#' #' @export
#' mcmc_violin <- function(object, ...) {
#'   UseMethod("mcmc_violin")
#' }
#'
#' #' @export
#' mcmc_violin.hmclearn <- function(object, burnin=NULL, ...) {
#'   data <- object$thetaDF
#'   if (!is.null(burnin)) {
#'     data <- data[-c(1:burnin), ]
#'   }
#'   bayesplot::mcmc_violin(data, ...)
#' }


#' #' @export
#' color_scheme_set <- function(...) {
#'   bayesplot::color_scheme_set(...)
#' }
#'
#' #' @export
#' color_scheme_get <- function(...) {
#'   bayesplot::color_scheme_get(...)
#' }



