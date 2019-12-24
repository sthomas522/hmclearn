
#' @export
mcmc_intervals <- function(object, ...) {
  UseMethod("mcmc_intervals")
}

#' @export
mcmc_intervals.hmclearn <- function(object, ...) {
  data <- object$thetaDF
  bayesplot::mcmc_intervals(data, ...)
}

#' @export
mcmc_areas <- function(object, ...) {
  UseMethod("mcmc_areas")
}

#' @export
mcmc_areas.hmclearn <- function(object, ...) {
  data <- object$thetaDF
  bayesplot::mcmc_areas(data, ...)
}

#' #' @export
#' color_scheme_set <- function(...) {
#'   bayesplot::color_scheme_set(...)
#' }
#'
#' #' @export
#' color_scheme_get <- function(...) {
#'   bayesplot::color_scheme_get(...)
#' }



