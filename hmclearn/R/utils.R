#' @importFrom stats aggregate
#' @importFrom stats quantile
#' @importFrom stats rbinom
#' @importFrom stats reshape
#' @importFrom stats rnorm
#' @importFrom stats rpois
#' @importFrom stats runif
#' @importFrom stats var
#' @importFrom utils tail

create_Uj <- function(uj, neg=TRUE) {
  q <- length(uj)
  if (q == 1) return(0)
  nr <- q
  nc <- q*(q-1)/2
  Uj <- matrix(0, nrow = nr, ncol = nc)
  for (kk in 2:nr)
  { Uj[kk, ((kk-1)*(kk-2)/2 + 1):(kk*(kk-1)/2)] <- uj[1:(kk-1)] }
  if (neg) {
    return(-Uj)
  } else {
    return(Uj)
  }
}


