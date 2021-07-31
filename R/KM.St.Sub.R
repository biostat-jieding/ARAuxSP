#' @title Get subgroup survival rates using KM
#'
#' @description This is a function used to calculate survival rates
#' in each subgroup using KM.
#'
#' @param yobs observed failure time.
#' @param delta censoring indicator.
#' @param X covariates.
#' @param G indicator matrix - one row indicates one group.
#'
#' @return return the subgroup survival rates \code{tSP}, which is
#' a vector of length \code{nrow(G)}.
#'
#' @examples
#' 1
#' @export KM.St.Sub
KM.St.Sub <- function(tstar,yobs,delta,G){
  # get the estimate of subgroup survival rates using KM
  K <- nrow(G)
  tSP <- rep(0, K)
  for(k in 1:K){
    idx <- (G[k,]==1)
    fit.k <- summary(survfit(Surv(yobs, delta) ~ 1, subset=idx))
    tSP[k] <- min( c(1,fit.k$surv)[ c(0,fit.k$time) <= tstar ] )
  }
  return( tSP )
}








