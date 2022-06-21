
#============== Auxiliary Information's Estimating Equations and derivatives  ==============#
#' @title subgroup survival rate estimating equation
#'
#' @description This is a function used to calculate the estimating equation
#' for the subgroup survival rates at t* under the additive risk model (
#' its derivation is based on 'Double Expectation')
#'
#' @aliases AR.AuxSP.EE AR.AuxSP.EE.Sigma
#'
#' @param bet betas that will be computed at
#' @param alp alpha that will be computed at
#' @inheritParams AR
#' @param auxinfo auxiliary ifnormation
#'
#' @return A vector \code{Psin}.
#'
#' @examples
#' 1
#' @export AR.AuxSP.EE
AR.AuxSP.EE <- function(bet,alp,yobs,delta,X,auxinfo){
  # Auxiliary Information's Estimating Equations at different time points

  # prepare
  p <- ncol(X) # number of covariates in internal data
  N <- length(yobs)   # sample size in internal data
  # the estimating equations in individual levels
  sumK <- nrow(auxinfo) # number of Groups
  Psi <- array(0, dim=c(N, sumK))
  for(iG in 1:sumK){ # iG <- 1
    sur <- exp(-alp[auxinfo[iG,'tstaridx']]-X%*%bet*auxinfo[iG,'tstar'])
    G <- auxinfo[iG,-c(1:4)]
    sprob <- auxinfo[iG,'sprob']
    Psi[,iG] <- (sur-sprob) * G
  }
  # output
  return(Psi)
}

AR.AuxSP.EE.D1 <- function(bet,alp,yobs,delta,X,auxinfo){
  # Auxiliary Information's Estimating Equations's first derivatives at different time points

  # prepare
  p <- ncol(X) # number of covariates in internal data
  N <- length(yobs)   # sample size in internal data
  # calculate psi.bet, psi.alp, psi.psi
  sur.allG <- t(sapply(1:nrow(auxinfo),function(iG){
    exp(-alp[auxinfo[iG,'tstaridx']]-X%*%bet*auxinfo[iG,'tstar'])}))
  psi.bet <- - (sur.allG*auxinfo[,-c(1:4)]*auxinfo[,'tstar']) %*% X / N
  psi.alp0 <- - as.matrix( apply( (sur.allG*auxinfo[,-c(1:4)]) ,1,mean) )
  psi.alp <- sapply(1:length(alp),function(itstaridx){psi.alp0*(auxinfo[,'tstaridx']==itstaridx)} )
  Psi.D1 <- cbind(psi.bet,psi.alp)
  # output
  return( Psi.D1 )
}



#==========================================================================#
# KM method for subgroup survival rates - [it is needed when considering uncertainty!]
#==========================================================================#

#============== KM method for subgroup survival rates ==============#
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
    fit.k <- summary(survival::survfit(survival::Surv(yobs, delta) ~ 1, subset=idx))
    tSP[k] <- min( c(1,fit.k$surv)[ c(0,fit.k$time) <= tstar ] )
  }
  return( tSP )
}

#' @export
KM <- function(tm,yobs,delta,type="right"){

  ### preparation ###
  N <- length(yobs)
  y.order <- order(yobs)
  y.sort <- yobs[y.order]
  delta.sort <- delta[y.order]
  yobs.1max <- max(yobs[delta==1])

  # # way 1
  # y.rank <- rank(yobs,ties.method="min")
  # w.sort <- delta.sort/(N:1) * c(1,cumprod(1-delta.sort/(N-(1:N)+1))[-N])
  # if(type=="right"){
  #   KMt <- 1-cumsum(w.sort)[y.rank] # right-continuous
  # }else{
  #   KMt <- 1 - c(0,cumsum(w.sort)[-N])[y.rank] # left-continuous
  # }
  # w <- w.sort[y.rank] # w*N = delta / KM(yobs,yobs,1-delta,type="left")
  # way 2
  prods <- 1-delta.sort/(N-(1:N)+1)
  if(type=="right"){
    KMt <- sapply(tm,function(tmi){prod(prods[y.sort<=tmi])}) # right-continuous
    # KMt[tm>=yobs.1max] <- 0
  }else{
    KMt <- sapply(tm,function(tmi){prod(prods[y.sort<tmi])}) # left-continuous
    # KMt[tm>yobs.1max] <- 0
  }

  # output
  return(KMt)

}
