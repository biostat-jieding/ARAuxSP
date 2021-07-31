
#==========================================================================#
# This file contains the elements relating to Lin and Ying (1994)
#==========================================================================#
#' @title Fit Additive Risk Model
#'
#' @description This is a function used to fit the additive risk model of Aalen (1989) based on
#' the estimation procedure proposed by Lin and Ying (1994).
#'
#' @aliases AR AR.fit
#'
#' @param formula formula object, with the response on the left of a ~ operator,
#' and the terms on the right. The response must be a survival object as
#' returned by the Surv function.
#' @param data a data.frame in which to interpret the variables named in the formula.
#' @param na.action a missing-data filter function.
#' This is applied to the model.frame after any subset argument has been used.
#'
#' @return A list \code{out} representing the fit,
#' which contains:
#' \code{coef} for estimated coefficients,
#' \code{convergence} for the optimizatoin is converged or not,
#' et al.
#'
#' @references Aalen, O. O. (1989).
#' A linear regression model for the analysis of life times.
#' Statistics in Medicine, 8(8):907–925.
#' @references Lin, D. Y. and Ying, Z. (1994).
#' Semiparametric analysis of the additive risk model.
#' Biometrika, 81(1):61–71.
#'
#' @examples
#' data("toydata1")
#' AR(Surv(yobs,delta)~X1+X2,data=toydata1)
#'
#' @importFrom survival Surv
#' @export AR
AR <- function(formula,data,na.action=na.omit){

  ## basic
  call <- match.call()
  indx <- match(c("formula", "data", "na.action"), names(call), nomatch = 0)
  if (indx[1] == 0)
    stop("A formula argument is required")

  ## prepare data
  sdata <- data.frame(data)
  mf <- model.frame(formula, sdata) # [dataframe: for latency part]
  mf <- na.action(mf) # dealting with nans
  N <- nrow(mf)
  Y <- model.response(mf)
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
  vars.name <- all.vars(formula)[-c(1,2)] # covariates
  X <-  as.matrix(sdata[,vars.name])
  yobs <- Y[,1]   # observed survival time
  delta <- Y[,2] # censoring indicator
  vars.num <- ncol(X) # num of covariates

  ## fit the model
  arfit <- AR.fit(yobs,delta,X)
  rownames( arfit$coef ) <- vars.name

  ## tidy the results
  fit <- list()
  class(fit) <- c("AR")
  fit$sdata <- arfit$sdata
  fit$call <- call
  fit$formula <- formula
  fit$vars.num <- vars.num
  fit$coef <- arfit$coef

  ## print the fitting results
  if (!is.null(cl <- fit$call)) {
    cat("Call:\n")
    dput(cl)
  }
  cat("\nAdditive Risk Model:\n")
  print(fit$coef)
  cat("\n")

  invisible(fit)

}


#' @export
AR.fit <- function(yobs,delta,X){
  # Additive Risk Model
  # From: Lin and Ying (1994)
  # Arguments:
  #   yobs: observed failure time
  #   delta:censoring indicator
  #   X:    covariates (should be a matrix and have colnames)

  # some preparation
  N <- length(yobs)
  p <- ncol(X)
  y.sort <- sort( yobs )
  y.rank <- rank( yobs, ties.method = 'min')
  X.bar.sort <- array(0, dim=c(N, p))
  for( j in 1:N){
    Y <- (yobs >= y.sort[j])
    X.bar.sort[j,] <- apply( X * Y, 2, sum ) / sum(Y)
  }
  y.sort.diff <- diff( c(0,y.sort) )
  # calculate A, B and d
  A0 <- B0 <- array(0, dim=c(p, p))
  d0 <- rep(0,p)
  for( i in 1:N ){ # i
    Ki <- y.rank[i]
    Xi.aug <- matrix(rep(X[i,], Ki), nrow=Ki, byrow = T)
    Ri <- Xi.aug - X.bar.sort[1:Ki,]
    di <- y.sort.diff[1:Ki]
    I2i <- t(Ri) %*% di %*% X[i,]
    I1i <- ( X[i,] - X.bar.sort[y.rank[i],] ) * delta[i]
    A0 <- A0 + I2i
    d0 <- d0 + I1i
    B0 <- B0 + delta[i]*( Ri[Ki,]%*%t(Ri[Ki,]) )
  }
  A <- A0/N; B <- B0/N; d <- d0/N
  # calculate the estimate of beta and SE
  Est    <- solve(A,d)
  Sigma  <- solve(A) %*% B %*% solve(A) # asymptotic var-cov matrix
  SE     <- sqrt( diag(Sigma)/N )
  zvalue <- Est/SE
  pvalue <- 2*(1-pnorm(abs(zvalue)))
  coef <- data.frame(Est=Est, SE=SE, zvalue=zvalue, pvalue=pvalue,
                     row.names=colnames(X))
  # output
  out <- list(
    sdata=list(yobs=yobs,delta=delta,X=X),
    coef=coef
  )
  return(out)
}

#==========================================================================#
# baseline cumulative hazard function in Lin and Ying (1994)
#==========================================================================#
#' @title baseline cumulative hazard function in Additive Risk Model
#'
#' @description This is a function used to get baseline cumulative hazard
#'   function in Lin and Ying (1994).
#'
#' @param object The result of an AR fit. Like: AR(), AR.AuxSP.GMM(),
#'   AR.AuxSP.AGMM() ....
#' @param tm time points that will be calculated at (length m)
#'   If 'NULL', calculate at yobs. (with sort!!!)
#' @return A list \code{out}.
#' @examples
#' AR(1,1)
#'
#' @export AR.Ht
AR.Ht <- function(object, tm=NULL){

  ### some preparation ###
  yobs  <- object$sdata$yobs
  delta <- object$sdata$delta
  X <- object$sdata$X
  bet <- object$coef[,1]
  N <- length(yobs); p <- ncol(X)

  ### calculate Lambda0(t) at time points yobs ###
  num.Risk  <- sapply( yobs, function(yi){sum( yobs >= yi )} )
  y.sort <- sort( yobs )
  y.sort.diff <- diff( c(0,y.sort) )
  X.bar.sort <- array(0, dim=c(N, p))
  cumhaz0.yobs <- rep(NA,N)
  for( j in 1:N){
    # calculate X.bar.sort till this time point yobsj
    Y <- (yobs >= y.sort[j])
    X.bar.sort[j,] <- apply( X * Y, 2, sum ) / sum(Y)
    # calculate the C(yobsj) at this time point
    Ct <- t(X.bar.sort[1:j,,drop=FALSE]) %*% y.sort.diff[1:j]
    # calculate the first part
    Ft <- sum(  (yobs<=y.sort[j]) * delta / num.Risk  )
    # Final Result
    cumhaz0.yobs[j] <- Ft - as.numeric( t(bet) %*% Ct )
  }

  ### calculate Lambda0(t) at time points tm ###
  if(is.null(tm)){
    tm <- y.sort
    cumhaz0 <- cumhaz0.yobs
  }else{
    cumhaz0 <- c(0,cumhaz0.yobs)[
      sapply(tm,function(tmi){sum( c(0, y.sort) <= tmi )}) ]
  }

  ### Output ###
  out <- data.frame(tm=tm,cumhaz0=cumhaz0)
  return(out)

}

#==========================================================================#
# S(t|X=x) function in Lin and Ying (1994)
#==========================================================================#
#' @title baseline cumulative hazard function in Additive Risk Model
#'
#' @description This is a function used to get baseline cumulative hazard
#'   function in Lin and Ying (1994).
#'
#' @param object The result of an AR fit. Like: AR(), AR.AuxSP.GMM(),
#'   AR.AuxSP.AGMM() ....
#' @param tm time points that will be calculated at (length m)
#'   If 'NULL', calculate at yobs. (without sort!!!)
#' @param x covariates that will be calculated at (length l)
#'   If 'NULL', calculate at X
#' @return A list \code{out}.
#' @examples
#' AR(1,1)
#'
#' @export AR.Stx
AR.Stx <- function(object,tm=NULL,x=NULL){

  ### some preparation ###
  yobs  <- object$sdata$yobs
  delta <- object$sdata$delta
  X <- object$sdata$X
  bet <- object$coef[,1]
  N <- length(yobs)
  if(is.null(tm)){tm <- yobs}
  if(is.null(x)){x <- X}
  if(is.null(dim(x))){x <- t(x)}
  l <- nrow(x); m <- length(tm)

  ### Obtain the Baseline Cumhaz by using AR.Ht ###
  cumhaz0 <- AR.Ht(object, tm=tm)[,2]

  ### Give estimates for S(t|x) - l*m ###
  Stx <- matrix(nrow=l,ncol=m,dimnames=list(paste("x",1:l,sep=""),paste("time",1:m,sep="")))
  for(im in 1:m){
    for(il in 1:l){
      Stx[il,im] <- exp(-cumhaz0[im]-sum(x[il,]*bet)*tm[im])
    }
  }

  ### Output ###
  out <- Stx #list(tm=tm, x=x, Stx=Stx)
  return(out)

}


#==========================================================================#
# Ling and Ying (1994)'s Estimating Equations and its Asymptotic Covariance Matrix
#==========================================================================#
#' @title Lin and Ying (1994)'s estimating equation
#'
#' @description This is a function used to calculate Lin and Ying (1994)'s
#' estimating equation. From: Lin and Ying (1994)'s formula:
#' \code{(1/n) (2.7)  = (1/n)Sum{ I1i-I2i bet }}
#' @aliases AR.EE AR.EE.Sigma
#'
#' @param bet betas that will be computed at
#' @inheritParams AR
#'
#' @return A vector \code{Phin}.
#'
#' @examples
#' 1
#' @export AR.EE
AR.EE <- function(bet, yobs, delta, X){

  # some preparation
  N <- length(yobs)
  p <- ncol(X)
  y.sort <- sort( yobs )
  y.rank <- rank( yobs, ties.method = 'min')
  X.bar.sort <- array(0, dim=c(N, p))
  for( j in 1:N){
    Y <- (yobs >= y.sort[j])
    X.bar.sort[j,] <- apply( X * Y, 2, sum ) / sum(Y)
  }
  y.sort.diff <- diff( c(0,y.sort) )
  # calculate A and d
  A0 <- array(0, dim=c(p, p))
  d0 <- rep(0,p)
  for( i in 1:N ){ # i
    Ki <- y.rank[i]
    Xi.aug <- matrix(rep(X[i,], Ki), nrow=Ki, byrow = T)
    Ri <- Xi.aug - X.bar.sort[1:Ki,]
    di <- y.sort.diff[1:Ki]
    I2i <- t(Ri) %*% di %*% X[i,]
    I1i <- ( X[i,] - X.bar.sort[y.rank[i],] ) * delta[i]
    A0 <- A0 + I2i
    d0 <- d0 + I1i
  }
  A <- A0/N; d <- d0/N
  # output the averaged estimating equations
  Phin <- d-A%*%bet
  return( as.vector(Phin) )

}

#' @export
AR.EE.Sigma <- function(yobs,delta,X){
  # the asymptotic covariance matrix of: sqrt(n)*AR.EE(beta0)
  # From: Lin and Ying (1994)'s formula [ B ]

  # prepare
  N <- length(yobs)
  p <- ncol(X)
  X.bar <- array(0, dim=c(N,p))  # X.bar matrix [ N×p ]
  for( j in 1:N){
    Y <- (yobs >= yobs[j])
    X.bar[j,] <- apply( X * Y, 2, sum ) / sum(Y)
  }
  xxx <- (X-X.bar)
  # output
  SigmaPhi <- t( xxx * delta ) %*% xxx / N
  return(SigmaPhi)

}


#==========================================================================#
# Other
#==========================================================================#
AR.St.Sub <- function(tstar,yobs,delta,X,G){
  # get the estimate of subgroup survival rates using AR

  N <- length(yobs)

  arfit <- AR(yobs,delta,X) # Lin and Ying's estimator
  bet.LY <- arfit$coef[,1]
  alp.LY <- AR.Ht(arfit,tm=tstar)[,2]

  sur.init <- exp(-alp.LY-(X%*%bet.LY)*tstar)
  a <- ( G%*%sur.init ) / N
  b <- apply(G,1,mean)
  tSP <- as.vector(a/b)

  return( tSP )
}




