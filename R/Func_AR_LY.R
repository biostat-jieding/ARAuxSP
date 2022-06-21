


#==========================================================================#
# Fit Additive Risk Model: All the elements relating to Lin and Ying (1994)'s AR Model ####
#==========================================================================#

#============== The main function that fit the AR model ==============####
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
#' Statistics in Medicine, 8(8):907-925.
#' @references Lin, D. Y. and Ying, Z. (1994).
#' Semiparametric analysis of the additive risk model.
#' Biometrika, 81(1):61-71.
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
  ## print fit for quick viewing
  invisible(fit)
}


#' @export
AR.fit <- function(yobs,delta,X){
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
    Rdi <- Ri * sqrt(di)
    I2i <- t(Rdi) %*% Rdi # I2i <- t(Ri) %*% di %*% X[i,]
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

#============== baseline cumulative hazard function in Lin and Ying (1994) ==============#
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

  ### to make it momtone
  cumhaz0.yobs.mon <- sapply(1:length(cumhaz0.yobs),function(icum){
    max(c(0,cumhaz0.yobs[1:icum]))
  })
  # sfit <- survfit(Surv(yobs, delta)~1)
  # cumhaz0.yobs.mon <- sfit$cumhaz
  # y.sort <- sfit$time

  ### calculate Lambda0(t) at time points tm ###
  if(is.null(tm)){
    tm <- y.sort
    cumhaz0 <- cumhaz0.yobs.mon
  }else{
    cumhaz0 <- c(0,cumhaz0.yobs.mon)[
      sapply(tm,function(tmi){sum( c(0, y.sort) <= tmi )}) ]
  }
  ### Output ###
  out <- data.frame(tm=tm,cumhaz0=cumhaz0)
  return(out)
}

#============== S(t|X=x) function in Lin and Ying (1994) ==============#
#' @title survival function in Additive Risk Model
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


#============== Estimators for subgroup survival rates using AR ==============#
#' @export
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




#============== Ling and Ying (1994)'s Estimating Equations and derivatives ==============#
#' @title Lin and Ying (1994)'s estimating equation
#'
#' @description This is a function used to calculate Lin and Ying (1994)'s
#' estimating equation. From: Lin and Ying (1994)'s formula:
#' \code{(1/n) (2.7)  = (1/n)Sum{ I1i-I2i bet }} (in individual level details)
#' @aliases AR.EE
#'
#' @param bet betas that will be computed at
#' @inheritParams AR
#'
#' @return A vector \code{Phin}.
#'
#' @examples
#' 1
#' @export AR.EE
AR.EE <- function(bet,yobs,delta,X){
  # return etimating equations: in individual level details

  # some preparation
  N <- length(yobs)
  p <- ncol(X)
  # sort and rank(Ki)
  y.sort <- sort( yobs )
  y.rank <- rank( yobs, ties.method = 'min')
  y.sort.diff <- diff( c(0,y.sort) )
  ###  Prepare values relating to LY's X: original covariates ###
  # X.bar matrix [ N*pX ]
  X.bar.sort <- array(0, dim=c(N, p))
  for( j in 1:N){ # j <- 10; i <- 1
    Y <- (yobs >= y.sort[j])
    X.bar.sort[j,] <- apply( X * Y, 2, sum ) / sum(Y)
  }
  # calculate Phi
  Phi <- array(0, dim=c(N, p))
  for( i in 1:N ){ # i
    Ki <- y.rank[i]
    Xi.aug <- matrix(rep(X[i,], Ki), nrow=Ki, byrow = T)
    Ri <- Xi.aug - X.bar.sort[1:Ki,]
    di <- y.sort.diff[1:Ki]
    Rdi <- Ri * sqrt(di)
    I2i <- t(Rdi) %*% Rdi # I2i <- t(Ri) %*% di %*% X[i,]
    I1i <- ( X[i,] - X.bar.sort[y.rank[i],] ) * delta[i]
    Phi[i,] <- I1i - I2i %*% bet
  }
  # output
  return(Phi)

}

#' @export
AR.EE.D1 <- function(yobs,delta,X){
  # return the first derivative of etimating equations

  # some preparation
  N <- length(yobs)
  p <- ncol(X)
  # sort and rank(Ki)
  y.sort <- sort( yobs )
  y.rank <- rank( yobs, ties.method = 'min')
  y.sort.diff <- diff( c(0,y.sort) )
  ###  Prepare values relating to LY's X: original covariates ###
  # X.bar matrix [ N*pX ]
  X.bar.sort <- array(0, dim=c(N, p))
  for( j in 1:N){ # j <- 10; i <- 1
    Y <- (yobs >= y.sort[j])
    X.bar.sort[j,] <- apply( X * Y, 2, sum ) / sum(Y)
  }
  # calculate Phi
  A0 <- array(0, dim=c(p, p))
  for( i in 1:N ){ # i
    Ki <- y.rank[i]
    Xi.aug <- matrix(rep(X[i,], Ki), nrow=Ki, byrow = T)
    Ri <- Xi.aug - X.bar.sort[1:Ki,]
    di <- y.sort.diff[1:Ki]
    I2i <- t(Ri) %*% di %*% X[i,]
    A0 <- A0 + I2i
  }
  A <- -A0/N
  # output
  return(A)
}

#' @export
AR.EE2 <- function(bet, yobs, delta, X){
  # return etimating equations: in individual level details

  # some preparation
  N <- length(yobs)
  pX <- ncol(X)
  # sort and rank(Ki)
  y.sort <- sort( yobs )
  y.rank <- rank( yobs, ties.method = 'min')
  y.sort.diff <- diff( c(0,y.sort) )
  ###  Prepare values relating to LY's X: original covariates ###
  # X.bar matrix [ N*pX ]
  X.bar.sort <- array(0, dim=c(N, pX))
  for( j in 1:N){ # j <- 10; i <- 1
    Y <- (yobs >= y.sort[j])
    X.bar.sort[j,] <- apply( X * Y, 2, sum ) / sum(Y)
  }
  # calculate Phi
  Phi <- array(0, dim=c(N, pX))
  for( i in 1:N ){ # i
    Ki <- y.rank[i]
    Xi.aug <- matrix(rep(X[i,], Ki), nrow=Ki, byrow = T)
    Ri <- Xi.aug - X.bar.sort[1:Ki,]
    di <- y.sort.diff[1:Ki]
    Di <- array(0,dim=c(Ki,Ki)); diag(Di) <- di;
    I2i <- t(Ri) %*% Di %*% Ri # t(Ri) %*% di %*% X[i,] #
    I1i <- ( X[i,] - X.bar.sort[y.rank[i],] ) * delta[i]
    Phi[i,] <- I1i - I2i %*% bet
  }
  # output
  return(Phi)

}

#' @export
AR.bet.Sigma <- function(yobs,delta,X){

  # some preparation
  N <- length(yobs)
  pX <- ncol(X)
  y.sort <- sort( yobs )
  y.rank <- rank( yobs, ties.method = 'min')
  X.bar.sort <- array(0, dim=c(N, pX))
  for( j in 1:N){
    Y <- (yobs >= y.sort[j])
    X.bar.sort[j,] <- apply( X * Y, 2, sum ) / sum(Y)
  }
  y.sort.diff <- diff( c(0,y.sort) )
  # calculate A, B and d
  A <- B <- array(0, dim=c(pX, pX))
  for( i in 1:N ){ # i <- 1
    Ki <- y.rank[i]
    Xi.aug <- matrix(rep(X[i,], Ki), nrow=Ki, byrow = T)
    Ri <- Xi.aug - X.bar.sort[1:Ki,]
    di <- y.sort.diff[1:Ki]
    I2i <- t(Ri) %*% di %*% X[i,]
    I1i <- ( X[i,] - X.bar.sort[y.rank[i],] ) * delta[i]
    A <- A + I2i
    B <- B + delta[i]*( Ri[Ki,]%*%t(Ri[Ki,]) )
  }
  A <- A/N; B <- B/N

  # calculate empirical V
  V <- solve(A) %*% B %*% solve(A)
  return(V)

}

#============== Influence function for Lin and Ying (1994)'s AR estimator ==============#
#' @export
AR.Influence <- function(yobs,delta,X){
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
  # calculate beta, I2
  A0 <- B0 <- array(0, dim=c(p, p)); d0 <- rep(0,p)
  dd1 <- list()
  dd2 <- array(NA, dim=c(p, N))
  for( i in 1:N ){ # i <- 1
    Ki <- y.rank[i]
    Xi.aug <- matrix(rep(X[i,], Ki), nrow=Ki, byrow = T)
    Ri <- Xi.aug - X.bar.sort[1:Ki,]
    di <- y.sort.diff[1:Ki]
    Rdi <- Ri * sqrt(di)
    I2i <- t(Rdi) %*% Rdi # I2i <- t(Ri) %*% di %*% X[i,]
    I1i <-  ( X[i,] - X.bar.sort[y.rank[i],] ) * delta[i]
    A0 <- A0 + I2i
    d0 <- d0 + I1i
    B0 <- B0 + delta[i]*( Ri[Ki,]%*%t(Ri[Ki,]) )
    # collect dd1 and dd2
    dd1 <- c(dd1,list(I2i)) # all Ai
    dd2[,i] <-  I1i
  }
  A <- A0/N; B <- B0/N; d <- d0/N
  A.inv <- solve(A)
  bet <- solve(A,d)
  # calculate influences
  sumYj <-  sapply(yobs,function(yobsj){sum(yobs >= yobsj)})
  Influs <- array(NA,dim=c(N,p))
  for( i in 1:N ){ # i <- 1
    Xi.aug <- matrix(rep(X[i,], N), nrow=N, byrow = T)
    XXbar <- Xi.aug - X.bar.sort[y.rank,]
    dd3i <- apply(delta * XXbar * (yobs[i] >= yobs) / sumYj,2,sum)
    dd1i <- as.vector(dd1[[i]]%*%bet)
    dd2i <- dd2[,i]
    dd <- - dd3i + dd2i - dd1i
    Influs[i,] <- A.inv%*%dd
  }
  # output
  return(Influs)
}







