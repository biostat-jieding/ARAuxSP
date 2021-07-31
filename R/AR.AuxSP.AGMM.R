

#' @title Fit additive risk model with auxiliary subgroup survival rates using AGMM
#'
#' @description This is a function used to fit the additive risk model
#' in the presence of potentially incomparable subgroup survival rates
#' as auxiliary information based on the GMM framework and penalization.
#' The method is motivated by Chen et al. (2020).
#'
#' @aliases AR.AuxSP.AGMM AR.AuxSP.AGMM.fit
#'
#'
#' @inheritParams AR
#' @param aux a list that should contain the auxiliary subgroup survival information.
#' It has three elements:
#' \code{tstar} auxiliary time point (only 1 is allowed);
#' \code{phi} auxiliary survival rates;
#' \code{G} indicator matrix whose one row indicates one group.
#' Note that we should let length(phi)=nrow(G).
#' @param lambdas user-specified tuning paramters.
#' default is NULL.
#' @param lambdas.num number tuning parameters that will be evaluated at.
#' default is 10.
#'
#' @return A list \code{out} representing the fit,
#' which contains the following elements:
#' \code{coef} for estimated coefficients,
#' \code{convergence} for the optimizatoin is converged or not.
#'
#' @references Chen, Z., Ning, J., Shen, Y., and Qin, J. (2020).
#' Combining primary cohort data with external aggregate information
#' without assuming comparability.
#' Biometrics. doi: https://doi.org/10.1111/biom.13356.
#'
#' @examples
#' ## import the data:
#' data("toydata1")
#' ## prepare the auxiliary information:
#' gXcut <- function(X){  # group function
#' rbind(  ( X[,1] >= 0.5 & X[,2] == 0),
#'         ( X[,1] <  0.5 & X[,2] == 0),
#'         ( X[,1] >= 0.5 & X[,2] == 1),
#'         ( X[,1] <  0.5 & X[,2] == 1)
#' ) * 1
#' }
#' phi.true <- c(0.646,0.732,0.433,0.491)
#' tau <- c(0.2,0,0,0)
#' phi <- phi.true - tau
#' aux <- list(tstar=0.5,phi=phi,
#'             G=gXcut(toydata1[,c("X1","X2")]))
#' ## fit the model:
#' AR.AuxSP.AGMM(Surv(yobs,delta)~X1+X2,data=toydata1,aux=aux,trace.cv=TRUE)
#'
#' @importFrom survival Surv
#' @export AR.AuxSP.AGMM
AR.AuxSP.AGMM <- function(formula,data,aux,na.action=na.omit,
                          lambdas=NULL,lambdas.num=10,
                          trace.cv=FALSE){

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
  tstar <- aux$tstar
  phi <- aux$phi
  G <- aux$G

  ## fit the model
  arfit <- AR.AuxSP.AGMM.fit(yobs,delta,X,tstar,phi,G,
                             lambdas=lambdas,lambdas.num=lambdas.num,
                             trace.cv=trace.cv)
  rownames( arfit$coef ) <- vars.name


  ## tidy the results
  fit <- list()
  class(fit) <- c("AR")
  fit$sdata <- arfit$sdata
  fit$call <- call
  fit$formula <- formula
  fit$vars.num <- vars.num
  fit$coef <- arfit$coef
  fit$coef.tau <- arfit$coef.tau
  fit$convergence <- arfit$convergence
  fit$IC.info <- arfit$IC.info

  ## print the fitting results
  if (!is.null(cl <- fit$call)) {
    cat("Call:\n")
    dput(cl)
  }
  cat("\nAdditive Risk Model with Auxiliary Subgroup Survival Information",
      "(Existing Potential Incomparability):\n")
  print(fit$coef)
  cat("\n")
  cat("Sparse Estimation Results:\n")
  print(fit$coef.tau)

  invisible(fit)

}

#' @export
AR.AuxSP.AGMM.fit <- function(yobs,delta,X,tstar,phi,G,
                          lambdas=NULL,lambdas.num=10,
                          trace.cv=FALSE){
  # Additive Risk Model with auxiliary subgroup survival rate
  #   using GMM and Adaptive lasso (profile NM algorithm)
  # From: DJ - Motivated by Chen et al. (2020) and Cheng (2013)
  # Arguments:
  #   tstar:  auxiliary time point (only 1)
  #   phi:    auxiliary survival rates [length(phi)=nrow(G)]
  #   G:      indicator matrix - one row indicates one group
  # Output:
  #   coef:
  #   coef.tau:
  #   IC.info:
  # How to implment:
  #   1. AR.AuxSP.AGMM(yobs,delta,X,tstar,phi,G)
  #   2. AR.AuxSP.AGMM(yobs,delta,X,tstar,phi,G,lambdas=0.05)
  #   3. AR.AuxSP.AGMM(yobs,delta,X,tstar,phi,G,lambdas.num=10)

  # Basic elements
  p <- ncol(X); K <- length(phi)
  N <- length(yobs); N0 <- sum(delta)

  ### Preparations before GMM procudure ###
  # Get the initial value for beta and alpha
  arfit <- AR.fit(yobs,delta,X) # Lin and Ying's estimator
  bet.LY <- arfit$coef[,1]
  alp.LY <- AR.Ht(arfit,tm=tstar)[,2]
  # get the initial estimate of tau
  # if(tauinit=="AB"){
  sur.init <- exp(-alp.LY-(X%*%bet.LY)*tstar)
  a <- ( G%*%sur.init - apply(G*phi,1,sum) ) / N
  b <- apply(G,1,mean)
  tau.init <- as.vector(a/b)
  # }else if(tauinit=="KM"){
  #   tau.init <- KM.St.Sub(tstar,yobs,delta,G) - phi
  # }
  # tau.init <- c(0.1,0,0,0)

  # set the weights for adaptive lasso
  w <- 1/abs(tau.init)
  # The Inverse of the Asymptotic Covariance Matrix
  Wc <- solve(AR.AuxSP.EE.Sigma(bet.LY,alp.LY,yobs,delta,X,tstar,phi+tau.init,G))
  W <- as.matrix(Matrix::bdiag(solve(AR.EE.Sigma(yobs,delta,X)),Wc))
  # Define the object function, initial value
  obj <- AR.AuxSP.AGMM.LossFunc
  rot <- c(bet.LY,alp.LY)

  ### GMM: Tuning parameter lam Using BIC from Andrews and Lu (2001) ###
  # define my lambdas points that will be estimated at
  if(is.null(lambdas)){ # if lambdas have not been defined before
    lambdas <- 2*sqrt(log(K))*N^(-1/2-1/4) * exp(seq(-3,3,length.out=lambdas.num))
  }else{
    lambdas.num <- length(lambdas)
  }
  # repeat to try different lams
  Res.all <- list()
  for( ilam in 1:lambdas.num){ # ilam <- 1
    if(trace.cv){cat("CV: Round",ilam,"...\n")}
    # Define the current lam
    lam <- lambdas[ilam]
    # optimize the object function using current lam
    res <- stats::optim(par = rot, fn = obj,method = "Nelder-Mead",
                        control = list(maxit = 500, fnscale=1),
                        lam=lam,w=w,W=W,Wc=Wc,yobs=yobs,
                        delta=delta,X=X,tstar=tstar,phi=phi,G=G)
    tau.c <- AR.AuxSP.AGMM.Tau.Profile(res$par[1:p],res$par[p+1],lam,w,Wc,X,tstar,phi,G)
    # Combine the estimates
    res$par <- c(res$par,tau.c)
    Res.all[[ilam]] <- res
  }
  # evaluate among all these results
  convergence.all <- sapply(Res.all, function(temp){ifelse(temp$convergence==0,T,F)})
  IC.all <- sapply(Res.all, function(temp){
    bet_alp.i <- temp$par[1:(p+1)]; tau.i <- temp$par[(p+2):(p+1+K)]
    AR.AuxSP.GMM.LossFunc(bet_alp.i,W,yobs,delta,X,tstar,phi+tau.i,G) +
      sum( abs(tau.i)>1e-6 ) * log(N)/N})
  # Select the minimal BIC
  IC.min.idx <- which.min(IC.all)
  res <- Res.all[[IC.min.idx]]; convergence <- convergence.all[IC.min.idx]
  bet <- res$par[1:p]; alp <- res$par[p+1]; tau <- res$par[(p+2):(p+1+K)]

  ### Variance Estimation: based on asymptotic var-cov matrix ###
  SE <- AR.AuxSP.AGMM.SE(bet,alp,tau,yobs,delta,X,tstar,phi,G)
  SE.bet <- SE$SE.bet
  SE.tau <- SE$SE.tau

  ### summary the final results ###
  zvalue.bet <- bet/SE.bet
  pvalue.bet <- 2*(1-pnorm(abs(zvalue.bet)))
  coef <- data.frame(Est=bet, SE=SE.bet, zvalue=zvalue.bet, pvalue=pvalue.bet,
                     row.names=colnames(X))
  zvalue.tau <- tau/SE.tau
  pvalue.tau <- 2*(1-pnorm(abs(zvalue.tau)))
  coef.tau <- data.frame(Est=tau, SE=SE.tau, zvalue=zvalue.tau, pvalue=pvalue.tau,
                         row.names=paste("tau",1:K,sep=""))

  ### Output the Results ###
  out <- list(
    sdata=list(yobs=yobs,delta=delta,X=X), # collect my data info
    coef=coef,
    coef.tau=coef.tau,
    IC.info=list(
      IC.all = data.frame(lam=lambdas, IC=IC.all),
      IC.min = c(lam=lambdas[IC.min.idx], IC=IC.all[IC.min.idx]) ),
    convergence=convergence # converge or not
  )
  return(out)

}

# Profiled ALasso GMM loss function with auxiliary subgroup survival rates
AR.AuxSP.AGMM.LossFunc <- function(beta_alpha,lam,w,W,Wc,yobs,delta,X,tstar,phi,G){
  # Profiled ALasso GMM loss function with auxiliary subgroup survival rates
  # From: DJ
  # Arguments:
  #   lam:  tuning parameter
  #   w:    adaptive weights
  #   W:    the Weight matrix (optimal one is the inverse of cov(U))
  #   Wc:   inverse of cov(AR.AuxSP.EE)

  # prepare
  N <- length(yobs)
  n_para <- length(beta_alpha)
  K <- length(phi)
  bet <- beta_alpha[-n_para]
  alp <- beta_alpha[n_para]
  # get the profiled tau
  tau <- AR.AuxSP.AGMM.Tau.Profile(bet,alp,lam,w,Wc,X,tstar,phi,G)

  # get U and Q
  U <- c( AR.EE(bet,yobs,delta,X),
          AR.AuxSP.EE(bet,alp,yobs,delta,X,tstar,phi+tau,G) )
  Q <- as.numeric( t(U)%*%W%*% U ) + lam * sum( abs(tau)*w )

  return( Q )
}

AR.AuxSP.AGMM.Tau.Profile <- function(bet,alp,lam,w,Wc,X,tstar,phi,G){
  # Profile Tau given beta and alpha (based on adaptive lasso)
  N <- nrow(X)
  K <- length(phi)
  sur <- exp(-alp-(X%*%bet)*tstar)
  a <- as.numeric( ( G%*%sur - apply(G*phi,1,sum) ) / N )
  b <- apply(G,1,mean)
  tau <- rep(0,K)
  for(k in 1:K){
    tau[k] <- Soft.Threshold( 2*a[k]*b[k]*Wc[k,k],lam*w[k] ) / (2*Wc[k,k]*b[k]^2)
  }
  return(tau)
}

AR.AuxSP.AGMM.SE <-  function(bet,alp,tau,yobs,delta,X,tstar,phi,G){
  # Variance Estimation for Shang and Wang's Method with Adaptibe lasso penalty

  # Prepare
  N <- length(yobs); p <- ncol(X); K <- nrow(G)
  y.sort <- sort( yobs )
  y.rank <- rank( yobs, ties.method = 'min') # Ki
  X.bar.sort <- array(0, dim=c(N, p))
  for( j in 1:N){
    Y <- (yobs >= y.sort[j])
    X.bar.sort[j,] <- apply( X * Y, 2, sum ) / sum(Y)
  }
  y.sort.diff <- c(y.sort[1], diff( y.sort ) )
  # calculate A, B and Sigma.LY.inv
  A0 <- B0 <- array(0, dim=c(p, p))
  for( i in 1:N ){ # i
    Ki <- y.rank[i]
    Xi.aug <- matrix(rep(X[i,], Ki), nrow=Ki, byrow = T)
    Ri <- Xi.aug - X.bar.sort[1:Ki,]
    di <- y.sort.diff[1:Ki]
    I2i <- t(Ri) %*% di %*% X[i,]
    A0 <- A0 + I2i
    B0 <- B0 + delta[i]*( Ri[Ki,]%*%t(Ri[Ki,]) )
  }
  A <- A0/N; B <- B0/N # A and B comes from Lin and Ying (1994)
  Sigma.LY.inv <- A %*% solve(B) %*% A

  # calculate psi.bet, psi.alp, psi.psi and H1
  sur <- as.vector( exp(-alp-(X%*%bet*tstar)) )
  psi.bet <- - t(t(G)*sur) %*% X * tstar / N
  psi.alp <- as.matrix( - colMeans(t(G)*sur) )
  tau.nonzero.idx <- ( abs(tau) > 1e-6 )
  psi.psi <- AR.AuxSP.EE.Sigma(bet,alp,yobs,delta,X,tstar,phi+tau,G)
  psi.psi.inv <- solve(psi.psi)
  H1 <- psi.psi.inv -
    ( psi.psi.inv %*% psi.alp %*% t(psi.alp) %*% psi.psi.inv ) /
    as.numeric( t(psi.alp) %*% psi.psi.inv %*%  psi.alp )

  # calculate matrix H0
  SE.tau <- rep(NA, K)
  if(any(tau.nonzero.idx)){
    psi.tau <- - diag(apply(G,1,mean))[,tau.nonzero.idx,drop=F]
    H0 <- H1 - H1 %*% psi.tau %*% solve( t(psi.tau) %*%H1%*%psi.tau ) %*% t(psi.tau) %*% H1
    # SE.tau
    M11 <- rbind( cbind( Sigma.LY.inv+t(psi.bet)%*%psi.psi.inv%*%psi.bet,
                         t(psi.bet)%*%psi.psi.inv%*%psi.alp ),
                  cbind( t(psi.alp)%*%psi.psi.inv%*%psi.bet,
                         t(psi.alp)%*%psi.psi.inv%*%psi.alp ) )
    M12 <- rbind( t(psi.bet)%*%psi.psi.inv%*%psi.tau, t(psi.alp)%*%psi.psi.inv%*%psi.tau )
    M22 <- t(psi.tau)%*%psi.psi.inv%*%psi.tau
    Sigma.tau <- solve( M22 - t(M12) %*% solve(M11) %*% M12 )
    SE.tau[tau.nonzero.idx] <- sqrt( diag(Sigma.tau)/N )
  }else{
    H0 <- H1
  }

  # calculate Sigma1, SE.bet
  Sigma0 <-  t(psi.bet) %*% H0 %*% psi.bet
  Sigma.bet <- solve( Sigma.LY.inv + Sigma0 )
  SE.bet <- sqrt( diag(Sigma.bet)/N )

  # return
  return( list(SE.bet=SE.bet, SE.tau=SE.tau) )

}



