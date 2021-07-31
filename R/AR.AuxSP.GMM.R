

#==========================================================================#
# Aux-AR using GMM
#==========================================================================#
#' @title Fit additive risk model with auxiliary subgroup survival rates using GMM
#'
#' @description This is a function used to fit the additive risk model
#' in the presence of subgroup survival rates as auxiliary information
#' based on the GMM framework.
#'
#' @aliases AR.AuxSP.GMM AR.AuxSP.GMM.fit AR.AuxSP.GMM.LossFunc AR.AuxSP.GMM.SE
#'
#' @inheritParams AR
#' @param aux a list that should contain the auxiliary subgroup survival information.
#' It has three elements:
#' \code{tstar} auxiliary time point (only 1 is allowed);
#' \code{phi} auxiliary survival rates;
#' \code{G} indicator matrix whose one row indicates one group.
#' Note that we should let length(phi)=nrow(G).
#'
#' @return A list \code{out} representing the fit,
#' which contains:
#' \code{coef} for estimated coefficients,
#' \code{convergence} for the optimizatoin is converged or not,
#' et al..
#'
#' @references Shang, W. and Wang, X. (2017).
#' The generalized moment estimation of the additive–multiplicative hazard model
#' with auxiliary survival information.
#'  \emph{Computational Statistics & Data Analysis},  \bold{112}:154–169.
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
#' aux <- list(tstar=0.5,phi=c(0.646,0.732,0.433,0.491),
#'             G=gXcut(toydata1[,c("X1","X2")]))
#' ## fit the model:
#' AR.AuxSP.GMM(Surv(yobs,delta)~X1+X2,data=toydata1,aux=aux)
#'
#' @importFrom survival Surv
#' @export AR.AuxSP.GMM
AR.AuxSP.GMM <- function(formula,data,aux,na.action=na.omit){

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
  arfit <- AR.AuxSP.GMM.fit(yobs,delta,X,tstar,phi,G)
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
  cat("\nAdditive Risk Model with Auxiliary Subgroup Survival Information:\n")
  print(fit$coef)
  cat("\n")

  invisible(fit)

}


#' @export
AR.AuxSP.GMM.fit <- function(yobs,delta,X,tstar,phi,G){

  # Basic elements
  p <- ncol(X); K <- length(phi); N <- length(yobs); N0 <- sum(delta)

  ### GMM procudure ###
  # Get the initial value for beta and alpha
  arfit <- AR.fit(yobs,delta,X) # Lin and Ying's estimator
  bet.LY <- arfit$coef[,1]
  alp.LY <- AR.Ht(arfit,tm=tstar)[,2]
  # # get the initial estimate of tau
  # phi.W <- phi
  sur.init <- exp(-alp.LY-(X%*%bet.LY)*tstar)
  a <- G%*%sur.init / N; b <- apply(G,1,mean)
  phi.W <- as.vector(a/b)   # tSP.AR(tstar,yobs,delta,X,G)
  # The Inverse of the Asymptotic Covariance Matrix
  Wc <- solve(AR.AuxSP.EE.Sigma(bet.LY,alp.LY,yobs,delta,X,tstar,phi.W,G))
  W <- as.matrix(Matrix::bdiag(solve(AR.EE.Sigma(yobs,delta,X)),Wc))
  # Define the object function, initial value
  obj <- AR.AuxSP.GMM.LossFunc
  rot <- c(bet.LY,alp.LY)
  # optimize the object function
  res <- stats::optim(par = rot, fn = obj,method = "BFGS",
                      control = list(maxit = 200, fnscale=1),
                      W=W,yobs=yobs,delta=delta,X=X,tstar=tstar,phi=phi,G=G)
  convergence <- ifelse(res$convergence==0,T,F)
  bet <- res$par[1:p]; alp <- res$par[p+1]

  ### Variance Estimation: Shang and Wang (2017) ###
  SE <- AR.AuxSP.GMM.SE(bet,alp,yobs,delta,X,tstar,phi,G)

  # summary the results
  zvalue <- bet/SE
  pvalue <- 2*(1-pnorm(abs(zvalue)))
  coef <- data.frame(Est=bet, SE=SE, zvalue=zvalue, pvalue=pvalue,
                     row.names=colnames(X))

  # Output the Results
  out <- list(
    sdata=list(yobs=yobs,delta=delta,X=X), # collect my data info
    coef=coef,
    convergence=convergence # converge or not
  )
  return(out)

}

#' @export
AR.AuxSP.GMM.LossFunc <- function(beta_alpha,W,yobs,delta,X,tstar,phi,G){
  # GMM loss function with auxiliary subgroup survival rates
  # From: Shang and Wang (2017)
  # Arguments:
  #   W:  the Weight matrix (optimal one is the inverse of cov(U))

  n_para <- length(beta_alpha)
  bet <- beta_alpha[-n_para]
  alp <- beta_alpha[n_para]
  U <- c( AR.EE(bet, yobs, delta, X),
          AR.AuxSP.EE(bet,alp,yobs,delta,X,tstar,phi,G) )
  return( as.numeric( t(U)%*%W%*% U ) )

}



#' @export
AR.AuxSP.GMM.SE <-  function(bet,alp,yobs,delta,X,tstar,phi,G){
  # Variance Estimation for Shang and Wang's Method

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

  # calculate psi.bet, psi.alp, psi.psi
  sur <- as.vector( exp(-alp-(X%*%bet*tstar)) )
  psi.bet <- - t(t(G)*sur) %*% X * tstar / N
  psi.alp <- as.matrix( - colMeans(t(G)*sur) )
  psi.psi <- AR.AuxSP.EE.Sigma(bet,alp,yobs,delta,X,tstar,phi,G)
  psi.psi.inv <- solve(psi.psi)

  # calculate matrix H1
  H1 <- psi.psi.inv -
    ( psi.psi.inv %*% psi.alp %*% t(psi.alp) %*% psi.psi.inv ) /
    as.numeric( t(psi.alp) %*% psi.psi.inv %*%  psi.alp )

  # calculate Sigma1 and SE
  Sigma1 <-  t(psi.bet) %*% H1 %*% psi.bet
  Sigma <- solve( Sigma.LY.inv + Sigma1 )

  # return
  return( sqrt( diag(Sigma)/N ) )

}


#==========================================================================#
# Auxiliary Information's Estimating Equations and its Asymptotic Covariance Matrix
#==========================================================================#
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
#' @param tstar auxiliary time point (only 1 is allowed).
#' @param phi auxiliary survival rates (length(phi)=nrow(G)).
#' @param G indicator matrix - one row indicates one group.
#'
#' @return A vector \code{Psin}.
#'
#' @examples
#' 1
#' @export AR.AuxSP.EE
AR.AuxSP.EE <- function(bet,alp,yobs,delta,X,tstar,phi,G){
  # prepare
  N <- length(yobs)
  K <- nrow(G)
  sur <- exp(-alp-X%*%bet*tstar)
  # output the averaged estimating equations
  Psin <- ( G%*%sur - apply(G*phi,1,sum) ) / N
  return( as.vector(Psin) )

}

#' @export
AR.AuxSP.EE.Sigma <- function(bet,alp,yobs,delta,X,tstar,phi,G){
  # the asymptotic covariance matrix of: sqrt(n)*Psi(beta, alpha)
  # From: Shang and Wang (2017)

  # prepare
  N <- length(yobs)
  K <- nrow(G)
  sur <- exp(-alp-X%*%bet*tstar)
  Psin2 <- array(0, dim=c(N, K))
  for(k in 1:K){
    Psin2[,k]   <- (sur-phi[k])^2 * G[k,]
  }
  SigmaPsi <- diag( colMeans(Psin2), nrow=K, ncol=K )
  # output
  return(SigmaPsi)

}






