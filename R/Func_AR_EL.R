#==========================================================================#
# AuxSP-AR using EL: Empirical likelihood method
#==========================================================================#

#============== The main function that fits the GMM model ==============#
#' @title Fit additive risk model with auxiliary subgroup survival rates using EL
#'
#' @description This is a function used to fit the additive risk model
#' in the presence of subgroup survival rates as auxiliary information
#' based on the empirical likelihood (EL) framework.
#'
#' @aliases AR.AuxSP.EL AR.AuxSP.EL.fit AR.AuxSP.EL.ProLik.EE AR.AuxSP.EL.ProLik.dEE
#'
#' @inheritParams AR
#' @param aux indicates the historical summary data. It is a list of lists, and each sub-list represents auxiliary information from a study.
#'   In each study, we have several time points and each time point contains the following four elements
#'   \code{tstar} the time point that the auxiliary information was calculated at;
#'   \code{gfun} a function used to identify the subgroup;
#'   \code{sprob} auxiliary subgroup survival rates for each subgroup at the current time point.
#'
#' @return A list \code{out} representing the fit,
#' which contains:
#' \code{coef} for estimated coefficients,
#' \code{convergence} for the optimizatoin is converged or not
#'
#' @references He J, Li H, Zhang S, Duan X. (2019).
#' Additive hazards model with auxiliary subgroup survival information.
#' \emph{Lifetime Data Analysis}, \bold{25}(1): 128-149.
#'
#' @examples
#' ## import the data:
#' data("toydata1")
#' ## prepare the auxiliary information:
#' aux <- list(
#'   study1 = list(
#'     time1 = list(
#'       tstar=0.5,
#'       sprob=c(0.646,0.732,0.433,0.491),
#'       gfunc=function(X){  # group function
#'         rbind(  ( X[,1] >= 0.5 & X[,2] == 0),
#'                 ( X[,1] <  0.5 & X[,2] == 0),
#'                 ( X[,1] >= 0.5 & X[,2] == 1),
#'                 ( X[,1] <  0.5 & X[,2] == 1)
#'         ) * 1
#'       }
#'     )
#'   )
#' )
#' ## fit the model:
#' sol.EL <- AR.AuxSP.EL(Surv(yobs,delta)~X1+X2,data=toydata1,aux=aux)
#'
#' @importFrom survival Surv
#' @export AR.AuxSP.EL
AR.AuxSP.EL <- function(formula,data,aux,na.action=na.omit){

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
  arfit <- AR.AuxSP.EL.fit(yobs,delta,X,aux)
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

  return(fit)

}


#============== The main function ==============#
AR.AuxSP.EL.fit <- function(yobs,delta,X,aux,maxit=30){

  ### Specify the dimension
  p <- ncol(X) # number of covariates in internal data
  N <- length(yobs)   # sample size in internal data

  ### prepare the auxinfo matrix
  tstar.unique <- unique(unlist(lapply(aux,function(studyi){sapply(studyi,function(timej){timej$tstar})})))
  auxinfo <- as.numeric()
  for(istudy in 1:length(aux)){ # istudy <- 2; itime <- 1
    for(itime in 1:length(aux[[istudy]])){
      aux.c <- aux[[istudy]][[itime]]
      K.c <- length(aux[[istudy]][[itime]]$sprob)
      auxinfo <- rbind(
        auxinfo,
        cbind(rep(istudy,K.c),
              rep(aux.c$tstar,K.c),
              rep(which(tstar.unique==aux.c$tstar),K.c),
              aux.c$sprob,
              aux.c$gfunc(X)))
    }
  }
  colnames(auxinfo) <- c('study','tstar','tstaridx','sprob',paste('ind',1:N,sep="")) # rename this matrix

  ### Initial value for beta and alpha using LY ###
  arfit <- AR.fit(yobs, delta, X)
  bet.init <- arfit$coef[,1]
  alp.init <- pmax(AR.Ht(arfit,tm=tstar.unique)[,2],0)
  v.init <- rep(0, p)  # Initial value for v
  xi.init <- rep(0, nrow(auxinfo))   # Initial value for xi

  # solve the estimating equations directly
  cytry <- try({
    numit <- 1 # number of iteration
    rot.old <- rot.all <- c(bet.init, alp.init, v.init, xi.init)
    repeat{ # Newton-Rapson
      EE  <- AR.AuxSP.EL.ProLik.EE(rot.old,yobs,delta,X,auxinfo)
      dEE <- AR.AuxSP.EL.ProLik.dEE(rot.old,yobs,delta,X,auxinfo)
      tm <- as.vector(MASS::ginv(dEE)%*%EE) # solve(dEE,EE) #
      rot.new <- rot.old - tm
      if(sqrt(mean((rot.new-rot.old)^2)) > 1e-6 & numit < maxit){
        rot.old <- rot.new
        rot.all <- rbind(rot.all,rot.new)
        numit <- numit + 1
      }else{
        rot.all <- rbind(rot.all,rot.new)
        break
      }
    }
  } ,
  silent=T);  rot.new
  convergence <- ifelse(class(cytry) == "try-error" || numit > maxit,F,T)
  bet <- rot.new[1:p]; alp <- rot.new[(p+1):(p+length(tstar.unique))]

  ### Variance Estimation###
  SE <- sqrt(diag(AR.AuxSP.EL.VCOV(bet,alp,yobs,delta,X,auxinfo))/N)

  ### summary the results ###
  zvalue <- bet/SE
  pvalue <- 2*(1-pnorm(abs(zvalue)))
  coef <- data.frame(Est=bet, SE=SE, zvalue=zvalue, pvalue=pvalue,
                     row.names=colnames(X))

  ### Output the Results ###
  out <- list(
    sdata=list(yobs=yobs,delta=delta,X=X), # collect my data info
    coef=coef,
    convergence=convergence # converge or not
  )
  return(out)

}


#==== score equations for the profiled log-empirical likelihood function ====#
#' @export
AR.AuxSP.EL.ProLik.EE <- function(par,yobs,delta,X,auxinfo)
{

  ### Prepare ###
  p <- ncol(X) # number of covariates in internal data
  N <- length(yobs)   # sample size in internal data
  sumK <- nrow(auxinfo) # number of Groups
  ntime <- length(unique(auxinfo[,'tstar']))
  bet <- par[1:p]
  alp <- par[(p+1):(p+ntime)]
  v   <- par[(p+ntime+1):(2*p+ntime)]
  xi  <- par[-(1:(2*p+ntime))]

  ### Prepare values relating to auxiliary information ###
  sur.allG <- t(sapply(1:nrow(auxinfo),function(iG){
    exp(-alp[auxinfo[iG,'tstaridx']]-X%*%bet*auxinfo[iG,'tstar'])}))
  Psi <- t((sur.allG-auxinfo[,'sprob'])*auxinfo[,-c(1:4)])
  Psi.xi <- as.vector( Psi %*% xi )
  Psidb.xi <- - X * as.vector(t(sur.allG*auxinfo[,-c(1:4)]*auxinfo[,'tstar']) %*% xi)
  Psida.xi <- - sapply(1:length(alp),function(itstaridx){
    t(sur.allG*auxinfo[,-c(1:4)]*(auxinfo[,'tstaridx']==itstaridx))%*%xi})

  ###  Prepare values relating to  LY ###
  # sort and rank(Ki)
  y.sort <- sort( yobs )
  y.rank <- rank( yobs, ties.method = 'min')
  # X.bar matrix [ N×p ]
  X.bar.sort <- array(0, dim=c(N, p))
  for( j in 1:N){ # j <- 10; i <- 1
    Y <- (yobs >= y.sort[j])
    X.bar.sort[j,] <- apply( X * Y, 2, sum ) / sum(Y)
  }
  y.sort.diff <- diff( c(0,y.sort) )
  # calculate A, B and d
  W <- Wdb.v <- array(0, dim=c(N, p))
  for( i in 1:N ){ # i
    Ki <- y.rank[i]
    Xi.aug <- matrix(rep(X[i,], Ki), nrow=Ki, byrow = T)
    Ri <- Xi.aug - X.bar.sort[1:Ki,]
    di <- y.sort.diff[1:Ki]
    I2i <- t(Ri) %*% di %*% X[i,]
    I1i <- ( X[i,] - X.bar.sort[y.rank[i],] ) * delta[i]
    W[i,] <- I1i - I2i %*% bet
    Wdb.v[i,] <- - t(I2i) %*% v # W-dbeta * v for each i
  }
  W.v <- as.vector( W %*% v ) # W * v

  ### Estimating Equation U ###
  denom <- 1 + W.v + Psi.xi
  U1_v <- apply( W / denom, 2, sum )     # prepare U1 for v
  U2_xi <- apply( Psi / denom , 2, sum )   # prepare U2 for xi
  U3_beta  <- apply( ( Wdb.v + Psidb.xi ) / denom , 2, sum )   # prepare U3 for beta
  U4_alpha <- apply( Psida.xi / denom, 2, sum )   # prepare U4 for alpha
  U <- c(U3_beta, U4_alpha, U1_v/N, U2_xi/N)   # prepare U
  return(U)

}

#==== numerical second derivative for the profiled log-empirical likelihood function ====#
#' @export
AR.AuxSP.EL.ProLik.dEE <- function(dpar,yobs,delta,X,auxinfo){
  dh <- 0.00000001
  dl <- length(dpar)
  te<-matrix(0, dl, dl)
  for(i in 1: dl)
  {
    s1 <- s2 <- dpar
    s1[i] <- s1[i] + dh
    s2[i] <- s2[i] - dh
    te[,i]<- ( AR.AuxSP.EL.ProLik.EE(s1,yobs,delta,X,auxinfo)-
                 AR.AuxSP.EL.ProLik.EE(s2,yobs,delta,X,auxinfo) ) / (2*dh)
  }
  return(te)
}





#==== variance-covariance matrix for the method ====#
#' @export
AR.AuxSP.EL.VCOV <- function(bet,alp,yobs,delta,X,auxinfo)
{

  ### Prepare ###
  p <- ncol(X) # number of covariates in internal data
  N <- length(yobs)   # sample size in internal data
  sumK <- nrow(auxinfo) # number of Groups
  ntime <- length(unique(auxinfo[,'tstar']))

  ### Prepare values relating to auxiliary information ###
  sur.allG <- t(sapply(1:nrow(auxinfo),function(iG){
    exp(-alp[auxinfo[iG,'tstaridx']]-X%*%bet*auxinfo[iG,'tstar'])}))
  Psi <- t((sur.allG-auxinfo[,'sprob'])*auxinfo[,-c(1:4)])
  B2 <- - (sur.allG*auxinfo[,-c(1:4)]*auxinfo[,'tstar']) %*% X / N
  B3 <- - sapply(1:length(alp),function(itstaridx){
    apply(sur.allG*auxinfo[,-c(1:4)]*(auxinfo[,'tstaridx']==itstaridx),1,sum) }) / N
  J <- t(Psi) %*% Psi / N

  ###  Prepare values relating to  LY ###
  # sort and rank(Ki)
  y.sort <- sort( yobs )
  y.rank <- rank( yobs, ties.method = 'min')
  # X.bar matrix [ N×p ]
  X.bar.sort <- array(0, dim=c(N, p))
  for( j in 1:N){ # j <- 10; i <- 1
    Y <- (yobs >= y.sort[j])
    X.bar.sort[j,] <- apply( X * Y, 2, sum ) / sum(Y)
  }
  y.sort.diff <- diff( c(0,y.sort) )
  # calculate A, B and d
  nS <- nB1 <- array(0, dim=c(p, p))
  for( i in 1:N ){ # i
    Ki <- y.rank[i]
    Xi.aug <- matrix(rep(X[i,], Ki), nrow=Ki, byrow = T)
    Ri <- Xi.aug - X.bar.sort[1:Ki,]
    di <- y.sort.diff[1:Ki]
    I2i <- t(Ri) %*% di %*% X[i,]
    nB1 <- nB1 - I2i
    nS <- nS + delta[i]*( Ri[Ki,]%*%t(Ri[Ki,]) )
  }
  S <- nS/N
  B1 <- nB1/N
  H <- t(B1)%*%MASS::ginv(S)%*%B1 + t(B2)%*%MASS::ginv(J)%*%B2

  ### calculate VCOV ###
  VCOV <- MASS::ginv( H-t(B2)%*%MASS::ginv(J)%*%B3 %*% MASS::ginv(t(B3)%*%MASS::ginv(J)%*%B3) %*% t(B3)%*%MASS::ginv(J)%*%B2 )

  ### output ###
  return(VCOV)

}

