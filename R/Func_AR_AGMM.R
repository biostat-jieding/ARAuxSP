


#==========================================================================#
# Aux-AR using AGMM: Fit additive risk model with auxiliary subgroup survival rates using AGMM
#==========================================================================#

#============== The main function that fit the AGMM model ==============#
#' @title Fit additive risk model with auxiliary subgroup survival rates using AGMM
#'
#' @description This is a function used to fit the additive risk model
#' in the presence of potentially incomparable subgroup survival rates
#' as auxiliary information based on the GMM framework and penalization.
#' The method is motivated by Chen et al. (2020).
#'
#' @aliases AR.AuxSP.AGMM AR.AuxSP.AGMM.fit AR.AuxSP.AGMM.Loss AR.AuxSP.AGMM.Profile AR.AuxSP.AGMM.Initial
#'
#'
#' @inheritParams AR
#' @param aux indicates the historical summary data. It is a list of lists, and each sub-list represents auxiliary information from a study.
#'   In each study, we have several time points and each time point contains the following four elements
#'   \code{tstar} the time point that the auxiliary information was calculated at;
#'   \code{gfun} a function used to identify the subgroup;
#'   \code{sprob} auxiliary subgroup survival rates for each subgroup at the current time point.
#' @param lambdas user-specified tuning paramters.
#' default is NULL.
#' @param lambdas.num number tuning parameters that will be evaluated at.
#' default is 10.
#' @param threshold a threshold. The absolute values of estimates of tau smaller than the threshold are treated as 0.
#' @param Profile a logical value. Whether use the profiled algorithm or not.
#' @param trace.cv a logical value. Print the calculating process or not
#'
#' @return A list \code{out} representing the fit,
#' which contains the following elements:
#' \code{coef} for estimated coefficients,
#' \code{convergence} for the optimizatoin is converged or not.
#'
#' @references Chen, Z., Ning, J., Shen, Y., and Qin, J. (2020).
#' Combining primary cohort data with external aggregate information
#' without assuming comparability.
#' Biometrics. 77(3): 1024-1036.
#'
#' @examples
#' ## import the data:
#' data("toydata1")
#' ## prepare the auxiliary information:
#' aux <- list(
#'   study1 = list(
#'     time1 = list(
#'       tstar=0.5,
#'       sprob=c(0.646,0.732,0.433,0.491)-c(0.2,0,0,0),
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
#' sol.AGMM <- AR.AuxSP.AGMM(Surv(yobs,delta)~X1+X2,data=toydata1,aux=aux,lambdas.num=5,trace.cv=TRUE)
#'
#' @importFrom survival Surv
#' @export AR.AuxSP.AGMM
AR.AuxSP.AGMM <- function(formula,data,aux,na.action=na.omit,
                          lambdas=NULL,lambdas.num=10,
                          threshold=1e-4,Profile=FALSE,
                          trace.cv=TRUE
){

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
  arfit <- AR.AuxSP.AGMM.fit(yobs,delta,X,aux,
                             lambdas=lambdas,lambdas.num=lambdas.num,
                             threshold=threshold,
                             Profile=Profile,
                             trace=trace.cv)
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

  return(fit)

}

#' @export
AR.AuxSP.AGMM.fit <- function(yobs,delta,X,aux,
                              lambdas=NULL,lambdas.num=7,
                              threshold=1e-4,
                              Profile=FALSE,
                              trace=TRUE
){

  ### Specify the dimension
  p <- ncol(X) # number of covariates in internal data
  N <- length(yobs)   # sample size in internal data
  M <- length(aux) # number of external studies
  L <- sapply(aux,function(studyi){length(studyi)}) # number of time points in each study
  K <- lapply(aux,function(studyi){sapply(studyi,function(timej){length(timej$sprob)})}) # number of groups in each time point and each study

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

  ### prepare initial values for bet and alp
  if(trace){cat("Preparing ...\n")}
  betalp.init <- AR.AuxSP.AGMM.Initial(yobs,delta,X,auxinfo,tstar.unique,
                                       nboot=300,alpha.tau=0.05,maxit=3)
  bet.init <- betalp.init$bet.init
  alp.init <- betalp.init$alp.init

  ### get the initial estimate of each of these sprobs (and taus)
  sur.allG.init <- sapply(1:nrow(auxinfo),function(iG){
    exp(-alp.init[auxinfo[iG,'tstaridx']]-X%*%bet.init*auxinfo[iG,'tstar'])})
  a <- apply(sur.allG.init * t(auxinfo[,-c(1:4)]),2,mean)
  b <- apply(auxinfo[,-c(1:4)],1,mean)
  sprob.init <- as.vector(a/b)
  tau.init <- sprob.init - auxinfo[,'sprob']
  w <- 1/abs((tau.init)^(1))  # omega <- 2

  ### The Inverse of the Asymptotic Covariance Matrix
  auxinfo.alter <- auxinfo; auxinfo.alter[,'sprob'] <- tau.init + auxinfo[,'sprob']
  Phin <- AR.EE(bet.init,yobs,delta,X) # internal data's estimating equation
  Psin <- AR.AuxSP.EE(bet.init,alp.init,yobs,delta,X,auxinfo.alter) # external data's estimating equation
  WPhi <- MASS::ginv(cov(Phin))
  WPsi <- MASS::ginv(cov(Psin))
  W <- as.matrix(Matrix::bdiag(WPhi,WPsi))

  # constraints matrices for maxLik (only need in direct method)
  if(Profile==FALSE){
    tau.up <- rep(0,nrow(auxinfo))
    tau.dw <- rep(0,nrow(auxinfo))
    tau.up[tau.init>=0] <- 1
    tau.dw[tau.init< 0] <- -1
    ineqA <- cbind(
      matrix(0,nrow=length(tstar.unique)+2*nrow(auxinfo),ncol=p),
      rbind(
        cbind(diag(length(tstar.unique)),matrix(0,nrow=length(tstar.unique),ncol=nrow(auxinfo))),
        cbind(matrix(0,nrow=2*nrow(auxinfo),ncol=length(tstar.unique)),rbind(diag(nrow(auxinfo)),-diag(nrow(auxinfo))))
      )
    )
    ineqB <- c(
      rep(0,length(tstar.unique)),
      -pmax(-auxinfo[,'sprob'],tau.dw),
      pmin(1-auxinfo[,'sprob'],tau.up)
    )
  }

  ### GMM: Tuning parameter lam Using BIC from Andrews and Lu (2001)
  # define my lambdas points that will be estimated at
  if(is.null(lambdas)){ # if lambdas have not been defined before
    if(lambdas.num==1){ timeset <- 1 }else{
      timeset <- exp(seq(-2,1,length.out=lambdas.num))
    }
    lambdas <- 2*sqrt(log(nrow(auxinfo)))*N^(-1/2-1/4) * timeset
  }else{ lambdas.num <- length(lambdas) }
  # repeat to try different lams
  bet.all <- array(NA,dim=c(lambdas.num,p))
  alp.all <- array(NA,dim=c(lambdas.num,length(tstar.unique)))
  tau.all <- array(NA,dim=c(lambdas.num,nrow(auxinfo)))
  obj.value.all <- convergence.all <- rep(NA,lambdas.num)
  for( ilam in 1:lambdas.num){ # ilam <- 3
    if(trace){cat("Tuning: Round",ilam,"...\n")}
    # Define the current lam
    lam <- lambdas[ilam]
    # optimize the object function using current lam
    if(Profile==FALSE){
      # solve the problem directly
      res <- stats::optim(
        par=c(bet.init,alp.init,tau.init),fn=AR.AuxSP.AGMM.Loss,
        method="L-BFGS-B",
        lower=c(rep(-Inf,p),rep(0,length(tstar.unique)),pmax(-auxinfo[,'sprob'],tau.dw)),
        upper=c(rep(Inf,p+length(tstar.unique)),pmin(1-auxinfo[,'sprob'],tau.up)),
        control=list(maxit=10000, fnscale=-1),
        lam=lam,w=w,WPhi=WPhi,WPsi=WPsi,yob=yobs,delta=delta,X=X,auxinfo=auxinfo
      );res
      obj.value.all[ilam] <- -res$value
      convergence.all[ilam] <- ifelse(res$convergence==0,T,F)
      bet.all[ilam,] <- res$par[c(1:p)]
      alp.all[ilam,] <- res$par[c((p+1):(p+length(tstar.unique)))]
      tau.all[ilam,] <- res$par[-c(1:(p+length(tstar.unique)))]
    }else if(Profile==TRUE){
      # solve the problem using profiled likelihood
      res <- stats::optim(
        par=c(bet.init,alp.init),fn=AR.AuxSP.AGMM.Profile,
        method="L-BFGS-B",
        lower=c(rep(-Inf,p),rep(0,length(tstar.unique))),
        upper=c(rep(Inf,p+length(tstar.unique))),
        control=list(maxit=10000, fnscale=-1),
        lam=lam,w=w,WPhi=WPhi,WPsi=WPsi,yob=yobs,delta=delta,X=X,
        auxinfo=auxinfo
      );res
      # Combine the estimates
      obj.value.all[ilam] <- -res$value
      convergence.all[ilam] <- ifelse(res$convergence==0,T,F)
      bet.all[ilam,] <- res$par[c(1:p)]
      alp.all[ilam,] <- res$par[-c(1:p)]
      tau.all[ilam,] <- AR.AuxSP.AGMM.Profile(bet_alp=res$par,lam,w,WPhi,WPsi,yobs,delta,X,auxinfo,tau.output=T)
    }
    # Combine the estimates
  }
  # evaluate among all these results
  tau.all.threshold <- ifelse(abs(tau.all)>=threshold,tau.all,0); tau.all.threshold
  IC.all <- sapply(1:lambdas.num, function(ilam){
    obj.value.all[ilam]-lambdas[ilam]*sum(abs(tau.all[ilam,])*w) + # loss part
      sum(tau.all.threshold[ilam,]!=0)*log(N)/N # BIC part
  })
  # Select the minimal BIC
  IC.min.idx <- which.min(IC.all)
  convergence <- convergence.all[IC.min.idx]
  bet <- bet.all[IC.min.idx,]; alp <- alp.all[IC.min.idx,];bet
  tau <- tau.all.threshold[IC.min.idx,]

  ### Variance Estimation ###
  tau.nonzero.idx <- ( tau != 0 ); tau.nonzero.num <- sum(tau.nonzero.idx)
  if(tau.nonzero.num>=(nrow(auxinfo)-1)){
    SE.bet <- arfit$coef[,2]
    SE.tau <- rep(NA,nrow(auxinfo))
  }else{
    auxinfo.alter <- auxinfo; auxinfo.alter[,'sprob'] <- tau + auxinfo[,'sprob']
    Phi.D1 <- cbind(AR.EE.D1(yobs,delta,X),matrix(0,nrow=p,ncol=length(alp)+tau.nonzero.num))
    Psi.D1 <- cbind(AR.AuxSP.EE.D1(bet,alp,yobs,delta,X,auxinfo.alter),
                    -diag(apply(auxinfo.alter[,-c(1:4)],1,mean))[,tau.nonzero.idx,drop=F])
    U.D1 <- rbind(Phi.D1,Psi.D1)
    Un <- cbind(AR.EE(bet,yobs,delta,X),AR.AuxSP.EE(bet,alp,yobs,delta,X,auxinfo.alter))
    Sigma0 <- MASS::ginv(t(U.D1)%*%W%*%U.D1)
    Sigma.Est <- Sigma0%*%t(U.D1)%*%W%*%cov(Un)%*%W%*%U.D1%*%Sigma0
    SE <- sqrt(diag(Sigma.Est)/N)
    SE.bet <- SE[1:p]
    SE.tau <- rep(NA,nrow(auxinfo))
    SE.tau[tau.nonzero.idx] <- SE[-c(1:(p+length(alp)))]
  }

  ### summary the final results ###
  zvalue.bet <- bet/SE.bet
  pvalue.bet <- 2*(1-pnorm(abs(zvalue.bet)))
  coef <- data.frame(Est=bet, SE=SE.bet, zvalue=zvalue.bet, pvalue=pvalue.bet,
                     row.names=colnames(X))
  zvalue.tau <- tau/SE.tau
  pvalue.tau <- 2*(1-pnorm(abs(zvalue.tau)))
  coef.tau <- data.frame(Est=tau, SE=SE.tau, zvalue=zvalue.tau, pvalue=pvalue.tau,
                         row.names=paste("tau",1:nrow(auxinfo),sep=""))

  # Output the Results
  out <- list(
    sdata=list(yobs=yobs,delta=delta,X=X), # collect my data info
    coef=coef,
    coef.tau=coef.tau,
    IC.info=list(
      IC.all = data.frame(lam=lambdas, IC=IC.all),
      IC.min = c(lam=lambdas[IC.min.idx], IC=IC.all[IC.min.idx]) ),
    max_fn = obj.value.all[IC.min.idx],
    convergence=convergence # converge or not
  )
  return(out)
}



### the object function (Define the object function)
#' @export
AR.AuxSP.AGMM.Loss <- function(bet_alp_tau,lam,w,WPhi,WPsi,yobs,delta,X,auxinfo){
  # for maximization
  # prepare
  p <- ncol(X); tau.idx <- c((length(bet_alp_tau)-nrow(auxinfo)+1):length(bet_alp_tau))
  bet <- bet_alp_tau[c(1:p)]; alp <- bet_alp_tau[-c(1:p,tau.idx)]; tau <- bet_alp_tau[tau.idx]
  # calculate final value
  U.Phi <- apply(AR.EE(bet,yobs,delta,X),2,mean)
  auxinfo.alter <- auxinfo; auxinfo.alter[,'sprob'] <- tau + auxinfo[,'sprob']
  U.Psi <- apply(AR.AuxSP.EE(bet,alp,yobs,delta,X,auxinfo.alter),2,mean) # external data's estimating equation
  Q <- as.numeric(t(U.Phi)%*%WPhi%*%U.Phi + t(U.Psi)%*%WPsi%*%U.Psi) + lam*sum(abs(tau)*w)
  return(-Q)
}


### the object function (Define the object function)
#' @export
AR.AuxSP.AGMM.Profile <- function(bet_alp,lam,w,WPhi,WPsi,yobs,delta,X,auxinfo,tau.output=FALSE){
  # for maximization
  # Profiled ALasso GMM loss function with auxiliary subgroup survival rates
  # prepare
  p <- ncol(X)
  bet <- bet_alp[c(1:p)]; alp <- bet_alp[-c(1:p)]
  # calculate the current profiled tau
  sur.allG <- t(sapply(1:nrow(auxinfo),function(iG){
    exp(-alp[auxinfo[iG,'tstaridx']]-X%*%bet*auxinfo[iG,'tstar'])}))
  a <- apply((sur.allG-auxinfo[,'sprob'])*auxinfo[,-c(1:4)],1,mean)
  b <- apply(auxinfo[,-c(1:4)],1,mean)
  tau.sol <- stats::optim(par=rep(0,length(a)),
                          fn=function(tau,lam,a,b,WPsi){as.vector(t(a-b*tau)%*%WPsi%*%(a-b*tau))+lam*sum(abs(tau)*w)},
                          control=list(maxit=5000,fnscale=1),
                          lam=lam,a=a,b=b,WPsi=WPsi)
  if(tau.output){tau <- tau.sol$par; return(tau)}
  # calculate final value
  U.Phi <- apply(AR.EE(bet,yobs,delta,X),2,mean)
  Q <- as.numeric(t(U.Phi)%*%WPhi%*%U.Phi) + tau.sol$value
  return(-Q)
}

#' @export
AR.AuxSP.AGMM.Initial <- function(yobs,delta,X,auxinfo,tstar.unique,
                                  nboot=300,alpha.tau=0.05,maxit=3){

  ### Prepare
  p <- ncol(X) # number of covariates in internal data
  N <- length(yobs)   # sample size in internal data

  ### boot to get the variance of the initial estimators and significant ones
  bet.init.all <- array(NA,dim=c(nboot,p))
  tau.init.all <- array(NA,dim=c(nboot,nrow(auxinfo)))
  for(iboot in 1:nboot){
    idx <- sample(1:N,N,replace=TRUE)

    ### Get the initial value for beta and alpha using Lin and Ying's estimator
    arfit <- AR.fit(yobs[idx],delta[idx],X[idx,])
    bet.init <- arfit$coef[,1]
    alp.init <- AR.Ht(arfit,tm=tstar.unique)[,2]

    ### get the initial estimate of each of these sprobs (and taus)
    sur.allG.init <- sapply(1:nrow(auxinfo),function(iG){
      exp(-alp.init[auxinfo[iG,'tstaridx']]-X[idx,]%*%bet.init*auxinfo[iG,'tstar'])})
    a <- apply(sur.allG.init * t(auxinfo[,-c(1:4)][,idx]),2,mean)
    b <- apply(auxinfo[,-c(1:4)][,idx],1,mean)
    sprob.init <- as.vector(a/b)
    tau.init <- sprob.init - auxinfo[,'sprob']

    ### combine
    bet.init.all[iboot,] <- bet.init
    tau.init.all[iboot,] <- tau.init
  }
  zvalue.tau.init <- apply(tau.init.all,2,mean)/apply(tau.init.all,2,sd)
  pvalue.tau <- round(2*(1-pnorm(abs(zvalue.tau.init))),5); pvalue.tau

  ## Get the initial value for beta and alpha using Lin and Ying's estimator
  arfit <- AR.fit(yobs,delta,X)
  bet.init <- arfit$coef[,1]
  alp.init <- pmax(AR.Ht(arfit,tm=tstar.unique)[,2],0)

  ## Refine the initial value using GMM based on significant taus
  tau.choose <- (pvalue.tau>alpha.tau) # alpha.tau <- 0.1
  if(any(tau.choose)){
    # Iterative GMM for more efficient estimates (maxit=1-One step; maxit=2-Two step)
    bet.old <- bet.init; alp.old <- alp.init
    numit <- 1
    repeat{
      ### The Inverse of the Asymptotic Covariance Matrix ###
      Phin <- AR.EE(bet.old,yobs,delta,X) # internal data's estimating equation
      Psin <- AR.AuxSP.EE(bet.old,alp.old,yobs,delta,X,auxinfo[tau.choose,,drop=F]) # external data's estimating equation
      WPhi <- solve(cov(Phin)); WPsi <- solve(cov(Psin)) # diag(nrow(auxinfo))#
      W <- as.matrix(Matrix::bdiag(WPhi,WPsi))
      ### Optimize the object function (Define initial value first)
      res <- stats::optim(
        par=c(bet.old,alp.old) ,fn=AR.AuxSP.GMM.Loss,
        method="BFGS",
        # method="L-BFGS-B",
        # lower=c(rep(-Inf,p),rep(0,length(tstar.unique))),
        # upper=c(rep(Inf,p+length(tstar.unique))),
        control=list(maxit=1e4,fnscale=-1),
        WPhi=WPhi,WPsi=WPsi,yobs=yobs,delta=delta,X=X,auxinfo=auxinfo[tau.choose,,drop=F]
      )
      convergence <- ifelse(res$convergence==0,T,F)
      bet <- res$par[1:p]; alp <- res$par[-c(1:p)]
      ### update or stop
      if(max(abs(c(bet,alp)-c(bet.old,alp.old))) >= 1e-6 & numit < maxit){
        bet.old <- bet; alp.old <- alp
        numit <- numit + 1
      }else{
        break
      }
    }
    bet.init <- bet
    alp.init <- pmax(AR.Ht(list(sdata=list(yobs=yobs,delta=delta,X=X),coef=data.frame(Est=bet)),tm=tstar.unique)[,2],0)
  }

  ### output
  out <- list(
    bet.init = bet.init,
    alp.init = alp.init
  )

}


