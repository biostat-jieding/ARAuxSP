

#==========================================================================#
# Aux-AR using GMM: Fit additive risk model with auxiliary subgroup survival rates using GMM
#==========================================================================#

#============== The main function that fits the GMM model ==============#
#' @title Fit additive risk model with auxiliary subgroup survival rates using GMM
#'
#' @description This is a function used to fit the additive risk model
#' in the presence of subgroup survival rates as auxiliary information
#' based on the GMM framework.
#'
#' @aliases AR.AuxSP.GMM AR.AuxSP.GMM.fit AR.AuxSP.GMM.Loss
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
#' @references Shang, W. and Wang, X. (2017).
#' The generalized moment estimation of the additive–multiplicative hazard model
#' with auxiliary survival information.
#'  \emph{Computational Statistics & Data Analysis},  \bold{112}:154???169.
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
#' sol.GMM <- AR.AuxSP.GMM(Surv(yobs,delta)~X1+X2,data=toydata1,aux=aux)
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

  ## fit the model
  arfit <- AR.AuxSP.GMM.fit(yobs,delta,X,aux)
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

#' @export
AR.AuxSP.GMM.fit <- function(yobs,delta,X,aux,maxit=3){

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

  ### Get the initial value for beta and alpha using Lin and Ying's estimator
  arfit <- AR.fit(yobs,delta,X)
  bet.init <- arfit$coef[,1]
  alp.init <- pmax(AR.Ht(arfit,tm=tstar.unique)[,2],0)

  ### Iterative GMM for more efficient estimates (maxit=1-One step; maxit=2-Two step)
  bet.old <- bet.init; alp.old <- alp.init
  numit <- 1
  repeat{
    ### The Inverse of the Asymptotic Covariance Matrix ###
    Phin <- AR.EE(bet.old,yobs,delta,X) # internal data's estimating equation
    Psin <- AR.AuxSP.EE(bet.old,alp.old,yobs,delta,X,auxinfo) # external data's estimating equation
    WPhi <- MASS::ginv(cov(Phin)); WPsi <- MASS::ginv(cov(Psin)) # diag(nrow(auxinfo))#
    W <- as.matrix(Matrix::bdiag(WPhi,WPsi))
    ### Optimize the object function (Define initial value first)
    res <- stats::optim(par=c(
      bet.old,alp.old) ,fn=AR.AuxSP.GMM.Loss,
      method="BFGS",
      # method="L-BFGS-B",
      # lower=c(rep(-Inf,p),rep(0,length(tstar.unique))),
      # upper=c(rep(Inf,p+length(tstar.unique))),
      control=list(maxit=10000,fnscale=-1),
      WPhi=WPhi,WPsi=WPsi,yobs=yobs,delta=delta,X=X,auxinfo=auxinfo
    ); res
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

  ### Variance Estimation - using standard formula derived in Hansen (1982) ###
  Phi.D1 <- cbind(AR.EE.D1(yobs,delta,X),matrix(0,nrow=p,ncol=length(alp)))
  Psi.D1 <- AR.AuxSP.EE.D1(bet,alp,yobs,delta,X,auxinfo)
  U.D1 <- rbind(Phi.D1,Psi.D1)
  Un <- cbind(AR.EE(bet,yobs,delta,X),AR.AuxSP.EE(bet,alp,yobs,delta,X,auxinfo))
  Sigma0 <- MASS::ginv(t(U.D1)%*%W%*%U.D1)
  Sigma.Est <- Sigma0%*%t(U.D1)%*%W%*%cov(Un)%*%W%*%U.D1%*%Sigma0
  SE <- sqrt(diag(Sigma.Est)[1:p]/N)

  ### summary the results
  zvalue <- bet/SE
  pvalue <- 2*(1-pnorm(abs(zvalue)))
  coef <- data.frame(Est=bet, SE=SE, zvalue=zvalue, pvalue=pvalue,
                     row.names=colnames(X))

  ### J test: Sargan Hansen J-test to check the validity of overidentifying restrictions
  Jtest.value <- -N*res$value
  dfreedom <- p+nrow(auxinfo)
  pvalue.Jtest <- 1-pchisq(Jtest.value,df=dfreedom)
  res.Jtest <- data.frame(Jvalue=Jtest.value,df=dfreedom,pvalue=pvalue.Jtest)

  ### Output the Results
  out <- list(
    sdata=list(yobs=yobs,delta=delta,X=X), # collect my data info
    coef=coef,
    alp=alp,
    Jtest=res.Jtest,
    convergence=convergence # converge or not
  )
  return(out)

}

#' @export
AR.AuxSP.GMM.Loss <- function(bet_alp,WPhi,WPsi,yobs,delta,X,auxinfo){ # AR.AuxSP.GMM.LossFunc
  # this loss is default for maximization

  # prepare
  p <- ncol(X)
  bet <- bet_alp[c(1:p)]
  alp <- bet_alp[-c(1:p)]
  # value of estimating equations
  U.Phi <- apply(AR.EE(bet,yobs,delta,X),2,mean) # internal data's estimating equation
  U.Psi <- apply(AR.AuxSP.EE(bet,alp,yobs,delta,X,auxinfo),2,mean)
  U <- c(U.Phi,U.Psi)
  # output
  return( - as.numeric( t(U.Phi)%*%WPhi%*% U.Phi + t(U.Psi)%*%WPsi%*% U.Psi ) )

}
