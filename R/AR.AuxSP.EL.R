
#==========================================================================#
# AR with tSP as aux using empirical likelihood
#==========================================================================#
#' @title Fit additive risk model with auxiliary subgroup survival rates using EL
#'
#' @description This is a function used to fit the additive risk model
#' in the presence of subgroup survival rates as auxiliary information
#' based on the empirical likelihood framework.
#'
#' @inheritParams AR
#' @param tstar auxiliary time point (only 1 is allowed).
#' @param phi auxiliary survival rates (length(phi)=nrow(G)).
#' @param G indicator matrix - one row indicates one group.
#'
#' @return A list \code{out} representing the fit,
#' which contains the following elements:
#' \code{coef} for estimated coefficients,
#' \code{convergence} for the optimizatoin is converged or not.
#'
#' @references He, J., Li, H., Zhang, S., and Duan, X. (2019).
#' Additive hazards model with auxiliary subgroup survival information.
#' Lifetime Data Analysis, 25(1):128–149.
#'
#' @examples
#' 1
#' @export AR.AuxSP.EL
AR.AuxSP.EL <- function(yobs,delta,X,tstar,phi,G,maxit=30,tol=1e-5){

  ### Specify the dimension ###
  p <- ncol(X);  K <- length(phi)

  ### Initial value for beta and alpha using LY ###
  # fit a add risk model: LY Estimator
  arfit <- AR(yobs, delta, X)
  bet.init <- arfit$coef[,1]
  alp.init <- AR.Ht(arfit,tm=tstar)[,2]
  v.init <- rep(0, p)  # Initial value for v
  xi.init <- rep(0, K)   # Initial value for xi
  # initial values
  rot.init <- c(bet.init, alp.init, v.init, xi.init)

  ### circle to find the solution ###
  maxit <- 30
  tol <- 1e-5
  cytry <- try({
    numit <- 1 # number of iteration
    rot.old <- rot.all <- rot.init
    repeat{ # quasi-Newton-Rapson
      EE  <- AR.AuxSP.EL.ProLik.EE(rot.old,yobs,delta,X,tstar,phi,G)
      dEE <- AR.AuxSP.EL.ProLik.dEE(rot.old,yobs,delta,X,tstar,phi,G)
      tm <- solve(dEE,EE)
      rot.new <- rot.old - tm          # 更新rot
      devv <-  #
        if(sqrt(mean((rot.new-rot.old)^2)) > tol & numit < maxit){
          rot.old <- rot.new
          rot.all <- rbind(rot.all,rot.new)
          numit <- numit + 1
        }else{
          rot.est <- rot.new
          rot.all <- rbind(rot.all,rot.est)
          break
        }
    }
  } ,
  silent=T);  rot.est
  # converge or not
  if( class(cytry) == "try-error" || numit > maxit){
    convergence <- FALSE
  }else{
    convergence <- TRUE
  }
  bet <- rot.est[1:p]; alp <- rot.est[p+1]

  ### Variance Estimation###
  SE <- arfit$coef[,2] # AR.AuxSP.EL.SE(bet,alp,yobs,delta,X,tstar,phi,G)

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

#----------------------------------------------------------------
# score equations for the profiled log-empirical likelihood function
#----------------------------------------------------------------
#' @export
AR.AuxSP.EL.ProLik.EE <- function(par,yobs,delta,X,tstar,phi,G)
{

  ### Prepare ###
  p <- ncol(X)
  K <- length(phi)
  N <- length(yobs)
  bet  <- par[1:p]
  alp <- par[p+1]
  v     <- par[(p+2):(2*p+1)]
  xi    <- par[(2*p+2):(2*p+1+K)]

  ### Prepare values relating to auxiliary information ###
  Psi <- Psidb.noX <- Psida <- array(0, dim=c(N, K))
  sur <- exp(-alp-X%*%bet*tstar)
  for(k in 1:K){
    Psi[,k]   <- (sur- phi[k]) * G[k,] # Psi itself
    Psidb.noX[,k] <- -sur * tstar * G[k,]	 # Psi-dbeta but with out X
    Psida[,k] <- -sur * G[k,]	         # Psi-dalpha
  }
  # times Psi, Psidb, Psida with xi
  Psi.xi <- as.vector( Psi %*% xi ) # Psi-xi for U1, U2 and U4
  Psidb.xi <- X * as.vector(Psidb.noX %*% xi)     # Psidb-xi
  Psida.xi <- Psida %*% xi     # Psida-xi

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
  U4_alpha <- sum( Psida.xi / denom )   # prepare U4 for alpha
  U <- c(U3_beta, U4_alpha, U1_v, U2_xi)   # prepare U
  return(U)

}

#' @export
AR.AuxSP.EL.ProLik.dEE <- function(dpar,yobs,delta,X,tstar,phi,G){
  dh <- 0.00000001
  dl <- length(dpar)
  te<-matrix(0, dl, dl)
  for(i in 1: dl)
  {
    s1<-s2<-dpar
    s1[i]<-s1[i] + dh
    s2[i]<-s2[i] - dh
    te[,i]<- ( AR.AuxSP.EL.ProLik.EE(s1,yobs,delta,X,tstar,phi,G)-
                 AR.AuxSP.EL.ProLik.EE(s2,yobs,delta,X,tstar,phi,G) )/(2*dh)
  }
  return(te)
}




# if(SE.return){
#
#   B1 <- - I2 / N
#   B2 <- t(Psidb) %*% X / N
#   B3 <- as.matrix( apply(Psida,2,mean) )
#   Sigma <- t( ( X - X.bar.sort[y.rank,] ) * delta ) %*%
#     ( X - X.bar.sort[y.rank,] * delta ) / N
#   Sigma.inv <- solve(Sigma)
#   J <- t(Psi) %*% Psi / N
#   J.inv <- solve(J)
#   H <- t(B1) %*% Sigma.inv %*% B1 + t(B2) %*% J.inv %*% B2
#   B2JB3 <- t(B2)%*%J.inv%*%B3
#   beta.VC.matrix <-
#     solve( H - B2JB3 %*% solve(t(B3)%*%J.inv%*%B3) %*%t(B2JB3) )
#   return(sqrt(beta.VC.matrix))
#
# }else{}

