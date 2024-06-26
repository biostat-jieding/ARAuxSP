% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Func_AR_GMM.R
\name{AR.AuxSP.GMM}
\alias{AR.AuxSP.GMM}
\alias{AR.AuxSP.GMM.fit}
\alias{AR.AuxSP.GMM.Loss}
\title{Fit additive risk model with auxiliary subgroup survival rates using GMM}
\usage{
AR.AuxSP.GMM(formula, data, aux, na.action = na.omit)
}
\arguments{
\item{formula}{formula object, with the response on the left of a ~ operator,
and the terms on the right. The response must be a survival object as
returned by the Surv function.}

\item{data}{a data.frame in which to interpret the variables named in the formula.}

\item{aux}{indicates the historical summary data. It is a list of lists, and each sub-list represents auxiliary information from a study.
In each study, we have several time points and each time point contains the following four elements
\code{tstar} the time point that the auxiliary information was calculated at;
\code{gfun} a function used to identify the subgroup;
\code{sprob} auxiliary subgroup survival rates for each subgroup at the current time point.}

\item{na.action}{a missing-data filter function.
This is applied to the model.frame after any subset argument has been used.}
}
\value{
A list \code{out} representing the fit,
which contains:
\code{coef} for estimated coefficients,
\code{convergence} for the optimizatoin is converged or not
}
\description{
This is a function used to fit the additive risk model
in the presence of subgroup survival rates as auxiliary information
based on the GMM framework.
}
\examples{
## import the data:
data("toydata1")
## prepare the auxiliary information:
aux <- list(
  study1 = list(
    time1 = list(
      tstar=0.5,
      sprob=c(0.646,0.732,0.433,0.491),
      gfunc=function(X){  # group function
        rbind(  ( X[,1] >= 0.5 & X[,2] == 0),
                ( X[,1] <  0.5 & X[,2] == 0),
                ( X[,1] >= 0.5 & X[,2] == 1),
                ( X[,1] <  0.5 & X[,2] == 1)
        ) * 1
      }
    )
  )
)
## fit the model:
sol.GMM <- AR.AuxSP.GMM(Surv(yobs,delta)~X1+X2,data=toydata1,aux=aux)

}
\references{
Shang, W. and Wang, X. (2017).
The generalized moment estimation of the additive–multiplicative hazard model
with auxiliary survival information.
\emph{Computational Statistics & Data Analysis},  \bold{112}:154???169.
}
