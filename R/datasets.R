
### This file contains all the documentation of datasets ###


#' toydata1 - a toy simulated dataset from additive risk model
#'
#' @description toydata1 is a toy simulated dataset from additive risk model.
#' The censoring rate is approximately 30%. There are two covariates and the
#' true corresponding covariate effects are 0.5 and 0.8, respectively.
#'
#' @format A data frame with 300 rows and 4 variables, including:
#' \describe{
#' \item{yobs}{the observed failure time}
#' \item{delta}{the censoring indicator}
#' \item{X1}{covariate 1 comes from U(0,1) distribution}
#' \item{X2}{covariate 2 comes from b(1,0.5) distribution}
#' }
#'
#' @examples
#' data("toydata1")
#' summary(toydata1$yobs)
#' 1-mean(toydata1$delta)  # censoring rate
#'
#'
#' @keywords datasets
#'
"toydata1"



#' ACTG175 - AIDS Clinical Trials Group Study 175
#'
#' @description ACTG 175 was a randomized clinical trial to compare monotherapy
#' with zidovudine or didanosine with combination therapy with zidovudine and
#' didanosine or zidovudine and zalcitabine in adults infected with the human
#' immunodeficiency virus type I whose CD4 T cell counts were between 200 and
#' 500 per cubic millimeter.
#'
#' @details The variable days contains right-censored time-to-event observations.
#' The data set includes the following post-randomization covariates: CD4 and CD8
#' T cell count at 20\\pm5 weeks and the indicator of whether or not the patient
#' was taken off-treatment before 96\\pm5 weeks.
#'
#'
#' @format This data frame contains the following columns:
#'
#' \describe{
#' \item{pidnum}{patient's ID number}
#' \item{age}{age in years at baseline}
#' }
#'
#' @source Hammer SM, et al. (1996), "A trial comparing nucleoside monotherapy
#' with combination therapy in HIV-infected adults with CD4 cell counts from
#' 200 to 500 per cubic millimeter.", New England Journal of Medicine, 335:1081–1090.
#'
#' @keywords datasets
#'
#' @examples
#'
#' data(ACTG175)
#' yobs  <- ACTG175$days / 365
#' delta <- ACTG175$cens
#' X     <- as.matrix( cbind(ACTG175['cd40']/100,
#'                     ACTG175['treat']) ) # summary(X0)
#'
"ACTG175"




