##' SAFE Feature Screening for Lasso and Adaptive Lasso Models
##' @title SurvELM ELMCox
##' @param x  The covariates(predictor variables) of training data.
##' @param y  Survival time and censored status of training data. First column be "time" and the second column be "status"
##' @param lambda Penalty Parameters.
##' @param beta Initial values for beta. set all to 1 as default for the Lasso PSH model. 
##' 
##' @return A list of elements
##'   \tabular{ll}{
##'       \code{delindex}    \tab  A vector of index in x that to be discarded \cr
##'       \code{numdel} \tab  The number of discarded features. \cr
##'   }
##' @author Hong Wang
##' @examples
##' set.seed(123)
##' require(survival)
##' @export
SFEcmprsk <- function(x,y, lambda, beta=dim(x)[2]) {
  result=SFEc(x,y,lambda,beta)
  res <- list()
  res$delindex=result$delinx
  res$numdel=result$m
}