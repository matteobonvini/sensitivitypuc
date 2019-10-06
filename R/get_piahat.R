#' get_piahat
#' 
#' get_piahat estimates the regression function P(A=1|X=x) and evaluate it
#' at newx. 
#'
#' @param a nx1 vector of treatments
#' @param x nxp matrix of covariates
#' @param newx qxp matrix of values that the regression estimates are evaluated at
#' @param family Currently allows gaussian or binomial to describe the error 
#' distribution. It is passed as argument to SL fun, should not contain link.
#' @param trunc_tol scalar specifying the amount that the predicted values need 
#' to be truncated by. pihat is in [trunc_tol, 1-trunc_tol]. Default is 0.05.
#' @export

get_piahat <- function(a, x, newx, family=binomial(), trunc_tol=0.05,
                       sl.lib=c("SL.earth","SL.gam","SL.glm",
                                "SL.glm.interaction","SL.mean", "SL.ranger")) {
  fit <- SuperLearner::SuperLearner(Y=a, X = x, newX=as.data.frame(newx), 
                                    SL.library=sl.lib, family=family)
  fitvals <- truncate_prob(fit$SL.predict, tol=trunc_tol)
  return(fitvals)
}