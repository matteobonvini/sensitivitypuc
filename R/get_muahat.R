#' get_muahat
#' 
#' get_muahat estimates the regression function E(Y|A=aval,X=x) and evaluate it
#' at newx. 
#' 
#' @param y nx1 vector of outcomes in [0, 1]
#' @param a nx1 vector of treatments
#' @param x nxp matrix of covariates
#' @param newx qxp matrix of values that the regression estimates are evaluated at
#' @param aval Scalar specifying which regression to estimate E(Y|A=aval,X=x)
#' @param family Currently allows gaussian or binomial to describe the error 
#' distribution. It is passed as argument to SL fun, should not contain link.
#' @param trunc_tol Scalar specifying the amount that the predicted values need 
#' to be truncated by. Yhat is in [trunc_tol, 1-trunc_tol]. Default is 0.
#' @export

get_muahat <- function(y, a, x, newx, aval, trunc_tol=0, family=gaussian(),
                       sl.lib=c("SL.earth","SL.gam","SL.glm",
                                "SL.glm.interaction","SL.mean", "SL.ranger")) {
  fit <- SuperLearner::SuperLearner(Y=y[a==aval], X=x[a==aval, ,drop=FALSE], 
                                    newX = as.data.frame(newx), 
                                    family=family, SL.library=sl.lib)
  fitvals <- truncate_prob(fit$SL.predict, tol=trunc_tol)
  return(fitvals)
}