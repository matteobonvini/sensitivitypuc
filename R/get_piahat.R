#' get_piahat
#' 
#' get_piahat estimates the regression function P(A=1|X=x) and evaluate it
#' at newx. 
#'
#' @param a nx1 vector of treatments
#' @param x nxp matrix of covariates
#' @param newx qxp matrix of values that the regression estimates are evaluated at
#' @param method a string specifying the method to use for estimation. Currently
#' supported is "logistic" and "ranger" (for random forest with ranger package)
#' @param trunc_tol scalar specifying the amount that the predicted values need 
#' to be truncated by. pihat is in [trunc_tol, 1-trunc_tol]. Default is 0.05.
#' @export

get_piahat <- function(a, x, newx, method="logistic", trunc_tol=0.05) {
  dat <- cbind(data.frame(a=as.factor(a)), x)
  if(method=="logistic"){
    fit <- glm(a ~ . , data=dat, family=binomial)
    fitvals <- predict(fit, newdata=as.data.frame(newx), type = "response")
  }
  if(method=="ranger"){
    fit <- ranger::ranger(as.formula(a ~ .), data=dat, probability=TRUE)
    fitvals <- predict(fit, data=newx)$predictions[, "1"]
  }
  fitvals <- truncate_prob(fitvals, tol=trunc_tol)
  return(fitvals)
}