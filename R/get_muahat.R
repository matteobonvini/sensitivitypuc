#' get_muahat
#' 
#' get_muahat estimates the regression function E(Y|A=aval,X=x) and evaluate it
#' at newx. 
#' 
#' @param y nx1 vector of outcomes in [0, 1]
#' @param a nx1 vector of treatments
#' @param x nxp matrix of covariates
#' @param newx qxp matrix of values that the regression estimates are evaluated at
#' @param aval scalar specifying which regression to estimate E(Y|A=aval,X=x)
#' @param method a string specifying the method to use for estimation. Currently
#' supported is "logistic" and "ranger" (for random forest with ranger package)
#' @param trunc_tol scalar specifying the amount that the predicted values need 
#' to be truncated by. Yhat is in [trunc_tol, 1-trunc_tol]. Default is 0.
#' @export

get_muahat <- function(y, a, x, newx, aval, method="logistic", trunc_tol=0) {
  
  dat <- cbind(data.frame(y=y[a==aval]), x[a==aval, , drop=FALSE])
  if(method=="logistic") {
    fit <- glm(as.factor(y) ~ . , data=dat, family=binomial)
    fitvals <- predict(fit, newdata=as.data.frame(newx), type = "response")
  }
  if(method=="ranger") {
    yvals <- unique(y)
    if(length(yvals)==2 & (all(yvals %in% c(0, 1)))) {
      dat$y <- as.factor(dat$y)
      fit <- ranger::ranger(as.formula(y ~ .), data=dat, probability=TRUE)
      fitvals <- predict(fit, data=as.data.frame(newx))$predictions[, "1"]
    } else {
      fit <- ranger::ranger(as.formula(y ~ .), data=dat)
      fitvals <- predict(fit, data=as.data.frame(newx))$predictions
    }
  }
  fitvals <- truncate_prob(fitvals, tol=trunc_tol)
  return(fitvals)
}