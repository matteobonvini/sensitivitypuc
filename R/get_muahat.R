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
#' @param ymin infimum of the support of y (predictions will be >= ymin)
#' @param ymax supremum of the support of y (predictions will be <= ymax)
#' @param family Currently allows gaussian or binomial to describe the error 
#' distribution. It is passed as argument to SL fun, should not contain link.
#' @export

get_muahat <- function(y, a, x, newx, aval, ymin, ymax, family = gaussian(), 
                       sl.lib = c("SL.earth","SL.gam","SL.glm", "SL.ranger",
                                  "SL.glm.interaction","SL.mean")) {
  y_aval <- y[a == aval]
  x_aval <- x[a == aval, , drop = FALSE]
  
  tried <- try(SuperLearner::SuperLearner(Y = y_aval, X = x_aval, family = family, 
                                          newX = as.data.frame(newx), SL.library = sl.lib)$SL.predict,
               silent = TRUE)
  if(inherits(tried, "try-error")) {
    message("SuperLearner throws this error:")
    message(tried[1])
    message(paste("\n I will use a GLM instead! Consider fitting nuisance", 
                  "regression functions separately and then pass them as",
                  "arg to get_bound()"))
    dat <- cbind(data.frame(y = y_aval), x_aval)
    fitvals <- predict.glm(glm(y ~ ., data = dat, family = family), 
                           newdata = as.data.frame(newx))
  } else {
    fitvals <- tried
  }
  
  fitvals[which(fitvals < ymin)] <- ymin
  fitvals[which(fitvals > ymax)] <- ymax
  
  return(fitvals)
}
