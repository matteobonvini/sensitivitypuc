#' @title Estimation of regression function E(Y|A = a, X)
#' 
#' @description \code{get_muahat} estimates the regression function 
#' E(Y|A = a, X) using the SuperLearner. 
#' 
#' @param y nx1 vector of outcomes in [0, 1]
#' @param a nx1 vector of treatments
#' @param x nxp \code{data.frame} of covariates
#' @param newx qxp \code{data.frame} of values that the regression estimates are 
#' evaluated at
#' @param aval Scalar specifying which regression to estimate E(Y|A=aval,X=x)
#' @param ymin Infimum of the support of y, so that predictions will be forced 
#' to be >= ymin.
#' @param ymax Supremum of the support of y, so that predictions will be forced 
#' to be <= ymax.
#' @param family family specifying the error distribution for outcome 
#' regression, currently \code{gaussian()} or \code{binomial()} supported. 
#' Link should not be specified. Default is \code{gaussian()}. 
#' @param sl.lib character vector specifying which libraries to use for the SL.
#' 
#' @section Details: 
#' If the SuperLearner returns an error, a GLM is fitted instead. In this
#' case the user is suggested to choose some other method of estimation and then
#' pass the estimates as arguments to other functions. Sometimes SuperLearner
#' returns warnings messages.
#' 
#' @return A n-dimensional vectors containing estimates of E(Y|A = aval, X) 
#' evaluated at newx.
#' 
#' @examples 
#' n <- 500
#' x <- data.frame(x1 = rnorm(n), x2 = runif(n))
#' newx <- data.frame(x1 = rnorm(20), x2 = runif(20))
#' a <- rbinom(n, 1, pnorm(x$x1))
#' y <- 2 + x$x1 - x$x2 + rnorm(n)
#' fits <- get_muahat(y, a, x, newx, 1, min(y), max(y), family = gaussian(), 
#'                    sl.lib = c("SL.mean", "SL.glm", "SL.gam"))
#' head(fits)
#' 
#' @seealso \code{\link{do_crossfit}} and \code{\link{get_piahat}}.
#' 
#' @references Van der Laan, M. J., Polley, E. C., & Hubbard, A. E. (2007). 
#' Super learner. 
#' \emph{Statistical applications in genetics and molecular biology}, 6(1).
#' 
#' @export

get_muahat <- function(y, a, x, newx, aval, ymin, ymax, family = gaussian(), 
                       sl.lib = c("SL.earth", "SL.gam", "SL.glm", "SL.ranger",
                                  "SL.glm.interaction", "SL.mean")) {
  
  y_aval <- y[a == aval]
  x_aval <- x[a == aval, , drop = FALSE]
  
  tried <- try(SuperLearner::SuperLearner(Y = y_aval, X = x_aval, 
                                          family = family, 
                                          newX = as.data.frame(newx), 
                                          SL.library = sl.lib)$SL.predict,
               silent = TRUE)
  
  if(inherits(tried, "try-error")) {
    
    message("SuperLearner throws this error:")
    message(tried[1])
    message(paste("I will use a GLM instead! Consider fitting nuisance", 
                  "regression functions separately and then pass them as",
                  "argument to, for instance, get_bound()"))
    
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

#' @title Estimation of regression function E(A|X)
#' 
#' @description \code{get_piahat} estimates the regression function E(A|X)
#' using the SuperLearner. 
#' 
#' @param a nx1 vector of treatments
#' @param x nxp \code{data.frame} of covariates
#' @param newx qxp \code{data.frame} of values that the regression estimates are 
#' evaluated at
#' @param family family specifying the error distribution for treatment
#' regression, currently only \code{binomial()} supported. 
#' Link should not be specified.
#' @param trunc_tol amount of tolerance allowed for truncating the propensity 
#' score in [tol, 1 - tol]. Default is 0.05.
#' @param sl.lib character vector specifying which libraries to use for the SL.
#' 
#' @section Details: 
#' If the SuperLearner returns an error, a GLM is fitted instead. In this
#' case the user is suggested to choose some other method of estimation and then
#' pass the estimates as arguments to other functions. Sometimes SuperLearner
#' returns warnings messages.
#' 
#' @return A n-dimensional vectors containing estimates of E(Y|A = aval, X) 
#' evaluated at newx.
#' 
#' @examples 
#' n <- 500
#' x <- data.frame(x1 = rnorm(n), x2 = runif(n))
#' newx <- data.frame(x1 = rnorm(20), x2 = runif(20))
#' a <- rbinom(n, 1, pnorm(x$x1))
#' fits <- get_piahat(a, x, newx, family = binomial(), 
#'                    sl.lib = c("SL.mean", "SL.glm", "SL.gam"))
#' head(fits)
#' 
#' @seealso \code{\link{do_crossfit}} and \code{\link{get_muahat}}.
#' 
#' @references Van der Laan, M. J., Polley, E. C., & Hubbard, A. E. (2007). 
#' Super learner. 
#' \emph{Statistical applications in genetics and molecular biology}, 6(1).
#' 
#' @export

get_piahat <- function(a, x, newx, family = binomial(), trunc_tol = 0.05,
                       sl.lib = c("SL.earth", "SL.gam", "SL.glm", "SL.mean", 
                                  "SL.glm.interaction", "SL.ranger")) {
  
  tried <- try(SuperLearner::SuperLearner(Y = a, X = x, family = family,
                                          newX = as.data.frame(newx),
                                          SL.library = sl.lib)$SL.predict,
               silent = TRUE)
  
  if(inherits(tried, "try-error")) {
    
    message("SuperLearner throws this error:")
    message(tried[1])
    message(paste("I will use a GLM instead! Consider fitting nuisance", 
                  "regression functions separately and then pass them as",
                  "arg to get_bound()"))
    
    dat <- cbind(data.frame(a = a), x)
    fitvals <- predict.glm(glm(a ~ ., data = dat, family = family), 
                           newdata = as.data.frame(newx))
    
  } else {
    fitvals <- tried
  }
  
  fitvals <- truncate_prob(fitvals, tol = trunc_tol)
  
  return(fitvals)
}
