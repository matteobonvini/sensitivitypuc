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
#' @return A list containing estimates of E(Y | A = aval, X):
#' \item{\code{testvals}}{nrow(\code{newx})-dimensional vector containing 
#' estimates of regression functions (computed using points \code{x}) evaluated
#' at test points \code{newx}.}
#' \item{\code{trainvals}}{nrow(\code{x})-dimensional vector containing 
#' estimates of regression functions (computed using points \code{x}) evaluated 
#' at train points \code{x}.}
#' 
#' @examples 
#' n <- 500
#' x <- data.frame(x1 = rnorm(n), x2 = runif(n))
#' newx <- data.frame(x1 = rnorm(20), x2 = runif(20))
#' a <- rbinom(n, 1, pnorm(x$x1))
#' y <- 2 + x$x1 - x$x2 + rnorm(n)
#' fits <- get_muahat(y, a, x, newx, 1, min(y), max(y), family = gaussian(), 
#'                    sl.lib = c("SL.mean", "SL.glm", "SL.gam"))
#' head(fits$testvals)
#' 
#' @seealso \code{\link{do_crossfit}} and \code{\link{get_piahat}}.
#' 
#' @references Van der Laan, M. J., Polley, E. C., & Hubbard, A. E. (2007). 
#' Super learner. 
#' \emph{Statistical applications in genetics and molecular biology}, 6(1).
#' 
#' @export

get_muahat <- function(y, a, x, newx, aval, ymin, ymax, family = gaussian(), 
                       sl.lib) {
  
  y_aval <- y[a == aval]
  x_aval <- x[a == aval, , drop = FALSE]
  
  tryfit <- try(SuperLearner::SuperLearner(Y = y_aval, X = x_aval, 
                                           family = family, SL.library = sl.lib), 
                silent = TRUE)
  
  if(inherits(tryfit, "try-error")) {
    
    message(paste0("SuperLearner throws this error for estimating mu", 
                   aval, ":"))
    message(tryfit[1])
    message(paste("I will use a GLM instead! Consider fitting nuisance", 
                  "regression functions separately and then pass them as",
                  "argument to, for instance, get_bound()"))
    
    dat <- cbind(data.frame(y = y_aval), x_aval)
    fit <- glm(y ~ ., data = dat, family = family)
    
    typepred <- ifelse(family$family == "binomial", "response", "link")
    
    testvals <- predict.glm(fit, newdata = as.data.frame(newx), type = typepred)
    trainvals <- predict.glm(fit, newdata = as.data.frame(x), type = typepred)
    
  } else {
    testvals <- predict(tryfit, as.data.frame(newx), onlySL = TRUE)$pred
    trainvals <- predict(tryfit, as.data.frame(x), onlySL = TRUE)$pred
  }
  
  testvals[which(testvals < ymin)] <- ymin
  testvals[which(testvals > ymax)] <- ymax
  trainvals[which(trainvals < ymin)] <- ymin
  trainvals[which(trainvals > ymax)] <- ymax
  
  out <- list(testvals = testvals, trainvals = trainvals)
  return(out)
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
#' regression, currently only \code{binomial()} supported. Default is 
#' \code{binomial()}. Link should not be specified.
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
#' @return A list containing estimates of E(A|X): 
#' \item{\code{testvals}}{nrow(\code{newx})-dimensional vector containing 
#' estimates of regression functions (computed using points \code{x}) evaluated
#' at test points \code{newx}.}
#' \item{\code{trainvals}}{nrow(\code{x})-dimensional vector containing 
#' estimates of regression functions (computed using points \code{x}) evaluated 
#' at train points \code{x}.}
#' 
#' @examples 
#' n <- 500
#' x <- data.frame(x1 = rnorm(n), x2 = runif(n))
#' newx <- data.frame(x1 = rnorm(20), x2 = runif(20))
#' a <- rbinom(n, 1, pnorm(x$x1))
#' fits <- get_piahat(a, x, newx, family = binomial(), 
#'                    sl.lib = c("SL.mean", "SL.glm", "SL.gam"))
#' head(fits$testvals)
#' 
#' @seealso \code{\link{do_crossfit}} and \code{\link{get_muahat}}.
#' 
#' @references Van der Laan, M. J., Polley, E. C., & Hubbard, A. E. (2007). 
#' Super learner. 
#' \emph{Statistical applications in genetics and molecular biology}, 6(1).
#' 
#' @export

get_piahat <- function(a, x, newx, family = binomial(), trunc_tol = 0.05, 
                       sl.lib) {
  
  tryfit <- try(SuperLearner::SuperLearner(Y = a, X = x, family = family,
                                           SL.library = sl.lib), 
                silent = TRUE)
  
  if(inherits(tryfit, "try-error")) {
    
    message("SuperLearner throws this error for estimating pi:")
    message(tryfit[1])
    message(paste("I will use a GLM instead! Consider fitting nuisance", 
                  "regression functions separately and then pass them as",
                  "arg to get_bound()"))
    
    dat <- cbind(data.frame(a = a), x)
    fit <- glm(a ~ ., data = dat, family = family)
    testvals <- predict.glm(fit, newdata = as.data.frame(newx), 
                            type = "response")
    trainvals <- predict.glm(fit, newdata = as.data.frame(x), type = "response")
    
  } else {
    testvals <- predict(tryfit, as.data.frame(newx), onlySL = TRUE)$pred
    trainvals <- predict(tryfit, as.data.frame(x), onlySL = TRUE)$pred
  }
  
  testvals <- truncate_prob(testvals, tol = trunc_tol)
  trainvals <- truncate_prob(trainvals, tol = trunc_tol)
  out <- list(testvals = testvals, trainvals = trainvals)
  
  return(out)
}
