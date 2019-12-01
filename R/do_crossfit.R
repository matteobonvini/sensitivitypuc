#' @title Estimation of regression functions using cross-fitting
#' 
#' @description \code{do_crossfit} estimates the nuisance regression functions 
#' using the SuperLearner and via cross-fitting. Cross-fitting allows the user 
#' to avoid imposing empirical process conditions on these functions, 
#' while still attaining, when possible, fast rates of convergence.
#' @importFrom doRNG %dorng%
#' @param y nx1 outcome vector in [0, 1]
#' @param a nx1 treatment received vector
#' @param x nxp \code{data.frame} of covariates. Variable must be named.
#' @param outfam family specifying the error distribution for outcome 
#' regression, currently \code{gaussian()} or \code{binomial()} supported. 
#' Link should not be specified. Default is \code{gaussian()}. 
#' @param treatfam family specifying the error distribution for treatment 
#' regression, currently \code{binomial()} supported.
#' Link should not be specified.
#' @param nsplits number of splits for the cross-fitting.
#' @param do_parallel boolean for whether parallel computing should be used.
#' Default is FALSE.
#' @param ncluster number of clusters used if parallel computing is used.
#' @param sl.lib character vector specifying which libraries to use for the SL.
#' 
#' @section Details: 
#' If the SuperLearner returns an error, a GLM is fitted instead. In this
#' case the user is suggested to choose some other method of estimation and then
#' pass the estimates as arguments to other functions. 
#' 
#' @return A nx4 matrix containing estimates of E(Y|A = 0, X), E(Y|A = 1, X),
#' P(A = 0|X), and P(A = a|X) evaluated at the observed values of X. 
#' 
#' @examples 
#' n <- 500
#' x <- data.frame(x1 = rnorm(n), x2 = runif(n))
#' a <- rbinom(n, 1, pnorm(x$x1))
#' y <- 2 + x$x1 - x$x2 + rnorm(n)
#' fits <- do_crossfit(y, a, x, min(y), max(y), outfam = gaussian(), 
#'                     treatfam = binomial(), nsplits = 5, 
#'                     sl.lib = c("SL.mean", "SL.glm", "SL.gam"))
#' head(fits)
#' 
#' @seealso \code{\link{get_muahat}} and \code{\link{get_piahat}}.
#'
#' @references Van der Laan, M. J., Polley, E. C., & Hubbard, A. E. (2007). 
#' Super learner. 
#' \emph{Statistical applications in genetics and molecular biology}, 6(1).
#' @references Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, 
#' E., Hansen, C., & Newey, W. K. (2016). Double machine learning for treatment 
#' and causal parameters (No. CWP49/16). \emph{cemmap working paper}.
#' 
#' @export

do_crossfit <- function(y, a, x, ymin, ymax, outfam = gaussian(), 
                        treatfam = binomial(), nsplits = 5, do_parallel = FALSE,
                        ncluster = NULL,
                        sl.lib = c("SL.earth","SL.gam","SL.glm", "SL.mean", 
                                   "SL.ranger","SL.glm.interaction")) {
  
  if(treatfam$family != "binomial") {
    stop("Currently only family = binomial() is supported for treatment model.")
  }
  
  n <- length(y)
  s <- sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
  
  train_idx <- lapply(1:nsplits, function(x) which(x != s))
  test_idx <- lapply(1:nsplits, function(x) which(x == s))
  
  if(do_parallel & nsplits > 1) {

    cl <- makeCluster(ncluster)
    registerDoParallel(cl)
    
    out <- foreach(i = 1:nsplits, .combine = rbind, 
                   .packages='sensitivitypuc') %dorng% {
                     
                     train <- train_idx[[i]]
                     test <- test_idx[[i]]
                     ytrain <- y[train]
                     atrain <- a[train]
                     xtrain <- x[train, , drop = FALSE]
                     xtest <- x[test, , drop = FALSE]
                     nt <- length(test)
                     
                     mu0hat <- get_muahat(y = ytrain, a = atrain, x = xtrain, 
                                          newx = xtest, aval = 0, 
                                          sl.lib = sl.lib, family = outfam, 
                                          ymin = ymin, ymax = ymax)
                     
                     mu1hat <- get_muahat(y = ytrain, a = atrain, x = xtrain, 
                                          newx = xtest, aval = 1, 
                                          sl.lib = sl.lib, family = outfam, 
                                          ymin = ymin, ymax = ymax)
                     
                     pi1hat <- get_piahat(a = atrain, x = xtrain, newx = xtest, 
                                          sl.lib = sl.lib, family = treatfam)
                     
                     pi0hat <- 1 - pi1hat
                     
                     res <- c(mu0hat, mu1hat, pi0hat, pi1hat, test)
                     matrix(res, ncol = 5, nrow = nt, byrow = FALSE,
                            dimnames=list(NULL, c("mu0", "mu1", "pi0", "pi1", 
                                                  "test_idx")))
                   }
    # 5th column contains indices for the test set, okay remove it after sorting
    out <- out[order(out[, 5]), -5]
    stopCluster(cl)
    
  } else {
    
    mu0hat <- mu1hat <- pi0hat <- pi1hat <- rep(NA, n)
    
    for (i in 1:nsplits) {
      
      train <- train_idx[[i]]
      test <- test_idx[[i]]
      
      if (nsplits == 1) {
        train <- test
      }
      
      ytrain <- y[train]
      atrain <- a[train]
      xtrain <- x[train, , drop = FALSE]
      xtest <- x[test, , drop = FALSE]
      
      mu0hat[test] <- get_muahat(y = ytrain, a = atrain, x = xtrain, 
                                 newx = xtest, aval = 0, sl.lib = sl.lib, 
                                 family = outfam, ymin = ymin, ymax = ymax)
      mu1hat[test] <- get_muahat(y = ytrain, a = atrain, x = xtrain, 
                                 newx = xtest, aval = 1, sl.lib = sl.lib, 
                                 family = outfam, ymin = ymin, ymax = ymax)
      pi1hat[test] <- get_piahat(a = atrain, x = xtrain, newx = xtest, 
                                 sl.lib = sl.lib, family = treatfam)
      pi0hat[test] <- 1 - pi1hat[test]
      
    } 
    
    res <- c(mu0hat, mu1hat, pi0hat, pi1hat)
    out <- matrix(res, ncol = 4, nrow = n,  byrow = FALSE, 
                  dimnames = list(NULL, c("mu0", "mu1", "pi0", "pi1")))
  }
  
  return(out)
}
