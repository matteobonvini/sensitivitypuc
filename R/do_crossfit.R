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
#' @param show_progress boolean for whether progress bar should be shown.
#' Default is FALSE. Currently, only available if do_parallel is FALSE.
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
                                   "SL.ranger","SL.glm.interaction"),
                        show_progress = FALSE) {
  
  if(do_parallel & show_progress) {
    stop("Currently, showing progress is only available if do_parallel = FALSE")
  }
  
  if(treatfam$family != "binomial") {
    stop("Currently only family = binomial() is supported for treatment model.")
  }
  
  n <- length(y)
  s <- sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
  
  train_idx <- lapply(1:nsplits, function(x) which(x != s))
  test_idx <- lapply(1:nsplits, function(x) which(x == s))
  
  tmp <- function(i, y, a, xmat, train_idx, test_idx, outfam, treatfam, ymin,
                  ymax, sl.lib) {
    ## Temporary function to avoid repeating code, i indicates the split and
    ## identifies the train set indices and the test set indices. 
    train <- train_idx[[i]]
    test <- test_idx[[i]]
    ytrain <- y[train]
    atrain <- a[train]
    xtrain <- xmat[train, , drop = FALSE]
    xtest <- xmat[test, , drop = FALSE]
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
  
  if(do_parallel) {
    ## Use dorng so that it easy to set the seed and reproduce the results.
    ## Issue is that I don't know how to easily show progress bar with dorng.
    cl <- makeCluster(ncluster)
    registerDoParallel(cl)
    
    out <- foreach(i = 1:nsplits, .combine = rbind, 
                   .packages='sensitivitypuc') %dorng% {
                     tmp(i, y = y, xmat = x, 
                         a = a, train_idx = train_idx, 
                         test_idx = test_idx, 
                         outfam = outfam, 
                         treatfam = treatfam, ymin = ymin, 
                         ymax = ymax, sl.lib = sl.lib)
                   }
    stopCluster(cl)
  } else {
    if(show_progress) {
      res <- pbapply::pblapply(1:nsplits, tmp, y = y, xmat = x, 
                               a = a, train_idx = train_idx, 
                               test_idx = test_idx, 
                               outfam = outfam, 
                               treatfam = treatfam, ymin = ymin, 
                               ymax = ymax, sl.lib = sl.lib)
      
    } else {
      res <- lapply(1:nsplits, tmp, y = y, xmat = x, 
                    a = a, train_idx = train_idx, 
                    test_idx = test_idx, 
                    outfam = outfam, treatfam = treatfam, 
                    ymin = ymin, ymax = ymax, sl.lib = sl.lib)
    }
    
    out <- do.call("rbind", res)
  }
  # "test_idx" contains indices for the test set, okay remove it after sorting
  # Need to keep track of test indices and sort so that order of observations
  # remains intact. 
  out <- out[order(out[, 5]), which(colnames(out) != "test_idx")]
  return(out)
} 
