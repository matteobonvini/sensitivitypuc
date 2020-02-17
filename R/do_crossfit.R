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
#' @param ymin scalar such that P(Y >= ymin) = 1.
#' @param ymax scalar such that P(Y <= ymax) = 1.
#' @param nsplits number of splits for the cross-fitting.
#' @param sl.lib character vector specifying which libraries to use for the SL.
#' @param outfam family specifying the error distribution for outcome 
#' regression, currently \code{gaussian()} or \code{binomial()} supported. 
#' Link should not be specified. Default is \code{gaussian()}. 
#' @param treatfam family specifying the error distribution for treatment 
#' regression, currently \code{binomial()} supported.
#' Link should not be specified.
#' @param show_progress boolean for whether progress bar should be shown.
#' Default is FALSE. Currently, only available if do_parallel is FALSE.
#' @param do_parallel boolean for whether parallel computing should be used.
#' Default is FALSE.
#' @param ncluster number of clusters used if parallel computing is used.
#' 
#' @section Details: 
#' If the SuperLearner returns an error, a GLM is fitted instead. In this
#' case, we suggest the user chooses some other method of estimation and then
#' pass the estimates as arguments to other functions. 
#' 
#' @return A list containing 
#' \item{\code{test}}{a nx4 matrix containing estimates of E(Y|A = 0, X), 
#' E(Y|A = 1, X), P(A = 0|X), and P(A = a|X) evaluated at the test points X. 
#' If the function is estimated using folds 1 and 2, the values return are the 
#' predictions corresponding to fold 3 (assuming nsplits = 3 in this case).}
#' \item{\code{train}}{a (n*(nsplits-1))x4 matrix containing estimates of 
#' E(Y|A = 0, X), E(Y|A = 1, X), P(A = 0|X), and P(A = a|X) evaluated at the 
#' train points X. If the function is estimated using folds 1 and 2, 
#' the values return are the predictions corresponding to folds 1 and 2.} 
#' \item{\code{order_obs}}{a n-dimensional vector specifying the order of the
#' observations after doing cross-fitting, where the order is given by fold num.
#' For instance, if unit 1 is in fold 3, unit 2 is in fold 1, unit 3 is in fold 
#' 1 and unit 4 is in fold 2, order_obs = c(2, 3, 4, 1).}
#' \item{\code{folds}}{a n-dimensional vectors specifying which fold each unit
#' falls into. For instance, if unit 1 is in fold 3, unit 2 is in fold 1, 
#' unit 3 is in fold 1 and unit 4 is in fold 2, folds = c(3, 1, 1, 2).}
#' 
#' @examples 
#' n <- 500
#' x <- data.frame(x1 = rnorm(n), x2 = runif(n))
#' a <- rbinom(n, 1, pnorm(x$x1))
#' y <- 2 + x$x1 - x$x2 + rnorm(n)
#' fits <- do_crossfit(y, a, x, min(y), max(y), outfam = gaussian(), 
#'                     treatfam = binomial(), nsplits = 5, 
#'                     sl.lib = c("SL.mean", "SL.glm", "SL.gam"))
#' head(fits$test)
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

do_crossfit <- function(y, a, x, ymin, ymax, nsplits, sl.lib, 
                        outfam = gaussian(), treatfam = binomial(), 
                        show_progress = FALSE, do_parallel = FALSE, 
                        ncluster = NULL){
  
  if(do_parallel & show_progress) {
    stop("Currently, showing progress is only available if do_parallel = FALSE")
  }
  
  if(treatfam$family != "binomial") {
    stop("Currently only family = binomial() is supported for treatment model.")
  }
  
  n <- length(y)
  s <- sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
  
  if(nsplits == 1) {
    train_idx <- test_idx <- list(1:n)
  } else {
    train_idx <- lapply(1:nsplits, function(x) which(x != s))
    test_idx <- lapply(1:nsplits, function(x) which(x == s))
  }
  
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
    nt <- length(test) + length(train)
    
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
    
    test_res <- cbind(mu0hat$testvals, mu1hat$testvals, 1-pi1hat$testvals, 
                      pi1hat$testvals, test, i, 1)
    train_res <- cbind(mu0hat$trainvals, mu1hat$trainvals, 1-pi1hat$trainvals, 
                       pi1hat$trainvals, train, i, 0)
    res <- rbind(test_res, train_res)
    out <- matrix(res, ncol = 7, nrow = nt, byrow = FALSE,
                  dimnames=list(NULL, c("mu0", "mu1", "pi0", "pi1", 
                                        "unit_num", "fold", "is_test")))
    return(out)
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

  test_out <- out[which(out[, "is_test"] == 1), 
                  which(colnames(out) != "is_test")]
  train_out <- out[which(out[, "is_test"] == 0), 
                   which(colnames(out) != "is_test")]
  out <- list(test = test_out, train = train_out, 
              order_obs = test_out[, "unit_num"],
              folds = s)
  return(out)
} 
