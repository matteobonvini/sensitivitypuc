#' do_crossfit
#' 
#' do_crossfit estimates the nuisance regression functions via cross-fitting.
#' This allows the user to avoid imposing empirical process conditions on these
#' functions, while still attaining fast rates of convergence.
#' 
#' @param y nx1 outcome vector in [0, 1]
#' @param a nx1 treatment received vector
#' @param x nxp matrix of covariates
#' @param outfam family specifying the error distribution for outcome 
#' regression, currently gaussian() or binomial(), should not specify link.
#' @param treatfam family specifying the error distribution for treatment 
#' regression, currently gaussian() or binomial(), should not specify link.
#' @param nsplits number of splits for the cross-fitting.
#' @param do_parallel boolean for whether parallel computing should be used
#' @param ncluster number of clusters used if parallel computing is used.
#' @param sl.lib vector specifying which libraries to use for the SL.
#' @export

do_crossfit <- function(y, a, x, ymin, ymax, outfam = gaussian(), 
                        treatfam = binomial(), nsplits = 5, do_parallel = FALSE,
                        ncluster = NULL,
                        sl.lib=c("SL.earth","SL.gam","SL.glm", "SL.mean", 
                                 "SL.ranger","SL.glm.interaction")) {
  n <- length(y)
  s <- sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
  
  train_idx <- lapply(1:nsplits, function(x) which(x != s))
  test_idx <- lapply(1:nsplits, function(x) which(x == s))
  
  if(do_parallel & nsplits > 1) {
    
    require(doRNG)
    
    cl <- makeCluster(ncluster)
    registerDoParallel(cl)
    
    out <- foreach(i = 1:nsplits, .combine = rbind, 
                   .packages='sensAteBounds') %dorng% {
                     
                     train <- train_idx[[i]]
                     test <- test_idx[[i]]
                     ytrain <- y[train]
                     atrain <- a[train]
                     xtrain <- x[train, , drop = FALSE]
                     xtest <- x[test, , drop = FALSE]
                     nt <- length(test)
                     
                     mu0hat <- get_muahat(y = ytrain, a = atrain, x = xtrain, newx = xtest, 
                                          aval = 0, sl.lib = sl.lib, family = outfam, 
                                          ymin = ymin, ymax = ymax)
                     
                     mu1hat <- get_muahat(y=ytrain, a=atrain, x=xtrain, newx=xtest, aval=1,
                                          sl.lib=sl.lib, family=outfam, ymin = ymin, 
                                          ymax = ymax)
                     
                     pi1hat <- get_piahat(a=atrain, x=xtrain, newx=xtest, sl.lib=sl.lib,
                                          family=treatfam)
                     
                     pi0hat <- 1 - pi1hat
                     
                     matrix(c(mu0hat, mu1hat, pi0hat, pi1hat, test), ncol=5, nrow = nt, byrow=FALSE,
                            dimnames=list(NULL, c("mu0", "mu1", "pi0", "pi1", "test_idx")))
                   }
    # 5th column contains the indices for the test set
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
      
      mu0hat[test] <- get_muahat(y = ytrain, a = atrain, x = xtrain, newx = xtest, 
                                 aval = 0, sl.lib = sl.lib, family = outfam, 
                                 ymin = ymin, ymax = ymax)
      mu1hat[test] <- get_muahat(y=ytrain, a=atrain, x=xtrain, newx=xtest, aval=1,
                                 sl.lib=sl.lib, family=outfam, ymin = ymin, 
                                 ymax = ymax)
      pi1hat[test] <- get_piahat(a=atrain, x=xtrain, newx=xtest, sl.lib=sl.lib,
                                 family=treatfam)
      pi0hat[test] <- 1 - pi1hat[test]
      
    } 
    out <- matrix(c(mu0hat, mu1hat, pi0hat, pi1hat), ncol=4, nrow=n, byrow=FALSE,
                  dimnames=list(NULL, c("mu0", "mu1", "pi0", "pi1")))
  }
  
  return(out)
}
