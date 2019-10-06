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
#' @param sl.lib vector specifying which libraries to use for the SL.
#' @param nsplits number of splits for the cross-fitting.
#' @export

do_crossfit <- function(y, a, x, outfam=gaussian(), treatfam=binomial(),
                        sl.lib=c("SL.earth","SL.gam","SL.glm", "SL.mean", 
                                 "SL.ranger","SL.glm.interaction"), nsplits=5) {
  n <- length(y)
  s <- sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
  mu0hat <- mu1hat <- pi0hat <- pi1hat <- rep(NA, n)
  
  for (vfold in 1:nsplits) {
    
    train <- s != vfold
    test <- s == vfold
    if (nsplits == 1) {
      train <- test
    }
    
    ytrain <- y[train]
    atrain <- a[train]
    xtrain <- x[train, , drop=FALSE]
    xtest <- x[test, , drop=FALSE]
    
    mu0hat[test] <- get_muahat(y=ytrain, a=atrain, x=xtrain, newx=xtest, aval=0,
                               sl.lib=sl.lib, family=outfam)
    mu1hat[test] <- get_muahat(y=ytrain, a=atrain, x=xtrain, newx=xtest, aval=1,
                               sl.lib=sl.lib, family=outfam)
    pi1hat[test] <- get_piahat(a=atrain, x=xtrain, newx=xtest, sl.lib=sl.lib,
                               family=treatfam)
    pi0hat[test] <- 1-pi1hat[test]
  } 
  out <- matrix(c(mu0hat, mu1hat, pi0hat, pi1hat), ncol=4, nrow=n, byrow=FALSE,
                dimnames=list(NULL, c("mu0", "mu1", "pi0", "pi1")))
  return(out)
}
