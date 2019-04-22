#' get_ifvals
#' 
#' get_ifvals is the main function for getting the influence function values
#' used to construct the estimator of either the lower or the upper bound curve.
#' 
#' @param eps vector of arbitrary length specifying the values for the 
#' proportion of confounding where the lower and upper bounds curves are 
#' evaluated.
#' @param upper boolean for whether the upper or the lower bound curve needs to
#' be computed.
#' @param nu nx1 vector of influence function values for the parameter 
#' E(E(Y|A = 1, X) - E(Y|A = 0, X)) 
#' @param tau nx1 vector of influence function values for the parameter
#' E(g_l(eta)) is upper=FALSE or E(g_u(eta)) if upper=TRUE. See manuscript.
#' @param ghatmat nxd matrix containing the values of g_l(eta)
#' if upper=FALSE or g_u(eta) if upper=TRUE, evaluated at each values of X. 
#' Each column of the matrix represents the vector of values for a specific 
#' value of delta, the maximum bias allowed for the confounded units. See 
#' manuscript. 
get_ifvals <- function(eps, upper, nu, tau, ghatmat) {
  
  n <- nrow(nu)
  neps <- ceiling(n*eps)
  ghatorder <- apply(ghatmat, 2, order)
  ndelta <- ncol(ghatmat)
  ghat <- sapply(1:ndelta, function(x) { ghatmat[ghatorder[, x], x] } )
  qhats <- matrix(NA, nrow=length(eps), ncol=ndelta)
  if(upper) { 
    qhats[which(n!=neps), ] <- ghat[n-neps, ]
    qhats[which(n==neps), ] <- ghat[1, ] - 0.01
    lambdahat <- do.call(cbind, lapply(neps, function(x) { c(rep(0, n-x), rep(1, x)) }))
    
  } else {
    qhats[which(neps!=0), ] <- ghat[neps, ]
    qhats[which(neps==0), ] <- ghat[1, ] - 0.01
    lambdahat <- do.call(cbind, lapply(neps, function(u) { c(rep(1, u), rep(0, n-u)) }))
  }
  nuorder <- unlist(lapply(1:ndelta, function(x) { nu[ghatorder[, x]] } ))
  if_gorder <- sapply(1:ndelta, function(x) { tau[ghatorder[, x], x] } )
  ifvals <- sweep(sweep(lambdahat, 1, if_gorder[, 1], "*"), 1, nuorder, "+")
  ifvals <- cbind(unlist(lapply(1:ndelta, function(x) { rep(x, n) } )), ifvals)
  colnames(ifvals) <- c("delta", eps)
  
  out <- list(ifvals=ifvals, lambda=lambdahat, quant=as.matrix(qhats))
  return(out)
}