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
#' @export
get_ifvals <- function(eps, upper, nu, tau, ghatmat) {
  
  # sample size
  n <- nrow(nu)
  neps <- ceiling(n*eps)
  
  # order the value of g(eta) increasingly so that it's easier to get term like
  # E(g I( g <= q_tau)) for q_tau being tau-quantile of g(eta).
  ghatorder <- apply(ghatmat, 2, order)
  ndelta <- ncol(ghatmat)
  ghat <- sapply(1:ndelta, function(x) { ghatmat[ghatorder[, x], x] } )
  .get_q <- Vectorize(function(delta) {
    if(upper){
      quant_eps <- 1-eps
    } else{
      quant_eps <- eps
    }
    out <- quantile(ghatmat[, which(colnames(ghatmat) == delta)], p = quant_eps)
    return(out)
  })
  qhats <- .get_q(colnames(ghatmat))
  
  if(upper) { 
    lambdahat <- do.call(cbind, lapply(neps, function(x) { c(rep(0, n-x), rep(1, x)) }))
    
  } else {
    lambdahat <- do.call(cbind, lapply(neps, function(x) { c(rep(1, x), rep(0, n-x)) }))
  }
  nuorder <- sapply(1:ndelta, function(x) { nu[ghatorder[, x]] } )
  if_gorder <- sapply(1:ndelta, function(x) { tau[ghatorder[, x], x] } )
  ifvals <- sapply(1:ndelta, function(x) { 
    sweep(lambdahat, 1, if_gorder[, x], "*") + nuorder[, x] })
  colnames(ifvals) <- colnames(ghatmat)
  ifvals <- cbind(unlist(lapply(1:length(eps), function(x) { rep(eps[x], n) } )), 
                  ifvals)
  colnames(ifvals) <- c("eps", colnames(ghatmat))
  
  out <- list(ifvals=ifvals, lambda=lambdahat, quant=as.matrix(qhats))
  return(out)
}
