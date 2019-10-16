#' get_ifvals
#' 
#' get_ifvals is the main function for getting the influence function values
#' used to construct the estimator of either the lower or the upper bound curve.
#' 
#' @param n number of observations in the dataset (n = length(nu) = nrow(ghatmat))
#' @param eps vector of arbitrary length specifying the values for the 
#' proportion of confounding where the lower and upper bounds curves are 
#' evaluated.
#' @param delta vector of delta values specifying the extent of unmeasured
#' confounding among S = 0 units, see manuscript. 
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

get_ifvals <- function(n, eps, delta, upper, nu, tau, ghatmat) {
  
  ndelta <- length(delta)
  neps <- length(eps)
  
  if(upper) {
    seq_eps <- 1-eps
  } else {
    seq_eps <- eps
  }
  
  qhats <- apply(ghatmat, 2, quantile, p = seq_eps, names = FALSE)
  
  # The following is because R selects the min as 0-quantile, matters only in
  # finite samples
  qhats[which(seq_eps == 0), ] <- apply(ghatmat, 2, min) - 1e-10
  
  ineq_sign <- ifelse(upper, ">", "<=")
  .get_indicator <- function(x) {
    out <- 1 * sweep(ghatmat, 2, qhats[x, ], ineq_sign)
    return(out)
  }
  lambda <- aperm(sapply(1:length(eps), .get_indicator, simplify = "array"),
                  c(1, 3, 2))
  
  # rqhats <- aperm(array(qhats, dim = c(neps, ndelta, 1)), c(3, 1, 2))
  # dimnames(rqhats) <- list(1, eps, delta)
  # .get_indicator <- function(x) {
  #   out <- 1 * sweep(ghatmat, 2, x, ineq_sign)
  #   return(out)
  # }
  
  # lambda <- apply(qhats, 1, .get_indicator)
  # lambda <- aperm(array(lambda, dim = c(ndelta, neps, n)), c(3, 2, 1))
  
  .get_lambdaq <- function(x) {
    out <- sweep(lambda[, , x], 2, qhats[, x], "*")
    return(out)
  }
  lambdaq <- sapply(1:ndelta, .get_lambdaq, simplify = "array")
  
  .get_ifs <- function(x) {
    out <- sweep(sweep(lambda[, , x], 1, tau[, x], "*"), 1, nu, "+")
    return(out)
  }
  ifs <- sapply(1:ndelta, .get_ifs, simplify = "array")
  
  out <- list(ifvals = ifs, lambda = lambda, quant = qhats, lambdaq = lambdaq)
  
  return(out)
}
