#' @title Multiplier bootstrap to construct uniform confidence bands around
#' the curves tracing the bounds on ATE.
#' 
#' @description \code{do_multboot} is the main function to implement 
#' the multiplier bootstrap to estimate uniformly valid bands around the bounds.
#' 
#' @param n sample size. 
#' @param ndelta number of delta parameters used. 
#' @param psihat Estimate of either lower or upper bound curve minus 
#' (epsilon x an estimate of the appropriate quantile). It is a 
#' (num eps points) x (ndelta) matrix.
#' @param sigmahat Estimate of the variance of psihat at each value of eps grid.
#'  It is a (num eps points) x (ndelta) matrix.
#' @param ifvals n x (num eps points) x ndelta array. It represents the 
#' "influence function" values associated with psihat.
#' @param alpha confidence level. Default is 0.05.
#' @param B number of rademacher rvs sampled. Default is 10000.
#' 
#' @return a ndelta-dimensional vector \code{calpha} equal to the z-score used 
#' to construct bands of the form psi(eps) \eqn{\pm} \code{calpha} * sigma(eps). 
#' 
#' @references Kennedy, E. H. (2019). Nonparametric causal effects based on 
#' incremental propensity score interventions. \emph{Journal of the American 
#' Statistical Association}, 114(526), 645-656.
#' @export

do_multboot <- function(n, ndelta, psihat, sigmahat, ifvals, alpha = 0.05, 
                        B = 10000) {

  ifvals2 <- sweep(sweep(ifvals, c(2, 3), psihat, "-"), c(2, 3), 
                   sqrt(sigmahat), "/")
  mult <- array(2 * rbinom(n * B * ndelta, 1, 0.5) - 1, dim = c(B, n, ndelta))
  
  maxvals <- sapply(1:ndelta, function(x) { 
    
    rad <- mult[, , x]
    ifs <- ifvals2[, , x]
    return(matrixStats::rowMaxs(rad %*% ifs))
    
    })

  calpha <- apply(maxvals / sqrt(n), 2, quantile, p = 1 - alpha)
  
  return(calpha)
  
}
