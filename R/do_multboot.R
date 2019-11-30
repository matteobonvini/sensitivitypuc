#' @title Multiplier bootstrap to construct uniform confidence bands around
#' the curves tracing the bounds on ATE.
#' 
#' @description \code{do_multboot} is the main function to implement 
#' the multiplier bootstrap to estimate uniformly valid bands around the bounds.
#' 
#' @param n sample size. 
#' @param psihat Estimate of either lower or upper bound curve minus 
#' (epsilon x an estimate of the appropriate quantile). It is a vector 
#' of length equal to the number of epsilon points at which the lower and upper 
#' bounds curves are evaluated.
#' @param sigmahat Estimate of the variance of psihat at each value of eps grid.
#' It is a vector of length equal to the number of epsilon points at which the 
#' lower and upper bounds curves are evaluated.
#' @param ifvals n times neps matrix, where n is number of observations, and 
#' neps is the number of epsilon points at which the lower and upper bounds 
#' curves are evaluated. It represents the "influence function" values 
#' associated with psihat.
#' @param alpha confidence level. Default is 0.05.
#' @param B number of rademacher rvs sampled. Default is 10000.
#' 
#' @return a scalar calpha equal to the multiplier used to construct bands of 
#' the form psi(eps) \pm calpha * sigma(eps). 
#' 
#' @references Kennedy, E. H. (2019). Nonparametric causal effects based on 
#' incremental propensity score interventions. \emph{Journal of the American 
#' Statistical Association}, 114(526), 645-656.
#' @export

do_multboot <- function(n, psihat, sigmahat, ifvals, alpha = 0.05, B = 10000) {
  
  ifvals2 <- sweep(sweep(ifvals, 2, psihat, "-"), 2, sqrt(sigmahat), "/")
  mult <- matrix(2 * rbinom(n * B, 1, 0.5) - 1, nrow = B, ncol = n)

  maxvals <- matrixStats::rowMaxs(mult %*% ifvals2 / sqrt(n))
  calpha <- quantile(maxvals, 1 - alpha)
  
  return(calpha)
  
}
