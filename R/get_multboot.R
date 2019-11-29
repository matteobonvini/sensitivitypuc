#' get_multboot
#' 
#' get_multboot is the main function to implement the multiplier bootstrap
#' 
#' @param n the number of observations
#' @param psihat Estimate of either lower or upper bound curve minus 
#' (epsilon times an estimate of the appropriate quantile). It is a vector 
#' of length equal to the number of epsilon points at which the lower and upper 
#' bounds curves are evaluated.
#' @param sigmahat Estimate of the variance of psihat at each value of eps grid.
#' It is a vector of length equal to the number of epsilon points at which the 
#' lower and upper bounds curves are evaluated.
#' @param ifvals n times neps matrix, where n is number of observations, and 
#' neps is the number of epsilon points at which the lower and upper bounds 
#' curves are evaluated. It represents the "influence function" values 
#' associated with psihat. See manuscript (Lemma 2)
#' @param alpha confidence level (default is 0.05)
#' @param B number of rademacher rvs sampled
#' @export

get_multboot <- function(n, psihat, sigmahat, ifvals, alpha = 0.05, B = 10000) {
  
  ifvals2 <- sweep(sweep(ifvals, 2, psihat, "-"), 2, sqrt(sigmahat), "/")
  mult <- matrix(2 * rbinom(n * B, 1, 0.5) - 1, nrow = B, ncol = n)

  maxvals <- matrixStats::rowMaxs(mult %*% ifvals2 / sqrt(n))
  calpha <- quantile(maxvals, 1 - alpha)
  
  return(calpha)
  
}
