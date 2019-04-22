#' get_multboot
#' 
#' get_multboot is the main function to implement the multiplier bootstrap
#' 
#' @param n the number of observations
#' @param psihat ADD
#' @param sigmahat ADD
#' @param ifvals ADD
#' @param alpha confidence level (default is 0.05)
#' @export

get_multboot <- function(n, psihat, sigmahat, ifvals, alpha=0.05) {
  
  psihat_mat <- matrix(rep(psihat, n), nrow=n, byrow=TRUE)
  sigmahat_mat <- matrix(rep(sigmahat, n), nrow=n ,byrow=TRUE)
  
  ifvals2 <- (ifvals-psihat_mat)/sqrt(sigmahat_mat)
  B <- 1000
  radem <- matrix(2*rbinom(n*B, 1, 0.5)-1, nrow=n, ncol=B)
  maxvals <- sapply(1:B, function(col){
    max(abs(apply(radem[,col]*ifvals2,2,sum)/sqrt(n))) } )
  calpha <- quantile(maxvals, 1-alpha)
  return(calpha)
  
}