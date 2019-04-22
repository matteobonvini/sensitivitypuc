#' if_tau
#' 
#' if_tau returns the influence function values for the parameter 
#' E(g_l(eta)) (upper=FALSE) or E(g_u(eta)) (upper=TRUE) evaluated at each X.
#' 
#' @param y nx1 vector of outcomes in [0, 1]
#' @param a nx1 vector of treatments
#' @param pi0 nx1 vector with estimates for P(A=0 | X=x) evaluated at each x
#' @param pi1 nx1 vector with estimates for P(A=1 | X=x) evaluated at each x
#' @param mu0 nx1 vector with estimates for E(Y |A=0, X=x) evaluated at each x
#' @param mu1 nx1 vector with estimates for E(Y |A=1, X=x) evaluated at each x
#' @param upper boolean (default is upper=FALSE) if TRUE if values for 
#' E(g_u(eta)) are returned, otherwise those for E(g_l(eta)) are returned

if_tau <- function(y, a, pi0, pi1, mu0, mu1, upper=FALSE) {
  
  term1 <- pi1/pi0*(1-a)*(y-mu0)+a*mu0
  term2 <- -pi0/pi1*a*(y-mu1)-(1-a)*mu1
  term3 <- -a
  
  out <- term1 + term2 + term3
  if(upper) {
    out <- out + 1
  }
  
  return(as.matrix(out))
  
}