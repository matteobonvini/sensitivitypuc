#' get_g
#' 
#' get_g is the function that computes g_l(upper=FALSE) or g_u (upper=TRUE) 
#' given estimates of the nuisance regression functions evaluated at each value 
#' of X and a value for max_bias (default=1) specifyin the max_bias allowed for 
#' the confounded units. 
#' 
#' @param pi0 nx1 vector with estimates for P(A=0 | X=x) evaluated at each x
#' @param pi1 nx1 vector with estimates for P(A=1 | X=x) evaluated at each x
#' @param mu0 nx1 vector with estimates for E(Y |A=0, X=x) evaluated at each x
#' @param mu1 nx1 vector with estimates for E(Y |A=1, X=x) evaluated at each x
#' @param upper boolean (default is upper=FALSE) if TRUE g_u(eta) is returned
#' @param max_bias scalar specifying maximum bias for confounded units 
#' (default is max_bias=1)
#' @export
get_g <- function(pi0, pi1, mu0, mu1, upper=FALSE, max_bias=1) {
  
  if(upper) {
    return(pi1*(mu0-pmax(mu0-max_bias, 0)) + pi0*(pmin(mu1+max_bias, 1)-mu1))
  }
  else {
    return(pi1*(mu0-pmin(mu0+max_bias, 1)) + pi0*(pmax(mu1-max_bias, 0)-mu1))
  }
}