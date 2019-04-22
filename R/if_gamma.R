#' if_gamma
#' 
#' if_gamma returns the influence function values for the parameter 
#' E(E(Y|A=a, X))) evaluated at each X.
#' @param y nx1 vector of outcomes in [0, 1]
#' @param a nx1 vector of treatments
#' @param aval scalar, if aval=a then it returns the if vals for E(E(Y|A=a, X)))
#' @param pia nx1 vector with estimates for P(A=aval | X) evaluated at each X. 
#' @param mua nx1 vector with estimates for E(Y | A=aval, X) evaluated at each 
#' X.
#' @export
if_gamma <- function(y, a, aval, pia, mua) {
  
  return(I(a==aval)*(y-mua)/pia + mua)
  
}