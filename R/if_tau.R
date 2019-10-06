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
#' @export
if_tau <- Vectorize(function(y, a, pi0, pi1, mu0, mu1, delta=1, upper=FALSE) {
  
  ymax <- max(y)
  ymin <- min(y)
  
  if(!upper) {
    e1 <- I(delta <= (ymax - mu0)) # deltal = pmin(delta, ymax - mu0hat)
    e2 <- I(-delta >= (ymin - mu1)) # deltau = pmax(-delta, ymin - mu1hat)
    
    # Recall that g = pi0g*deltau - pi1g*deltal. Write uncentered IFs:
    
    # IF for E(pi(x)*(ymax - mu0(x)))
    if1 <- a*ymax - (pi1/pi0*(1-a)*(y-mu0) + a*mu0) 
    # IF for E((1-pi(x))*(ymin - mu1(x)))
    if2 <- (1-a)*ymin - (pi0/pi1*a*(y-mu1) + (1-a)*mu1)
    # IF for E(pi(x)*delta)
    if3 <- a*delta
    # IF for E((1-pi(x))*(-delta))
    if4 <- -(1-a)*delta
  } else {
    e1 <- I(-delta >= (ymin - mu0)) # deltal = pmax(-delta, ymin - mu0hat)
    e2 <- I(delta <= (ymax - mu1)) # deltau = pmin(delta, ymax - mu1hat)
    # Recall that g = pi0g*deltau - pi1g*deltal
    # write uncentered IFs
    # IF for E(pi(x)*(ymin - mu0(x)))
    if1 <- a*ymin - (pi1/pi0*(1-a)*(y-mu0) + a*mu0) 
    # IF for E((1-pi(x))*(ymax - mu1(x)))
    if2 <- (1-a)*ymax - (pi0/pi1*a*(y-mu1) + (1-a)*mu1)
    # IF for E(pi(x)*(-delta))
    if3 <- -a*delta
    # IF for E((1-pi(x))*delta)
    if4 <- (1-a)*delta
  }
  
  term1 <- I(e1 & e2) * (if4 - if3)
  term2 <- I(!e1 & e2) * (if4 - if1)
  term3 <- I(e1 & !e2) * (if2 - if3)
  term4 <- I(!e1 & !e2) * (if2 - if1)
  out <- term1 + term2 + term3 + term4
  return(as.matrix(out))
}, vectorize.args = "delta")
