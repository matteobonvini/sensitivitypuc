#' get_eps_zero
#' 
#' get_eps_zero computes an estimate of esp_zero and an associated 
#' confidence interval based on Theorem 4 in manuscript. 
#' 
#' @param eps vector of arbitrary length specifying the values for the 
#' proportion of confounding where the lower and upper bounds curves are 
#' evaluated.
#' @param lb vector of length(eps) containing the values for the lower bound
#' @param ub vector of length(eps) containing the values for the upper bound
#' @param ql vector of length(eps) containing the eps quantiles of g_l(eta)
#' @param qu vector of length(eps) containing the eps quantiles of g_u(eta)
#' @param ifvals_lb matrix nxlength(eps) containg the if values for lb:
#' lb - ql * lambdal, where lambdal is the indicator for g_l(eta) <= ql
#' @param ifvals_ub matrix nxlength(eps) containg the if values for ub:
#' ub - qu * lambdau, where lambdau is the indicator for g_u(eta) > qu
#' @param alpha scalar specifying the confidence level (default is alpha=0.05)
#' @export


get_eps_zero <- function(n, eps, lb, ub, ql, qu, ifvals_lb, ifvals_ub, 
                         delta, alpha=0.05){
  
  idx <- which.min(abs(lb*ub))
  est <- eps[idx]
  lb_zero <- lb[idx]
  ub_zero <- ub[idx]
  ql_zero <- ql[idx]
  qu_zero <- qu[idx]
  phi_zero <- ub_zero*ifvals_lb[which(ifvals_lb$eps==est), "ifvals"] +
              lb_zero*ifvals_ub[which(ifvals_ub$eps==est), "ifvals"]
  se_eps_zero <- sqrt((ub_zero*ql_zero + lb_zero*qu_zero)^(-2)*var(phi_zero)/n)
  ci <- get_ci(est, se_eps_zero, qnorm(1-alpha/2))
  out <- data.frame(delta = delta)
  out$est <- est
  out$ci_lo <- ci[, 1]
  out$ci_hi <- ci[, 2]
  return(out)
}
