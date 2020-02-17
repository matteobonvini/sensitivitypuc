#' @title Estimation of the minimum proportion of confounded units such that
#' the bounds on the average treatment effect contain zero (eps0) 
#' 
#' @description \code{get_eps_zero} computes an estimate of eps0 and Wald-type
#' confidence interval. 
#' 
#' @param n sample size.
#' @param eps vector of arbitrary length specifying the values for the 
#' proportion of confounding where the lower and upper bounds curves are 
#' evaluated.
#' @param pt_est point-estimate under no unmeasured confounding. Generally,
#' \code{mean(nuhat)}, where \code{nuhat} is a vector of influence fn values
#' for the parameter E{E(Y|A = 1, X) - E(Y|A = 0, X)}. 
#' @param lb length(eps) x (num of delta values) matrix containing the values 
#' for the lower bounds.
#' @param ub length(eps) x (num of delta values) matrix containing the values 
#' for the upper bounds.
#' @param ql a list of size \code{nsplits}, where element j is a
#' (num eps) x (num delta) matrix containing estimates of eps-quantile of 
#' g(etab) for lower bound computed using all obs except those in fold j.
#' @param qu a list of size \code{nsplits}, where element j is a
#' (num eps) x (num delta) matrix containing estimates of eps-quantile of 
#' g(etab) for upper bound computed using all obs except those in fold j.
#' @param ifvals_lb a list of size \code{nsplits}, where element j
#' is n x length(eps) x length(delta) array containing
#' values for \code{ifvals_lb} - \code{lambda_lb} * \code{q_lb} computed using
#' regression functions estimated from all obs except those in fold j and 
#' evaluated at obs in fold j.
#' @param ifvals_ub a list of size \code{nsplits}, where element j
#' is n x length(eps) x length(delta) array containing
#' values for \code{ifvals_ub} - \code{lambda_ub} * \code{q_ub} computed using
#' regression functions estimated from all obs except those in fold j and 
#' evaluated at obs in fold j.
#' @param delta vector of delta values specifying the values of delta used
#' to bound maximal confounding among S = 0 units. Length(\code{delta}) should 
#' match the 2nd dimension of \code{lb}, \code{ub}, \code{ql}, \code{qu} and the
#' 3rd dimension of \code{ifvals_lb} and \code{ifvals_ub}. 
#' @param alpha scalar specifying the confidence level. Default is 0.05.
#' 
#' @return A length(delta)x5 \code{data.frame} with values of delta, estimate
#' of eps0, \code{max}(0, ci_lo), \code{min}(1, ci_hi), variance of estimate of 
#' eps0. 
#' 
#' @seealso \code{\link{if_gamma}}, \code{\link{if_tau}},  
#' \code{\link{get_ifvals}}.
#' @export


get_eps_zero <- function(n, eps, pt_est, lb, ub, ql, qu, ifvals_lb, 
                         ifvals_ub, delta = 1, alpha = 0.05) {
  
  ndelta <- length(delta)
  
  find_zero <- function(x) { which.min(abs(x)) }
  
  is_ub <- I(pt_est < 0)
  if(is_ub){
    idx <- apply(ub, 2, find_zero)
    quant <- qu
    ifvals <- ifvals_ub
  } else {
    idx <- apply(lb, 2, find_zero)
    quant <- ql
    ifvals <- ifvals_lb
  }
  
  est <- eps[idx]
  
  ifvals_delta <- function(dat1, dat2) {
    rr <- dim(dat1)
    vapply(1:ndelta, function(x) { dat1[, idx[x], x] / dat2[idx[x], x] },
           FUN.VALUE = array(0, dim = rr[1]))
  }

  ifvals <- mapply(ifvals_delta, ifvals, quant, SIMPLIFY = FALSE)
  
  var_eps0 <- apply(abind::abind(ifvals, along = 1), 2, var) / n
  se_eps_zero <- sqrt(var_eps0)
  
  ci_lo <- est - qnorm(1 - alpha / 2) * se_eps_zero
  ci_hi <- est + qnorm(1 - alpha / 2) * se_eps_zero
  ci_lo[ci_lo < 0] <- 0
  ci_hi[ci_hi > 1] <- 1
  
  out <- data.frame(delta = delta, est = est, ci_lo = ci_lo, ci_hi = ci_hi,
                    var_eps0 = var_eps0)
  rownames(out) <- NULL
  
  return(out)
}
