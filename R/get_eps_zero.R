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
#' @param lb vector of length(eps) containing the values for the lower bound.
#' @param ub vector of length(eps) containing the values for the upper bound.
#' @param ql vector of length(eps) containing the eps quantiles of g_l(eta).
#' @param qu vector of length(eps) containing the eps quantiles of g_u(eta).
#' @param ifvals_lb matrix nxlength(eps) containg the if values for lb.
#' @param ifvals_ub matrix nxlength(eps) containg the if values for ub.
#' @param delta vector of arbitrary length specifyin the values of delta used
#' to bound maximal confounding among S = 0 units. Default is delta = 1, which 
#' imposes no assumption if the outcome Y is bounded. 
#' @param alpha scalar specifying the confidence level. Default is 0.05.
#' 
#' @return A length(delta)x5 \code{data.frame} with values of delta, estimate
#' of eps0, max(0, ci_lo), min(1, ci_hi), variance of estimate of eps0. 
#' 
#' @examples 
#' eps <- seq(0, 0.1, 0.001)
#' delta <- c(0.5, 1)
#' n <- 500
#' x <- data.frame(x1 = rnorm(n), x2 = runif(n))
#' a <- rbinom(n, 1, pnorm(x$x1))
#' y <- 2 + x$x1 - x$x2 + rnorm(n)
#' ymin <- min(y)
#' ymax <- max(y)
#' nuis_fns <- do_crossfit(y, a, x, min(y), max(y), outfam = gaussian(), 
#'                         treatfam = binomial(), nsplits = 5, 
#'                         sl.lib = c("SL.mean", "SL.glm", "SL.gam"))
#' pi0hat <- nuis_fns[, "pi0"]
#' pi1hat <- nuis_fns[, "pi1"]
#' mu0hat <- nuis_fns[, "mu0"]
#' mu1hat <- nuis_fns[, "mu1"]
#' 
#' psi0hat <- if_gamma(y = y, a = a, aval = 0, pia = pi0hat, mua = mu0hat)
#' psi1hat <- if_gamma(y = y, a = a, aval = 1, pia = pi1hat, mua = mu1hat)
#' nuhat <- as.matrix(psi1hat - psi0hat)
#' guhat <- pi0hat * (ymax - mu1hat) - pi1hat * (ymin - mu0hat)
#' glhat <- pi0hat * (ymin - mu1hat) - pi1hat * (ymax - mu0hat)
#' glhat <- glhat %*% t(delta)
#' guhat <- guhat %*% t(delta)
#' 
#' tauhat_lb <- if_tau(y = y, a = a, ymin = ymin, ymax = ymax, pi0 = pi0hat, 
#'                     pi1 = pi1hat, mu0 = mu0hat, mu1 = mu1hat, upper = FALSE)
#' tauhat_ub <- if_tau(y = y, a = a, ymin = ymin, ymax = ymax, pi0 = pi0hat, 
#'                     pi1 = pi1hat, mu0 = mu0hat, mu1 = mu1hat, upper = TRUE)
#' tauhat_lb <- tauhat_lb %*% t(delta)
#' tauhat_ub <- tauhat_ub %*% t(delta)
#'                     
#' list_lb <- get_ifvals(n = n, eps = eps, delta = delta, upper = FALSE, 
#'                       nu = nuhat, tau = tauhat_lb, ghatmat = glhat)
#' list_ub <- get_ifvals(n = n, eps = eps, delta = delta, upper = TRUE, 
#'                       nu = nuhat, tau = tauhat_ub, ghatmat = guhat)
#'                       
#' ql <- list_lb$quant
#' qu <- list_ub$quant
#' 
#' phibar_lb <- list_lb$ifvals - list_lb$lambdaq
#' phibar_ub <- list_ub$ifvals - list_ub$lambdaq
#' 
#' lb <- apply(list_lb$ifvals, c(2, 3), mean)
#' ub <- apply(list_ub$ifvals, c(2, 3), mean)
#' 
#' get_eps_zero(n, eps, lb, ub, ql, qu, phibar_lb, phibar_ub, delta = delta)
#' 
#' @seealso \code{\link{if_gamma}}, \code{\link{if_tau}}, 
#' \code{\link{get_ifvals}}.
#' @export


get_eps_zero <- function(n, eps, lb, ub, ql, qu, ifvals_lb, ifvals_ub, 
                         delta = 1, alpha = 0.05) {
  
  find_zero <- function(x) { which.min(abs(x)) }
  idx <- apply(lb * ub, 2, find_zero)
  ndelta <- length(delta)
  
  est <- eps[idx]
  lb_zero <- diag(as.matrix(lb[idx, ]))
  ub_zero <- diag(as.matrix(ub[idx, ]))
  ql_zero <- diag(as.matrix(ql[idx, ]))
  qu_zero <- diag(as.matrix(qu[idx, ]))

  if_l <-  sapply(1:ndelta, function(x) { ifvals_lb[, idx[x], x] } )
  if_u <-  sapply(1:ndelta, function(x) { ifvals_ub[, idx[x], x] } )
  
  phi_zero <- sweep(if_l, 2, ub_zero, "*") + sweep(if_u, 2, lb_zero, "*")
  der <- (ub_zero * ql_zero + lb_zero * qu_zero) ^ (-2)
  vars <- apply(phi_zero, 2, var) / n
  var_eps0 <- der * vars
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
