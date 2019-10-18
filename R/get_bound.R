#' get_bound
#' 
#' get_bound is the main function to estimate the lower and upper bounds curves.
#' @param y nx1 outcome vector in [0, 1]
#' @param a nx1 treatment received vector
#' @param x nxp matrix of covariates
#' @param ymin infimum of the support of y
#' @param ymax supremum of the support of y
#' @param outfam family specifying the error distribution for outcome 
#' regression, currently gaussian() or binomial(), should not specify link.
#' @param treatfam family specifying the error distribution for treatment 
#' regression, currently gaussian() or binomial(), should not specify link.
#' @param model a string specifying the assumption placed on S when 
#' computing the bounds. Currently only "x" and "xa".
#' @param eps vector of arbitrary length specifying the values for the 
#' proportion of confounding where the lower and upper bounds curves are 
#' evaluated.
#' @param delta vector of arbitrary length specifying the maximum bias allowed
#' for the confounded units. If delta is not equal to 1, then we use a plug-in
#' estimator. 
#' @param nsplits number of splits for the cross-fitting (default=5)
#' @param do_mult_boot boolean for whether uniform bands via the multiplier 
#' bootstrap need to computed. Default is do_mult_boot=TRUE if delta = 1. If 
#' delta is not equal to 1, then pointwise coverage of the curves is computed. 
#' @param do_eps_zero boolean for whether estimate of espilon_zero shoul be
#' computed (default is do_eps_zero=TRUE)
#' @param alpha confidence level, default is 0.05
#' @param B number of simulated radamacher random variables for unif bands.
#' @param nuis_fns Optional. A nx4 matrix specifying the estimated regression
#' functions evaluated at the observed x, columns should be named: pi0, pi1, 
#' mu1, mu0. Default is NULL so that regressions are estimated using the
#' SuperLearner via cross-fitting. 
#' @param plugin boolean for whether the estimator for the bounds is of plug-in
#' type (ensuring monotonicity in epsilon): uses g(etahat) rather than its
#' estimator based on influence functions when computing term that multiplies 
#' the indicator. See manuscript (g(etahat) instead of tauhat). 
#' @param do_parallel boolean for whether parallel computing should be used
#' @param ncluster number of clusters used if parallel computing is used.
#' @export

get_bound <- function(y, a, x, ymin, ymax, outfam, treatfam, model = "x", 
                      eps = 0, delta = 0, nsplits = 5, do_mult_boot = TRUE, 
                      do_eps_zero = TRUE, alpha = 0.05, B = 10000, 
                      nuis_fns = NULL, plugin = FALSE, 
                      sl.lib = c("SL.earth","SL.gam","SL.glm", "SL.mean", 
                                 "SL.ranger", "SL.glm.interaction"),
                      do_parallel = FALSE, ncluster = NULL) {
  
  if(is.null(nuis_fns)) {
    nuis_fns <- do_crossfit(y = y, a = a, x = x, outfam = outfam, 
                            treatfam = treatfam, nsplits = nsplits, 
                            sl.lib = sl.lib, ymin = ymin, ymax = ymax, 
                            do_parallel = do_parallel, ncluster = ncluster)
  }
  pi0hat <- nuis_fns[, "pi0"]
  pi1hat <- nuis_fns[, "pi1"]
  mu0hat <- nuis_fns[, "mu0"]
  mu1hat <- nuis_fns[, "mu1"]
  n <- length(y)
  ndelta <- length(delta)
  neps <- length(eps)
  
  # Estimate term E(mu1(X) - mu0(X)) (ATE under no unmeasured confounding)
  psi0 <- if_gamma(y = y, a = a, aval = 0, pia = pi0hat, mua = mu0hat)
  psi1 <- if_gamma(y = y, a = a, aval = 1, pia = pi1hat, mua = mu1hat)
  nuhat <- as.matrix(psi1 - psi0)
  
  if(model == "x") {
    
    pi0g <- pi0hat
    pi1g <- pi1hat
    
  } else if(model == "xa") {
    
    pi0g <- 1 - a
    pi1g <- a
    
  } else {
    
    stop("model not supported!")
    
  }
  
  glhat <- pi0g * (ymin - mu1hat) - pi1g * (ymax - mu0hat)
  guhat <- pi0g * (ymax - mu1hat) - pi1g * (ymin - mu0hat)
  glhat <- glhat %*% t(delta)
  guhat <- guhat %*% t(delta)
  
  colnames(glhat) <- colnames(guhat) <- delta
  
  if(!plugin) {
    
    tauhat_lb <- if_tau(y = y, a = a, ymin = ymin, ymax = ymax, pi0 = pi0hat, 
                        pi1 = pi1hat, mu0 = mu0hat, mu1 = mu1hat, upper = FALSE)
    tauhat_ub <- if_tau(y = y, a = a, ymin = ymin, ymax = ymax, pi0 = pi0hat, 
                        pi1 = pi1hat, mu0 = mu0hat, mu1 = mu1hat, upper = TRUE)
    tauhat_lb <- tauhat_lb %*% t(delta)
    tauhat_ub <- tauhat_ub %*% t(delta)
    
  } else {
    
    tauhat_lb <- glhat
    tauhat_ub <- guhat
    
  }
  colnames(tauhat_lb) <- colnames(tauhat_ub) <- delta
  
  list_lb <- get_ifvals(n = n, eps = eps, delta = delta, upper = FALSE, 
                        nu = nuhat, tau = tauhat_lb, ghatmat = glhat)
  list_ub <- get_ifvals(n = n, eps = eps, delta = delta, upper = TRUE, 
                        nu = nuhat, tau = tauhat_ub, ghatmat = guhat)
  
  lambda_l <- list_lb$lambda
  lambda_u <- list_ub$lambda
  lambdaq_l <- list_lb$lambdaq
  lambdaq_u <- list_ub$lambdaq
  ifvals_l <- list_lb$ifvals
  ifvals_u <- list_ub$ifvals
  q_l <- list_lb$quant
  q_u <- list_ub$quant
  
  phibar_l <- ifvals_l - lambdaq_l
  phibar_u <- ifvals_u - lambdaq_u
  
  est_l <- apply(ifvals_l, c(2, 3), mean)
  est_u <- apply(ifvals_u, c(2, 3), mean)
  var_l <- apply(phibar_l, c(2, 3), var)
  var_u <- apply(phibar_u, c(2, 3), var)
  
  if(do_mult_boot) {
    lb_estbar <- apply(phibar_l, c(2, 3), mean)
    ub_estbar <- apply(phibar_u, c(2, 3), mean)
    
    temp_fn <- function(x, psihat, sigmahat, ifvals) {
      out <- get_multboot(n = n, psihat = psihat[, x], 
                          sigmahat = sigmahat[, x], ifvals = ifvals[, , x], 
                          alpha = alpha / 2, B = B)
      return(out)
    }
    calpha_lb <- sapply(1:ndelta, temp_fn, psihat = lb_estbar, sigmahat = var_l,
                        ifvals = phibar_l)
    calpha_ub <- sapply(1:ndelta, temp_fn, psihat = -ub_estbar, 
                        sigmahat = var_u, ifvals = -phibar_u)
    calpha_lb <- matrix(rep(calpha_lb, neps), ncol = ndelta, nrow = neps, 
                        byrow = TRUE)
    calpha_ub <- matrix(rep(calpha_lb, neps), ncol = ndelta, nrow = neps, 
                        byrow = TRUE)
  } else {
    calpha_lb <- calpha_ub <- qnorm(1-alpha/2)
  }
  
  if(do_eps_zero) {
    eps_zero <- get_eps_zero(n = n, eps = eps, lb = est_l, ub = est_u, ql = q_l,
                             qu = q_u, ifvals_lb = phibar_l, 
                             ifvals_ub = phibar_u, delta = delta, alpha = alpha)
  } else {
    eps_zero <- NULL
  }
  
  get_ci <- function(muhat, sigmahat, c) {
    lo <- muhat - sigmahat * c
    hi <- muhat + sigmahat * c
    out <- aperm(array(cbind(lo, hi), dim = c(neps, ndelta, 2)), c(1, 3, 2))
    return(out)
  }
  
  ci_lb <- get_ci(est_l, sqrt(var_l/n), calpha_lb)
  ci_ub <- get_ci(est_u, sqrt(var_u/n), calpha_ub)
  ci_lb_pt <- get_ci(est_l, sqrt(var_l/n), qnorm(1-alpha/2))
  ci_ub_pt <- get_ci(est_u, sqrt(var_u/n), qnorm(1-alpha/2))
  
  temp <- array(cbind(est_l, est_u, sqrt(var_l), sqrt(var_u)), 
                dim = c(neps, ndelta, 4))
  cim04 <- apply(temp, c(1, 2), get_im04, n = n, alpha = alpha)
  
  ci_im04_l <- get_ci(est_l, sqrt(var_l/n), cim04)[, 1, ]
  ci_im04_u <- get_ci(est_u, sqrt(var_u/n), cim04)[, 2, ]      
  ci_im04 <- aperm(array(cbind(ci_im04_l, ci_im04_u), dim = c(neps, ndelta, 2)), 
                   c(1, 3, 2))
  
  # Return results in a user-friendly format
  temp_fn <- function(x) {
    out <- cbind(est_l[, x], est_u[, x], ci_lb[, , x], ci_ub[, , x],
                 ci_lb_pt[, , x], ci_ub_pt[, , x], ci_im04[, , x])
    return(out)
  }
  
  res <- sapply(1:ndelta, temp_fn, simplify = "array")
  
  dim2names <- c("lb", "ub", "ci_lb_lo_unif", "ci_lb_hi_unif", 
                 "ci_ub_lo_unif", "ci_ub_hi_unif", "ci_lb_lo_pt", "ci_lb_hi_pt",
                 "ci_ub_lo_pt", "ci_ub_hi_pt", "ci_im04_lo", "ci_im04_hi")
  dimnames(res) <- list(eps, dim2names, delta)
  
  out <- list(bounds = res, var_ub = var_u, var_lb = var_l,
              eps_zero = eps_zero,  q_l = q_l, q_u = q_u,
              lambda_l = lambda_l, lambda_u = lambda_u,
              ifvals_lb = list_lb$ifvals, ifvals_ub = list_ub$ifvals,
              nuis_fns = nuis_fns, nuhat = nuhat, glhat = glhat, guhat = guhat, 
              tauhat_l = tauhat_lb, tauhat_u = tauhat_ub, phibar_lb = phibar_l, 
              phibar_ub = phibar_u, mult_calpha_lb = calpha_lb,
              mult_calpha_ub = calpha_ub, im04_calpha = cim04)
  
  return(out)
}
