#' get_bound
#' 
#' get_bound is the main function to estimate the lower and upper bounds curves.
#' @param y nx1 outcome vector in [0, 1]
#' @param a nx1 treatment received vector
#' @param x nxp matrix of covariates
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
#' for the confounded units.
#' @param nsplits number of splits for the cross-fitting (default=5)
#' @param do_mult_boot boolean for whether uniform bands via the multiplier 
#' bootstrap need to computed (default is do_mult_boot=TRUE)
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
#' @importFrom magrittr %>%
#' @export

get_bound <- function(y, a, x, outfam, treatfam, model="x", eps=0, delta=0, 
                      nsplits=5, do_mult_boot=TRUE, do_eps_zero=TRUE, 
                      alpha=0.05, B=10000, nuis_fns=NULL, 
                      sl.lib=c("SL.earth","SL.gam","SL.glm", "SL.mean", 
                               "SL.ranger","SL.glm.interaction"),
                      plugin=FALSE) {
  .reshape_fn <- function(x, id, var, val) {
    out <- reshape2::melt(as.data.frame(x), id.vars=id, variable.name=var, 
                          value.name=val)
    rownames(out) <- c()
    return(out)
  }
  if(is.null(nuis_fns)) {
    nuis_fns <- do_crossfit(y=y, a=a, x=x, outfam=outfam, treatfam=treatfam,
                            nsplits=nsplits, sl.lib = sl.lib)
  }
  pi0hat <- nuis_fns[, "pi0"]
  pi1hat <- nuis_fns[, "pi1"]
  mu0hat <- nuis_fns[, "mu0"]
  mu1hat <- nuis_fns[, "mu1"]
  n <- length(y)
  ymax <- max(y)
  ymin <- min(y)
  ndelta <- length(delta)
  neps <- length(eps)
  
  # Estimate term E(mu1(X) - mu0(X)) (ATE under no unmeasured confounding)
  psi0 <- if_gamma(y=y, a=a, aval=0, pia=pi0hat, mua=mu0hat)
  psi1 <- if_gamma(y=y, a=a, aval=1, pia=pi1hat, mua=mu1hat)
  nuhat <- as.matrix(psi1-psi0)
  
  if(model == "x") {
    pi0g <- pi0hat
    pi1g <- pi1hat
  } else if(model == "xa") {
    pi0g <- 1 - a
    pi1g <- a
  } else {
    stop("model not supported!")
  }
  .get_g <- Vectorize(function(delta, upper) {
    if(upper){
      deltal <- pmax(-delta, ymin - mu0hat)
      deltau <- pmin(delta, ymax - mu1hat)
    } else {
      deltal <- pmin(delta, ymax - mu0hat)
      deltau <- pmax(-delta, ymin - mu1hat)
    }
    return(pi0g*deltau - pi1g*deltal)
  })
  glhat <- .get_g(delta, upper = FALSE)
  guhat <- .get_g(delta, upper = TRUE)
  colnames(glhat) <- colnames(guhat) <- delta
  
  if(!plugin) {
  tauhat_lb <- if_tau(y=y, a=a, pi0=pi0hat, pi1=pi1hat, mu0=mu0hat,
                      mu1=mu1hat, delta=delta, upper=FALSE)
  tauhat_ub <- if_tau(y=y, a=a, pi0=pi0hat, pi1=pi1hat, mu0=mu0hat,
                      mu1=mu1hat, delta=delta, upper=TRUE)
  } else {
    tauhat_lb <- glhat
    tauhat_ub <- guhat
  }
  colnames(tauhat_lb) <- colnames(tauhat_ub) <- delta
  
  list_lb <- get_ifvals(eps=eps, upper=FALSE, nu=nuhat, tau=tauhat_lb,
  ghatmat=glhat)
  list_ub <- get_ifvals(eps=eps, upper=TRUE, nu=nuhat, tau=tauhat_ub,
  ghatmat=guhat)
  
  lambdaq_lb <- do.call(rbind, lapply(1:ndelta, function(x) { 
    cbind(sweep(list_lb$lambda, 2, list_lb$quant[, x], "*"), delta[x]) }))
  colnames(lambdaq_lb) <- c(eps, "delta")
  lambdaq_lb <- .reshape_fn(lambdaq_lb, "delta", "eps", "lambdaq_lb")
  lambdaq_lb <- lambdaq_lb[order(lambdaq_lb$delta, lambdaq_lb$eps), ]
  
  lambdaq_ub <- do.call(rbind, lapply(1:ndelta, function(x) { 
    cbind(sweep(list_ub$lambda, 2, list_ub$quant[, x], "*"), delta[x]) }))
  colnames(lambdaq_ub) <- c(eps, "delta")
  lambdaq_ub <- .reshape_fn(lambdaq_ub, "delta", "eps", "lambdaq_ub")
  lambdaq_ub <- lambdaq_ub[order(lambdaq_ub$delta, lambdaq_ub$eps), ]
  
  ifvals_lb <- .reshape_fn(list_lb$ifvals, "eps", "delta", "ifvals_lb")
  ifvals_lb <- ifvals_lb[order(ifvals_lb$delta, ifvals_lb$eps), ]
  ifvals_ub <- .reshape_fn(list_ub$ifvals, "eps", "delta", "ifvals_ub")
  ifvals_ub <- ifvals_ub[order(ifvals_ub$delta, ifvals_ub$eps), ]
  
  phibar_lb <- ifvals_lb[, c("delta", "eps")]
  phibar_lb$ifvals <- ifvals_lb$ifvals_lb - lambdaq_lb$lambdaq_lb
  
  phibar_ub <- ifvals_ub[, c("delta", "eps")]
  phibar_ub$ifvals <- ifvals_ub$ifvals_ub - lambdaq_ub$lambdaq_ub
  
  lb_est <- summarize_by(ifvals_lb, "mean", delta, eps)
  ub_est <- summarize_by(ifvals_ub, "mean", delta, eps)
  lb_var <- summarize_by(phibar_lb, "var", delta, eps)
  ub_var <- summarize_by(phibar_ub, "var", delta, eps)
  
  if((length(delta)==1) && (delta==1)){
    if(do_mult_boot) {
      lb_estbar <- summarize_by(phibar_lb, "mean", delta, eps)
      ub_estbar <- summarize_by(phibar_ub, "mean", delta, eps)
      calpha_lb <- get_multboot(n=n, psihat=lb_estbar$ifvals, 
                                sigmahat=lb_var$ifvals,
                                ifvals=phibar_lb$ifvals, alpha=alpha/2, B=B)
      calpha_ub <- get_multboot(n=n, psihat=-ub_estbar$ifvals, 
                                sigmahat=ub_var$ifvals,
                                ifvals=-phibar_ub$ifvals, alpha=alpha/2, B=B)
    } else {
      calpha_lb <- calpha_ub <- qnorm(1-alpha/2) 
    }

  } else {
    calpha_lb <- calpha_ub <- qnorm(1-alpha/2)
  }
  if(do_eps_zero) {
    .get_eps_zero_wrap <- function(dd) {
      out <- get_eps_zero(n=n, eps=eps, lb=lb_est$ifvals_lb[lb_est$delta==dd], 
                          ub=ub_est$ifvals_ub[ub_est$delta==dd], 
                          ql=list_lb$quant[, colnames(list_lb$quant) == dd], 
                          qu=list_ub$quant[, colnames(list_ub$quant) == dd], 
                          ifvals_lb=phibar_lb[phibar_lb$delta == dd, ], 
                          ifvals_ub=phibar_ub[phibar_ub$delta == dd, ], 
                          delta=dd, alpha=alpha)
      return(out)
    }
    eps_list <- do.call(rbind, lapply(delta, .get_eps_zero_wrap))
  } else {
    eps_list <- list(eps_zero=NA, eps_zero_lo=NA, eps_zero_hi=NA)
  }
  ci_lb <- get_ci(lb_est$ifvals_lb, sqrt(lb_var$ifvals/n), calpha_lb)
  ci_ub <- get_ci(ub_est$ifvals_ub, sqrt(ub_var$ifvals/n), calpha_ub)
  ci_lb_ptwise <- get_ci(lb_est$ifvals_lb, sqrt(lb_var$ifvals/n), 
                         qnorm(1-alpha/2))
  ci_ub_ptwise <- get_ci(ub_est$ifvals_ub, sqrt(ub_var$ifvals/n), 
                          qnorm(1-alpha/2))
  cim04 <- get_im04(n=length(y), lb=lb_est$ifvals_lb, ub=ub_est$ifvals_ub,
                    sigma_l=sqrt(lb_var$ifvals), sigma_u=sqrt(ub_var$ifvals),
                    alpha=alpha)
  ci_ptwise_im04 <- matrix(c(lb_est$ifvals_lb - cim04*sqrt(lb_var$ifvals/n),
                             ub_est$ifvals_ub + cim04*sqrt(ub_var$ifvals/n)),
                           ncol=2, nrow=length(lb_est$eps))

  res <- data.frame(eps=lb_est$eps, delta=lb_est$delta,
                    lb=lb_est$ifvals_lb, ub=ub_est$ifvals_ub,
                    no_zero=1*(lb_est$ifvals_lb*ub_est$ifvals_ub > 0),
                    ci_lo=ci_lb[, 1], ci_hi=ci_ub[, 2],
                    ci_lb_hi=ci_lb[, 2], ci_ub_lo=ci_ub[, 1],
                    ci_lo_ptwise=ci_lb_ptwise[, 1],
                    ci_hi_ptwise=ci_ub_ptwise[, 2],
                    ci_lo_ptwise_im04=ci_ptwise_im04[, 1],
                    ci_hi_ptwise_im04=ci_ptwise_im04[, 2],
                    eps_zero=eps_list$eps_zero,
                    eps_zero_lo=eps_list$eps_zero_lo,
                    eps_zero_hi=eps_list$eps_zero_hi)
  
  out <- list(bounds=res, var_ub=ub_var, var_lb=lb_var,
              ifvals_lb=list_lb$ifvals, ifvals_ub=list_ub$ifvals,
              mu0hat=mu0hat, mu1hat=mu1hat, pi0hat=pi0hat, pi1hat=pi1hat,
              nuhat=nuhat, glhat=glhat, guhat=guhat, tauhat_l=tauhat_lb,
              tauhat_ub=tauhat_ub, phibar_lb=phibar_lb, phibar_ub=phibar_ub,
              lb_var=lb_var, ub_var=ub_var, mult_calpha_lb=calpha_lb,
              mult_calpha_ub=calpha_ub, im04_calpha=cim04)
  
  return(out)
}
