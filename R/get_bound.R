#' get_bound
#' 
#' get_bound is the main function to estimate the lower and upper bounds curves.
#' @param y nx1 outcome vector in [0, 1]
#' @param a nx1 treatment received vector
#' @param x nxp matrix of covariates
#' @param out_model a string specifying the outcome model for E(Y|A=a, X)
#' currently only "logistic" and "ranger" supported.
#' @param treat_model a string specifying the outcome model for P(A=1| X)
#' currently only "logistic" and "ranger" supported.
#' @param eps vector of arbitrary length specifying the values for the 
#' proportion of confounding where the lower and upper bounds curves are 
#' evaluated.
#' @param delta vector of arbitrary length specifying the maximum bias allowed
#' for the confounded units.
#' @param nsplits number of splits for the cross-fitting (default=5)
#' @param do_mult_boot boolean for whether uniform bands via the multiplier 
#' bootstrap need to computed (default is do_mult_boot=TRUE)
#' @importFrom magrittr %>%
#' @export

get_bound <- function(y, a, x, out_model, treat_model, eps=0, delta=0, 
                      nsplits=5, do_mult_boot=TRUE) {
  
  n <- length(y)
  s <- sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
  mu0hat <- mu1hat <- pi0hat <- pi1hat <- rep(NA, n)
  
  for (vfold in 1:nsplits) {
    
    train <- s != vfold
    test <- s == vfold
    if (nsplits == 1) {
      train <- test
    }
    
    ytrain <- y[train]
    atrain <- a[train]
    xtrain <- x[train, , drop=FALSE]
    xtest <- x[test, , drop=FALSE]
    
    mu0hat[test] <- get_muahat(y=ytrain, a=atrain, x=xtrain, newx=xtest, aval=0,
                               method=out_model)
    mu1hat[test] <- get_muahat(y=ytrain, a=atrain, x=xtrain, newx=xtest, aval=1,
                               method=out_model)
    pi1hat[test] <- get_piahat(a=atrain, x=xtrain, newx=xtest,
                               method=treat_model)
    pi0hat[test] <- 1-pi1hat[test]
    
  }
  
  psi0 <- if_gamma(y=y, a=a, aval=0, pia=pi0hat, mua=mu0hat)
  psi1 <- if_gamma(y=y, a=a, aval=1, pia=pi1hat, mua=mu1hat)
  nuhat <- as.matrix(psi1-psi0)
  
  tauhat_lb <- if_tau(y=y, a=a, pi0=pi0hat, pi1=pi1hat, mu0=mu0hat, mu1=mu1hat,
                      upper=FALSE)
  tauhat_ub <- tauhat_lb + 1
  
  glhats <- sapply(delta, get_g, pi0=pi0hat, pi1=pi1hat, mu0=mu0hat, mu1=mu1hat, 
                   upper=FALSE)
  guhats <- glhats + 1
  
  list_lb <- get_ifvals(eps=eps, upper=FALSE, nu=nuhat, tau=tauhat_lb, 
                        ghatmat=glhats)
  list_ub <- get_ifvals(eps=eps, upper=TRUE, nu=nuhat, tau=tauhat_ub, 
                        ghatmat=guhats)
  
  lambdaq_lb <- do.call(rbind, lapply(1:ncol(list_lb$quant), function(x) { 
    sweep(list_lb$lambda, 2, t(list_lb$quant)[x, ], "*") }))
  lambdaq_ub <- do.call(rbind, lapply(1:ncol(list_ub$quant), function(x) { 
    sweep(list_ub$lambda, 2, t(list_ub$quant)[x, ], "*") }))
  
  phibar_lb <- list_lb$ifvals-cbind(0, lambdaq_lb)
  phibar_ub <- list_ub$ifvals-cbind(0, lambdaq_ub)
  lb_estbar <- apply(phibar_lb, 2, mean)[-1]
  ub_estbar <- apply(phibar_ub, 2, mean)[-1]
  
  summarize_by <- function(x, fun) {
    out <- as.data.frame(x) %>% 
      dplyr::group_by(delta) %>% 
      dplyr::summarize_all(fun) %>%
      as.matrix()
    return(out)
  }
  
  lb_est <- summarize_by(list_lb$ifvals, "mean")
  ub_est <- summarize_by(list_ub$ifvals, "mean")
  lb_var <- summarize_by(phibar_lb, "var")
  ub_var <- summarize_by(phibar_ub, "var")
  
  # if(lb_est[which(eps==0 | delta==0)][1] > 0){
  #   no_zero <- 1*(lb_est > 0 & ub_est > 0)
  # } else{
  #   no_zero <- 1*(lb_est < 0 & ub_est < 0)
  # }
  no_zero <- 0
  if(length(delta)==1 & do_mult_boot){
    calpha_lb <- get_multboot(n=n, psihat=lb_estbar, 
                              sigmahat=lb_var[,-1],
                              ifvals=phibar_lb[,-1], alpha=0.05/2)
    calpha_ub <- get_multboot(n=n, psihat=-ub_estbar, 
                              sigmahat=ub_var[,-1],
                              ifvals=-phibar_ub[,-1], alpha=0.05/2)
    # stopifnot(calpha_ub > 1.96 && calpha_lb > 1.96)
    # print(calpha_ub); print(calpha_lb)
    eps_tilde_idx <- which.min(abs(sapply(2:(length(eps)+1), function(x) { 
      lb_est[, x]*ub_est[, x] } )))
    eps_tilde <- eps[eps_tilde_idx]
    lb_tilde <- lb_est[eps_tilde_idx]
    ub_tilde <- ub_est[eps_tilde_idx]
    ql_tilde <- list_lb$quant[eps_tilde_idx]
    qu_tilde <- list_ub$quant[eps_tilde_idx]
    phi_tilde <- phibar_lb[,eps_tilde_idx]*lb_tilde + phibar_ub[,eps_tilde_idx]*ub_tilde
    se_eps_tilde <- sqrt((lb_tilde*ql_tilde - ub_tilde*qu_tilde)^(-2)*var(phi_tilde)/n)
    eps_tilde_lo <- eps_tilde - 1.96*se_eps_tilde
    eps_tilde_hi <- eps_tilde + 1.96*se_eps_tilde
  } else {
    calpha_lb <- calpha_ub <- 1.96 
    eps_tilde <- eps_tilde_lo <- eps_tilde_hi <- NA
  }
  ci_lo <- lb_est - calpha_lb*sqrt(lb_var/n)
  ci_hi <- ub_est + calpha_ub*sqrt(ub_var/n)
  ci_lb_hi <- lb_est + calpha_lb*sqrt(lb_var/n)
  ci_ub_lo <- ub_est - calpha_ub*sqrt(ub_var/n)
  ci_lo_ptwise <- lb_est - 1.96*sqrt(lb_var/n)
  ci_hi_ptwise <- ub_est + 1.96*sqrt(ub_var/n)
  
  lb_est <- reshape2::melt(as.data.frame(lb_est), id.vars="delta", 
                 variable.name = "eps", value.name="lb")
  ub_est <- reshape2::melt(as.data.frame(ub_est), id.vars="delta", 
                 variable.name = "eps", value.name="ub")
  ci_lo <- reshape2::melt(as.data.frame(ci_lo), id.vars="delta", 
                variable.name = "eps", value.name="ci_lo")
  ci_hi <- reshape2::melt(as.data.frame(ci_hi), id.vars="delta", 
                variable.name = "eps", value.name="ci_hi")
  ci_lb_hi <- reshape2::melt(as.data.frame(ci_lb_hi), id.vars="delta", 
                   variable.name = "eps", value.name="ci_lb_hi")
  ci_ub_lo <- reshape2::melt(as.data.frame(ci_ub_lo), id.vars="delta", 
                   variable.name = "eps", value.name="ci_ub_lo")
  ci_lo_ptwise <- reshape2::melt(as.data.frame(ci_lo_ptwise), id.vars="delta", 
                       variable.name = "eps", value.name="ci_lo_ptwise")
  ci_hi_ptwise <- reshape2::melt(as.data.frame(ci_hi_ptwise), id.vars="delta", 
                       variable.name = "eps", value.name="ci_hi_ptwise")
  
  res <- data.frame(eps=lb_est$eps, delta=lb_est$delta,
                    lb=lb_est$lb, ub=ub_est$ub,
                    # no_zero=no_zero,
                    ci_lo=ci_lo$ci_lo, ci_hi=ci_hi$ci_hi,
                    ci_lb_hi=ci_lb_hi$ci_lb_hi, ci_ub_lo=ci_ub_lo$ci_ub_lo,
                    ci_lo_ptwise=ci_lo_ptwise$ci_lo_ptwise,
                    ci_hi_ptwise=ci_hi_ptwise$ci_hi_ptwise,
                    eps_tilde=eps_tilde,
                    eps_tilde_lo=eps_tilde_lo,
                    eps_tilde_hi=eps_tilde_hi)
  
  # res <- data.frame(tau=tauhat_lb[, 1], ghat=glhats[, 1])
  
  out <- list(bounds=res, var_ub=ub_var, var_lb=lb_var,
              ifvals_lb=list_lb$ifvals, ifvals_ub=list_ub$ifvals,
              mu0hat=mu0hat, mu1hat=mu1hat, pi0hat=pi0hat, pi1hat=pi1hat,
              nuhat=nuhat)
  return(out)
}