#' @title Estimation of bounds on ATE as a function of the proportion of 
#' unmeasured confounding.
#' 
#' @description \code{get_bound} is the main function to estimate the lower and 
#' upper bounds curves as a function of eps, the proportion of 
#' unmeasured confounding.
#' @param y nx1 outcome vector in [0, 1].
#' @param a nx1 treatment received vector.
#' @param x nxp \code{data.frame} of covariates.
#' @param ymin infimum of the support of y.
#' @param ymax supremum of the support of y.
#' @param outfam family specifying the error distribution for outcome 
#' regression, currently \code{gaussian()} or \code{binomial()} supported. 
#' Link should not be specified.
#' @param treatfam family specifying the error distribution for treatment 
#' regression, currently \code{binomial()} supported.
#' Link should not be specified.
#' @param sl.lib character vector specifying which libraries to use for the SL.
#' @param model a string specifying the assumption placed on S when 
#' computing the bounds. Currently only "x" (\eqn{S \perp (Y, A) | X)} and "xa"
#' (\eqn{S \perp Y | A, X}). Default is model = "x".
#' @param eps vector of arbitrary length specifying the values for the 
#' proportion of confounding where the lower and upper bounds curves are 
#' evaluated. Default is 0 (no unmeasured confounding).
#' @param delta vector of arbitrary length specifyin the values of delta used
#' to bound maximal confounding among S = 0 units. Default is delta = 1, which 
#' imposes no assumption if the outcome Y is bounded. 
#' @param nsplits number of splits for the cross-fitting. Default is 5.
#' @param do_mult_boot boolean for whether uniform bands via the multiplier 
#' bootstrap need to computed. Default is do_mult_boot=TRUE.
#' @param do_eps_zero boolean for whether estimate of espilon_zero shoul be
#' computed. Default is do_eps_zero=TRUE.
#' @param alpha confidence level. Default is 0.05.
#' @param B number of rademacher rvs sampled. Default is 10000.
#' @param nuis_fns Optional. A nx4 matrix specifying the estimated regression
#' functions evaluated at the observed x, columns should be named: pi0, pi1, 
#' mu1, mu0. Default is NULL so that regressions are estimated using the
#' SuperLearner via cross-fitting. 
#' @param plugin boolean for whether the estimator for the bounds is of plug-in
#' type: uses g(etahat) rather than its estimator based on influence functions 
#' when computing term that multiplies the indicator.  So g(etahat) instead of 
#' tauhat. 
#' @param do_rearrange bollean for whether the precedure by Chernozhukov et al
#' (2008) should be applied to the estimators of the bounds and the CIs.
#' @param do_parallel boolean for whether parallel computing should be used
#' @param ncluster number of clusters used if parallel computing is used.
#' @param show_progress boolean for whether progress bar in estimating 
#' regression functions should be shown. Default is FALSE. Currently, only 
#' available if do_parallel is FALSE.
#' 
#' @return A list containing
#' \item{\code{bounds}}{a length(eps)x12xlength(delta) array, where, for each 
#' eps and delta, it has the estimates of lower bound (\code{lb}), upper bound
#' (\code{ub}), lower uniform band for lower bound (\code{ci_lb_lo_unif}), 
#' upper uniform band for upper bound  (\code{ci_ub_hi_unif}),
#' lower pointwise band for lower bound (\code{ci_lb_lo_pt}), 
#' upper pointwise band for lower bound  (\code{ci_lb_hi_pt}), 
#' lower pointwise band for upper bound (\code{ci_ub_lo_pt}), 
#' upper pointwise band for upper bound  (\code{ci_ub_hi_pt}), 
#' lower confidence band using Imbens & Manski (2004) procedure 
#' (\code{ci_im04_lo}), upper confidence band using Imbens & Manski (2004) 
#' procedure(\code{ci_im04_hi}).}
#' \item{\code{var_ub}}{estimate of the variance of the upper bound curve as fn 
#' of eps.}
#' \item{\code{var_lb}}{estimate of the variance of the lower bound curve as fn 
#' of eps.}
#' \item{\code{eps_zero}}{a length(delta)x5 \code{data.frame} with values of 
#' delta, estimate of eps0, max(0, ci_lo), min(1, ci_hi), variance of estimate 
#' of eps0.}
#' \item{\code{q_lb}}{a list of size \code{nsplits}, where element j is a
#' (num eps) x (num delta) matrix containing estimates of eps-quantile of 
#' g(etab) for lower bound computed using all obs except those in fold j}.
#' \item{\code{q_ub}}{a list of size \code{nsplits}, where element j is a
#' (num eps) x (num delta) matrix containing estimates of eps-quantile of 
#' g(etab) for upper bound computed using all obs except those in fold j}.
#' \item{\code{lambda_lb}}{a list of size \code{nsplits}, where element j
#' contains a nxlength(eps)xlength(delta) array containing 
#' the indicator ghatmat <= q, where q is eps-quantile of ghatmat and 
#' ghatmat is g(eta) for lower bound, computed using regression functions
#' estimated using all folds but j and evaluated using obs at fold j.}
#' \item{\code{lambda_ub}}{a list of size \code{nsplits}, where element j
#' contains a nxlength(eps)xlength(delta) array containing 
#' the indicator ghatmat > q, where q is (1-eps)-quantile of ghatmat and 
#' ghatmat is g(eta) for upper bound, computed using regression functions
#' estimated using all folds but j and evaluated using obs at fold j.}
#' \item{\code{ifvals_lb}}{a n x length(eps) x length(delta) array containing 
#' the influence functions for lower bound evaluated at the observed X as a 
#' function of epsilon and delta.}
#' \item{\code{ifvals_ub}}{a n x length(eps) x length(delta) array containing 
#' the influence functions for upper bound evaluated at the observed X as a 
#' function of epsilon and delta.}
#' \item{\code{nuis_fns}}{a list containing estimates of regression functions
#' evaluated at test obs and train obs. See \code{\link{do_crossfit}}.}
#' \item{\code{nuhat}}{a nx1 matrix containing the influence function values 
#' for E{E(Y|A = 1, X) - E(Y|A = 0, X)} computed using regressions fns estimated 
#' using all obs except those in fold j and evaluated at obs in fold j.}
#' \item{\code{glhat}}{a list of size \code{nsplits}, where element j
#' is (num obs in split j) x length(delta) matrix containing estimates of g(eta) 
#' for lower bound computed using regressions fns estimated using all obs except 
#' those in fold j and evaluated at obs in fold j.}
#' \item{\code{guhat}}{a list of size \code{nsplits}, where element j
#' is (num obs in split j) x length(delta) matrix containing estimates of g(eta) 
#' for upper bound computed using regressions fns estimated using all obs except 
#' those in fold j and evaluated at obs in fold j.}
#' \item{\code{glhat_train}}{a list of size \code{nsplits}, where element j
#' is (n - num obs in split j) x length(delta) matrix containing estimates of 
#' g(eta) for lower bound computed using regressions fns estimated using all obs 
#' except those in fold j and evaluated at obs \emph{not} in fold j.}
#' \item{\code{guhat_train}}{a list of size \code{nsplits}, where element j
#' is (n - num obs in split j) x length(delta) matrix containing estimates of 
#' g(eta) for upper bound computed using regressions fns estimated using all obs 
#' except those in fold j and evaluated at obs \emph{not} in fold j.}
#' \item{\code{tauhat_lb}}{a list of size \code{nsplits}, where element j 
#' is (num obs in split j) x 1 matrix containing the influence function
#'  values of the parameter E(g(eta)) with g(eta) for the lower bound computed
#'  using regression fns estimated from all obs except those in fold j and 
#'  evaluated at obs in fold j.}
#' \item{\code{tauhat_ub}}{a list of size \code{nsplits}, where element j 
#' is (num obs in split j) x 1 matrix containing the influence function
#'  values of the parameter E(g(eta)) with g(eta) for the upper bound computed
#'  using regression fns estimated from all obs except those in fold j and 
#'  evaluated at obs in fold j.}
#' \item{\code{phibar_lb}}{a n x length(eps) x length(delta) array containing
#' values for \code{ifvals_lb} - \code{lambda_lb} * \code{q_lb}.}
#' \item{\code{phibar_ub}}{a n x length(eps) x length(delta) array containing
#' values for \code{ifvals_ub} - \code{lambda_ub} * \code{q_ub}.}
#' \item{\code{phibar_lb_fold}}{a list of size \code{nsplits}, where element j
#' is n x length(eps) x length(delta) array containing
#' values for \code{ifvals_lb} - \code{lambda_lb} * \code{q_lb} computed using
#' regression functions estimated from all obs except those in fold j and 
#' evaluated at obs in fold j.}
#' \item{\code{phibar_ub_fold}}{a list of size \code{nsplits}, where element j
#' is n x length(eps) x length(delta) array containing
#' values for \code{ifvals_ub} - \code{lambda_ub} * \code{q_ub} computed using
#' regression functions estimated from all obs except those in fold j and 
#' evaluated at obs in fold j.}
#' \item{\code{mult_calpha_lb}}{a ndelta-dimensional vector containing 
#' \code{calpha} equal to the z-score used to construct uniform bands for the 
#' lower bound of the form psi(eps) \eqn{\pm} \code{calpha} * sigma(eps).}
#' \item{\code{mult_calpha_ub}}{a ndelta-dimensional vector containing 
#' \code{calpha} equal to the z-score used to construct uniform bands for the 
#' upper bound of the form psi(eps) \eqn{\pm} \code{calpha} * sigma(eps).}
#' \item{\code{im04_calpha}}{a ndelta-dimensional vector containing to the 
#' z-score used to construct the confidence interval for partially identified 
#' ATE as in Imbens & Manski (2004).}
#' 
#' @section Details:
#' As done in the paper, one can see that g(eta) for the lower bound is equal to
#' g(eta) for the upper bound minus delta * (ymax - ymin). Therefore the IFs for
#' E(g(eta)) follows the same relation. They are keep separated just for code
#' clarity. 
#' 
#' @examples 
#' n <- 1000
#' eps <- seq(0, 1, 0.001)
#' delta <- c(0.25, 0.5, 1)
#' a <- rbinom(n, 1, 0.5)
#' x <- as.data.frame(matrix(rnorm(2*n), ncol = 2, nrow = n))
#' ymin <- 0
#' ymax <- 1
#' y <- runif(n, ymin, ymax)
#' res <- get_bound(y = y, a = a, x = x, ymin = ymin, ymax = ymax, 
#'                  outfam = gaussian(),  treatfam = binomial(), 
#'                  model = "x", eps = eps, delta = delta, 
#'                  do_mult_boot = TRUE, do_eps_zero = TRUE, nsplits = 5, 
#'                  alpha = 0.05, B = 1000, sl.lib = "SL.glm")
#' print(res$eps_zero)
#' print(head(res$bounds[, , 1]))
#' 
#' @references Imbens, G. W., & Manski, C. F. (2004). Confidence intervals for 
#' partially identified parameters. \emph{Econometrica}, 72(6), 1845-1857.
#' @references Van der Laan, M. J., Polley, E. C., & Hubbard, A. E. (2007). 
#' Super learner. 
#' \emph{Statistical applications in genetics and molecular biology}, 6(1).
#' @references Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, 
#' E., Hansen, C., & Newey, W. K. (2016). Double machine learning for treatment 
#' and causal parameters (No. CWP49/16). \emph{cemmap working paper}.
#' @references Kennedy, E. H. (2019). Nonparametric causal effects based on 
#' incremental propensity score interventions. \emph{Journal of the American 
#' Statistical Association}, 114(526), 645-656.
#' @references Chernozhukov, V., Fernandez-Val, I., & Galichon, A. (2009). 
#' Improving point and interval estimators of monotone functions by 
#' rearrangement. \emph{Biometrika}, 96(3), 559-575.
#' 
#' @export

get_bound <- function(y, a, x, ymin, ymax, outfam, treatfam, sl.lib, 
                      model = "x", eps = 0, delta = 0, nsplits = 5, 
                      do_mult_boot = TRUE, do_eps_zero = TRUE, alpha = 0.05, 
                      B = 10000, nuis_fns = NULL, plugin = FALSE, 
                      do_rearrange = FALSE, do_parallel = FALSE, 
                      ncluster = NULL, show_progress = FALSE) {
  
  if(is.null(nuis_fns)) {
    nuis_fns <- do_crossfit(y = y, a = a, x = x, outfam = outfam, 
                            treatfam = treatfam, nsplits = nsplits, 
                            sl.lib = sl.lib, ymin = ymin, ymax = ymax, 
                            do_parallel = do_parallel, ncluster = ncluster,
                            show_progress = show_progress)
  }
  
  pi0hat <- nuis_fns$test[, "pi0"]
  pi1hat <- nuis_fns$test[, "pi1"]
  mu0hat <- nuis_fns$test[, "mu0"]
  mu1hat <- nuis_fns$test[, "mu1"]
  
  pi0hat_train <- nuis_fns$train[, "pi0"]
  pi1hat_train <- nuis_fns$train[, "pi1"]
  mu0hat_train <- nuis_fns$train[, "mu0"]
  mu1hat_train <- nuis_fns$train[, "mu1"]
  
  folds_train <- as.factor(nuis_fns$train[, "fold"])
  folds_test <- as.factor(nuis_fns$test[, "fold"])
  nobs_fold <- table(folds_test)
  nsplits <- length(unique(folds_train))
  
  n <- length(y)
  ndelta <- length(delta)
  neps <- length(eps)
  order_obs <- nuis_fns$order_obs
  y <- y[order_obs]
  a <- a[order_obs]
  x <- x[order_obs, , drop = FALSE]
  
  # Estimate term E(mu1(X) - mu0(X)) (ATE under no unmeasured confounding)
  psi0 <- if_gamma(y = y, a = a, aval = 0, pia = pi0hat, mua = mu0hat)
  psi1 <- if_gamma(y = y, a = a, aval = 1, pia = pi1hat, mua = mu1hat)
  nuhat <- psi1 - psi0
  nuhat_fold <- split(nuhat, folds_test)
  
  if(model == "x") {
    
    pi0g <- pi0hat
    pi1g <- pi1hat
    
    pi0g_train <- pi0hat_train
    pi1g_train <- pi1hat_train
    
  } else if(model == "xa") {
    
    pi0g <- pi0g_train <- 1 - a
    pi1g <- pi1g_train <- a
    
  } else {
    
    stop("model not supported!")
    
  }
  
  datg <- data.frame(ymin = ymin, ymax = ymax, pi0g = pi0g, pi1g = pi1g, 
                     mu0 = mu0hat, mu1 = mu1hat)
  
  gl <- get_g(data = datg, delta = delta, upper = FALSE)
  gu <- get_g(data = datg, delta = delta, upper = TRUE)
  
  min_gl <- apply(gl, 2, min)
  min_gu <- apply(gu, 2, min)
  max_gl <- apply(gl, 2, max)
  max_gu <- apply(gu, 2, max)
  
  glhat <- by(gl, folds_test, identity, simplify = FALSE)
  guhat <- by(gu, folds_test, identity, simplify = FALSE)
  
  datg_train <- data.frame(ymin = ymin, ymax = ymax, pi0g = pi0g_train, 
                           pi1g = pi1g_train, mu0 = mu0hat_train, 
                           mu1 = mu1hat_train)
  glhat_train <- by(datg_train, folds_train, get_g, delta = delta, upper = FALSE,
                    simplify = FALSE)
  guhat_train <- by(datg_train, folds_train, get_g, delta = delta, upper = TRUE,
                    simplify = FALSE)
  
  qhats_lb <- lapply(glhat_train, get_quant_g, eps = eps, min_g = min_gl,
                     max_g = max_gl)
  qhats_ub <- lapply(guhat_train, get_quant_g, eps = 1 - eps, min_g = min_gu,
                     max_g = max_gu)
  
  if(!plugin) {
    
    dat_tau <- data.frame(y = y, a = a, ymin = ymin, ymax = ymax, pi0 = pi0hat,
                          pi1 = pi1hat, mu0 = mu0hat, mu1 = mu1hat)
    tauhat_lb <- by(dat_tau, folds_test, if_tau, delta = delta, upper = FALSE, 
                    simplify = FALSE)
    tauhat_ub <- by(dat_tau, folds_test, if_tau, delta = delta, upper = TRUE, 
                    simplify = FALSE)
    
  } else {
    
    tauhat_lb <- glhat
    tauhat_ub <- guhat
    
  }
  
  list_lb <- get_ifvals(n = n, neps = neps, ndelta = ndelta, upper = FALSE, 
                        nu = nuhat_fold, tau = tauhat_lb, g = glhat, 
                        quant = qhats_lb, nsplits = nsplits, 
                        nobs_fold = nobs_fold)
  
  list_ub <- get_ifvals(n = n, neps = neps, ndelta = ndelta, upper = TRUE, 
                        nu = nuhat_fold, tau = tauhat_ub, g = guhat,
                        quant = qhats_ub, nsplits = nsplits, 
                        nobs_fold = nobs_fold)
  
  lambda_l <- list_lb$lambda
  lambda_u <- list_ub$lambda
  lambdaq_l <- list_lb$lambdaq
  lambdaq_u <- list_ub$lambdaq
  ifvals_l <- list_lb$ifvals
  ifvals_u <- list_ub$ifvals
  
  phibar_l_fold <- Map("-", ifvals_l, lambdaq_l)
  phibar_u_fold <- Map("-", ifvals_u, lambdaq_u)
  
  .get_est_fold <- function(x, fun) { 
    array( apply(x, c(2, 3), fun), dim = c(neps, ndelta) ) 
    }
  
  est_l_fold <- array(vapply(ifvals_l, .get_est_fold, fun = mean,
                             FUN.VALUE = array(0, dim = c(neps, ndelta))),
                      dim = c(neps, ndelta, nsplits))
  est_u_fold <- array(vapply(ifvals_u, .get_est_fold, fun = mean,
                             FUN.VALUE = array(0, dim = c(neps, ndelta))),
                      dim = c(neps, ndelta, nsplits))
    
  est_l <- array(apply(est_l_fold, c(1, 2), mean), dim = c(neps, ndelta))
  est_u <- array(apply(est_u_fold, c(1, 2), mean), dim = c(neps, ndelta))
  
  phibar_l <- abind::abind(phibar_l_fold, along = 1)
  phibar_u <- abind::abind(phibar_u_fold, along = 1)
  
  var_l <- array(apply(phibar_l, c(2, 3), var), dim = c(neps, ndelta))
  var_u <- array(apply(phibar_u, c(2, 3), var), dim = c(neps, ndelta))
  
  if(do_mult_boot) {

    estbar_l_fold <- array(vapply(phibar_l_fold, .get_est_fold, fun = mean, 
                                  FUN.VALUE = array(0, dim = c(neps, ndelta))),
                           dim = c(neps, ndelta, nsplits))
    estbar_u_fold <- array(vapply(phibar_u_fold, .get_est_fold, fun = mean, 
                                  FUN.VALUE = array(0, dim = c(neps, ndelta))),
                           dim = c(neps, ndelta, nsplits))
    
    lb_estbar <- array(apply(estbar_l_fold, c(1, 2), mean), 
                       dim = c(neps, ndelta))
    ub_estbar <- array(apply(estbar_u_fold, c(1, 2), mean), 
                       dim = c(neps, ndelta))
    
    calpha_lb <- do_multboot(n = n, ndelta = ndelta, psihat = lb_estbar, 
                             sigmahat = var_l, ifvals = phibar_l,
                             alpha = alpha / 2, B = B)
    calpha_ub <- do_multboot(n = n, ndelta = ndelta, psihat = -ub_estbar, 
                             sigmahat = var_u, ifvals = -phibar_u, 
                             alpha = alpha / 2, B = B)
    calpha_lb <- matrix(rep(calpha_lb, neps), ncol = ndelta, nrow = neps, 
                        byrow = TRUE)
    calpha_ub <- matrix(rep(calpha_ub, neps), ncol = ndelta, nrow = neps, 
                        byrow = TRUE)
  } else {
    calpha_lb <- calpha_ub <- qnorm(1 - alpha / 2)
  }
  
  if(do_eps_zero) {

    eps_zero <- get_eps_zero(n = n, pt_est = mean(nuhat), eps = eps, lb = est_l, 
                             ub = est_u, ql = qhats_lb, qu = qhats_ub, 
                             ifvals_lb = phibar_l_fold, 
                             ifvals_ub = phibar_u_fold, 
                             delta = delta, alpha = alpha)
    
  } else {
    eps_zero <- NULL
  }
  
  get_ci <- function(muhat, sigmahat, c) {
    lo <- muhat - sigmahat * c
    hi <- muhat + sigmahat * c
    out <- aperm(array(cbind(lo, hi), dim = c(neps, ndelta, 2)), c(1, 3, 2))
    return(out)
  }
  
  ci_lb <- get_ci(est_l, sqrt(var_l / n), calpha_lb)
  ci_ub <- get_ci(est_u, sqrt(var_u / n), calpha_ub)
  
  ci_lb_pt <- get_ci(est_l, sqrt(var_l / n), qnorm(1 - alpha / 2))
  ci_ub_pt <- get_ci(est_u, sqrt(var_u / n), qnorm(1 - alpha / 2))
  
  temp <- array(cbind(est_l, est_u, sqrt(var_l), sqrt(var_u)), 
                dim = c(neps, ndelta, 4))
  cim04 <- apply(temp, c(1, 2), get_im04, n = n, alpha = alpha)
  
  ci_im04_l <- get_ci(est_l, sqrt(var_l/n), cim04)[, 1, , drop = FALSE]
  ci_im04_u <- get_ci(est_u, sqrt(var_u/n), cim04)[, 2, , drop = FALSE] 
  
  if(do_rearrange) {
    
    ci_lb <- array(apply(ci_lb, c(2, 3), sort, decreasing = TRUE),
                   dim = c(neps, 2, ndelta))
    ci_ub <- array(apply(ci_ub, c(2, 3), sort, decreasing = FALSE),
                   dim = c(neps, 2, ndelta))
    est_l <- array(apply(est_l, 2, sort, decreasing = TRUE),
                   dim = c(neps, ndelta))
    est_u <- array(apply(est_u, 2, sort, decreasing = FALSE),
                   dim = c(neps, ndelta))
    ci_im04_l <- array(apply(ci_im04_l, 2, sort, decreasing = TRUE),
                       dim = c(neps, ndelta))
    ci_im04_u <- array(apply(ci_im04_u, 2, sort, decreasing = FALSE),
                       dim = c(neps, ndelta))
    
  }
  
  ci_im04 <- aperm(array(cbind(ci_im04_l, ci_im04_u), dim = c(neps, ndelta, 2)), 
                   c(1, 3, 2))

  # Return results in a user-friendly format
  temp_fn <- function(x) {
    out1 <- array(c(est_l[, x], est_u[, x]), dim = c(neps, 2))
    out2 <- array(c(ci_lb[, 1 , x], ci_ub[, 2 , x], ci_lb_pt[, , x], 
                    ci_ub_pt[, , x], ci_im04[, , x]), dim = c(neps, 10))
    out <- array(c(out1, out2), dim = c(neps, 12))
    return(out)
  }
  
  res <- array(vapply(1:ndelta, temp_fn, 
                      FUN.VALUE = array(0, dim = c(neps, 12))),
               dim = c(neps, 10, ndelta))
  
  dim2names <- c("lb", "ub", "ci_lb_lo_unif", "ci_ub_hi_unif", "ci_lb_lo_pt", 
                 "ci_lb_hi_pt", "ci_ub_lo_pt", "ci_ub_hi_pt", "ci_im04_lo", 
                 "ci_im04_hi")
  dimnames(res) <- list(eps, dim2names, delta)
  
  out <- list(bounds = res, var_ub = var_u, var_lb = var_l,
              eps_zero = eps_zero,  q_lb = qhats_lb, q_ub = qhats_ub,
              lambda_lb = lambda_l, lambda_ub = lambda_u,
              ifvals_lb = ifvals_l, ifvals_ub = ifvals_u,
              nuis_fns = nuis_fns, nuhat = nuhat, glhat = glhat, guhat = guhat,
              glhat_train = glhat_train, guhat_train = guhat_train, 
              tauhat_l = tauhat_lb, tauhat_u = tauhat_ub, phibar_lb = phibar_l, 
              phibar_ub = phibar_u, phibar_lb_fold = phibar_l_fold,
              phibar_ub_fold = phibar_u_fold, mult_calpha_lb = calpha_lb,
              mult_calpha_ub = calpha_ub, im04_calpha = cim04)
  
  return(out)
}
