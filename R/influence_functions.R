#' @title Estimation of influence functions for the bounds on ATE
#' 
#' @description \code{get_ifvals} is the main function for getting the influence 
#' function values used to construct the estimator of either the lower and upper 
#' bound curves.
#' 
#' @param n sample size.
#' @param eps vector of arbitrary length specifying the values for the 
#' proportion of confounding where the lower and upper bounds curves are 
#' evaluated.
#' @param upper boolean for whether the upper or the lower bound curve needs to
#' be computed.
#' @param nu nx1 vector of influence function values for the parameter 
#' E(E(Y|A = 1, X) - E(Y|A = 0, X)).
#' @param tau nx1 vector of influence function values for the parameter
#' E(g_l(eta)) is upper=FALSE or E(g_u(eta)) if upper=TRUE.
#' @param ghatmat nxd matrix containing the values of g_l(eta)
#' if upper=FALSE or g_u(eta) if upper=TRUE, evaluated at each values of X. 
#' Each column of the matrix represents the vector of values for a specific 
#' value of delta, the maximum bias allowed for the confounded units.
#' @param delta vector of arbitrary length specifyin the values of delta used
#' to bound maximal confounding among S = 0 units. Default is delta = 1, which 
#' imposes no assumption if the outcome Y is bounded. 
#' 
#' @return A list containing:
#' \item{ifvals} a n x length(eps) x length(delta) array containing 
#' the influence functions (varphi_l or varphi_u in manuscript) evaluated at 
#' the observed X as a function of epsilon and delta; 
#' \item{lambda} a n x length(eps) x length(delta) array containing 
#' the indicator ghatmat <= q if upper = FALSE or ghatmat > q, where q is 
#' eps-quantile or (1-eps)-quantile of ghatmat.
#' \item{quant} a length(eps) x length(delta) matrix containing the estimated
#' quantiles of ghatmat as a function of eps and delta. If upper = FALSE, 
#' eps-quantiles are returned, if upper = TRUE, (1-eps)-quantiles are returned.
#' \item{lambdaq} the product lambda * quant as a function of eps and delta.
#' @examples 
#' eps <- seq(0, 0.1, 0.01)
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
#' head(list_lb$quant)
#' head(list_ub$quant)
#' head(list_lb$ifvals[, , 1])
#' head(list_ub$ifvals[, , 1])
#' head(list_lb$lambdaq[, , 1])
#' head(list_ub$lambdaq[, , 1])
#' 
#' @seealso \code{\link{if_gamma}}, \code{\link{if_tau}}.
#' @export

get_ifvals <- function(n, eps, upper, nu, tau, ghatmat, delta = 1) {
  
  ndelta <- length(delta)
  neps <- length(eps)
  
  if(upper) {
    seq_eps <- 1-eps
  } else {
    seq_eps <- eps
  }
  
  qhats <- apply(ghatmat, 2, quantile, p = seq_eps, names = FALSE)
  
  # The following is because R selects the min as 0-quantile, matters only in
  # finite samples
  qhats[which(seq_eps == 0), ] <- apply(ghatmat, 2, min) - 1e-10
  
  ineq_sign <- ifelse(upper, ">", "<=")
  .get_indicator <- function(x) {
    out <- 1 * sweep(ghatmat, 2, qhats[x, ], ineq_sign)
    return(out)
  }
  
  lambda <- aperm(sapply(1:length(eps), .get_indicator, simplify = "array"),
                  c(1, 3, 2))
  
  .get_lambdaq <- function(x) {
    out <- sweep(lambda[, , x], 2, qhats[, x], "*")
    return(out)
  }
  lambdaq <- sapply(1:ndelta, .get_lambdaq, simplify = "array")
  
  .get_ifs <- function(x) {
    out <- sweep(sweep(lambda[, , x], 1, tau[, x], "*"), 1, nu, "+")
    return(out)
  }
  ifs <- sapply(1:ndelta, .get_ifs, simplify = "array")
  
  out <- list(ifvals = ifs, lambda = lambda, quant = qhats, lambdaq = lambdaq)
  
  return(out)
}

#' @title Estimation of influence function for E(E(Y|A=a, X))
#' 
#' @description \code{if_gamma} computes the influence function values for the 
#' parameter E(E(Y|A=a, X))) evaluated at each X.
#' @param y nx1 vector of outcomes in [0, 1]
#' @param a nx1 vector of treatments
#' @param aval scalar, if aval=a then it returns the if vals for E(E(Y|A=a, X)))
#' @param pia nx1 vector with estimates for P(A=aval | X) evaluated at each X. 
#' @param mua nx1 vector with estimates for E(Y|A=aval, X) evaluated at each X.
#' 
#' @return a length(y)x1 matrix containing the influence function values at 
#' each observed X of the parameter E(E(Y|A=a, X)).
#' 
#' @seealso \code{\link{get_ifvals}}, \code{\link{if_tau}}.
#' @export
if_gamma <- function(y, a, aval, pia, mua) {
  
  return(as.matrix(I(a == aval) * (y - mua) / pia + mua))
  
}

#' @title Estimation of influence function for E(g(eta))
#' 
#' @description \code{if_tau} returns the influence function values for the 
#' parameter E(g_l(eta)) (upper=FALSE) or E(g_u(eta)) (upper=TRUE) at each X.
#' 
#' @param y nx1 vector of outcomes in [0, 1]
#' @param a nx1 vector of treatments
#' @param pi0 nx1 vector with estimates for P(A=0 | X=x) evaluated at each x
#' @param pi1 nx1 vector with estimates for P(A=1 | X=x) evaluated at each x
#' @param mu0 nx1 vector with estimates for E(Y |A=0, X=x) evaluated at each x
#' @param mu1 nx1 vector with estimates for E(Y |A=1, X=x) evaluated at each x
#' @param upper boolean (default is upper=FALSE) if TRUE if values for 
#' E(g_u(eta)) are returned, otherwise those for E(g_l(eta)) are returned
#' 
#' @return a length(y)x1 matrix containing the influence function values at 
#' each observed X of the parameter E(g(eta)).
#' 
#' @section Details:
#' As done in the paper, one can see that g(eta) for the lower bound is equal to
#' g(eta) for the upper bound minus delta * (ymax - ymin). Therefore the IFs for
#' E(g(eta)) follows the same relation. They are keep separated just for code
#' clarity. 
#' @seealso \code{\link{if_gamma}}, \code{\link{get_ifvals}}.
#' @export
if_tau <- function(y, a, ymin, ymax, pi0, pi1, mu0, mu1, upper = FALSE) {
  
  if(!upper) {
    # IF for E{ pi(x) * (ymax - mu0(x)) }
    if1 <- a * ymax - (pi1 / pi0 * (1 - a) * (y - mu0) + a * mu0) 
    # IF for E{ (1-pi(x)) * (ymin - mu1(x)) }
    if2 <- (1 - a) * ymin - (pi0 / pi1 * a * (y - mu1) + (1 - a) * mu1)
  } else {
    # IF for E{ pi(x) * (ymin - mu0(x)) }
    if1 <- a * ymin - (pi1 / pi0 * (1 - a) * (y - mu0) + a * mu0) 
    # IF for E{ (1-pi(x)) * (ymax - mu1(x)) }
    if2 <- (1 - a) * ymax - (pi0 / pi1 * a * (y - mu1) + (1 - a) * mu1)
  }
  
  return(as.matrix(if2 - if1))
}
