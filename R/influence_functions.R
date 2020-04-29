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
#' @param nu a list of size \code{nsplits}, where element j 
#' is a (num obs in split j)x1 matrix containing the influence function values 
#' for E{E(Y|A = 1, X) - E(Y|A = 0, X)} computed using regressions fns estimated 
#' using all obs except those in fold j and evaluated at obs in fold j.
#' @param tau a list of size \code{nsplits}, where element j 
#' is a (num obs in split j)x1 matrix containing the influence function values 
#' for E{g(eta)} computed using regressions fns estimated 
#' using all obs except those in fold j and evaluated at obs in fold j.
#' @param ghatmat a list of size \code{nsplits}, where element j 
#' is a (num obs in split j)xlength(delta) matrix containing the 
#' values of g(eta) for lower bound if \code{upper = FALSE}, for the upper
#' bound otherwise, computed using regressions fns estimated 
#' using all obs except those in fold j and evaluated at obs in fold j.
#' @param delta vector of delta values. \code{length}(delta) should match 
#' \code{ncol(ghatmat)}.  
#' 
#' @return A list containing:
#' \item{ifvals}{a n x length(eps) x length(delta) array containing 
#' the influence functions (varphi_l or varphi_u in manuscript) evaluated at 
#' the observed X as a function of epsilon and delta;}
#' \item{lambda}{a list of size \code{nsplits}, where element j
#' contains a nxlength(eps)xlength(delta) array containing 
#' the indicator ghatmat <= q (if \code{upper = FALSE}) or
#' ghatmat > q (if \code{upper = TRUE}), where q is eps-quantile or 1-eps 
#' quantile of ghatmat and, computed using regression functions
#' estimated using all folds but j and evaluated using obs at fold j.}
#' \item{quant}{a list of size \code{nsplits}, where element j is a
#' (num eps) x (num delta) matrix containing estimates of eps-quantile 
#' (if \code{upper = FALSE}) or (1-eps) quantiles (if \code{upper = TRUE}) of 
#' \code{ghatmat}, where the regression functions are computed and evaluated 
#' using all obs except those in fold j. See \code{\link{get_quant_g}}.}
#' \item{lambdaq}{a list of size \code{nsplits} containing the product 
#' \code{lambda} * \code{quant}.}

#' @seealso \code{\link{if_gamma}}, \code{\link{if_tau}}.
#' @export

get_ifvals <- function(n, neps, ndelta, upper, nu, tau, g, nsplits, quant,
                       nobs_fold) {
  
  indic_fn <- function(g, q, multiply_by_q, upper) {
    
    mult <- ifelse(multiply_by_q, q, 1)
    
    if(upper) {
      out <- mult * I(g > q) 
    } else {
      out <- mult * I(g <= q)
    }
    
    return(out)
    
  }
  
  .get_lambda <- function(delta_idx, g, q, upper, multiply_by_q, nobs_fold) {

    out <- vapply(q[, delta_idx], FUN = indic_fn, g = g[, delta_idx],
                  multiply_by_q =  multiply_by_q, upper = upper,
                  FUN.VALUE = array(0, dim = nobs_fold))

    return( out )
  }
  
  .get_lambda_wrapper <- function(g, q, upper, multiply_by_q, nobs_fold) {
    vapply(1:ndelta, FUN = .get_lambda, g = g, q = q, 
           multiply_by_q = multiply_by_q, upper = upper, nobs_fold = nobs_fold,
           FUN.VALUE = array(0, dim = c(nobs_fold, neps)))
  }
  
  lambda <- mapply(.get_lambda_wrapper, g, quant, nobs_fold,
                   MoreArgs = list(upper = upper, multiply_by_q = FALSE),
                   SIMPLIFY = FALSE)
  # print(sum(lambda[[9]]))
  lambdaq <- mapply(.get_lambda_wrapper, g, quant, nobs_fold,
                    MoreArgs = list(upper = upper, multiply_by_q = TRUE),
                    SIMPLIFY = FALSE)
  
  .get_ifs <- function(delta_idx, lambda, tau, nu, nobs_fold) {
    out <- sweep(sweep(lambda[, , delta_idx, drop = FALSE], 1, 
                       tau[, delta_idx, drop = FALSE], "*"), 1, nu, "+")
    out <- array(out, dim = c(nobs_fold, neps))
    return(out)
  }
  
  .get_ifs_wrapper <- function(lambda, tau, nu, nobs_fold){
    vapply(1:ndelta, FUN = .get_ifs, lambda = lambda, tau = tau, nu = nu, 
           nobs_fold = nobs_fold,
           FUN.VALUE = array(0, dim = c(nobs_fold, neps)))
  }
  
  ifs <- mapply(.get_ifs_wrapper, lambda, tau, nu, nobs_fold, SIMPLIFY = FALSE)
  
  out <- list(ifvals = ifs, lambda = lambda, lambdaq = lambdaq)
  
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
#' @param data data.frame of dimensions n x 8, with columns: y, a, ymin, ymax, 
#' pi0g (= P(A = 0 | X = x) if model "x" or 1-A if model "xa"),
#' pi1g (= P(A = 1 | X = x) if model "x" or A if model "xa")
#' mu1 (E(Y | A = 1, X = x)) and mu0 (E(Y | A = 0, X = x)).
#' @param delta vector of values for the delta parameter.
#' @param upper boolean (default is upper=FALSE) if TRUE if values for 
#' E(g_u(eta)) are returned, otherwise those for E(g_l(eta)) are returned
#' 
#' @return a length(y)xlength(delta) matrix containing the influence function 
#' values at each observed X of the parameter E(g(eta)).
#' 
#' @section Details:
#' As done in the paper, one can see that g(eta) for the lower bound is equal to
#' g(eta) for the upper bound minus delta * (ymax - ymin). Therefore the IFs for
#' E(g(eta)) follows the same relation. They are keep separated just for code
#' clarity. 
#' @seealso \code{\link{if_gamma}}, \code{\link{get_ifvals}}.
#' @export
if_tau <- function(data, delta, upper = FALSE) {
  
  y <- data$y
  a <- data$a
  ymin <- data$ymin
  ymax <- data$ymax
  pi0 <- data$pi0
  pi1 <- data$pi1
  mu0 <- data$mu0
  mu1 <- data$mu1
  
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
  tau <- (if2 - if1) %*% t(delta)
  colnames(tau) <- delta
  return( tau )
}

#' @title Estimation of the function g(eta)
#' 
#' @description \code{get_g} returns the values for the g_l(eta) 
#' (upper=FALSE) or g_u(eta) (upper=TRUE) at each X.
#' 
#' @param data data.frame of dimensions n x 6, with columns: ymin, ymax, 
#' pi0g (= P(A = 0 | X = x) if model "x" or 1-A if model "xa"),
#' pi1g (= P(A = 1 | X = x) if model "x" or A if model "xa")
#' mu1 (E(Y | A = 1, X = x)) and mu0 (E(Y | A = 0, X = x)).
#' @param delta vector of values for the delta parameter.
#' @param upper boolean (default is upper=FALSE) if TRUE values for 
#' g_u(eta) are returned, otherwise those for g_l(eta) are returned.
#' 
#' @return a nrow(data) x ndelta matrix containing the values of g(eta) at 
#' each observed X.
#' 
#' @section Details:
#' As done in the paper, one can see that g(eta) for the lower bound is equal to
#' g(eta) for the upper bound minus delta * (ymax - ymin). 
#' They are keep separated just for code clarity. 
#' @seealso \code{\link{if_gamma}}, \code{\link{get_ifvals}}, 
#' \code{\link{if_tau}}.
#' @export

get_g <- function(data, delta, upper) {
  # Computes the function g in the paper (Theorem 1)
  ymin <- data$ymin
  ymax <- data$ymax
  pi0g <- data$pi0g
  pi1g <- data$pi1g
  mu0hat <- data$mu0
  mu1hat <- data$mu1
  
  if(upper) {
    
    gfn <- pi0g * (ymax - mu1hat) - pi1g * (ymin - mu0hat)
    
  } else {
    
    gfn <- pi0g * (ymin - mu1hat) - pi1g * (ymax - mu0hat)
    
  }
  
  out <- gfn %*% t(delta)
  colnames(out) <- delta
  
  return(out)
}

#' @title Estimation of the quantiles of g(eta)
#' 
#' @description \code{get_quant_g} returns the quantiles for the function g(eta)
#' 
#' @param g a matrix n x delta g (generally evaluated at the X in the train set),
#' where each column represents the function g(eta) for a given value of delta.
#' @param eps vector of epsilon values (prop. unmeasured confounding)
#' @param min_g the min value of g. The zero-quantile is set to be the 
#' min(g) - 1e-10, so that for eps = 0, the ATE is point-identified. This is
#' because R takes the minimum to be the zero-quantile.
#' @param max_g the max value of g. The one-quantile is set to be the max(g).
#' 
#' @return a length(eps) x ncol(g) (= ndelta) matrix, where each entry  (i,j)
#' equals quantile(g[, delta[j]], p = eps[i])

#' @seealso \code{\link{if_gamma}}, \code{\link{get_ifvals}}, 
#' \code{\link{if_tau}}, \code{\link{get_g}}.
#' @export

get_quant_g <- function(g, eps, min_g, max_g) {
  # Given a function g of dim n x length(delta), it computes the eps-quantiles
  # for each column (i.e. each value of delta). 
  out <- matrix(apply(g[, , drop = FALSE], 2, quantile, p = eps, 
                      names = FALSE), nrow = length(eps), ncol = ncol(g))
  out[which(eps == 0), ] <- min_g - 1e-10
  out[which(eps == 1), ] <- max_g
  return(out)
}
