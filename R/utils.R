# Some helper functions for the experiments contained in the paper

#' summarize_by
#' Wrapper for group_by + summarize_all from dplyr
#' 
#' @param x a matrix
#' @param fun a function (as a string) that we want to apply to the values in 
#' the matrix
#' @param ... some columns in x which we want to apply fun conditional on
#' @export
summarize_by <- function(x, fun, ...) {
  out <- as.data.frame(x) %>% 
    dplyr::group_by(...) %>% 
    dplyr::summarize_all(fun) %>%
    as.data.frame()
  return(out)
}

#' get_ci
#' 
#' get_ci constructs a Wald-type confidence interval
#' 
#' @param muhat the n x 1 vector of means
#' @param sigmahat the variance
#' @param c the multiplier (e.g. 1.96)
#' @export
get_ci <- function(muhat, sigmahat, c) {
  lo <- muhat - c*sigmahat
  hi <- muhat + c*sigmahat
  return(cbind(lo, hi))
}

#' rmse
#' 
#' rmse computes RMSE, see Section 5 in paper
#' @param simmat a n x nsim matrix, where n is the number of epsilon points used
#' in evaluation and nsim is the number of simulations
#' @param truth a n x 1 vector containing the true values for the curve
#' @param n is the sample size
#' @export 
rmse <- function(simmat, truth, n) {
  out <- apply(sweep(simmat, 1, truth, "-"), 2, function(x) sqrt(sum(x^2)))
  return(mean(out) * sqrt(n))
}

#' coverage
#' 
#' coverage computes the frequentist coverage
#' 
#' @param lo a n x nsim matrix, where n is the number of epsilon points used
#' in evaluation and nsim is the number of simulations (generally lb)
#' If it is a nsim-dim vector, the length is interpreted as the number of 
#' simulations of a scalar quantity. 
#' @param hi a n x nsim matrix, where n is the number of epsilon points used
#' in evaluation and nsim is the number of simulations (generally ub)
#' If it is a nsim-dim vector, the length is interpreted as the number of 
#' simulations of a scalar quantity. 
#' @param truth_lo a n x 1 vector containing the true values for the curve (lb)
#' If lo is a vector, truth_lo must be scalar.
#' @param truth_hi a n x 1 vector containing the true values for the curve (ub)
#' If hi is a vector, truth_hi must be scalar.
#' @export
coverage <- function(lo, hi, truth_lo, truth_hi) {
  if(is.matrix(lo)) {
    lb_cont <- sweep(lo, 1, truth_lo, "<=") 
    ub_cont <- sweep(hi, 1, truth_hi, ">=")
    coverage <- mean(apply(lb_cont * ub_cont, 2, function(x) all(x == 1))) 
  }
  else {
    coverage <- mean(lo <= truth_lo & truth_hi <= hi)
  }
  return(100 * coverage)
}

#' bias
#' 
#' bias computes the bias, see Section 5 for definition
#' 
#' @param simmat a n x nsim matrix, where n is the number of epsilon points used
#' in evaluation and nsim is the number of simulations. If it is a nsim-dim 
#' vector, the length is interpreted as the number of simulations of a scalar
#' quantity. 
#' @param truth a n x 1 vector containing the true values for the curve.
#' If simmat is a vector, truth must be scalar.
#' @export
bias <- function(simmat, truth) {
  if(is.matrix(simmat)) {
    out <- apply(sweep(simmat, 1, truth, "-"), 2, mean)
  } else {
    out <- simmat - truth
  }
  return(100 * mean(abs(out)))
}
