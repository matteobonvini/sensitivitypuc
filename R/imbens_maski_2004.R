#' @title Estimation of multiplier for confidence interval for ATE as described
#' in Imbens and Manski (2004)
#' 
#' @description \code{get_im04} computes c_alpha as in Imbens and Manski (2004).
#' @param dat vector containing estimate of lower bound, estimate of upper bound,
#' estimate of standard deviation of lb, estimate of standard deviation of ub in
#' this order.
#' @param n sample size.
#' @param alpha desired significance level (capped at \code{pnorm(5)}). 
#' 
#' @return a scalar equal to the multiplier used to construct the confidence 
#' interval for partially identified ATE. 
#' 
#' @examples 
#' n <- 100
#' dat <- c(-0.4, 0.3, 0.2, 0.1)
#' get_im04(dat, n, 0.05) # large identification interval width
#' dat2 <- c(-0.4, -0.395, 0.2, 0.1)
#' get_im04(dat2, n, 0.05) # small identification interval width
#' 
#' @references Imbens, G. W., & Manski, C. F. (2004). Confidence intervals for 
#' partially identified parameters. \emph{Econometrica}, 72(6), 1845-1857.
#' @export 

get_im04 <- function(dat, n, alpha) {
  
  lb <- dat[1]
  ub <- dat[2] 
  sigma_l <- dat[3] 
  sigma_u <- dat[4] 
  
  imfn <- function(x) {
    out <- pnorm(x + sqrt(n) * (ub - lb) / max(sigma_l, sigma_u)) - pnorm(-x) - 
      (1 - alpha)
    return(out)
  }
  
  calpha <- uniroot(f = imfn, interval = c(0, 5), 
                    tol = .Machine$double.eps)$root
  
  return(calpha)
}
