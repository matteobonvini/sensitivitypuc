#' get_im04
#' 
#' get_im04 computes c_alpha as in Imbens and Manski (2004)
#' @param n number of observations
#' @param lb estimate of lower bound
#' @param ub estimate of upper bound
#' @param sigma_l estimate of standard deviation of lb
#' @param sigma_u estimate of standard deviation of ub
#' @param alpha desired significance level
#' @export 

get_im04 <- function(dat, n, alpha) {
  lb <- dat[1]
  ub <- dat[2] 
  sigma_l <- dat[3] 
  sigma_u <- dat[4] 
  imfn <- function(x) {
    out <- pnorm(x + sqrt(n)*(ub-lb)/max(sigma_l, sigma_u)) - pnorm(-x) - 
      (1-alpha)
    return(out)
  }
  calpha <- uniroot(f = imfn, interval = c(0, 5), tol = .Machine$double.eps)$root
  return(calpha)
}
