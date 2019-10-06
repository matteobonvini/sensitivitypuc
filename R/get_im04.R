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

get_im04 <- Vectorize(function(n, lb, ub, sigma_l, sigma_u, alpha) {
  imfn <- function(x) {
    out <- pnorm(x + sqrt(n)*(ub-lb)/max(sigma_l, sigma_u)) - pnorm(-x) - 
      (1-alpha)
    return(out)
  }
  calpha <- uniroot(f = imfn, interval = c(0, 5))$root
  return(calpha)
}, c("lb", "ub", "sigma_l", "sigma_u"))
