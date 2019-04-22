#' truncate_prob
#' 
#' truncate_prob truncates a vector of observations in [0, 1] so that each 
#' observation is in [tol, 1-tol]
#' 
#' @param x a vector
#' @param tol scalar used to make each observation in x be in [tol, 1-tol]
#' @export
truncate_prob <- function(x, tol=0.05) {
  newx <- x
  newx[newx>=1-tol] <- 1-tol
  newx[newx<=tol] <- tol
  return(newx)
}