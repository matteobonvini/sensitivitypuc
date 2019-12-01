###################################################
# True regression functions for DGP of Simulation #
###################################################
require(truncnorm)
require(cubature)

## Global parameters ##
effect <- 0.05
xlb <- -2
xub <- 2

# P(U = 1 | X)
pu_x <- function(x1) {
  return(0.5)
}

# P(S = 1 | X)
ps_x <- function(x1) {
  return(pnorm(x1))
}

# P(A = 1 | U, S , X)
pi_xus <- function(x1, u, s){
  return(0.5*pnorm(x1) + s*0.25 + (1-s)*0.5*u)
}

# P(A = 1 | U, X)
pi_xu <- function(x1, u) {
  ps1 <- ps_x(x1)
  return(pi_xus(x1, u, 0)*(1-ps1) + pi_xus(x1, u, 1)*ps1)
}

# P(A = 1 | X)
pix <- function(x1) {
  return(pi_xu(x1, 0)*(1-pu_x(x1)) + pi_xu(x1, 1)*pu_x(x1))
}

# P(Y^0 = 1 | X, U, S, A)
expect_y0 <- function(x1, x2, u){
  return(0.25 + 0.5 * pnorm(x1 + x2) - 0.5 * effect - 0.1 * u)
}

# P(Y^1 = 1 | X, U, S, A)
expect_y1 <- function(x1, x2, u){
  return(0.25 + 0.5 * pnorm(x1 + x2) + 0.5 * effect - 0.1 * u)
}

# P(Y = 1 | X, A = 0)
mu0x <- function(x1, x2) {
  u0 <- expect_y0(x1, x2, 0)*(1-pi_xu(x1, 0))/(1-pix(x1))
  u1 <- expect_y0(x1, x2, 1)*(1-pi_xu(x1, 1))/(1-pix(x1))
  return((1-pu_x(x1))*u0 + pu_x(x1)*u1)
}

# P(Y = 1 | X, A = 1)
mu1x <- function(x1, x2) {
  u0 <- expect_y1(x1, x2, 0)*pi_xu(x1, 0)/pix(x1)
  u1 <- expect_y1(x1, x2, 1)*pi_xu(x1, 1)/pix(x1)
  return((1-pu_x(x1))*u0 + pu_x(x1)*u1)
}

# g(eta) = pi(x)*mu0(x) + (1-pi(x))*(1-mu1(x)), uses worst-case delta = 1
gx <- function(x1, x2) {
  return((1-pix(x1))*(1-mu1x(x1, x2)) + pix(x1)*mu0x(x1, x2))
}

# Function to generate data
gen_data <- function(n) {
  
  require(truncnorm)
  
  x1 <- rtruncnorm(n, xlb, xub)
  x2 <- rtruncnorm(n, xlb, xub)
  u <- rbinom(n, 1, pu_x(x1))
  ps <- ps_x(x1)
  s <- rbinom(n, 1, ps)
  a <- rbinom(n, 1, pi_xus(x1, u, s))
  y0 <- rbinom(n, 1, expect_y0(x1, x2, u))
  y1 <- rbinom(n, 1, expect_y1(x1, x2, u))
  y <- (1-a) * y0 + a * y1
  
  return(data.frame(y=y, a=a, x1=x1, x2=x2, u=u, s=s))
}

# Density of (X1, X2) = product of truncnorm densities since X1 \ind X2. 
densx <- function(x1, x2) {
  dtruncnorm(x1, xlb, xub)*dtruncnorm(x2, xlb, xub)
}

# E_X(P(A = 1 | X))
pi1 <- adaptIntegrate(f=function(x) { 
  pix(x[1])*densx(x[1], x[2]) 
}, lowerLimit=c(xlb, xlb), upperLimit=c(xub, xub), absError=1e-15)$integral

# E_X(mu0(X))
mu0 <- adaptIntegrate(f=function(x) { 
  mu0x(x[1], x[2])*densx(x[1], x[2]) 
}, lowerLimit=c(xlb, xlb), upperLimit=c(xub, xub), absError=1e-15)$integral

# E_X(mu1(X))
mu1 <- adaptIntegrate(f=function(x) { 
  mu1x(x[1], x[2])*densx(x[1], x[2]) 
}, lowerLimit=c(xlb, xlb), upperLimit=c(xub, xub), absError=1e-15)$integral

# estimate the integral E[g(x) * I{g(x) <= quants}] (upper = FALSE) or 
# E[g(x) * I{g(x) > quants}] (upper = TRUE) as a sample average
get_g_term <- function(g, quants, upper) {
  if(upper) {
    lambdas <- outer(g, quants, ">")
  } else {
    lambdas <- outer(g, quants, "<=")
  }
  out <- apply(sweep(lambdas, 1, g, "*"), 2, mean)
  return(out)
}