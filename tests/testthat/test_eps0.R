context("Epsilon 0")

test_that("coverage eps0 is correct", {
  
  require(pbapply)
  require(uniftest)
  # Data generating process:
  # X ~ Unif(0, 1); A ~ Bern(0.5); Y = A*X + (1-A)
  # so that pi(x) = 0.5; mu1(x) = x; mu0(x) = 1; g(x) = 1 - 0.5*x
  # eps-quantile of g(x) = 0.5*(1+eps)
  # E(mu1(x) - mu0(x)) = -0.5
  # upper bound is psiu(eps) = -0.25(eps^2 - 4*eps + 2)
  # hence eps0 is 0.5*(4 - sqrt(8)) = 0.5858
  # lower bound is psil(eps) = 0.25(eps^2 - 2*eps - 2)
  
  eps <- seq(0.57, 0.60, 0.0001)
  eps0_true <- 0.5*(4 - sqrt(8))
  
  # Start simulation to check coverage
  nsim <- 500
  n <- 10000
  alpha <- 0.2
  
  sim_fn <- function() {
    # Simulation function to estimate eps0
    # Generate data according to the model above
    x <- runif(n, 0, 1)
    a <- rbinom(n, 1, 0.5)
    y <- a*x + (1-a)
    
    # Regression functions are known exacty
    mu1x <- x
    mu0x <- rep(1, n)
    pi1x <- pi0x <- rep(0.5, n)
    nuis_fns <- matrix(NA, ncol=4, nrow=n,
                       dimnames=list(NULL, c("mu1", "mu0", "pi1", "pi0")))
    nuis_fns[, "mu1"] <- mu1x
    nuis_fns[, "mu0"] <- mu0x
    nuis_fns[, "pi1"] <- pi1x
    nuis_fns[, "pi0"] <- pi0x
    gx <- 1 - 0.5*x
    
    # P(A = a | X)
    piax <- a*pi1x + (1-a)*pi0x
    # E(Y|A = a, X)
    muax <- a*mu1x + (1-a)*mu0x
    nu <- (2*a-1) * (y-muax) / piax + mu1x - mu0x 
    tau <- (1-2*a) * (y-muax) / piax * (1-piax) + a*mu0x + (1-a)*(1-mu1x)
    phiu0 <- nu + tau * I( gx > (1-eps0_true/2))  
    psil0 <- 0.25*(eps0_true^2-2*eps0_true-2)
    tildephi0 <- (phiu0 - (1-eps0_true/2)*I( gx > (1-eps0_true/2))) * psil0  
    var_eps0 <- (psil0*(1-eps0_true/2))^(-2) * var(tildephi0)
    
    tmp <- get_bound(y=y, a=a, x=as.data.frame(x), 
                      outfam=NULL, treatfam=NULL, model="x", 
                      eps=eps, delta=1, nsplits=NULL, do_mult_boot=FALSE, 
                      do_eps_zero=TRUE, nuis_fns=nuis_fns, alpha=alpha)$eps_zero
    out <- c(tmp$est, tmp$ci_lo, tmp$ci_hi, var_eps0)
    names(out) <-  c("eps0", "eps0_lo", "eps0_hi", "var0")
    
    return(out)
  }
  # Simulation begins
  sims <- pbreplicate(nsim, sim_fn())
  cvg <- coverage(sims["eps0_lo", ], sims["eps0_hi", ], eps0_true, eps0_true)
  eps0 <- sims["eps0", ]
  bias0 <- bias(sims["eps0", ], eps0_true)
  var0 <- mean(sims["var0", ])
  pvs <- pnorm(sqrt(n)*(eps0 - eps0_true) / sqrt(var0))
  # Both coverage and bias are expressed in %
  expect_true(abs(cvg - 100*(1-alpha)) <= 5)
  expect_true(bias0 <= 1)
  # If distribution of eps0 is correct, then we would expect pvs to be unif
  expect_true(frosini.unif.test(pvs)$p.value > 0.05)
})

