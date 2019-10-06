context("Multiplier Bootstrap")

test_that("uniform coverage is correct", {
  
  require(pbapply)
  # Data generating process:
  # X ~ Unif(0, 1); A ~ Bern(0.5); Y = A*X + (1-A)
  # so that pi(x) = 0.5; mu1(x) = x; mu0(x) = 1; g(x) = 1 - 0.5*x
  # eps-quantile of g(x) = 0.5*(1+eps)
  # E(mu1(x) - mu0(x)) = -0.5
  # upper bound is psiu(eps) = -0.25(eps^2 - 4*eps + 2)
  # hence eps0 is 0.5*(4 - sqrt(8)) = 0.5858
  # lower bound is psil(eps) = 0.25(eps^2 - 2*eps - 2)
  
  eps <- seq(0, 1, 0.05)
  
  # Start simulation to check coverage
  nsim <- 500
  n <- 1000
  alpha <- 0.2
  
  psiu <- -0.25*(eps^2 - 4*eps + 2)
  psil <- 0.25*(eps^2 - 2*eps - 2)
  
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
    
    temp <- get_bound(y=y, a=a, x=as.data.frame(x), 
                      outfam=NULL, treatfam=NULL, model="x", 
                      eps=eps, delta=1, nsplits=NULL, do_mult_boot=TRUE,
                      B=10000, do_eps_zero=FALSE, nuis_fns=nuis_fns, 
                      alpha=alpha)
    
    out <- matrix(c(temp$bounds$ci_lo, temp$bounds$ci_hi, 
                    rep(temp$mult_calpha_lb, length(eps)),
                    rep(temp$mult_calpha_ub, length(eps))), 
                  ncol=4, nrow=length(eps),
                  dimnames=list(NULL, c("ci_lo", "ci_hi", "calpha_lb", 
                                        "calpha_ub")))
    return(out)
  }
  # Simulation begins
  sims <- pbreplicate(nsim, sim_fn())
  
  # Make sure the multiplier for CI coming out of bootstrap procedure is no
  # smaller than that for pointwise CI
  calpha_lb <- sims[1, "calpha_lb", ]
  expect_true(all(calpha_lb >= qnorm(1-alpha/2)))
  calpha_ub <- sims[1, "calpha_ub", ]
  expect_true(all(calpha_ub >= qnorm(1-alpha/2)))
  
  # Check uniform coverage (%) is okay
  cvg <- coverage(sims[, "ci_lo", ], sims[, "ci_hi", ], 
                  psil, psiu)
  expect_true(abs(cvg - 100*(1-alpha)) <= 5)
})

