context("Confidence Interval from Imbens & Manski (2004)")

test_that("Coverage of Imbens & Manski (2004) is correct", {
  
  require(pbapply)
  # Data generating process:
  # X ~ Unif(0, 1); A ~ Bern(0.5); Y = A*X + (1-A)
  # so that pi(x) = 0.5; mu1(x) = x; mu0(x) = 1; g(x) = 1 - 0.5*x
  # eps-quantile of g(x) = 0.5*(1+eps)
  # E(mu1(x) - mu0(x)) = -0.5
  # upper bound is psiu(eps) = -0.25(eps^2 - 4*eps + 2)
  # hence eps0 is 0.5*(4 - sqrt(8)) = 0.5858
  # lower bound is psil(eps) = 0.25(eps^2 - 2*eps - 2)
  
  eps <- seq(0, 0.1, 0.001)
  
  # Start simulation to check coverage
  nsim <- 1000
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
    nuis_fns <- matrix(NA, ncol = 4, nrow = n,
                       dimnames = list(NULL, c("mu1", "mu0", "pi1", "pi0")))
    nuis_fns[, "mu1"] <- mu1x
    nuis_fns[, "mu0"] <- mu0x
    nuis_fns[, "pi1"] <- pi1x
    nuis_fns[, "pi0"] <- pi0x
    gx <- 1 - 0.5*x
    
    temp <- get_bound(y = y, a = a, x = as.data.frame(x), ymin = 0, ymax = 1, 
                      outfam = NULL, treatfam = NULL, model = "x", eps = eps, 
                      delta = 1, nsplits = NULL, do_mult_boot = FALSE,
                      do_eps_zero = FALSE, nuis_fns = nuis_fns, alpha = alpha)
    # Machine precision 
    calpha <- round(temp$im04_calpha, 15)
    zalpha <- round(qnorm(1-alpha), 15)
    zalpha2 <- round(qnorm(1-alpha/2), 15)
    expect_true(all(zalpha <= calpha & calpha <= zalpha2))
    
    ci_lo_im04 <- round(temp$bounds[, "ci_im04_lo", 1], 15)
    ci_hi_im04 <- round(temp$bounds[, "ci_im04_hi", 1], 15)
    ci_lo_pt <- round(temp$bounds[, "ci_lb_lo_pt", 1], 15)
    ci_hi_pt <- round(temp$bounds[, "ci_ub_hi_pt", 1], 15)
    # IM 04 Ci should be contained in the pointwise bands covering the 
    # identification region
    expect_true(all(ci_lo_pt <= ci_lo_im04 & ci_hi_im04 <= ci_hi_pt))
    
    out <- matrix(c(ci_lo_pt, ci_hi_pt, ci_lo_im04, ci_hi_im04, calpha),
                  ncol = 5, nrow = length(eps),
                  dimnames = list(NULL, c("ci_lo_pt", "ci_lo_pt", "ci_lo_im04", 
                                          "ci_hi_im04", "calpha")))
    return(out)
  }
  # Simulation begins
  sims <- pbreplicate(nsim, sim_fn())

  # Check pointwise coverage is okay at, e.g., eps = eps[5]
  cvg <- Vectorize(function(x) {
    coverage(sims[x, "ci_lo_im04", ], sims[x, "ci_hi_im04", ], 
             psil[x], psil[x])
  })
  cvgs <- cvg(1:length(eps))
  # Both coverage and bias are expressed in %
  expect_true(all(abs(cvgs - 100*(1-alpha)) <= 5))
})

