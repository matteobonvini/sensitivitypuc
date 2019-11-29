context("Epsilon 0")


require(pbapply)
require(uniftest)

check_result <- function(sims, eps0_true, alpha, n) {
  # Check that coverage is correct and the the z-scores are correct
  cvg <- coverage(sims["eps0_lo", ], sims["eps0_hi", ], eps0_true, eps0_true)
  print(paste0("Coverage is ", cvg, "%"))
  eps0 <- sims["eps0", ]
  bias0 <- bias(sims["eps0", ], eps0_true)
  
  # pvalues using estimated variance
  pvs <- pnorm((eps0 - eps0_true) / sqrt(sims["var0", ]))
  # pvalues using true variance
  pvs_true_var <- pnorm((eps0 - eps0_true) / sqrt(sims["var0_true", ]))
  hist(pvs); hist(pvs_true_var)
  # Both coverage and bias are expressed in %
  expect_true(abs(cvg - 100*(1-alpha)) <= 5)
  # Bias should also be small (in % less than 1)
  expect_true(bias0 <= 1)
  # If distribution of eps0 is correct, then we would expect pvs to be unif,
  # regardless of whether we use true variance or sample variance
  expect_true(frosini.unif.test(pvs)$p.value > 0.05)
  expect_true(frosini.unif.test(pvs_true_var)$p.value > 0.05)
}

outnames <- c("eps0", "eps0_lo", "eps0_hi", "var0", "var0_true")

test_that("coverage eps0 is correct", {

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
    
    # This code computes the variance of eps_zero using truth
    gx <- 1 - 0.5*x
    piax <- a*pi1x + (1-a)*pi0x
    muax <- a*mu1x + (1-a)*mu0x
    nu <- (2*a-1) * (y-muax) / piax + mu1x - mu0x
    tau <- (1-2*a) * (y-muax) / piax * (1-piax) + a*mu0x + (1-a)*(1-mu1x)
    phiu0 <- nu + tau * I( gx > (1 - eps0_true/2))
    psil0 <- 0.25 * (eps0_true^2 - 2 * eps0_true - 2)
    tildephi0 <- (phiu0 - (1 - eps0_true/2)*I( gx > (1 - eps0_true/2)))
    var_eps0 <- (psil0 * (1 - eps0_true/2))^(-2) * var(tildephi0 * psil0) / n
    
    tmp <- get_bound(y = y, a = a, x = as.data.frame(x), ymin = 0, ymax = 1, 
                     outfam = NULL, treatfam = NULL, model = "x", 
                     eps = eps, delta = 1, nsplits = NULL, do_mult_boot = FALSE, 
                     do_eps_zero = TRUE, nuis_fns = nuis_fns, alpha = alpha)
    tmp <- tmp$eps_zero
    out <- c(tmp$est, tmp$ci_lo, tmp$ci_hi, tmp$var_eps0, var_eps0)
    names(out) <-  outnames
    
    return(out)
  }
  # Simulation begins
  sims <- pbreplicate(nsim, sim_fn())
  check_result(sims, eps0_true, alpha, n)
  
})

test_that("coverage eps0 is correct using simulation from paper", {
  
  source("simulation_true_regression_functions.R")
  ## Load true values ##
  truth <- readRDS("./data/truth_simulation.RData")
  n <- 5000
  eps <- attributes(truth)$eps0_seq
  alpha <- 0.10
  nsim <- 500
  
  sim_fn <- function() {
    
    df <- gen_data(n) 
    
    cnames <- c("mu1", "mu0", "pi1", "pi0")
    nuis_fns <- matrix(NA, ncol = 4, nrow = n, dimnames = list(NULL, cnames))
    nuis_fns[, "mu1"] <- mu1x(df$x1, df$x2)
    nuis_fns[, "mu0"] <- mu0x(df$x1, df$x2)
    nuis_fns[, "pi1"] <- pix(df$x1)
    nuis_fns[, "pi0"] <- 1 - nuis_fns[, "pi1"]
    
    y <- df$y; a <- df$a; x <- df[, c("x1", "x2")]
    
    piax <- a*nuis_fns[, "pi1"] + (1-a)*nuis_fns[, "pi0"]
    muax <- a*nuis_fns[, "mu1"] + (1-a)*nuis_fns[, "mu0"]
    nu <- (2*a-1) * (y-muax) / piax + nuis_fns[, "mu1"] - nuis_fns[, "mu0"]
    tau <- (1-2*a) * (y-muax) / piax * (1-piax) + a*nuis_fns[, "mu0"] + 
      (1-a)*(1-nuis_fns[, "mu1"]) 
    gux <- nuis_fns[, "pi0"] * (1-nuis_fns[, "mu1"]) + nuis_fns[, "pi1"] * nuis_fns[, "mu0"]
    glx <- gux - 1
    
    psiu0 <- truth[1, "ub0"]
    psil0 <- truth[1, "lb0"]
    ql0 <- truth[1, "ql0"]
    qu0 <- truth[1, "qu0"]
    eps0 <- truth[1, "eps_zero"]
    der <- psiu0 * ql0
    var_eps0 <- var(psiu0 * (nu + I(glx <= ql0) * (tau-ql0) - eps0)) * der^(-2)
    sl.lib <- c("SL.mean", "SL.glm", "SL.glm.interaction")
    tmp <- get_bound(y = y, a = a, x = x, ymin = 0, ymax = 1, outfam = binomial(),
                     treatfam = binomial(), model = "x", eps = eps, delta = 1,
                     nsplits = 5, do_mult_boot = FALSE, B = NULL,
                     do_eps_zero = TRUE, nuis_fns = NULL, alpha = alpha,
                     do_rearrange = FALSE, sl.lib = sl.lib)
    # tmp <- get_bound(y = y, a = a, x = x, ymin = 0, ymax = 1, outfam = NULL, 
    #                  treatfam = NULL, model = "x", eps = eps, delta = 1, 
    #                  nsplits = 5, do_mult_boot = FALSE, B = NULL, 
    #                  do_eps_zero = TRUE, nuis_fns = nuis_fns, alpha = alpha, 
    #                  do_rearrange = FALSE)
    
    tmp <- tmp$eps_zero
    out <- c(tmp$est, tmp$ci_lo, tmp$ci_hi, tmp$var_eps0, var_eps0 / n)
    names(out) <-  outnames
    
    return(out)
  }
  
  sims <- pbreplicate(nsim, sim_fn())
  check_result(sims, eps0_true = truth[1, "eps_zero"], alpha, n)
})

