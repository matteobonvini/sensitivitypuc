context("Multiplier Bootstrap")

require(pbapply)

# Make sure the multiplier for unif CI is no smaller than that for ptwise CI
check_multiplier <- function(sims, alpha) {
  calpha_lb <- sims[, "calpha_lb", ]
  expect_true(all(calpha_lb >= qnorm(1-alpha/2)))
  calpha_ub <- sims[, "calpha_ub", ]
  expect_true(all(calpha_ub >= qnorm(1-alpha/2)))
}

# Check uniform coverage (%) is okay
check_cvg <- function(sims, psil, psiu, alpha) {
  cvg <- coverage(sims[, "ci_lo", ], sims[, "ci_hi", ], psil, psiu)
  print(paste0("Coverage is ", cvg, "% when truth is ", 100 * (1-alpha), "%"))
  expect_true(abs(cvg - 100 * (1 - alpha)) <= 5)
}

test_that("uniform coverage is correct", {
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
  n <- 100
  alpha <- 0.50
  
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
    cnames <- c("mu1", "mu0", "pi1", "pi0")
    nuis_fns <- matrix(NA, ncol = 4, nrow = n, dimnames = list(NULL, cnames))
    nuis_fns[, "mu1"] <- mu1x
    nuis_fns[, "mu0"] <- mu0x
    nuis_fns[, "pi1"] <- pi1x
    nuis_fns[, "pi0"] <- pi0x
    gx <- 1 - 0.5*x
    
    cnames <- c("ci_lo", "ci_hi", "calpha_lb", "calpha_ub")
    
    get_res <- function(do_rearrange) {
      tmp <- get_bound(y = y, a = a, x = as.data.frame(x), ymin = 0, ymax = 1, 
                       outfam = NULL, treatfam = NULL, model = "x", 
                       eps = eps, delta = 1, nsplits = NULL, do_mult_boot = TRUE,
                       B = 10000, do_eps_zero = FALSE, nuis_fns = nuis_fns, 
                       alpha = alpha, do_rearrange = do_rearrange)
      bounds <- tmp$bounds[, , 1]
      
      out <- cbind(bounds[, "ci_lb_lo_unif"], bounds[, "ci_ub_hi_unif"],
                   tmp$mult_calpha_lb, tmp$mult_calpha_lb)
      colnames(out) <- cnames
      return(out)
    }
    
    res_r <- get_res(do_rearrange = TRUE)
    res_nr <- get_res(do_rearrange = FALSE)
    
    # To avoid randomness from sampling multiplier bootstrap, we do the
    # re-arranging on the non-rearranged bands for fair comparison
    covered_r <- all(I(sort(res_nr[, "ci_lo"], decreasing = TRUE) <= psil & 
                         psiu <= sort(res_nr[, "ci_hi"], decreasing = FALSE)))
    covered_nr <- all(I(res_nr[, "ci_lo"] <= psil & psiu <= res_nr[, "ci_hi"]))
    expect_true(covered_r >= covered_nr)
    
    out <- array(c(res_r, res_nr), dim = c(nrow(res_r), ncol(res_r), 2),
                 dimnames = list(NULL,  cnames, c("r", "nr")))
    return(out)
  }
  
  # Simulation begins
  sims <- pbreplicate(nsim, sim_fn())
  
  check_multiplier(sims[, , "r", ], alpha)
  check_multiplier(sims[, , "nr", ], alpha)
  check_cvg(sims_r, psil, psiu, alpha)
  check_cvg(sims_nr, psil, psiu, alpha)
})

test_that("uniform coverage is correct using truth in simulation paper", {
  
  source("simulation_true_regression_functions.R")
  ## Load true values ##
  truth <- readRDS("./data/truth_simulation.RData")
  n <- 5000
  eps <- attributes(truth)$eps_seq
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
    
    tmp <- get_bound(y = df$y, a = df$a, x = df[, c("x1", "x2")], ymin = 0,
                     ymax = 1, outfam = NULL, treatfam = NULL, model = "x",
                     eps = eps, delta = 1, nsplits = NULL, do_mult_boot = TRUE,
                     B = 10000, do_eps_zero = FALSE, nuis_fns = nuis_fns,
                     alpha = alpha, do_rearrange = FALSE)
    
    cnames <- c("ci_lo", "ci_hi", "calpha_lb", "calpha_ub")
    
    bounds <- tmp$bounds[, , 1]
    
    out <- cbind(bounds[, "ci_lb_lo_unif"], bounds[, "ci_ub_hi_unif"],
                 tmp$mult_calpha_lb, tmp$mult_calpha_ub)
    colnames(out) <- cnames
    return(out)
  }
  
  sims <- pbreplicate(nsim, sim_fn())
  
  check_multiplier(sims, alpha)
  check_multiplier(sims, alpha)
  cvg <- coverage(sims[, "ci_lo", ], sims[, "ci_hi", ], truth[, "lb"], truth[, "ub"])
  print(paste0("Coverage is ", cvg, "% when truth is ", 100 * (1-alpha), "%"))
  # tolerance for miscoverage
  tol <- 0.03
  # Because we allow alpha/2 error on each curve, by bonferroni coverage should
  # satisfy the following inequalities
  expect_true(100 *(1-alpha - tol) <= cvg & cvg <= 100 *(1-alpha/2 + tol))
})
