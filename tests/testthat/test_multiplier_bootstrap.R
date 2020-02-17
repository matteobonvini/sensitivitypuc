context("Multiplier Bootstrap")

require(pbapply)

# Make sure the multiplier for unif CI is no smaller than that for ptwise CI
check_multiplier <- function(sims, alpha) {
  calpha_lb <- sims[, "calpha_lb", ]
  expect_true(all(calpha_lb >= qnorm(1 - alpha / 2)))
  calpha_ub <- sims[, "calpha_ub", ]
  expect_true(all(calpha_ub >= qnorm(1 - alpha / 2)))
}

# Check uniform coverage is okay
check_cvg <- function(nsim, eps, truth_lo, truth_hi, ci_lo, ci_hi, alpha, tol,
                      point_ident = TRUE) {
  truth_mat_lo <- matrix(truth_lo, nrow = length(eps), ncol = nsim, 
                         byrow = FALSE)
  truth_mat_hi <- matrix(truth_hi, nrow = length(eps), ncol = nsim, 
                         byrow = FALSE)
  is_contained <- I(ci_lo <= truth_mat_lo & truth_mat_hi <= ci_hi)
  cvg <- mean(apply(is_contained, 2, all))
  print(paste0("Coverage is ", cvg, " when truth is ", (1 - alpha)))
  if(point_ident) {
    expect_true(abs(cvg - (1 - alpha)) <= tol)
  } else {
    expect_true(1 - alpha - tol <= cvg & cvg <= 1 - alpha / 2 + tol)
  }
}

nsim <- 500

test_that("uniform coverage is correct", {
  # Data generating process:
  # X ~ Unif(0, 1); A ~ Bern(0.5); Y = A*X + (1-A)
  # so that pi(x) = 0.5; mu1(x) = x; mu0(x) = 1; g(x) = 1 - 0.5*x
  # eps-quantile of g(x) = 0.5*(1+eps)
  # E(mu1(x) - mu0(x)) = -0.5
  # upper bound is psiu(eps) = -0.25(eps^2 - 4*eps + 2)
  # hence eps0 is 0.5*(4 - sqrt(8)) = 0.5858
  # lower bound is psil(eps) = 0.25(eps^2 - 2*eps - 2)
  
  eps <- seq(0, 1, 0.1)
  
  # Start simulation to check coverage
  n <- 2000
  alpha <- 0.2
  
  psiu <- -0.25*(eps^2 - 4*eps + 2)
  psil <- 0.25*(eps^2 - 2*eps - 2)
  
  sim_fn <- function() {
    # Simulation function to estimate eps0
    # Generate data according to the model above
    x <- runif(n, 0, 1)
    a <- rbinom(n, 1, 0.5)
    y <- a*x + (1-a)
    
    mu1x <- x
    mu0x <- rep(1, n)
    pi1x <- pi0x <- rep(0.5, n)
    cnames <- c("mu1", "mu0", "pi1", "pi0")
    nuis_fns <- matrix(NA, ncol = 7, nrow = 2 * n,
                       dimnames = list(NULL, c("mu1", "mu0", "pi1", "pi0",
                                               "unit_num", "fold", "is_test")))
    nuis_fns[, "mu1"] <- rep(mu1x, 2)
    nuis_fns[, "mu0"] <- rep(mu0x, 2)
    nuis_fns[, "pi1"] <- rep(pi1x, 2)
    nuis_fns[, "pi0"] <- rep(pi0x, 2)
    
    nuis_fns[, "unit_num"] <- rep(1:n, 2)
    nuis_fns[, "fold"] <- 1
    nuis_fns[, "is_test"] <- c(rep(1, n), rep(0, n))
    
    nuis_fns_list <- list(test = nuis_fns[1:n, ], 
                          train = nuis_fns[(n+1):(2*n), ],
                          order_obs = nuis_fns[1:n, "unit_num"])
    gx <- 1 - 0.5*x
    
    cnames <- c("ci_lb_lo", "ci_ub_hi", "calpha_lb", "calpha_ub")
    
    get_res <- function(do_rearrange) {
      tmp <- get_bound(y = y, a = a, x = as.data.frame(x), ymin = 0, ymax = 1, 
                       outfam = NULL, treatfam = NULL, model = "x", 
                       eps = eps, delta = 1, nsplits = 1, do_mult_boot = TRUE,
                       B = 10000, do_eps_zero = FALSE, nuis_fns = nuis_fns_list, 
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
    covered_r <- all(I(sort(res_nr[, "ci_lb_lo"], decreasing = TRUE) <= psil & 
                         psiu <= sort(res_nr[, "ci_ub_hi"], decreasing = FALSE)))
    covered_nr <- all(I(res_nr[, "ci_lb_lo"] <= psil & psiu <= res_nr[, "ci_ub_hi"]))
    expect_true(covered_r >= covered_nr)
    
    out <- array(c(res_r, res_nr), dim = c(nrow(res_r), ncol(res_r), 2),
                 dimnames = list(NULL,  cnames, c("r", "nr")))
    return(out)
  }
  
  # Simulation begins
  sims <- pbreplicate(nsim, sim_fn())
  
  tol <- 0.05
  
  check_multiplier(sims[, , "r", ], alpha)
  check_multiplier(sims[, , "nr", ], alpha)
  # because we compute the calpha z-scores via multiplier bootstrap to guarantee
  # one-sided coverage for each bound at level 1-alpha/2, we use -Inf and Inf as
  # needed to compute one-sided coverage.
  check_cvg(nsim, eps, psil, psil, sims[, "ci_lb_lo", "r", ], Inf, alpha / 2, tol)
  check_cvg(nsim, eps, psil, psil, sims[, "ci_lb_lo", "nr", ], Inf, alpha / 2, tol)
  check_cvg(nsim, eps, psiu, psiu, -Inf, sims[, "ci_ub_hi", "r", ], alpha / 2, tol)
  check_cvg(nsim, eps, psiu, psiu, -Inf, sims[, "ci_ub_hi", "nr", ], alpha / 2, tol)
  check_cvg(nsim, eps, psil, psiu, sims[, "ci_lb_lo", "r", ], 
            sims[, "ci_ub_hi", "r", ], alpha, tol, FALSE)
  check_cvg(nsim, eps, psil, psiu, sims[, "ci_lb_lo", "nr", ], 
            sims[, "ci_ub_hi", "nr", ], alpha, tol, FALSE)
  
})

test_that("uniform coverage is correct using truth in simulation paper", {
  
  source("simulation_true_regression_functions.R")
  ## Load true values ##
  truth <- readRDS("./data/truth_simulation.RData")
  n <- 5000
  eps <- attributes(truth)$eps_seq
  alpha <- 0.10
  nsplits <- 10
  
  sim_fn <- function() {
    
    df <- gen_data(n) 
    
    cnames <- c("mu1", "mu0", "pi1", "pi0",
                "unit_num", "fold", "is_test")
  
    if(nsplits == 1) {
      train_idx <- test_idx <- list(1:n)
    } else {
      s <- sample(rep(1:nsplits, ceiling(n/nsplits))[1:n])
      train_idx <- lapply(1:nsplits, function(x) which(x != s))
      test_idx <- lapply(1:nsplits, function(x) which(x == s))
    }
    test_set <- train_set <- c()
    
    for(i in 1:nsplits) {
      # Simulate cross-fitting even if nuisance functions are known exactly
      idx_test <- test_idx[[i]]
      idx_train <- train_idx[[i]]
      
      newtest <- cbind(mu1x(df$x1[idx_test], df$x2[idx_test]),
                       mu0x(df$x1[idx_test], df$x2[idx_test]),
                       pix(df$x1[idx_test]), 1 - pix(df$x1[idx_test]),
                       idx_test, i, 1)
      test_set <- rbind(test_set, newtest)
      
      newtrain <- cbind(mu1x(df$x1[idx_train], df$x2[idx_train]),
                        mu0x(df$x1[idx_train], df$x2[idx_train]),
                        pix(df$x1[idx_train]), 1 - pix(df$x1[idx_train]),
                        idx_train, i, 0)
      train_set <- rbind(train_set, newtrain)
    }
    
    colnames(test_set) <- cnames
    colnames(train_set) <- cnames
    
    nuis_fns_list <- list(test = test_set,
                          train = train_set,
                          order_obs = test_set[, "unit_num"])
    
    y <- df$y; a <- df$a; x <- df[, c("x1", "x2")]
    
    tmp <- get_bound(y = df$y, a = df$a, x = df[, c("x1", "x2")], ymin = 0,
                     ymax = 1, outfam = NULL, treatfam = NULL, model = "x",
                     eps = eps, delta = 1, nsplits = nsplits, do_mult_boot = TRUE,
                     B = 10000, do_eps_zero = FALSE, nuis_fns = nuis_fns_list,
                     alpha = alpha, do_rearrange = FALSE)
    
    cnames <- c("ci_lo", "ci_hi", "calpha_lb", "calpha_ub")
    
    bounds <- tmp$bounds[, , 1]
    
    out <- cbind(bounds[, "ci_lb_lo_unif"], bounds[, "ci_ub_hi_unif"],
                 tmp$mult_calpha_lb, tmp$mult_calpha_ub)
    colnames(out) <- cnames
    return(out)
  }
  
  sims <- pbreplicate(nsim, sim_fn())
  
  tol <- 0.05
  
  check_multiplier(sims, alpha)
  check_multiplier(sims, alpha)
  
  check_cvg(nsim, eps, truth[, "lb"], truth[, "lb"], sims[, "ci_lo", ], Inf, 
            alpha / 2, tol)
  check_cvg(nsim, eps, truth[, "ub"], truth[, "ub"], -Inf, sims[, "ci_hi", ], 
            alpha / 2 , tol)
  check_cvg(nsim, eps, truth[, "lb"], truth[, "ub"], sims[, "ci_lo", ], 
            sims[, "ci_hi", ], alpha, tol, FALSE)
  
})

test_that("Multiplier bootstrap coverage is correct with sample splitting", {
  
  # Data generating process:
  # X ~ Unif(0, 1); A ~ Bern(expit(x))
  # Y = t1 * A + t0 * (1 - A) + beta(aa, bb), t0 > t1
  # so that pi(x) = expit(x); mu1(x) = t1; mu0(x) = t0
  # E(mu1(x) - mu0(x)) = t1 - t0
  aa <- 3; bb <- 2; mean_beta <- (aa / (aa + bb))
  t0 <- 1; t1 <- 0.5
  bmin <- -mean_beta; bmax <- 1 + mean_beta; ymax <- bmax + t0; ymin <- bmin + t1
  
  expit <- function(x) { exp(x) / (1 + exp(x)) }
  logit <- function(x) { log( x / (1 - x) ) }
  
  gx_fn <- function(x) {
    (1 - expit(x)) * (ymax - t1) - expit(x) * (ymin - t0)
  } 
  
  ming <- min(gx_fn(0), gx_fn(1))
  maxg <- max(gx_fn(0), gx_fn(1))
  
  cdf <- Vectorize(function(u) {
    if(ming < u & u < maxg) {
      p <- ( - u + ymax - t1 ) / (ymax + ymin - t1 - t0)
      return(1 - logit(p))
    } else {
      if(u <= ming) return(0) else return(1)
    }
  })
  
  dens_gx <- function(u) {
    der <- (t0 + t1 - ymin - ymax) / ( (- t0 + ymin + u) * (t1 + u - ymax) )
    return(I(ming <= u & u <= maxg) * der)
  }
  
  get_quant <- Vectorize(function(eps) {
    
    if(eps == 1) {
      out <- maxg
    } else {
      out <- uniroot(f = function(q) { cdf(q) - eps }, lower = ming - 1e-5, 
                     upper = maxg + 1e-5)$root
    }
    
    return(out)
    
  })
  
  lb <- Vectorize(function(eps) {
    
    if(eps == 0) {
      quant <- get_quant(0) - 0.0001
    } else {
      quant <- get_quant(eps)
    }
    
    t1 - t0 + integrate(f = function(u) { 
      u * I(u <= quant) * dens_gx(u)
    }, lower = ming, upper = maxg, abs.tol = 1e-15)$value - eps * (ymax - ymin)
    
  })
  
  ub <- Vectorize(function(eps) {
    quant <- get_quant(1 - eps)
    t1 - t0 + integrate(f = function(u) { 
      u * I(u > quant) * dens_gx(u)
    }, lower = ming, upper = maxg, abs.tol = 1e-15)$value
  })
  
  eps_seq <- seq(0, 1, 0.2)
  lbvals <- lb(eps_seq)
  ubvals <- ub(eps_seq)
  
  # Start simulation to check coverage
  n <- 5000
  alpha <- 0.20
  nsplits <- 15
  B <- 10000
  
  sim_fn <- function() {
    # Simulation function to estimate eps0
    # Generate data according to the model above
    x <- runif(n, 0, 1); a <- rbinom(n, 1, expit(x))
    x <- as.data.frame(x)
    
    y <- t1 * a + t0 * (1 - a) + rbeta(n, aa, bb) - mean_beta
    
    outfam <- gaussian(); treatfam <- binomial(); model <- "x"
    eps <- eps_seq; delta <- 1; sl.lib <- c("SL.mean", "SL.glm")
    do_mult_boot <- TRUE; do_eps_zero <- FALSE; nuis_fns <- NULL 
    do_rearrange <- TRUE; show_progress <- FALSE; do_parallel <- TRUE
    ncluster <- 3; plugin <- FALSE
    
    tmp <- get_bound(y = y, a = a, x = x, ymin = ymin, ymax = ymax,
                     outfam = outfam, treatfam = treatfam, model = model,
                     eps = eps, delta = delta, nsplits = nsplits, 
                     sl.lib = sl.lib, alpha = alpha, do_mult_boot = do_mult_boot, 
                     do_eps_zero = do_eps_zero, nuis_fns = nuis_fns, 
                     do_rearrange = do_rearrange, B = B, 
                     show_progress = show_progress, do_parallel = do_parallel, 
                     ncluster = ncluster, plugin = plugin)
    
    cnames <- c("est_lb", "est_ub", "sigma_lb", "sigma_ub", 
                "ci_lb_lo", "ci_ub_hi", "calpha_lb", "calpha_ub")
    
    out <- cbind(tmp$bounds[, "lb", 1], tmp$bounds[, "ub", 1],
                 tmp$var_l[, 1], tmp$var_u[, 1],
                 tmp$bounds[, "ci_lb_lo_unif", 1], 
                 tmp$bounds[, "ci_ub_hi_unif", 1],
                 tmp$mult_calpha_lb, tmp$mult_calpha_ub)
    colnames(out) <- cnames
    
    return(out)
    
  }
  
  sims <- pbreplicate(nsim, sim_fn())
  
  if(length(eps_seq) > 1) {
    # if we are expecting uniform coverage over more than 1 eps val, the
    # z-score should be greater than z_{1 - alpha/2}. 
    check_multiplier(sims, alpha)
    check_multiplier(sims, alpha)
  }
  
  # tolerance for miscoverage
  tol <- 0.05
  
  lb_cvg <- sapply(1:length(eps_seq), function(x) { 
    mean(sims[x, "ci_lb_lo", ] <= lbvals[x]) })
  print(paste0("For LB, coverage is ", lb_cvg, " at eps = ", eps_seq, 
               " when ptwise is ", (1 - alpha / 2)))

  ub_cvg <- sapply(1:length(eps_seq), function(x) { 
    mean(sims[x, "ci_ub_hi", ] >= ubvals[x]) })
  print(paste0("For UB, coverage is ", ub_cvg, " at eps = ", eps_seq, 
               " when ptwise is ", (1 - alpha / 2)))
  
  check_cvg(nsim, eps_seq, lbvals, lbvals, sims[, "ci_lb_lo", ], Inf, alpha / 2, tol)
  check_cvg(nsim, eps_seq, ubvals, ubvals, -Inf, sims[, "ci_ub_hi", ], alpha / 2 , tol)
  check_cvg(nsim, eps_seq, lbvals, ubvals, sims[, "ci_lb_lo", ], sims[, "ci_ub_hi", ], 
            alpha, tol, FALSE)
  
})
