context("Epsilon 0")

require(pbapply)
require(uniftest)

check_result <- function(sims, eps0_true, alpha) {
  # Check that coverage is correct and the the z-scores are correct
  cvg <- coverage(sims["eps0_lo", ], sims["eps0_hi", ], eps0_true, eps0_true)
  print(paste0("Coverage is ", cvg, "% when truth is ", 100*(1-alpha), "%"))
  eps0 <- sims["eps0", ]
  bias0 <- bias(sims["eps0", ], eps0_true)
  print(paste0("Bias is ", bias0))
  
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
  nsim <- 1000
  n <- 5000
  alpha <- 0.5

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
    nuis_fns <- matrix(NA, ncol = 7, nrow = 2 * n,
                       dimnames = list(NULL, c("mu1", "mu0", "pi1", "pi0",
                                               "unit_num", "fold", "is_test")))
    nuis_fns[, "mu1"] <- rep(mu1x, 2)
    nuis_fns[, "mu0"] <- rep(mu0x, 2)
    nuis_fns[, "pi1"] <- rep(pi1x, 2)
    nuis_fns[, "pi0"] <- rep(pi0x, 2)

    nuis_fns[, "unit_num"] <- rep(1:n, 2)
    nuis_fns[, "fold"] <- 1

    nuis_fns_list <- list(test = nuis_fns[1:n, ],
                          train = nuis_fns[(n+1):(2*n), ],
                          order_obs = nuis_fns[1:n, "unit_num"])

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
                     eps = eps, delta = 1, nsplits = 1, do_mult_boot = FALSE,
                     do_eps_zero = TRUE, nuis_fns = nuis_fns_list, alpha = alpha)
    tmp <- tmp$eps_zero
    out <- c(tmp$est, tmp$ci_lo, tmp$ci_hi, tmp$var_eps0, var_eps0)
    names(out) <-  outnames

    return(out)
  }
  # Simulation begins
  sims <- pbreplicate(nsim, sim_fn())
  check_result(sims, eps0_true, alpha)

})

test_that("coverage eps0 is correct when using sample splitting", {

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
  
  eps_seq <- seq(0.3, 0.4, 0.00001)
  quants <- get_quant(eps_seq)
  lbvals <- lb(eps_seq)
  ubvals <- ub(eps_seq)
  eps0_true <- eps_seq[which.min(abs(lbvals * ubvals))]

  # Start simulation to check coverage
  nsim <- 500
  n <- 5000
  alpha <- 0.2
  nsplits <- 3

  sim_fn <- function() {
    # Simulation function to estimate eps0
    # Generate data according to the model above
    x <- runif(n, 0, 1)
    a <- rbinom(n, 1, expit(x))
    
    y <- t1 * a + t0 * (1 - a) + rbeta(n, aa, bb) - mean_beta

    # Regression functions are known exacty
    mu1x <- rep(t1, n)
    mu0x <- rep(t0, n)
    pi1x <- expit(x)
    pi0x <- 1 - expit(x)

    # This code computes the variance of eps_zero using truth
    gx <- gx_fn(x)
    piax <- a * pi1x + (1 - a) * pi0x
    muax <- a * mu1x + (1 - a) * mu0x
    nu <- (2 * a - 1) * (y - muax) / piax + mu1x - mu0x
    tau <- (1 - 2 * a) * (y - muax) / piax * (1 - piax) + a * (mu0x - ymin) + 
      (1 - a) * (ymax - mu1x)
    qeps0u <- get_quant(1 - eps0_true)
    psil0 <- lb(eps0_true)
    tildephi0 <- nu + (tau - qeps0u) * I(gx > qeps0u)
    var_eps0 <- var( (tildephi0 * psil0 ) / (psil0 * qeps0u)) / n
    
    tmp <- get_bound(y = y, a = a, x = as.data.frame(x), ymin = ymin, ymax = ymax,
                     outfam = gaussian(), treatfam = binomial(), model = "x",
                     eps = seq(0.3, 0.4, 0.001), delta = 1, nsplits = nsplits, 
                     sl.lib = c("SL.mean", "SL.glm"), alpha = alpha, 
                     do_mult_boot = FALSE, do_eps_zero = TRUE, nuis_fns = NULL, 
                     do_rearrange = FALSE)

    tmp0 <- tmp$eps_zero
    out <- c(tmp0$est, tmp0$ci_lo, tmp0$ci_hi, tmp0$var_eps0, var_eps0)
    names(out) <-  outnames

    return(out)
    
  }
  # Simulation begins
  sims <- pbreplicate(nsim, sim_fn())
  check_result(sims, eps0_true, alpha)

})
 
test_that("coverage eps0 is correct using simulation from paper", {

  source("simulation_true_regression_functions.R")
  ## Load true values ##
  truth <- readRDS("./data/truth_simulation.RData")
  
  psiu0 <- truth[1, "ub0"]
  psil0 <- truth[1, "lb0"]
  ql0 <- truth[1, "ql0"]
  qu0 <- truth[1, "qu0"]
  eps0_true <- truth[1, "eps_zero"]
  der <- psiu0 * ql0
  
  n <- 5000
  alpha <- 0.3
  nsim <- 2000
  nsplits <- 3
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
    
    nuis_fns_list <- list(test = test_set, train = train_set,
                          order_obs = test_set[, "unit_num"])
    
    mu1hat <- nuis_fns_list$test[, "mu1"]
    mu0hat <- nuis_fns_list$test[, "mu0"]
    pi1hat <- nuis_fns_list$test[, "pi1"]
    pi0hat <- nuis_fns_list$test[, "pi0"]
    order_obs <- nuis_fns_list$order_obs
    
    y <- df$y
    a <- df$a
    x <- df[, c("x1", "x2")]
    
    yord <- y[order_obs]
    aord <- a[order_obs]
    
    piax <- aord * pi1hat + (1 - aord) * pi0hat
    muax <- aord * mu1hat + (1 - aord) * mu0hat
    nu <- (2 * aord - 1) * (yord - muax) / piax + mu1hat - mu0hat
    tau <- (1 - 2 * aord) * (yord - muax) / piax * (1 - piax) + aord * mu0hat +
      (1 - aord) * (1 - mu1hat)
    gux <- pi0hat * (1 - mu1hat) + pi1hat * mu0hat
    glx <- gux - 1

    var_eps0 <- var(psiu0 * (nu + I(glx <= ql0) * (tau - ql0) - eps0_true) / der) / n
    
    eps <- attr(truth, "eps0_seq")
    tmp <- get_bound(y = y, a = a, x = x, ymin = 0, ymax = 1, outfam = NULL,
                     treatfam = NULL, model = "x", eps = eps, delta = 1,
                     nsplits = nsplits, do_mult_boot = FALSE, B = NULL,
                     do_eps_zero = TRUE, nuis_fns = nuis_fns_list, alpha = alpha,
                     do_rearrange = FALSE)
    tmp0 <- tmp$eps_zero
    out <- c(tmp0$est, tmp0$ci_lo, tmp0$ci_hi, tmp0$var_eps0, var_eps0)
    names(out) <-  outnames

    return(out)
  }

  sims <- pbreplicate(nsim, sim_fn())
  check_result(sims, eps0_true, alpha)
})

