context("Length of bound (plugin estimator used)")

test_that("length at delta/eps in {0, 1} is {0, 1}", {
  nsim  <- 500
  n <- 500 
  
  width_test <- function(ymin, ymax, model) {
    
    a <- rbinom(n, 1, 0.5)
    x <- as.data.frame(matrix(rnorm(2*n), ncol = 2, nrow = n))
    y <- runif(n, ymin, ymax)
    max_width <- ymax - ymin
    
    res <- get_bound(y = y, a = a, x = x, ymin = ymin, ymax = ymax, 
                     outfam = gaussian(),  treatfam = binomial(), 
                     model = model, eps = c(0, 1), delta = c(0, 1), 
                     do_mult_boot = FALSE, do_eps_zero = FALSE, nsplits = 1, 
                     alpha = 0.05, B = NULL, sl.lib = "SL.glm")$bounds
    
    
    max_widthc <- as.character(max_width)
    expect_equal(res["0", "lb" , "0"] - res["0", "ub" , "0"], 0)
    expect_equal(res["1", "lb" , "0"] - res["1", "ub" , "0"], 0)
    expect_equal(res["0", "lb" , "1"] - res["0", "ub" , "1"], 0)
    expect_equal(res["1", "ub" , "1"] - res["1", "lb" , "1"], max_width)
    
  }
  
  width_test_wrapper <- function() {
    
    width_test(0, 1, "x")
    width_test(0, 1, "xa")
    ymin <- runif(1, -100, 100)
    ymax <- ymin + runif(1, 0, 200)
    width_test(ymin, ymax, "x")
    width_test(ymin, ymax, "xa")
    
  }
  
  pbreplicate(nsim, width_test_wrapper())
})

test_that("width of bounds always in [0, 1]", {
  
  nsim <- 500
  n <- 500
  
  range_test <- function(ymin, ymax, model) {
    a <- rbinom(n, 1, 0.5)
    x <- as.data.frame(matrix(rnorm(2*n), ncol = 2, nrow = n))
    y <- runif(n, ymin, ymax)
    max_width <- ymax - ymin
    delta_seq <- seq(0, 1, length.out = 10)
    
    res <- get_bound(y = y, a = a, x = x, ymin = ymin, ymax = ymax, 
                     outfam = gaussian(),  treatfam = binomial(), 
                     model = model, eps = c(0, 1), delta = delta_seq, 
                     do_mult_boot = FALSE, do_eps_zero = FALSE, nsplits = 1, 
                     alpha = 0.05, B = NULL, plugin = TRUE,
                     sl.lib = c("SL.mean", "SL.glm"))$bounds
    
    length_bounds <- res[, "ub", ] - res[, "lb", ]
    tol <- 1e-10
    length_bounds <- round(length_bounds, tol)
    max_width <- round(max_width, tol)
    expect_true(all(0 <= length_bounds & length_bounds <= max_width))
    
  }
  
  range_test_wrapper <- function() {
    
    range_test(0, 1, "x")
    range_test(0, 1, "xa")
    ymin <- runif(1, -100, 100)
    ymax <- ymin + runif(1, 0, 200)
    range_test(ymin, ymax, "x")
    range_test(ymin, ymax, "xa")
  }
  
  pbreplicate(nsim, range_test_wrapper())
})

test_that("width is increasing in epsilon", {
  
  nsim  <- 500
  n <- 500
  
  mono_test <- function(ymin, ymax, model) {
    y <- runif(n, ymin, ymax)
    a <- rbinom(n, 1, 0.5)
    x <- as.data.frame(matrix(rnorm(2*n), ncol = 2, nrow = n))
    
    max_width <- ymax - ymin
    delta_seq <- seq(0, 1, length.out = 10)
    eps_seq <- seq(0, 1, length.out = 50)
    res <- get_bound(y = y, a = a, x = x, ymin = ymin, ymax = ymax, 
                     outfam = gaussian(),  treatfam = binomial(), 
                     model = model, eps = eps_seq, delta = delta_seq, 
                     do_mult_boot = FALSE, do_eps_zero = FALSE, nsplits = 1, 
                     alpha = 0.05, B = NULL, plugin = TRUE,
                     sl.lib = c("SL.mean", "SL.glm"))$bounds
    
    length_bounds <- res[, "ub", ] - res[, "lb", ]
    tol <- 1e-10
    diff_lengths <- round(apply(length_bounds, 2, diff, lag = 1), tol)
    expect_true(all(diff_lengths >= 0))
    
  }
  
  mono_test_wrapper <- function() {
    
    mono_test(0, 1, "x")
    mono_test(0, 1, "xa")
    
    ymin <- runif(1, -100, 100)
    ymax <- ymin + runif(1, 0, 200)
    
    mono_test(ymin, ymax, "x")
    mono_test(ymin, ymax, "xa")
    
  }
  
  pbreplicate(nsim, mono_test_wrapper())
  
})
