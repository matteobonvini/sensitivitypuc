context("Length of bound")

test_that("length at delta/eps in {0, 1} is {0, 1}", {
  
  for(j in 1:10) {
    n <- 50000
    y <- rbinom(n, 1, 0.5)
    a <- rbinom(n, 1, 0.5)
    x <- as.data.frame(matrix(rnorm(2*n), ncol = 2, nrow = n))
    
    res1 <- get_bound(y=y, a=a, x=x, outfam=binomial(), treatfam=binomial(), 
                      model="x", eps=c(0, 1), delta=c(0, 1), nsplits=1, 
                      do_mult_boot=FALSE, do_eps_zero=FALSE, alpha=0.05, 
                      B=NULL, nuis_fns=NULL,  sl.lib=c("SL.glm"))$bounds
    
    expect_equal(res1[res1$eps == 0 & res1$delta == 0, "lb"] - 
                   res1[res1$eps == 0 & res1$delta == 0, "ub"], 0)
    expect_equal(res1[res1$eps == 1 & res1$delta == 0, "lb"] - 
                   res1[res1$eps == 1 & res1$delta == 0, "ub"], 0)
    expect_equal(res1[res1$eps == 0 & res1$delta == 1, "lb"] - 
                   res1[res1$eps == 0 & res1$delta == 1, "ub"], 0)
    expect_true(abs(res1[res1$eps == 1 & res1$delta == 1, "ub"] - 
                      res1[res1$eps == 1 & res1$delta == 1, "lb"] - 1) <= 1e-5)
    
    # Check that the same holds if Y is continuous in [0, 1]
    y2 <- runif(n)
    res2 <- get_bound(y=y2, a=a, x=x, outfam=gaussian(), treatfam=binomial(), 
                      model="x", eps=c(0, 1), delta=c(0, 1), nsplits=1, 
                      do_mult_boot=FALSE, do_eps_zero=FALSE, alpha=0.05, 
                      B=NULL, nuis_fns=NULL,  sl.lib=c("SL.glm"))$bounds
    
    expect_equal(res2[res2$eps == 0 & res2$delta == 0, "lb"] - 
                   res2[res2$eps == 0 & res2$delta == 0, "ub"], 0)
    expect_equal(res2[res2$eps == 1 & res2$delta == 0, "lb"] - 
                   res2[res2$eps == 1 & res2$delta == 0, "ub"], 0)
    expect_equal(res2[res2$eps == 0 & res2$delta == 1, "lb"] - 
                   res2[res2$eps == 0 & res2$delta == 1, "ub"], 0)
    expect_true(abs(res2[res2$eps == 1 & res2$delta == 1, "ub"] - 
                     res2[res2$eps == 1 & res2$delta == 1, "lb"] - 1) <= 1e-4)
  }
})

test_that("width of bounds always in [0, 1]", {
  
  for(j in 1:10) {
    n <- 5000
    y <- rbinom(n, 1, 0.5)
    a <- rbinom(n, 1, 0.5)
    x <- as.data.frame(matrix(rnorm(2*n), ncol = 2, nrow = n))
    
    res1 <- get_bound(y=y, a=a, x=x, outfam=binomial(), treatfam=binomial(), 
                      model="x", eps=seq(0, 1, 0.01), delta=c(0, 1), nsplits=1, 
                      do_mult_boot=FALSE, do_eps_zero=FALSE, alpha=0.05, 
                      B=NULL, nuis_fns=NULL,  sl.lib=c("SL.glm"))$bounds
    length_bounds <- res1[, "ub"] - res1[, "lb"]
    expect_true(all(0 <= length_bounds & length_bounds <= 1))
    
    # Check that the same holds if Y is continuous in [0, 1]
    y2 <- runif(n)
    res2 <- get_bound(y=y2, a=a, x=x, outfam=gaussian(), treatfam=binomial(), 
                      model="x", eps=seq(0, 1, 0.01), delta=c(0, 1), nsplits=1, 
                      do_mult_boot=FALSE, do_eps_zero=FALSE, alpha=0.05, 
                      B=NULL, nuis_fns=NULL,  sl.lib=c("SL.glm"))$bounds
    length_bounds2 <- res2[, "ub"] - res2[, "lb"]
    expect_true(all(0 <= length_bounds2 & length_bounds2 <= 1))
  }
})

test_that("width is increasing in epsilon (using plug-in estimator)", {
  
  for(j in 1:10) {
    n <- 500
    y <- rbinom(n, 1, 0.5)
    a <- rbinom(n, 1, 0.5)
    x <- as.data.frame(matrix(rnorm(2*n), ncol = 2, nrow = n))
    
    res1 <- get_bound(y=y, a=a, x=x, outfam=binomial(), treatfam=binomial(), 
                      model="x", eps=seq(0, 1, 0.001), delta=1, nsplits=1, 
                      do_mult_boot=FALSE, do_eps_zero=FALSE, alpha=0.05, 
                      B=NULL, nuis_fns=NULL,  sl.lib=c("SL.glm"), 
                      plugin=TRUE)$bounds
    length_bounds <- res1[, "ub"] - res1[, "lb"]
    expect_true(all(diff(length_bounds, 1) >= 0))
    
    # Check that the same holds if Y is continuous in [0, 1]
    y2 <- runif(n)
    res2 <- get_bound(y=y2, a=a, x=x, outfam=gaussian(), treatfam=binomial(), 
                      model="x", eps=seq(0, 1, 0.001), delta=1, nsplits=1, 
                      do_mult_boot=FALSE, do_eps_zero=FALSE, alpha=0.05, 
                      B=NULL, nuis_fns=NULL,  sl.lib=c("SL.glm"),
                      plugin=TRUE)$bounds
    length_bounds2 <- res2[, "ub"] - res2[, "lb"]
    expect_true(all(diff(length_bounds2, 1) >= 0))
  }
})
