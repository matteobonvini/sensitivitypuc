context("Length of bound")
library(sensAteBounds)

test_that("length at eps=0 is 0", {
  for(j in 1:10) {
    n <- 500
    y <- rbinom(n,1,0.5)
    y2 <- runif(n)
    a <- rbinom(n,1,0.5)
    x <- as.data.frame(matrix(rnorm(2*n), ncol=2, nrow=n))
    res1 <- get_bound(y=y, a=a, x=x, out_model="logistic", treat_model="logistic",
                      eps=c(0,1), delta=0, nsplits=2, do_mult_boot=FALSE)$bounds
    res2 <- get_bound(y=y, a=a, x=x, out_model="ranger", treat_model="ranger",
                      eps=c(0,1), delta=0, nsplits=2, do_mult_boot=FALSE)$bounds
    res3 <- get_bound(y=y2, a=a, x=x, out_model="ranger", treat_model="ranger",
                      eps=c(0,1), delta=0, nsplits=2, do_mult_boot=FALSE)$bounds
    res4 <- get_bound(y=y, a=a, x=x, out_model="logistic", treat_model="logistic",
                      eps=seq(0, 1, length.out=100), delta=0, nsplits=2, 
                      do_mult_boot=FALSE)$bounds
    res5 <- get_bound(y=y2, a=a, x=x, out_model="ranger", treat_model="ranger",
                      eps=seq(0, 1, length.out=100), delta=0, nsplits=2, 
                      do_mult_boot=FALSE)$bounds
    expect_equal(res1[1, "lb"] - res1[1, "ub"], 0)
    expect_equal(res1[2, "ub"] - res1[2, "lb"], 1)
    expect_equal(res2[1, "lb"] - res2[1, "ub"], 0)
    expect_equal(res2[2, "ub"] - res2[2, "lb"], 1)
    expect_equal(res3[1, "lb"] - res3[1, "ub"], 0)
    expect_equal(res3[2, "ub"] - res3[2, "lb"], 1)
    expect_true(all(res4[, "lb"] <= res4[, "ub"]))
    expect_true(all(res5[, "lb"] <= res5[, "ub"]))
    expect_true(all(-1 <= res4[, "lb"] & res4[, "lb"] <= 1))
    expect_true(all(-1 <= res5[, "lb"] & res5[, "lb"] <= 1))
    expect_true(all(-1 <= res4[, "ub"] & res4[, "ub"] <= 1))
    expect_true(all(-1 <= res5[, "ub"] & res5[, "ub"] <= 1))
  }
})