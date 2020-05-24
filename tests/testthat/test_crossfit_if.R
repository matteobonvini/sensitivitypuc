context("Crossfitting and IFs")

require(pbapply)

test_that("Cross-fitting is correct", {
  
  nsplits <- c(1, 3, 5, 15, 20)
  nsim <- 50
  eps_seq <- seq(0, 1, 0.25)
  
  sim_fn <- function() {
    
    for(ee in eps_seq) {
      for(jj in nsplits) {
        
        n <- 500
        a <- rbinom(n, 1, 0.5)
        y <- rnorm(n)
        corx <- -0.5
        sigma <- matrix(c(2, corx, corx, 1), byrow = TRUE, ncol = 2, nrow = 2)
        x <- as.data.frame(MASS::mvrnorm(n, c(0, 1), sigma))
        ymin <- min(y); ymax <- max(y); sl.lib <- c("SL.glm")
        
        est <- get_bound(y = y, a = a, x = x, ymin = ymin, ymax = ymax, 
                         model = "x", eps = ee, delta = 1, sl.lib = sl.lib, 
                         outfam = gaussian(), treatfam = binomial(),
                         nsplits = jj, alpha = 0.10, do_mult_boot = FALSE,
                         do_eps_zero = TRUE, do_rearrange = FALSE,
                         do_parallel = FALSE, ncluster = 3)
        
        s <- est$nuis_fn$folds
        
        pi0 <- pi1 <- mu1 <- mu0 <- lambda_l <- lambda_u <- rep(NA, n)
        
        qu <- ql <- list(nsplits)
        
        gl <- gu <- nu <- tauu <- taul <- matrix(NA, ncol = 1, nrow = n)
        
        qu <- ql <- pi0t <- pi1t <- mu0t <- mu1t <- glt <- gut <- 
          vector("list", length = jj)

        for(k in 1:jj) {
          if(jj == 1) {
            train <- 1:n
            test <- 1:n
          } else {
            train <- which(s != k)
            test <- which(s == k)
          }
          
          dat <- cbind(data.frame(a = a[train]), x[train, ])
          fita <- glm(a ~ ., data = dat, family = binomial())
          pi1[test] <- predict.glm(fita, newdata = x[test, ], type = "response")
          pi1t[[k]] <- predict.glm(fita, newdata = x[train, ], type = "response")
          pi0 <- 1 - pi1
          pi0t[[k]] <- 1 - pi1t[[k]]
          
          train0 <- intersect(train, which(a == 0))
          daty0 <- cbind(data.frame(y = y[train0]), x[train0, ])
          fity0 <- glm(y ~ ., data = daty0, family = gaussian())
          mu0[test] <- predict.glm(fity0, newdata = x[test, ])
          mu0t[[k]] <- predict.glm(fity0, newdata = x[train, ])
          
          train1 <- intersect(train, which(a == 1))
          daty1 <- cbind(data.frame(y = y[train1]), x[train1, ])
          fity1 <- glm(y ~ ., data = daty1, family = gaussian())
          mu1[test] <- predict.glm(fity1, newdata = x[test, ])
          mu1t[[k]] <- predict.glm(fity1, newdata = x[train, ])
          
          gu[test, ] <- pi0[test] * (ymax - mu1[test]) - pi1[test] * (ymin - mu0[test])
          gl[test, ] <- pi0[test] * (ymin - mu1[test]) - pi1[test] * (ymax - mu0[test])
          
          gut[[k]] <- pi0t[[k]] * (ymax - mu1t[[k]]) - pi1t[[k]] * (ymin - mu0t[[k]])
          glt[[k]] <- pi0t[[k]] * (ymin - mu1t[[k]]) - pi1t[[k]] * (ymax - mu0t[[k]])
          
          qu[[k]] <- quantile(gu[test, ], p = 1 - ee)
          ql[[k]] <- quantile(gl[test, ], p = ee)
          
          nu[test, ] <- a[test] / pi1[test] * (y[test] - mu1[test]) - 
            (1 - a[test]) / pi0[test] * (y[test] - mu0[test]) +
            mu1[test] - mu0[test]
          
          mua <- a[test] * mu1[test] + (1 - a[test]) * mu0[test]
          pia <- a[test] * pi1[test] + (1 - a[test]) * pi0[test]
          
          tauu[test, ] <- (1 - 2 * a[test]) * (y[test] - mua) / pia * (1 - pia) - 
            a[test] * (ymin - mu0[test]) + (1 - a[test]) * (ymax - mu1[test])
          taul[test, ] <- tauu[test, ] - (ymax - ymin)
          
          lambda_l[test] <- 1 * I(gl[test] <= ql[[k]])
          lambda_u[test] <- 1 * I(gu[test] > qu[[k]])
          
        }
        
        if(ee == 0) {
          qu <- rep(list(max(gu)), jj)
          ql <- rep(list(min(gl) - 1e-10), jj)
          lambda_l <- lambda_u <- rep(0, n)
        }
        if(ee == 1) {
          qu <- rep(list(min(gu) - 1e-10), jj)
          ql <- rep(list(max(gl)), jj)
          lambda_l <- lambda_u <- rep(1, n)
        }
        
        mu1_pckg <- est$nuis_fns$test[, "mu1"]
        mu0_pckg <- est$nuis_fns$test[, "mu0"]
        pi1_pckg <- est$nuis_fns$test[, "pi1"]
        pi0_pckg <- est$nuis_fns$test[, "pi0"]
        
        expect_true(sum(abs(mu1[order(s)] - mu1_pckg)) < 1e-15)
        expect_true(sum(abs(mu0[order(s)] - mu0_pckg)) < 1e-15)
        expect_true(sum(abs(pi1[order(s)] - pi1_pckg)) < 1e-15)
        expect_true(sum(abs(pi0[order(s)] - pi0_pckg)) < 1e-15)
        
        mu1t_pckg <- est$nuis_fns$train[, "mu1"]
        mu0t_pckg <- est$nuis_fns$train[, "mu0"]
        pi1t_pckg <- est$nuis_fns$train[, "pi1"]
        pi0t_pckg <- est$nuis_fns$train[, "pi0"]
        
        expect_true(sum(abs(unlist(mu1t) - mu1t_pckg)) < 1e-15)
        expect_true(sum(abs(unlist(mu0t) - mu0t_pckg)) < 1e-15)
        expect_true(sum(abs(unlist(pi1t) - pi1t_pckg)) < 1e-15)
        expect_true(sum(abs(unlist(pi0t) - pi0t_pckg)) < 1e-15)
        
        datg <- data.frame(ymin = ymin, ymax = ymax, pi0g = pi0_pckg, 
                           pi1g = pi1_pckg, mu0 = mu0_pckg, mu1 = mu1_pckg)
        
        gl_pckg <- abind::abind(est$glhat, along = 1)
        gu_pckg <- abind::abind(est$guhat, along = 1)
        
        expect_true(sum(abs(gu[order(s), ] - gu_pckg)) < 1e-15)
        expect_true(sum(abs(gl[order(s), ] - gl_pckg)) < 1e-15)
        expect_true(sum(abs(gu[order(s), ] - gl[order(s)] - (ymax - ymin))) < 1e-10)
        
        datgt <- data.frame(ymin = ymin, ymax = ymax, pi0g = pi0t_pckg, 
                            pi1g = pi1t_pckg, mu0 = mu0t_pckg, mu1 = mu1t_pckg)
        
        glt_pckg <- abind::abind(est$glhat_train, along = 1)
        gut_pckg <- abind::abind(est$guhat_train, along = 1)
        
        expect_true(sum(abs(unlist(gut) - gut_pckg)) < 1e-15)
        expect_true(sum(abs(unlist(glt) - glt_pckg)) < 1e-15)
        expect_true(sum(abs(unlist(gut) - unlist(glt) - (ymax - ymin))) < 1e-10)
        
        nu_pckg <- unlist(est$nuhat)
        expect_true(sum(abs(nu[order(s), ] - nu_pckg)) < 1e-12)
        
        tauu_pckg <- abind::abind(est$tauhat_u, along = 1)
        taul_pckg <- abind::abind(est$tauhat_l, along = 1)
        expect_true(sum(abs(tauu[order(s), ] - tauu_pckg)) < 1e-12)
        expect_true(sum(abs(taul[order(s), ] - taul_pckg)) < 1e-12)
        
        ql_pckg <- abind::abind(est$q_l, along = 1)
        qu_pckg <- abind::abind(est$q_u, along = 1)
        expect_true(sum(abs(as.matrix(unlist(ql)) - ql_pckg)) < 1e-12)
        expect_true(sum(abs(as.matrix(unlist(qu)) - qu_pckg)) < 1e-12)
        
        lb_fold <- by(nu + lambda_l * taul, as.factor(s), 
                 function(x) apply(x, 2, mean))
        lb <- mean(unlist(lb_fold))
        
        ub_fold <- by(nu + lambda_u * tauu, as.factor(s), 
                      function(x) apply(x, 2, mean))
        ub <- mean(unlist(ub_fold))
        
        lb_pckg <- est$bounds[, "lb", 1]
        ub_pckg <- est$bounds[, "ub", 1]
        
        expect_true(abs(lb - lb_pckg) < 1e-15)
        expect_true(abs(ub - ub_pckg) < 1e-15)

        lambdaq_lb <- Map("*", split(lambda_l, as.factor(s)), ql)
        lambdaq_ub <- Map("*", split(lambda_u, as.factor(s)), qu)
        
        ifs_lb <- Map("-", split(nu + lambda_l * taul, as.factor(s)),
                      lambdaq_lb)
        ifs_ub <- Map("-", split(nu + lambda_u * tauu, as.factor(s)),
                      lambdaq_ub)
        
        var_lb <- var(unlist(ifs_lb))
        var_ub <- var(unlist(ifs_ub))
        
        var_lb_pckg <- est$var_lb
        var_ub_pckg <- est$var_ub
        
        expect_true(abs(var_lb - var_lb_pckg) < 1e-15)
        expect_true(abs(var_ub - var_ub_pckg) < 1e-15)
        
        ifs_lb_eps0 <- Map("/", ifs_lb, ql)
        ifs_ub_esp0 <- Map("/", ifs_ub, qu)
        
        var_eps0 <- ifelse(mean(nu) < 0, var(unlist(ifs_ub_esp0)) / n,
                           var(unlist(ifs_lb_eps0)) / n)
        
        var_eps0_pckg <- est$eps_zero[, "var_eps0"]
        
        expect_true(abs(var_eps0 - var_eps0_pckg) < 1e-15)
        
        # print(paste0("Done with testing (eps = ", ee, ", nsplits = ", jj, ")"))
      
      }
    }
  }
  
  sims <- pbreplicate(nsim, sim_fn())
  
})
