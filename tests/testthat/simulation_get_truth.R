###############################################
## Compute truth for simulation of Section 5 ##
###############################################
rm(list = ls())

set.seed(1000)

library(pbapply)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("simulation_true_regression_functions.R")


# For the simulation it's important that eps_seq is a subset of eps0_seq
# this saves computational time. See simulation.R file.
eps_seq <- seq(0, 0.2, 0.01)
eps0_seq <- seq(0, 0.2, 0.0001)

nsim <- 500
col_names <- c("eps", "lb", "ub", "eps_zero", "ql0", "qu0", "lb0", "ub0")

# Simulate data to estimate term E{g * I(g <= q)} as sample average
sim_fun <- function() {
  
  x1vals <- rtruncnorm(1e4, xlb, xub)
  x2vals <- rtruncnorm(1e4, xlb, xub)
  gub <- gx(x1vals, x2vals)
  glb <- gub - 1
  
  quants_lb <- quantile(glb, eps_seq)
  quants_ub <- quantile(gub, 1-eps_seq)
  quants_lb0 <- quantile(glb, eps0_seq)
  quants_ub0 <- quantile(gub, eps0_seq)
  
  g_term_lb <- get_g_term(g = glb, quants = quants_lb, upper = FALSE)
  g_term_ub <- get_g_term(g = gub, quants = quants_ub, upper = TRUE)
  g_term_lb0 <- get_g_term(g = glb, quants = quants_lb0, upper = FALSE)
  g_term_ub0 <- get_g_term(g = gub, quants = quants_ub0, upper = TRUE)
  
  true_lb <- mu1 - mu0 + g_term_lb
  true_ub <- mu1 - mu0 + g_term_ub
  true_lb0 <- mu1 - mu0 + g_term_lb0
  true_ub0 <- mu1 - mu0 + g_term_ub0
  
  idx_eps0 <- which.min(abs(true_lb0))
  true_eps0 <- eps0_seq[idx_eps0]
  true_ql_eps0 <- quants_lb0[idx_eps0]
  true_qu_eps0 <- quants_ub0[idx_eps0]
  true_lb_eps0 <- true_lb0[idx_eps0]
  true_ub_eps0 <- true_ub0[idx_eps0]
  
  out <- matrix(c(eps_seq, true_lb, true_ub, rep(true_eps0, length(eps_seq)),
                  rep(true_ql_eps0, length(eps_seq)), 
                  rep(true_qu_eps0, length(eps_seq)),
                  rep(true_lb_eps0, length(eps_seq)), 
                  rep(true_ub_eps0, length(eps_seq))),
                ncol = length(col_names), nrow = length(eps_seq), 
                dimnames = list(NULL, col_names))
  return(out)
}

res <- pbreplicate(nsim, sim_fun())
# The truth is the average of the results across simulations
out <- apply(res, c(1, 2), mean)
# Sanity check that the variance of the results across simulations is not large
sds <- apply(res, c(1, 2), sd) / sqrt(nsim)

attr(out, "tau") <- effect
attr(out, "pi") <- pi1
attr(out, "mu0") <- mu0
attr(out, "mu1") <- mu1
attr(out, "eps0_seq") <- seq(0, 0.2, 0.001)
attr(out, "eps_seq") <- eps_seq
saveRDS(out, file = "./data/truth_simulation.RData")
