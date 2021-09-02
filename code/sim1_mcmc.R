################################################################################
### Description: This code runs the first scenario of the simulation study 
###              presented in the paper. 
### Output: A data set called "res_sim1.RData" that contains posterior samples 
###         for the GMTD model in each case of the first scenario, running time 
###         for each case, and true_weight.
################################################################################
rm(list = ls())
path <- "path/to/save/results"
library(mtd)

################################################################################
### simulate data
################################################################################
set.seed(42)
J <- 5
expDecay <- function(x, la) exp(-la*x)
true_weight <- expDecay(1:J, 1)
true_weight <- true_weight / sum(true_weight)
param <- list(rho = c(0.7, 0.5, 0.3, 0.1, 0.1), mu = 10, sigma2 = 100)
init_data <- c(5, 12, 7, 8, 10)
ndata <- 2000
sim_data <- rMTD(J, true_weight, "gaussian", param, ndata, trun = TRUE, init_data = init_data)

################################################################################
### some general settings
################################################################################
mcmc_param <- list(niter = 165000, nburn = 5000, nthin = 20)

test_size <- length(sim_data$y)
test_dat <- sim_data$y[1:test_size]

weight_name <- c("dir", "sb", "cdp")
kk <- length(weight_name)
res_L5 <- res_L15 <- res_L25 <- vector("list", length = kk)
runtime_L5 <- runtime_L15 <- runtime_L25 <- array(NA, dim = c(kk, 3))

################################################################################
### fit the model with order 5
################################################################################
mtdorder <- 5
prior <- list(mu_0 = 0, sigma2_0 = 100,                  # hyperparam for mu 
              u_0 = 2, v_0 = 0.1,                        # hyperparam for sigma2 
              dir_shape = rep(1 / mtdorder, mtdorder),   # hyperparam for DIR
              alpha = 1, alpha_0 = 5, a_G0 = 1, b_G0 = 3)        # hyperparam for SB and CDP
starting <- list(rho = rep(0, mtdorder), sigma2 = var(test_dat))
tuning <- list(step_size = 0.01)

set.seed(32)
for (j in 1:kk) {
    runtime <- system.time(res <- tsMTD(obs = test_dat, 
                                        mtdorder = mtdorder, 
                                        family = "gaussian",
                                        weight = weight_name[j],
                                        prior = prior, 
                                        tuning = tuning, 
                                        starting = starting, 
                                        mcmc_param = mcmc_param))
    res_L5[[j]] <- res
    runtime_L5[j,] <- runtime[1:3]
}

################################################################################
### fit the model with order 15
################################################################################
mtdorder <- 15
prior$alpha <- 2
prior$b_G0 <- 6
prior$dir_shape <- rep(1 / mtdorder, mtdorder)
starting$rho <- rep(0, mtdorder)

set.seed(32)
for (j in 1:kk) {
    runtime <- system.time(res <- tsMTD(obs = test_dat, 
                                        mtdorder = mtdorder, 
                                        family = "gaussian",
                                        weight = weight_name[j],
                                        prior = prior, 
                                        tuning = tuning, 
                                        starting = starting, 
                                        mcmc_param = mcmc_param))
    res_L15[[j]] <- res
    runtime_L15[j,] <- runtime[1:3]
}

################################################################################
### fit the model with order 25
################################################################################
mtdorder <- 25
prior$alpha <- 3
prior$b_G0 <- 7
prior$dir_shape <- rep(1 / mtdorder, mtdorder)
starting$rho <- rep(0, mtdorder)

set.seed(32)
for (j in 1:kk) {
    runtime <- system.time(res <- tsMTD(obs = test_dat, 
                                        mtdorder = mtdorder, 
                                        family = "gaussian",
                                        weight = weight_name[j],
                                        prior = prior, 
                                        tuning = tuning, 
                                        starting = starting, 
                                        mcmc_param = mcmc_param))
    res_L25[[j]] <- res
    runtime_L25[j,] <- runtime[1:3]
}

################################################################################
### save results
################################################################################
save(res_L5, res_L15, res_L25, runtime_L5, runtime_L15, runtime_L25, true_weight, 
     file = paste0(path, 'res_sim1.RData'))




