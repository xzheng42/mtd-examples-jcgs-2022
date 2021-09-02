################################################################################
### Description: this code runs the second simulation study presented in the 
###              Supplementary Material. There are two models. They are
###              negative binomial MTD (NBMTD) and Poisson MTD (PMTD).
### Output: a data set called "res_sim3.RData" that contains posterior samples 
###         of the NBMTD and PMTD models with different orders, running time 
###         for each case, simulated data, prior, mcmc parameters, and
###         true parameters of the simulation.
################################################################################
rm(list = ls())
path <- "path/to/save/results"
library(mtd)

################################################################################
### simulate data from a negative binomial MTD
################################################################################
set.seed(42)
J <- 5
expDecay <- function(x, la) exp(-la*x)
true_weight <- expDecay(1:J, 1)
true_weight <- true_weight / sum(true_weight)
param <- list(la = 5, ga = 3, kap = 3, eta = 2)
init_data <- c(5, 12, 7, 8, 10)
ndata <- 2000
sim_data <- rMTD(J, true_weight, "negative binomial", param, ndata, trun = TRUE, init_data = init_data)

################################################################################
### general setting for MCMC
################################################################################
nn <- length(sim_data$y)
mcmc_param <- list(niter = 85000, nburn = 5000, nthin = 10)
prior <- list(u_th = 2, v_th = 2, u_psi = 6, v_psi = 2, u_kap = 2, v_kap = 1, 
              alpha = 1, alpha_0 = 5, a_G0 = 1, b_G0 = 3)
starting <- list(kap = 2)
tuning <- list(se_kap = 0.15)

res_sb_list <- vector("list", length = length(2))
res_cdp_list <- vector("list", length = length(2))
for (j in 1:2) res_sb_list[[j]] <- vector("list", length = 2)
for (j in 1:2) res_cdp_list[[j]] <- vector("list", length = 2)
runtime_list <- vector("list", length = 2)
runtime_L5 <- runtime_L15 <- array(NA, dim = c(4, 3))

################################################################################
### fit the model with order L = 5
################################################################################
mtdorder <- 5
set.seed(32)

# fit the NBMTD with two priors for the weights --------------------------------
runtime1 <- system.time(res1 <- tsMTD(obs = sim_data$y, 
                                     mtdorder = mtdorder, 
                                     family = "negative binomial", 
                                     weight = "sb", 
                                     prior = prior, 
                                     tuning = tuning, 
                                     starting = starting, 
                                     mcmc_param = mcmc_param))
res_sb_list[[1]][[1]] <- res1
runtime_L5[1, ] <- runtime1[1:3] 

runtime2 <- system.time(res2 <- tsMTD(obs = sim_data$y, 
                                      mtdorder = mtdorder, 
                                      family = "negative binomial", 
                                      weight = "cdp", 
                                      prior = prior, 
                                      tuning = tuning, 
                                      starting = starting, 
                                      mcmc_param = mcmc_param))
res_cdp_list[[1]][[1]] <- res2
runtime_L5[2, ] <- runtime2[1:3] 

# fit the PMTD with two priors for the weights ---------------------------------
runtime3 <- system.time(res3 <- tsMTD(obs = sim_data$y, 
                                      mtdorder = mtdorder, 
                                      family = "poisson", 
                                      weight = "sb", 
                                      prior = prior, 
                                      tuning = tuning, 
                                      starting = starting, 
                                      mcmc_param = mcmc_param))
res_sb_list[[1]][[2]] <- res3
runtime_L5[3, ] <- runtime3[1:3] 
runtime4 <- system.time(res4 <- tsMTD(obs = sim_data$y, 
                                      mtdorder = mtdorder, 
                                      family = "poisson", 
                                      weight = "cdp", 
                                      prior = prior, 
                                      tuning = tuning, 
                                      starting = starting, 
                                      mcmc_param = mcmc_param))
res_cdp_list[[1]][[2]] <- res4
runtime_L5[4, ] <- runtime4[1:3]

################################################################################
### fit the model with order L = 15
################################################################################
mtdorder <- 15
prior$alpha <- 2
prior$b_G0 <- 6
set.seed(32)

# fit the NBMTD with two priors for the weights --------------------------------
runtime1 <- system.time(res1 <- tsMTD(obs = sim_data$y, 
                                      mtdorder = mtdorder, 
                                      family = "negative binomial", 
                                      weight = "sb", 
                                      prior = prior, 
                                      tuning = tuning, 
                                      starting = starting, 
                                      mcmc_param = mcmc_param))
res_sb_list[[2]][[1]] <- res1
runtime_L15[1, ] <- runtime1[1:3] 

runtime2 <- system.time(res2 <- tsMTD(obs = sim_data$y, 
                                      mtdorder = mtdorder, 
                                      family = "negative binomial", 
                                      weight = "cdp", 
                                      prior = prior, 
                                      tuning = tuning, 
                                      starting = starting, 
                                      mcmc_param = mcmc_param))
res_cdp_list[[2]][[1]] <- res2
runtime_L15[2, ] <- runtime2[1:3] 

# fit the PMTD with two priors for the weights ---------------------------------
runtime3 <- system.time(res3 <- tsMTD(obs = sim_data$y, 
                                      mtdorder = mtdorder, 
                                      family = "poisson", 
                                      weight = "sb", 
                                      prior = prior, 
                                      tuning = tuning, 
                                      starting = starting, 
                                      mcmc_param = mcmc_param))
res_sb_list[[2]][[2]] <- res3
runtime_L15[3, ] <- runtime3[1:3] 

runtime4 <- system.time(res4 <- tsMTD(obs = sim_data$y, 
                                      mtdorder = mtdorder, 
                                      family = "poisson", 
                                      weight = "cdp",
                                      prior = prior, 
                                      tuning = tuning, 
                                      starting = starting, 
                                      mcmc_param = mcmc_param))
res_cdp_list[[2]][[2]] <- res4
runtime_L15[4, ] <- runtime4[1:3] 

################################################################################
### save results
################################################################################
runtime_list[[1]] <- runtime_L5
runtime_list[[2]] <- runtime_L15
true_param <- list(la = 5, ga = 3, kap = 3, eta = 2, weight = true_weight)
save(res_sb_list, res_cdp_list, runtime_list, sim_data, prior, mcmc_param, 
     starting, tuning, true_param, file = paste0(path, "res_sim3.RData"))



