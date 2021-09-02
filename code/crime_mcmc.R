################################################################################
### Description: this code runs the mcmc algorithms for the crime data example
###              and computes the randomized residuals.
### Output: a data set called "res_crime_L20.RData" that contains data, 
###         posterior samples for the PMTD model, running times, 
###         computed randomized residuals, priors and mcmc parameters. 
################################################################################
rm(list = ls())
path1 <- "path/to/import/data"
path2 <- "path/to/save/results"
library(mtd)
library(readr)

################################################################################
### load data
################################################################################
dat_dt_count <- read_csv("path/to/data/crime.csv", skip = 3)
dat_dt_count$Date <- as.Date(dat_dt_count$Date, "%m/%d/%y")
dat_dt_count$n <- as.integer(dat_dt_count$n)
obs <- unlist(dat_dt_count[,2])

################################################################################
### fit the model with order L = 20
################################################################################
# set up the model
mcmc_param <- list(niter = 85000, nburn = 5000, nthin = 10)
prior <- list(u_la = 2, v_la = 1, u_th = 2, v_th = 2,
              alpha = 2, alpha_0 = 5, a_G0 = 1, b_G0 = 8)
mtdorder <- 20

# fit the model
set.seed(32)
runtime1 <- system.time(res1 <- tsMTD(obs = obs, 
                                      mtdorder = mtdorder, 
                                      family = "poisson",
                                      weight = "sb",
                                      prior = prior, 
                                      tuning = NULL, 
                                      starting = NULL,
                                      mcmc_param = mcmc_param))

runtime2 <- system.time(res2 <- tsMTD(obs = obs, 
                                      mtdorder = mtdorder, 
                                      family = "poisson",
                                      weight = "cdp",
                                      prior = prior, 
                                      tuning = NULL,
                                      starting = NULL,
                                      mcmc_param = mcmc_param))

################################################################################
### compute randomized quantile residuals
################################################################################
runtime3 <- system.time(res1_rqr <- rqrMTD(res1, obs, "poisson"))
runtime4 <- system.time(res2_rqr <- rqrMTD(res2, obs, "poisson"))

################################################################################
### save results
################################################################################
runtime <- rbind(runtime1, runtime2, runtime3, runtime4)
save(res1, res2, res1_rqr, res2_rqr, runtime, mcmc_param, prior, dat_dt_count, 
     file = paste0(path2, 'res_crime_L20.RData'))




