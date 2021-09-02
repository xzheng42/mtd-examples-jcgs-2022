################################################################################
### Description: this code runs the mcmc algorithms for the precipitation data 
###             example and computes the randomized residuals.
### Output: a data set called "res_precip_L10.RData" that contains posterior
###         samples for the Lomax MTD regression model, running times, computed 
###         randomized residuals, data, multiplicative factor, prior and mcmc parameters.
################################################################################
rm(list = ls())
path <- "path/to/save/results"
library(mtd)
library(xts)
library(lubridate)
library(hddtools)

################################################################################
### load data
################################################################################
allMOPEX <- catalogueMOPEX()
raw_dat <- tsMOPEX(allMOPEX$USGS_ID[31])
day_precip <- xts(raw_dat$P, order.by = date(raw_dat))
week_precip <- apply.weekly(day_precip, sum)
week_precip_8203 <- week_precip['1982/2003']
obs <- as.numeric(week_precip_8203)

################################################################################
### compute multiplicative factor
################################################################################
cyclength <- 52
freq <- 2 * pi / cyclength
ndata <- length(obs)
xx <- cbind(cos(freq * (1:ndata)), sin(freq * (1:ndata)),
            cos(2 * freq * (1:ndata)), sin(2 * freq * (1:ndata)),
            cos(3 * freq * (1:ndata)), sin(3 * freq * (1:ndata)))

################################################################################
### fit the model with order L = 10
################################################################################
# set prior
mtdorder <- 10
mcmc_param <- list(niter = 85000, nburn = 5000, nthin = 10)
prior <- list(u_alpha = 6, v_alpha = 1, u_phi = 3, v_phi = 20,
              alpha = 1, alpha_0 = 5, a_G0 = 1, b_G0 = 6.5)
starting <- list(phi = 5, bb = rep(0, ncol(xx)))
tuning <- list(se_phi = 0.7, se_bb = rep(1, length = ncol(xx)) / 4)

# run algorithms
set.seed(32)
runtime1 <- system.time(res1 <- tsRegMTD(formula = obs ~ xx,
                                         mtdorder = mtdorder, 
                                         family = "lomax",
                                         weight = "sb",
                                         prior = prior, 
                                         tuning = tuning, 
                                         starting = starting, 
                                         mcmc_param = mcmc_param))

runtime2 <- system.time(res2 <- tsRegMTD(formula = obs ~ xx,
                                         mtdorder = mtdorder, 
                                         family = "lomax",
                                         weight = "cdp",
                                         prior = prior, 
                                         tuning = tuning, 
                                         starting = starting, 
                                         mcmc_param = mcmc_param))

################################################################################
### compute randomized quantile residuals
################################################################################
runtime3 <- system.time(res1_rqr <- rqrMTD(res1, obs, "Lomax", xx))
runtime4 <- system.time(res2_rqr <- rqrMTD(res2, obs, "Lomax", xx))

################################################################################
### save results
################################################################################
runtime <- rbind(runtime1, runtime2, runtime3, runtime4)
save(res1, res2, res1_rqr, res2_rqr, runtime, week_precip_8203, xx, prior, mcmc_param, 
     starting, tuning, file = paste0(path, 'res_precip_L10.RData'))

