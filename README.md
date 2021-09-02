# Description

This repository provides R code for reproducing the results in Section 5 of the paper

Zheng, X., Kottas, A., and Sansó, B. (2021), “On Construction and
Estimation of Stationary Mixture Transition Distribution Models,”
[*arXiv:2010.12696*](https://arxiv.org/abs/2010.12696)

The R Code are in `code/` and data files are in `data/`. 
To run the code, install R package `mtd` from https://github.com/xzheng42/mtd.
It takes a long time to run the simulation studies and is recommended to
implement them on a cluster. 

-   Simulation studies: `sim1_mcmc`, `sim2_mcmc`, `sim3_mcmc`,
    `plot_sim1_weights`, `plot_sim2_weights`, `plot_sim3_results`.

-   Crime data example: `crime_mcmc`, `plot_crime_results`.

-   Precipitation data example: `precip_mcmc`, `plot_precip_results`.

