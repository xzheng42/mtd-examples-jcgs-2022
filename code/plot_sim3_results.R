################################################################################
### Description: make plots of the results of the second simulation study 
###              presented in the Supplementary Material.
### Output: Figures 3, 4, 5 and Table 1 of the Supplementary Material.
################################################################################
rm(list = ls())
library(mtd)
library(ggplot2)
library(pammtools)

################################################################################
### load data
################################################################################
path1 <- "path/to/import/data"
path2 <- "path/to/save/results"
load(paste0(path1, 'res_sim3.RData'))

################################################################################
### Figure 3 of Supplementary Material
################################################################################
# plot function ----------------------------------------------------------------
plotWeight <- function(res, mtdorder, prior_mu, true_weight, 
                       xlab_name = '', ylab_name = '', weight_name = NULL) {
  
  if (mtdorder <= 5) {
    plot_break <- 0:length(true_weight)
    plot_label <- 0:length(true_weight)
  } else {
    plot_break <- seq(0, mtdorder, by = 5)
    plot_label <- seq(0, mtdorder, by = 5)
  }
  
  kk <- length(true_weight)
  post_mu <- rowMeans(res$weight)
  post_qq <- apply(res$weight, 1, function(x) quantile(x, probs = c(0.025, 0.975)))
  weight <- cbind(c(true_weight, rep(0, mtdorder - kk)), prior_mu, post_mu, t(post_qq))
  weight <- rbind(weight, tail(weight, 1))
  weight <- as.data.frame(cbind(0:mtdorder, weight))
  colnames(weight) <- c("idx", "true", "prior_mu", "post_mu", "q025", "q975")
  plot_weight_save <- ggplot(weight, aes(x = idx)) + 
    geom_stepribbon(aes(ymin = q025, ymax = q975, fill = "Posterior 95% CI"), alpha = 0.15) + 
    geom_step(aes(y = post_mu, color = "Posterior mean"), size = 1.5) + 
    geom_step(aes(y = prior_mu, color = "Prior mean"), linetype = "dotdash", size = 1.5) + 
    geom_step(aes(y = true, color = "True"), linetype = "dashed", size = 1.5) + 
    scale_fill_manual("", values = "blue") + 
    scale_colour_manual("", values = c("blue", "brown1","black")) + 
    scale_x_continuous(breaks = plot_break, labels = plot_label) + 
    theme_bw() + ylim(0, 0.8) + labs(fill = "", x = xlab_name, y = ylab_name) + 
    theme(legend.position = "none",
          plot.margin = unit(c(0.1,0.1,1,0.1), "cm"),
          plot.caption = element_text(hjust = 0.5, vjust = -3, size = 26),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          axis.text.x = element_text(margin = unit(c(0.3, 0, 0.3, 0), "cm")),
          axis.text.y = element_text(margin = unit(c(0, 0.3, 0, 0.3), "cm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  if (!is.null(weight_name)) {
    plot_weight_save <- plot_weight_save + labs(caption = weight_name)
  }
  
  plot_weight_save
  
}

# order L = 5 ------------------------------------------------------------------
mtdorder <- 5
alpha <- 1
a_G0 <- 1
b_G0 <- 3
cutoff <- seq(0, 1, length = mtdorder + 1)
CDP_prior_mu <- diff(pbeta(cutoff, a_G0, b_G0))
SB_mu <- 1 / (1 + alpha) * (alpha / (1 + alpha))^((1:(mtdorder - 1)) - 1) 
SB_prior_mu <- c(SB_mu, 1 - sum(SB_mu))

png(paste0(path2, 1, "_", mtdorder, "_sb.png"), units = "in", width = 6, height = 5, res = 300)
print(plotWeight(res_sb_list[[1]][[1]], mtdorder, SB_prior_mu, true_param$weight, "", "Probability"))
dev.off()

png(paste0(path2, 1, "_", mtdorder, "_cdp.png"), units = "in", width = 6, height = 5, res = 300)
print(plotWeight(res_cdp_list[[1]][[1]], mtdorder, SB_prior_mu, true_param$weight))
dev.off()

png(paste0(path2, 2, "_", mtdorder, "_sb.png"), units = "in", width = 6, height = 5, res = 300)
print(plotWeight(res_sb_list[[1]][[2]], mtdorder, SB_prior_mu, true_param$weight))
dev.off()

png(paste0(path2, 2, "_", mtdorder, "_cdp.png"), units = "in", width = 6, height = 5, res = 300)
print(plotWeight(res_cdp_list[[1]][[2]], mtdorder, SB_prior_mu, true_param$weight))
dev.off()

# order L = 15 -----------------------------------------------------------------
mtdorder <- 15
alpha <- 2
a_G0 <- 1
b_G0 <- 6
cutoff <- seq(0, 1, length = mtdorder + 1)
CDP_prior_mu <- diff(pbeta(cutoff, a_G0, b_G0))
SB_mu <- 1 / (1 + alpha) * (alpha / (1 + alpha))^((1:(mtdorder - 1)) - 1) 
SB_prior_mu <- c(SB_mu, 1 - sum(SB_mu))

png(paste0(path2, 1, "_", mtdorder, "_sb.png"), units = "in", width = 6, height = 5.3, res = 300)
print(plotWeight(res_sb_list[[2]][[1]], mtdorder, SB_prior_mu, true_param$weight, "Lag", "Probability", "(a) NBMTD-SB"))
dev.off()

png(paste0(path2, 1, "_", mtdorder, "_cdp.png"), units = "in", width = 6, height = 5.3, res = 300)
print(plotWeight(res_cdp_list[[2]][[1]], mtdorder, SB_prior_mu, true_param$weight, "Lag", "", "(b) NBMTD-CDP"))
dev.off()

png(paste0(path2, 2, "_", mtdorder, "_sb.png"), units = "in", width = 6, height = 5.3, res = 300)
print(plotWeight(res_sb_list[[2]][[2]], mtdorder, SB_prior_mu, true_param$weight, "Lag", "", "(c) PMTD-SB"))
dev.off()

png(paste0(path2, 2, "_", mtdorder, "_cdp.png"), units = "in", width = 6, height = 5.3, res = 300)
print(plotWeight(res_cdp_list[[2]][[2]], mtdorder, SB_prior_mu, true_param$weight, "Lag", "", "(d) PMTD-CDP"))
dev.off()

################################################################################
### Figure 4 of Supplementary Material
################################################################################
# plot function for NBMTD ------------------------------------------------------
plot_nbmtd_est_marg <- function(res_sb, res_cdp, gridval, true_pmf, obs, ylab_name, caption_name) {
  
  yy <- sapply(gridval, function(x) dnbinom(x, res_sb$kap, 1-(1/res_sb$psi-1)/(1-res_sb$th)))
  yy_mu <- colMeans(yy)
  yy_qq <- apply(yy, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
  zz <- sapply(gridval, function(x) dnbinom(x, res_cdp$kap, 1-(1/res_cdp$psi-1)/(1-res_cdp$th)))
  zz_mu <- colMeans(zz)
  zz_qq <- apply(zz, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
  dens <- data.frame(grid = gridval, true = true_pmf,
                     sb = yy_mu, sb025 = yy_qq[1,], sb975 = yy_qq[2,],
                     cdp = zz_mu, cdp025 = zz_qq[1,], cdp975 = zz_qq[2,])
  dat = data.frame(obs = obs)
  
  plot_hist_save <- ggplot(data = dens) + theme_bw() + ylim(0, 0.12) + 
    geom_histogram(data = dat, aes(x = obs, y=..density..), binwidth = 4, fill="white", color="darkgrey", size = 1) + 
    geom_point(aes(x = grid, y = true), shape = 1, size = 2, stroke = 0.8) + 
    # geom_line(aes(x = grid, y = true), size = 1) + 
    geom_line(aes(x = grid, y = sb), color = "brown1") + 
    geom_line(aes(x = grid, y = cdp), color = "blue") + 
    geom_line(aes(x = grid, y = sb025), color = "brown1", linetype = "dashed") + 
    geom_line(aes(x = grid, y = sb975), color = "brown1", linetype = "dashed") + 
    geom_line(aes(x = grid, y = cdp025), color = "blue", linetype = "dashed") + 
    geom_line(aes(x = grid, y = cdp975), color = "blue", linetype = "dashed") + 
    labs(fill = "", x = "Data", y = ylab_name, caption = caption_name) + 
    theme(legend.position = "none",
          plot.margin = unit(c(0.5,0.1,1,0.1), "cm"),
          plot.caption = element_text(hjust = 0.5, vjust = -3, size = 24),
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 18),
          axis.text.x = element_text(margin = unit(c(0.3, 0, 0.3, 0), "cm")),
          axis.text.y = element_text(margin = unit(c(0, 0.3, 0, 0.3), "cm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  plot_hist_save
  
}

# Plot function for PMTD -------------------------------------------------------
plot_pmtd_est_marg <- function(res_sb, res_cdp, gridval, true_pmf, obs, ylab_name, caption_name) {
  
  yy <- sapply(gridval, function(x) dpois(x, res_sb$la / (1 - res_sb$th)))
  yy_mu <- colMeans(yy)
  yy_qq <- apply(yy, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
  zz <- sapply(gridval, function(x) dpois(x, res_cdp$la / (1 - res_cdp$th)))
  zz_mu <- colMeans(zz)
  zz_qq <- apply(zz, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
  dens <- data.frame(grid = gridval, true = true_pmf,
                     sb = yy_mu, sb025 = yy_qq[1,], sb975 = yy_qq[2,],
                     cdp = zz_mu, cdp025 = zz_qq[1,], cdp975 = zz_qq[2,])
  dat = data.frame(obs = obs)
  
  plot_hist_save <- ggplot(data = dens) + theme_bw() + ylim(0, 0.12) + 
    geom_histogram(data = dat, aes(x = obs, y=..density..), binwidth = 4, fill="white", color="darkgrey", size = 1) + 
    geom_point(aes(x = grid, y = true), shape = 1, size = 2, stroke = 0.8) + 
    # geom_line(aes(x = grid, y = true), size = 1) + 
    geom_line(aes(x = grid, y = sb), color = "brown1") +
    geom_line(aes(x = grid, y = cdp), color = "blue") + 
    geom_line(aes(x = grid, y = sb025), color = "brown1", linetype = "dashed") + 
    geom_line(aes(x = grid, y = sb975), color = "brown1", linetype = "dashed") +
    geom_line(aes(x = grid, y = cdp025), color = "blue", linetype = "dashed") +
    geom_line(aes(x = grid, y = cdp975), color = "blue", linetype = "dashed") + 
    labs(fill = "", x = "Data", y = ylab_name, caption = caption_name) + 
    theme(legend.position = "none",
          plot.margin = unit(c(0.5,0.1,1,0.1), "cm"),
          plot.caption = element_text(hjust = 0.5, vjust = -3, size = 24),
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 18),
          axis.text.x = element_text(margin = unit(c(0.3, 0, 0.3, 0), "cm")),
          axis.text.y = element_text(margin = unit(c(0, 0.3, 0, 0.3), "cm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  plot_hist_save
  
}

# Plot the estimated marginals -------------------------------------------------
xx <- seq(0, max(sim_data$y) + 2, by = 1)
true_nb_pmf <- dnbinom(xx, true_param$kap, true_param$eta / (true_param$la + true_param$ga + true_param$eta))

png(paste0(path2, "est_marg_", 1, "_", 5, ".png"), units = "in", width = 6, height = 5, res = 300)
print(plot_nbmtd_est_marg(res_sb_list[[1]][[1]], res_cdp_list[[1]][[1]], 
                          xx, true_nb_pmf, sim_data$y, "Density", "(a) NBMTD (L = 5)"))
dev.off()

png(paste0(path2, "est_marg_", 1, "_", 15, ".png"), units = "in", width = 6, height = 5, res = 300)
print(plot_nbmtd_est_marg(res_sb_list[[2]][[1]], res_cdp_list[[2]][[1]], 
                          xx, true_nb_pmf, sim_data$y, "", "(b) NBMTD (L = 15)"))
dev.off()

png(paste0(path2, "est_marg_", 2, "_", 5, ".png"), units = "in", width = 6, height = 5, res = 300)
print(plot_pmtd_est_marg(res_sb_list[[1]][[2]], res_cdp_list[[1]][[2]], 
                         xx, true_nb_pmf, sim_data$y, "Density", "(c) PMTD (L = 5)"))
dev.off()

png(paste0(path2, "est_marg_", 2, "_", 15, ".png"), units = "in", width = 6, height = 5, res = 300)
print(plot_pmtd_est_marg(res_sb_list[[2]][[2]], res_cdp_list[[2]][[2]], 
                         xx, true_nb_pmf, sim_data$y, "", "(d) PMTD (L = 15)"))
dev.off()

################################################################################
### Figure 5 of Supplementary Material
################################################################################
# a function that produces the plot and empirical coverage probabilities for the NBMTD
get_nbmtd_pred_intl <- function(res_sb, res_cdp, mtdorder, lags, probs, trun_obs, caption_name=NULL) {
  
  res_sb_pred <- predMTD(res = res_sb, family = "negative binomial", lags = lags, probs = probs)
  res_cdp_pred <- predMTD(res = res_cdp, family = "negative binomial", lags = lags, probs = probs)
  
  kk <- length(trun_obs)
  pred <- as.data.frame(cbind(trun_obs, t(res_sb_pred), t(res_cdp_pred)))
  pred <- cbind(1:kk, pred)
  colnames(pred) <- c("date", "obs", "m1_q025", "m1_q975", "m2_q025", "m2_q975")
  
  cover_sb_boolen <- sapply(1:kk, function(x) pred[x,2] >= pred[x,3] & pred[x,2] <= pred[x,4])
  cover_cdp_boolen <- sapply(1:kk, function(x) pred[x,2] >= pred[x,5] & pred[x,2] <= pred[x,6])
  pred_sb_cover <- sum(cover_sb_boolen) / length(cover_sb_boolen)
  pred_cdp_cover <- sum(cover_cdp_boolen) / length(cover_cdp_boolen)
  
  plot_pred_save <- ggplot(data = pred[1:800,]) + theme_bw() + ylim(0, 62) + 
    geom_point(aes(x = date, y = obs), size = 2, shape = 1) + 
    geom_line(aes(x = date, y = m1_q025), color = "brown1") + 
    geom_line(aes(x = date, y = m1_q975), color = "brown1") + 
    geom_line(aes(x = date, y = m2_q025), color = "blue", linetype = "dashed") + 
    geom_line(aes(x = date, y = m2_q975), color = "blue", linetype = "dashed") + 
    labs(x = "", y = "Count", caption = caption_name) + 
    theme(legend.position = "none",
          plot.margin = unit(c(0, 0.1, 0.5, 0.1), "cm"),
          plot.caption = element_text(hjust = 0.5, vjust = -1.2, size = 28),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 24),
          axis.text.x = element_text(margin = unit(c(0.3, 0, 0.3, 0), "cm")),
          axis.text.y = element_text(margin = unit(c(0, 0.3, 0, 0.3), "cm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  list(pred_sb_cover = pred_sb_cover,
       pred_cdp_cover = pred_cdp_cover,
       plot_pred_save = plot_pred_save)  
  
}

# a function that produces the plot and empirical coverage probabilities for the PMTD
get_pmtd_pred_intl <- function(res_sb, res_cdp, mtdorder, lags, probs, trun_obs, caption_name=NULL) {
  
  res_sb_pred <- predMTD(res = res_sb, family = "poisson", lags = lags, probs = probs)
  res_cdp_pred <- predMTD(res = res_cdp, family = "poisson", lags = lags, probs = probs)
  
  kk <- length(trun_obs)
  pred <- as.data.frame(cbind(trun_obs, t(res_sb_pred), t(res_cdp_pred)))
  pred <- cbind(1:kk, pred)
  colnames(pred) <- c("date", "obs", "m1_q025", "m1_q975", "m2_q025", "m2_q975")
  
  cover_sb_boolen <- sapply(1:kk, function(x) pred[x,2] >= pred[x,3] & pred[x,2] <= pred[x,4])
  cover_cdp_boolen <- sapply(1:kk, function(x) pred[x,2] >= pred[x,5] & pred[x,2] <= pred[x,6])
  pred_sb_cover <- sum(cover_sb_boolen) / length(cover_sb_boolen)
  pred_cdp_cover <- sum(cover_cdp_boolen) / length(cover_cdp_boolen)
  
  plot_pred_save <- ggplot(data = pred[1:800,]) + theme_bw() + ylim(0, 62) + 
    geom_point(aes(x = date, y = obs), size = 2, shape = 1) + 
    geom_line(aes(x = date, y = m1_q025), color = "brown1") + 
    geom_line(aes(x = date, y = m1_q975), color = "brown1") + 
    geom_line(aes(x = date, y = m2_q025), color = "blue", linetype = "dashed") + 
    geom_line(aes(x = date, y = m2_q975), color = "blue", linetype = "dashed") + 
    labs(x = "", y = "Count", caption = caption_name) + 
    theme(legend.position = "none",
          plot.margin = unit(c(0, 0.1, 0.5, 0.1), "cm"),
          plot.caption = element_text(hjust = 0.5, vjust = -1.2, size = 28),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 24),
          axis.text.x = element_text(margin = unit(c(0.3, 0, 0.3, 0), "cm")),
          axis.text.y = element_text(margin = unit(c(0, 0.3, 0, 0.3), "cm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  list(pred_sb_cover = pred_sb_cover,
       pred_cdp_cover = pred_cdp_cover,
       plot_pred_save = plot_pred_save)  
  
}

# obtain predictive intervals and coverage probabilities -----------------------
set.seed(42)
lags_L5 <- lagmat(sim_data$y, 5)
lags_L15 <- lagmat(sim_data$y, 15)
trun_obs_L5 <- sim_data$y[-(1:5)]
trun_obs_L15 <- sim_data$y[-(1:15)]
probs <- c(0.025, 0.975)
nbmtd_pred_L5 <- get_nbmtd_pred_intl(res_sb_list[[1]][[1]], res_cdp_list[[1]][[1]], 5, lags_L5, probs, trun_obs_L5, "(a) NBMTD (L = 5)")
nbmtd_pred_L15 <- get_nbmtd_pred_intl(res_sb_list[[2]][[1]], res_cdp_list[[2]][[1]], 15, lags_L15, probs, trun_obs_L15, "(b) NBMTD (L = 15)")
pmtd_pred_L5 <- get_pmtd_pred_intl(res_sb_list[[1]][[2]], res_cdp_list[[1]][[2]], 5, lags_L5, probs, trun_obs_L5, "(c) PMTD (L = 5)")
pmtd_pred_L15 <- get_pmtd_pred_intl(res_sb_list[[2]][[2]], res_cdp_list[[2]][[2]], 15, lags_L15, probs, trun_obs_L15, "(d) PMTD (L = 15)")

# panels in Figure 5
png(paste0(path2, "pred", 1, "_", 5, ".png"), units = "in", width = 18, height = 6, res = 300)
print(nbmtd_pred_L5$plot_pred_save)
dev.off()

png(paste0(path2, "pred", 1, "_", 15, ".png"), units = "in", width = 18, height = 6, res = 300)
print(nbmtd_pred_L15$plot_pred_save)
dev.off()

png(paste0(path2, "pred", 2, "_", 5, ".png"), units = "in", width = 18, height = 6, res = 300)
print(pmtd_pred_L5$plot_pred_save)
dev.off()

png(paste0(path2, "pred", 2, "_", 15, ".png"), units = "in", width = 18, height = 6, res = 300)
print(pmtd_pred_L15$plot_pred_save)
dev.off()

################################################################################
### Table 1 of Supplementary Material
################################################################################
L5_cover <- c(nbmtd_pred_L5$pred_sb_cover, nbmtd_pred_L5$pred_cdp_cover, pmtd_pred_L5$pred_sb_cover, pmtd_pred_L5$pred_cdp_cover)
L15_cover <- c(nbmtd_pred_L15$pred_sb_cover, nbmtd_pred_L15$pred_cdp_cover, pmtd_pred_L15$pred_sb_cover, pmtd_pred_L15$pred_cdp_cover)
pred_tbl <- rbind(L5_cover, L15_cover)
colnames(pred_tbl) <- c("NBMTD-SB", "NBMTD-CDP", "PMTD-SB", "PMTD-CDP")
rownames(pred_tbl) <- c("L = 5", "L = 15")
knitr::kable(pred_tbl, digits = 3, format = "latex")

