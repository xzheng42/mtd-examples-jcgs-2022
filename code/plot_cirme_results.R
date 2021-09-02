################################################################################
### Description: make plots of the results of crime data example
### Output: 1) posterior summaries of model parameters;
###         2) Figure 2 of the paper;
###         3) Figure 6 of the Supplementary Material. 
################################################################################
rm(list = ls())
library(mtd)
library(ggplot2)
library(gridExtra)

################################################################################
### load data
################################################################################
path1 <- "path/to/import/data"
path2 <- "path/to/save/results"
load(paste0(path1, 'res_crime_L20.RData'))

################################################################################
### posterior estimate
################################################################################
post_la <- summarize(rbind(res1$la, res2$la))
post_th <- summarize(rbind(res1$th, res2$th))
post_phi <- summarize(rbind(res1$la / (1 - res1$th), res2$la / (1 - res2$th)))
tbl <- rbind(post_la, post_th, post_phi)
rownames(tbl) <- c("$\\lambda$", "$\\theta$", "$\\phi$")
colnames(tbl) <- c("PMTD(20)-SB", "PMTD(20)-CDP")
knitr::kable(t(tbl), align = 'c', format = "latex", escape = FALSE)

################################################################################
### Figure 3
################################################################################
# top panel --------------------------------------------------------------------
set.seed(32)
mtdorder <- nrow(res1$weight)
probs <- c(0.025, 0.975)
obs <- unlist(dat_dt_count[,2])
lags <- lagmat(obs, mtdorder)
res1_pred <- predPMTD(mtdorder = mtdorder, 
                      weight = res1$weight, 
                      lags = lags, 
                      la = res1$la,
                      th = res1$th, 
                      probs = c(0.025, 0.975))
res2_pred <- predPMTD(mtdorder = mtdorder, 
                      weight = res2$weight, 
                      lags = lags, 
                      la = res2$la,
                      th = res2$th, 
                      probs = c(0.025, 0.975))

pred <- as.data.frame(cbind(obs[-(1:mtdorder)], t(res1_pred), t(res2_pred)))
pred <- cbind(dat_dt_count$Date[-(1:mtdorder)], pred)
colnames(pred) <- c("date", "obs", "m1_q025", "m1_q975", "m2_q025", "m2_q975")
plot_pred_save <- ggplot(data = pred) + theme_bw() + 
  geom_point(aes(x = date, y = obs), shape = 1, size = 3) + 
  geom_line(aes(x = date, y = m1_q025), col = "brown1") + 
  geom_line(aes(x = date, y = m1_q975), col = "brown1") + 
  geom_line(aes(x = date, y = m2_q025), col = "blue", linetype = "dashed") + 
  geom_line(aes(x = date, y = m2_q975), col = "blue", linetype = "dashed") + 
  labs(x = "Year", y = "Number of incidents", 
       caption = "(a) 95% one-step posterior predictive intervals for the crime data") + 
  theme(plot.margin = unit(c(0.5,0.5,1,0.5), "cm"),
        plot.caption = element_text(hjust = 0.5, vjust = -3, size = 27),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 21),
        axis.text.x = element_text(margin = unit(c(0.3, 0, 0.3, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.3, 0, 0.3), "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
png(paste0(path2, "crime_pred", ".png"), units = "in", width = 18, height = 7, res = 300)
print(plot_pred_save)
dev.off()

# bottom panels ----------------------------------------------------------------
plotWeight <- function(res, mtdorder, prior_mu, ylab_name, weight_name) {
  
  plot_break <- seq(0, mtdorder, by = 5)
  plot_label <- seq(0, mtdorder, by = 5)
  post_mu <- rowMeans(res$weight)
  post_qq <- apply(res$weight, 1, function(x) quantile(x, probs = c(0.025, 0.975)))
  weight <- cbind(prior_mu, post_mu, t(post_qq))
  weight <- rbind(weight, tail(weight, 1))
  weight <- as.data.frame(cbind(0:mtdorder, weight))
  colnames(weight) <- c("idx", "prior_mu", "post_mu", "q025", "q975")
  
  plot_weight_save <- ggplot(weight, aes(x = idx)) + 
    pammtools::geom_stepribbon(aes(ymin = q025, ymax = q975, fill = "Posterior 95% CI"), alpha = 0.15) + 
    geom_step(aes(y = post_mu, color = "Posterior mean"), size = 1.5) + 
    geom_step(aes(y = prior_mu, color = "Prior mean"), linetype = "dashed", size = 1.5) + 
    scale_fill_manual("", values = "blue") + 
    scale_colour_manual("", values = c("blue", "brown1","black")) + 
    scale_x_continuous(breaks = plot_break, labels = plot_label) + 
    theme_bw() + ylim(0, 0.45) + 
    labs(fill = "", x = "Lag", y = ylab_name, caption = weight_name) + 
    theme(legend.position = "none",
          plot.margin = unit(c(0.5,0.5,1,0.5), "cm"),
          plot.caption = element_text(hjust = 0.5, vjust = -3, size = 26),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          axis.text.x = element_text(margin = unit(c(0.3, 0, 0.3, 0), "cm")),
          axis.text.y = element_text(margin = unit(c(0, 0.3, 0, 0.3), "cm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  plot_weight_save
  
}

cutoff = seq(0, 1, length = mtdorder + 1)
SB_mu <- 1 / (1 + prior$alpha) * (prior$alpha / (1 + prior$alpha))^((1:(mtdorder-1))-1)
SB_prior_mu <- c(SB_mu, 1 - sum(SB_mu))
CDP_prior_mu <- diff(pbeta(cutoff, prior$a_G0, prior$b_G0))

png(paste0(path2, "20_SB", ".png"), units = "in", width = 9, height = 5.5, res = 300)
plotWeight(res1, mtdorder, SB_prior_mu, "Probability", "(b) SB")
dev.off()

png(paste0(path2, "20_CDP", ".png"), units = "in", width = 9, height = 5.5, res = 300)
plotWeight(res2, mtdorder, CDP_prior_mu, "", "(c) CDP")
dev.off()

################################################################################
### Figure 6 of the Supplementary Material 
################################################################################
plotRQR <- function(rqr){
  
  q_mu <- rowMeans(rqr)
  q_025 <- apply(rqr, 1, function(x) quantile(x, 0.025))
  q_975 <- apply(rqr, 1, function(x) quantile(x, 0.975))
  qq <- qqnorm(q_mu[order(q_mu)], plot.it = FALSE)
  resid <- data.frame(z = qq$x, mu = qq$y, q025 = q_025, q975 = q_975)
  
  # qq-plot
  plot_gg_save <- ggplot(data = resid) + theme_bw() + 
    geom_abline(slope = 1, intercept = 0) + 
    geom_point(aes(x = z, y = mu), size = 0.5) + 
    geom_line(aes(x = z, y = q025[order(q025)]), linetype = "dashed", size = 0.8) + 
    geom_line(aes(x = z, y = q975[order(q975)]), linetype = "dashed", size = 0.8) +
    scale_x_continuous(breaks = seq(-3, 3, by = 1)) + 
    scale_y_continuous(breaks = seq(-3, 3, by = 1)) + 
    labs(x = "Theoretical quantiles", y = "Sample quantiles") + 
    theme(plot.margin = unit(c(1.5,0.5,1,0.5), "cm"),
          plot.caption = element_text(hjust = 0.5, vjust = -3, size = 26),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 20),
          axis.text.x = element_text(margin = unit(c(0.3, 0, 0.3, 0), "cm")),
          axis.text.y = element_text(margin = unit(c(0, 0.3, 0, 0.3), "cm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # histogram
  zz <- seq(-3.5, 3.5, length = 100)
  yy <- data.frame(zz = zz, yy = dnorm(zz))
  q_mu_kde <- density(q_mu, bw = 0.5)
  mu_kde <- data.frame(xx = q_mu_kde$x, yy = q_mu_kde$y)
  mu_kde <- mu_kde[mu_kde$xx >= -3.5 & mu_kde$xx <= 3.5,]
  plot_hist_save <- ggplot(data = resid, aes(x = mu)) + theme_bw() + 
    geom_histogram(aes(y=..density..), binwidth = 0.7, fill="white", color="darkgrey", size = 0.8) + 
    labs(x = "Posterior mean residuals", y = "Density") + 
    scale_x_continuous(breaks = seq(-3, 3, by = 1)) + ylim(0, 0.45) + 
    geom_line(data = yy, aes(x = zz, y = yy), size = 1) +
    geom_line(data = mu_kde, aes(x = xx, y = yy), linetype = "dashed", size = 0.8) +
    theme(plot.margin = unit(c(1.5,0.5,1,0.5), "cm"),
          plot.caption = element_text(hjust = 0.5, vjust = -3, size = 26),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 20),
          axis.text.x = element_text(margin = unit(c(0.3, 0, 0.3, 0), "cm")),
          axis.text.y = element_text(margin = unit(c(0, 0.3, 0, 0.3), "cm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # acf
  sacf <- acf(q_mu, plot = FALSE)
  band <- qnorm((1 + 0.95) / 2) / sqrt(sacf$n.used)
  sacf <- data.frame(lag = sacf$lag, acf = sacf$acf)
  plot_acf_save <- ggplot(data = sacf, aes(x = lag, y = acf)) + theme_bw() + 
    geom_hline(aes(yintercept = 0)) + 
    geom_segment(aes(xend = lag, yend = 0)) + 
    geom_hline(aes(yintercept = -band), linetype = "dashed") + 
    geom_hline(aes(yintercept = band), linetype = "dashed") + 
    labs(x = "Lag", y = "ACF") + 
    theme(plot.margin = unit(c(1.5,0.5,1,0.5), "cm"),
          plot.caption = element_text(hjust = 0.5, vjust = -3, size = 26),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 20),
          axis.text.x = element_text(margin = unit(c(0.3, 0, 0.3, 0), "cm")),
          axis.text.y = element_text(margin = unit(c(0, 0.3, 0, 0.3), "cm")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  list(plot_gg_save, plot_hist_save, plot_acf_save)
  
}

res1_rar_plot <- plotRQR(res1_rqr)
res2_rar_plot <- plotRQR(res2_rqr)

png(paste0(path2, "mc_sb", ".png"), units = "in", width = 18, height = 6, res = 300)
grid.arrange(res1_rar_plot[[1]], res1_rar_plot[[2]], res1_rar_plot[[3]], nrow = 1)
dev.off()

png(paste0(path2, "mc_cdp", ".png"), units = "in", width = 18, height = 6, res = 300)
grid.arrange(res2_rar_plot[[1]], res2_rar_plot[[2]], res2_rar_plot[[3]], nrow = 1)
dev.off()

