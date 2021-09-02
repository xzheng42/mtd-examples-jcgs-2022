################################################################################
### Description: make plots of the results of crime data example
### Output: 1) posterior summaries of model parameters;
###         2) Figure 3 of the paper;
###         3) Figure 3 of the paper (Figure 6 of the Supplementary Material).
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
load(paste0(path1, 'res_precip_L10.RData'))

################################################################################
### posterior estimate
################################################################################
post_alpha <- summarize(rbind(res1$alpha, res2$alpha), k = 2)
post_phi <- summarize(rbind(res1$phi, res2$phi), k = 2)
post_bb <- cbind(summarize(res1$bb, k = 2), summarize(res2$bb, k = 2))
tbl <- rbind(post_alpha, post_phi, post_bb)
rownames(tbl) <- c("$\\alpha$", "$\\phi$", "$\\beta_1$", "$\\beta_2$", 
                   "$\\beta_3$", "$\\beta_4$", "$\\beta_5$", "$\\beta_6$")
colnames(tbl) <- c("LMTD(5)-SB", "LMTD(5)-CDP")
print(tbl)
knitr::kable(tbl, align = 'c', format = "latex", escape = FALSE)

################################################################################
### Figure 3 
################################################################################
# left panels ------------------------------------------------------------------
plotWeight <- function(res, mtdorder, prior_mu, weight_name) {
  
  plot_break <- seq(0, mtdorder, by = 2)
  plot_label <- seq(0, mtdorder, by = 2)
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
    theme_bw() + ylim(0, 1) + 
    labs(fill = "", x = "Lag", y = "Probability", caption = weight_name) + 
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

mtdorder <- nrow(res1$weight)
cutoff = seq(0, 1, length = mtdorder + 1)
SB_mu <- 1 / (1 + prior$alpha) * (prior$alpha / (1 + prior$alpha))^((1:(mtdorder-1))-1)
SB_prior_mu <- c(SB_mu, 1 - sum(SB_mu))
CDP_prior_mu <- diff(pbeta(cutoff, prior$a_G0, prior$b_G0))

plot_weight1 <- plotWeight(res1, mtdorder, SB_prior_mu, "(a) SB")
plot_weight2 <- plotWeight(res2, mtdorder, CDP_prior_mu, "(b) CDP")
plot_weight <- grid.arrange(plot_weight1, plot_weight2, nrow = 2, ncol = 1)

# right panels -----------------------------------------------------------------
obs <- as.numeric(week_precip_8203)
precip_date <- date(week_precip_8203)
idx <- which(precip_date == "2000-01-02")
size <- length(obs) - idx + 1
selec_date <- tail(precip_date, size)
selec_dat <- cbind(selec_date, as.data.frame(tail(obs, size)))
colnames(selec_dat) <- c("date", "obs")

set.seed(2)
simu1 <- simRegLMTD(mtdorder, xx, rowMeans(res1$weight), mean(res1$alpha), mean(res1$phi),
                    rowMeans(res1$bb), head(obs, mtdorder), length(obs)-mtdorder)
simu2 <- simRegLMTD(mtdorder, xx, rowMeans(res2$weight), mean(res2$alpha), mean(res2$phi),
                    rowMeans(res2$bb), head(obs, mtdorder), length(obs)-mtdorder)
selec_simu1 <- cbind(selec_date, as.data.frame(tail(simu1$y, size)))
selec_simu2 <- cbind(selec_date, as.data.frame(tail(simu2$y, size)))
colnames(selec_simu1) <- c("date", "simu1")
colnames(selec_simu2) <- c("date", "simu2")

plot_obs <- ggplot(data = selec_dat) + theme_bw() + 
  geom_line(aes(x = date, y = obs), size = 0.8) + ylim(0, 100) + 
  labs(x = "", y = "Observation") + 
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20),
        axis.text.x = element_text(margin = unit(c(0.3, 0, 0.3, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.3, 0, 0.3), "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot_simu1 <- ggplot(data = selec_simu1) + theme_bw() + 
  geom_line(aes(x = date, y = simu1), size = 0.8) + ylim(0, 100) + 
  labs(x = "", y = "Simulated data (SB)") + 
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), "cm"),
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20),
        axis.text.x = element_text(margin = unit(c(0.3, 0, 0.3, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.3, 0, 0.3), "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

plot_simu2 <- ggplot(data = selec_simu2) + theme_bw() + 
  geom_line(aes(x = date, y = simu2), size = 0.8) + ylim(0, 100) + 
  labs(x = "Year", y = "Simulated data (CDP)", 
       caption = "(c) Observed and simulated precipitation data") + 
  theme(plot.margin = unit(c(0, 0.5,1,0.5), "cm"),
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 20),
        axis.text.x = element_text(margin = unit(c(0.3, 0, 0.3, 0), "cm")),
        axis.text.y = element_text(margin = unit(c(0, 0.3, 0, 0.3), "cm")),
        plot.caption = element_text(hjust = 0.5, vjust = -3, size = 26),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

png(paste0(path2, "all_in_one", ".png"), units = "in", width = 18, height = 11, res = 300)
grid.arrange(plot_weight, 
             arrangeGrob(plot_obs, plot_simu1, plot_simu2, nrow = 3, heights = c(1,1,1.2)), 
             nrow = 1, widths = c(1,2))
dev.off()

################################################################################
### Figure 3 of the paper (Figure 6 of the Supplementary Material)
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
    labs(x = "Theoretical quantiles", y = "Sample quantiles", caption = "(a) Quantile-quantile plot") + 
    theme(plot.margin = unit(c(0.5,0.5,1,0.5), "cm"),
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
    labs(x = "Posterior mean residuals", y = "Density", caption = "(b) Histogram") + 
    scale_x_continuous(breaks = seq(-3, 3, by = 1)) + ylim(0, 0.45) + 
    geom_line(data = yy, aes(x = zz, y = yy), size = 1) +
    geom_line(data = mu_kde, aes(x = xx, y = yy), linetype = "dashed", size = 0.8) +
    theme(plot.margin = unit(c(0.5,0.5,1,0.5), "cm"),
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
    labs(x = "Lag", y = "ACF", caption = "(c) Autocorrelation") + 
    theme(plot.margin = unit(c(0.5,0.5,1,0.5), "cm"),
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




