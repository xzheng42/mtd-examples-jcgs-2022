################################################################################
### Description: Make plots of the weights of the GMTD model ran in the second
###              scenario of the simulation study presented in the paper
### Output: Figure 1 of the paper, and Figures 1-2 of the Supplementary Material
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
load(paste0(path1, 'res_sim2.RData'))

################################################################################
### plot function
################################################################################
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
    theme_bw() + ylim(0, 0.6) + labs(fill = "", x = xlab_name, y = ylab_name) + 
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

################################################################################
### Figure 1 of the paper
################################################################################
mtdorder <- 15
alpha <- 2
a_G0 <- 1
b_G0 <- 6
dir_shape <- rep(1 / mtdorder, mtdorder)
cutoff_L15 <- seq(0, 1, length = mtdorder + 1)

DIR_prior_mu <- diff(pbeta(cutoff_L15, 1, 1))
SB_mu <- 1 / (1 + alpha) * (alpha / (1 + alpha))^((1:(mtdorder - 1)) - 1) 
SB_prior_mu <- c(SB_mu, 1 - sum(SB_mu))
CDP_prior_mu <- diff(pbeta(cutoff_L15, a_G0, b_G0))

png(paste0(path2, "15_DIR", ".png"), units = "in", width = 6, height = 5.3, res = 300)
plotWeight(res_L15[[1]], mtdorder, DIR_prior_mu, true_weight, "Lag",  "Probability", "(a) DIR")
dev.off()

png(paste0(path2, "15_SB", ".png"), units = "in", width = 6, height = 5.3, res = 300)
plotWeight(res_L15[[2]], mtdorder, SB_prior_mu, true_weight, "Lag", "", "(b) SB")
dev.off()

png(paste0(path2, "15_CDP", ".png"), units = "in", width = 6, height = 5.3, res = 300)
plotWeight(res_L15[[3]], mtdorder, CDP_prior_mu, true_weight, "Lag", "", "(c) CDP")
dev.off()

################################################################################
### Figures 1-2 of the Supplementary Material
################################################################################
# order L = 5 ------------------------------------------------------------------
mtdorder <- 5
alpha <- 1
a_G0 <- 1
b_G0 <- 3
dir_shape <- rep(1 / mtdorder, mtdorder)
cutoff_L5 <- seq(0, 1, length = mtdorder + 1)

DIR_prior_mu <- diff(pbeta(cutoff_L5, 1, 1))
SB_mu <- 1 / (1 + alpha) * (alpha / (1 + alpha))^((1:(mtdorder - 1)) - 1) 
SB_prior_mu <- c(SB_mu, 1 - sum(SB_mu))
CDP_prior_mu <- diff(pbeta(cutoff_L5, a_G0, b_G0))

png(paste0(path2, "5_DIR", ".png"), units = "in", width = 6, height = 5, res = 300)
plotWeight(res_L5[[1]], mtdorder, DIR_prior_mu, true_weight, "",  "Probability")
dev.off()

png(paste0(path2, "5_SB", ".png"), units = "in", width = 6, height = 5, res = 300)
plotWeight(res_L5[[2]], mtdorder, SB_prior_mu, true_weight)
dev.off()

png(paste0(path2, "5_CDP", ".png"), units = "in", width = 6, height = 5, res = 300)
plotWeight(res_L5[[3]], mtdorder, CDP_prior_mu, true_weight)
dev.off()

# order L = 15 -----------------------------------------------------------------
mtdorder <- 15
alpha <- 2
a_G0 <- 1
b_G0 <- 6
dir_shape <- rep(1 / mtdorder, mtdorder)
cutoff_L15 <- seq(0, 1, length = mtdorder + 1)

DIR_prior_mu <- diff(pbeta(cutoff_L15, 1, 1))
SB_mu <- 1 / (1 + alpha) * (alpha / (1 + alpha))^((1:(mtdorder - 1)) - 1) 
SB_prior_mu <- c(SB_mu, 1 - sum(SB_mu))
CDP_prior_mu <- diff(pbeta(cutoff_L15, a_G0, b_G0))

png(paste0(path2, "15_DIR_SM", ".png"), units = "in", width = 6, height = 5, res = 300)
plotWeight(res_L15[[1]], mtdorder, DIR_prior_mu, true_weight, "",  "Probability")
dev.off()

png(paste0(path2, "15_SB_SM", ".png"), units = "in", width = 6, height = 5, res = 300)
plotWeight(res_L15[[2]], mtdorder, SB_prior_mu, true_weight)
dev.off()

png(paste0(path2, "15_CDP_SM", ".png"), units = "in", width = 6, height = 5, res = 300)
plotWeight(res_L15[[3]], mtdorder, CDP_prior_mu, true_weight)
dev.off()

# order L = 25 -----------------------------------------------------------------
mtdorder <- 25
alpha <- 3
a_G0 <- 1
b_G0 <- 7
dir_shape <- rep(1 / mtdorder, mtdorder)
cutoff_L25 <- seq(0, 1, length = mtdorder + 1)

DIR_prior_mu <- diff(pbeta(cutoff_L25, 1, 1))
SB_mu <- 1 / (1 + alpha) * (alpha / (1 + alpha))^((1:(mtdorder - 1)) - 1) 
SB_prior_mu <- c(SB_mu, 1 - sum(SB_mu))
CDP_prior_mu <- diff(pbeta(cutoff_L25, a_G0, b_G0))

png(paste0(path2, "25_DIR", ".png"), units = "in", width = 6, height = 5.3, res = 300)
plotWeight(res_L25[[1]], mtdorder, DIR_prior_mu, true_weight, "Lag",  "Probability", "(a) DIR")
dev.off()

png(paste0(path2, "25_SB", ".png"), units = "in", width = 6, height = 5.3, res = 300)
plotWeight(res_L25[[2]], mtdorder, SB_prior_mu, true_weight, "Lag", "", "(b) SB")
dev.off()

png(paste0(path2, "25_CDP", ".png"), units = "in", width = 6, height = 5.3, res = 300)
plotWeight(res_L25[[3]], mtdorder, CDP_prior_mu, true_weight, "Lag", "", "(c) CDP")
dev.off()

