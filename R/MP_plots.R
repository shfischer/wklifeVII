### ------------------------------------------------------------------------ ###
### create plots for publication ####
### ------------------------------------------------------------------------ ###
### uses data prepared in MP_analysis.R

library(FLCore)
library(ggplotFL)
library(tidyverse)
library(foreach)
library(cowplot)
library(Cairo)
library(doParallel)

### path to final scenario
path_default <- "output/observation_error/new_baseline/"

### load lhist
stocks <- read.csv("input/stock_list_full2.csv", stringsAsFactors = FALSE)
### reference points
refpts <- readRDS("input/refpts.rds")

### theme for plotting
theme_paper <- theme_bw(base_size = 8, base_family = "serif")

### ------------------------------------------------------------------------ ###
### default catch rule & fishing history ####
### ------------------------------------------------------------------------ ###

### load MSE results
files <- paste0(path_default, "corrected/", 
                rep(c("one-way", "roller-coaster"), each = 29), 
                "/", stocks$stock, ".rds")
names(files) <- rep(stocks$stock_short, 2)
res_corrected <- lapply(files, readRDS)

### extract quants
res_corrected <- foreach(x = seq_along(res_corrected),
                         fhist = rep(c("one-way", "roller-coaster"), each = 29),
                         stock = names(files)) %do% {
  tmp <- res_corrected[[x]][c("fbar_rel", "ssb_rel", "catch_rel")]
  tmp <- lapply(tmp, iterMedians)
  tmp <- as.data.frame(FLQuants(tmp))
  tmp$fhist = fhist
  tmp$stock = stock
  return(tmp)
}
df_res <- do.call(rbind, res_corrected)

### format for plotting
df_res <- df_res %>% spread(qname, data) %>%
  mutate(year = year - 100,
         group = ifelse(year < 1, "history", "projection"))
df_res <- df_res %>% bind_rows(df_res %>% filter(year == 0) %>%
                             mutate(group = "projection"))
df_res <- df_res %>% 
  mutate(group = factor(group, levels = c("history", "projection")))

### plot elements
plot_ssb_rel <- df_res %>%
  ggplot(aes(x = year, y = ssb_rel, group = stock)) +
  geom_line(size = 0.15) +
  facet_grid(fhist ~ group, scales = "free_x", space = "free_x") +
  theme_paper +
  scale_x_continuous(breaks = c(-25, 0, 25, 50, 75, 100), expand = c(0, 0)) +
  labs(y = expression(italic(SSB/B[MSY]))) +
  theme(strip.text.y = element_blank(), panel.spacing.x = unit(0, units = "cm"),
        plot.margin = unit(x = c(4, 6, 4, 4), units = "pt"))
plot_fbar_rel <- df_res %>%
  ggplot(aes(x = year, y = fbar_rel, group = stock)) +
  geom_line(size = 0.15) +
  facet_grid(fhist ~ group, scales = "free_x", space = "free_x") +
  theme_paper +
  scale_x_continuous(breaks = c(-25, 0, 25, 50, 75, 100), expand = c(0, 0)) +
  labs(y = expression(italic(F/F[MSY]))) +
  theme(strip.text.y = element_blank(), panel.spacing.x = unit(0, units = "cm"),
        plot.margin = unit(x = c(4, 6, 4, 4), units = "pt"))
plot_catch_rel <- df_res %>%
  ggplot(aes(x = year, y = catch_rel, group = stock)) +
  geom_line(size = 0.15) +
  facet_grid(fhist ~ group, scales = "free_x", space = "free_x") +
  theme_paper +
  scale_x_continuous(breaks = c(-25, 0, 25, 50, 75, 100), expand = c(0, 0)) +
  labs(y = "catch/MSY") +
  theme(panel.spacing.x = unit(0, units = "cm"))

### combine the three plots
plot_grid(plot_ssb_rel, plot_fbar_rel, plot_catch_rel, 
          rel_widths = c(1, 1, 1.1), nrow = 1)
### and save
ggsave(filename = "output/plots/paper_revision/trajectories_rel.png",
       width = 17, height = 7, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper_revision/trajectories_rel.jpeg", 
       quality = 100,
       width = 17, height = 7, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper_revision/trajectories_rel.pdf",
       width = 17, height = 7, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### penalized regression ####
### ------------------------------------------------------------------------ ###
library(glmnet)

### load data
glmnet_res <- readRDS(paste0(path_default, "corrected/one-way/glmnet.rds"))

### plots
p1 <- ggplot(data = glmnet_res$data,
             aes(x = k, y = f_rel)) +
  geom_point(size = 0.3) + theme_paper + 
  lims(x = c(0, NA), y = c(0, NA)) +
  geom_abline(data = data.frame(
    intercept = c(glmnet_res$glmnet$f_rel["(Intercept)",],
                  glmnet_res$lms[glmnet_res$lms$var == "f_rel", "Intercept"]),
    slope = c(glmnet_res$glmnet$f_rel_new["k",],
              glmnet_res$lms[glmnet_res$lms$var == "f_rel", "k"]),
    model = c("lasso regression", "linear regression")),
    aes(intercept = intercept, slope = slope, linetype = model), size = 0.3, 
    show.legend = FALSE) +
  labs(y = expression(italic(F/F[MSY])), x = "") +
  theme(axis.title.x = element_blank())
p2 <- ggplot(data = glmnet_res$data,
             aes(x = k, y = ssb_rel)) +
  geom_point(size = 0.3) + theme_paper + 
  lims(x = c(0, NA), y = c(0, NA)) +
  geom_abline(data = data.frame(
    intercept = c(glmnet_res$glmnet$ssb_rel["(Intercept)",],
                  glmnet_res$lms[glmnet_res$lms$var == "ssb_rel", "Intercept"]),
    slope = c(glmnet_res$glmnet$ssb_rel["k",],
              glmnet_res$lms[glmnet_res$lms$var == "ssb_rel", "k"]),
    model = c("lasso regression", "linear regression")),
    aes(intercept = intercept, slope = slope, linetype = model), size = 0.3, 
    show.legend = FALSE) +
  labs(y = expression(italic(SSB/B[MSY])), x = "") +
  theme(axis.title.x = element_blank())
p3 <- ggplot(data = glmnet_res$data,
             aes(x = k, y = risk_collapse)) +
  geom_point(size = 0.3) + theme_paper + 
  lims(x = c(0, NA), y = c(0, NA)) +
  geom_abline(data = data.frame(
    intercept = c(glmnet_res$glmnet$risk_collapse["(Intercept)",],
                  glmnet_res$lms[glmnet_res$lms$var == "risk_collapse", 
                                 "Intercept"]),
    slope = c(glmnet_res$glmnet$risk_collapse["k",],
              glmnet_res$lms[glmnet_res$lms$var == "risk_collapse", "k"]),
    model = c("lasso regression", "linear regression")),
    aes(intercept = intercept, slope = slope, linetype = model), size = 0.3, 
    show.legend = FALSE) +
  labs(y = expression(collapse~risk), x = "") +
  theme(axis.title.x = element_blank())
p4 <- ggplot(data = glmnet_res$data,
             aes(x = k, y = risk_blim)) +
  geom_point(size = 0.3) + theme_paper + 
  lims(x = c(0, NA), y = c(0, NA)) +
  geom_abline(data = data.frame(
    intercept = c(glmnet_res$glmnet$risk_blim["(Intercept)",],
                  glmnet_res$lms[glmnet_res$lms$var == "risk_blim", 
                                 "Intercept"]),
    slope = c(glmnet_res$glmnet$risk_blim["k",],
              glmnet_res$lms[glmnet_res$lms$var == "risk_blim", "k"]),
    model = c("lasso regression", "linear regression")),
    aes(intercept = intercept, slope = slope, linetype = model), size = 0.3, 
    show.legend = FALSE) +
  labs(y = expression(italic(B[lim])~risk), x = "") +
  theme(axis.title.x = element_blank())
p5 <- ggplot(data = glmnet_res$data[, ],
             aes(x = k, y = yield_rel)) +
  geom_point(size = 0.3) + theme_paper + 
  lims(x = c(0, NA), y = c(0, NA)) +
  geom_abline(data = data.frame(
    intercept = c(glmnet_res$glmnet$yield_rel["(Intercept)",],
                  glmnet_res$lms[glmnet_res$lms$var == "yield_rel", 
                                 "Intercept"]),
    slope = c(glmnet_res$glmnet$yield_rel["k",],
              glmnet_res$lms[glmnet_res$lms$var == "yield_rel", "k"]),
    model = c("lasso regression", "linear regression")),
    aes(intercept = intercept, slope = slope, linetype = model), size = 0.3, 
    show.legend = FALSE) +
  labs(y = expression(catch/MSY), x = expression(italic(k)))
p6 <- ggplot(data = glmnet_res$data[, ],
             aes(x = k, y = iav)) +
  geom_point(size = 0.3) + theme_paper + 
  lims(x = c(0, NA), y = c(0, NA)) +
  geom_abline(data = 
    data.frame(intercept = c(glmnet_res$glmnet$iav["(Intercept)",],
                             glmnet_res$lms[glmnet_res$lms$var == "iav",
                                            "Intercept"]),
               slope = c(glmnet_res$glmnet$iav["k",],
                         glmnet_res$lms[glmnet_res$lms$var == "iav", "k"]),
               model = c("lasso\nregression", 
                              "linear\nregression")),
              aes(intercept = intercept, slope = slope, linetype = model),
              size = 0.3) +
  labs(y = "ICV", x = expression(italic(k))) +
  theme(legend.position = "bottom", 
        legend.key.size = unit(0.7, units = "lines"))
### extract legend
legend <- get_legend(p6)
### combine plots
plot_grid(plot_grid(p1, p2, p3, p4, p5, p6 + theme(legend.position = "none"),
                    align = "vh", nrow = 3),
          legend, ncol = 1, rel_heights = c(3, 0.2))

### save plot
ggsave(filename = "output/plots/paper_revision/glmnet.png",
       width = 8.5, height = 12, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper_revision/glmnet.jpeg", quality = 100,
       width = 8.5, height = 12, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper_revision/glmnet.pdf",
       width = 8.5, height = 12, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### clustering ####
### ------------------------------------------------------------------------ ###
library(dtwclust)
library(ggdendro)

cluster <- readRDS(paste0(path_default, 
                          "corrected/one-way/plot_data_cluster.rds"))

### plot
p1 <- ggplot(cluster$curves[cluster$curves$k %in% 1:4, ],
             aes(x = year, y = value, group = stock, 
                 linetype = stock == "cluster", alpha = stock == "cluster",
                 size = stock == "cluster", colour = as.factor(cluster))) +
  geom_line() + 
  facet_grid(ifelse(k > 1, 
                    paste(k, "clusters"), 
                    paste(k, "cluster")) ~ paste("cluster", cluster)) +
  theme_paper +
  scale_linetype_manual("", values = c("longdash", "solid"), 
                        labels = c("stock", "cluster")) +
  scale_alpha_manual("", values = c(0.7, 1), 
                     labels = c("stock", "cluster")) +
  scale_size_manual("", values = c(0.15, 0.5), 
                    labels = c("stock", "cluster")) +
  scale_colour_discrete(guide = FALSE) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
  labs(y = expression(italic(SSB/B[MSY]))) +
  theme(legend.position = c(0.87, 0.87),
        legend.background = element_blank(),
        legend.key = element_blank(),
        strip.text.y = element_blank())

p2 <- ggplot(data = cluster$bars[cluster$bars$n_cluster %in% 1:4, ],
             aes(x = k_pos, y = k, fill = as.factor(cluster))) +
  geom_bar(stat = "identity", colour = "black", size = 0.1) +
  scale_fill_discrete("cluster") +
  facet_grid(ifelse(n_cluster > 1, 
                    paste(n_cluster, "clusters"), 
                    paste(n_cluster, "cluster")) ~ label) + 
  theme_paper + labs(x = "stocks") +
  theme(axis.text.x = element_text(colour = "white"),
        axis.ticks.x = element_blank(),
        strip.background.x = element_blank()) +
  labs(y = expression(italic(k)))

p3 <- ggdendrogram(cluster$dend) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text.y = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_text(colour = "black", angle = 90, hjust = 0.5),
        axis.text.x = element_text(size = 8),
        text = element_text(size = 8, family = "serif")) +
  labs(x = "stocks", y = "DTW distance")

p4 <- ggplot(segment(cluster$dend) %>%
         mutate(yend = ifelse(yend == 0, -2, yend))) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), size = 0.4) +
  geom_text(data = label(cluster$dend),
               aes(label = label, x = x, y = 0, colour = as.factor(`4`)),
            show.legend = FALSE, hjust = 1, angle = 90, nudge_y = -5,
            family = "serif", size = 2.5) +
  coord_cartesian(ylim = c(-20, max(cluster$dend$segments$y))) + 
  #coord_trans(y = log) +
  theme_paper +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.title.x = element_blank(),#element_text(colour = NA),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()) +
  labs(y = "DTW distance")

plot_grid(p4, plot_grid(p1, p2, ncol = 2, rel_widths = c(1.7, 1), 
                        labels = c("B", "C"), label_fontfamily = "serif"),
          nrow = 2, rel_heights = c(1, 2), labels = c("A", ""), 
          label_fontfamily = "serif")

ggsave(filename = "output/plots/paper_revision/clustering.png",
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper_revision/clustering.jpeg", quality = 100,
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper_revision/clustering.pdf",
       width = 17, height = 13, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### multiplier ####
### ------------------------------------------------------------------------ ###

data_mult <- readRDS(paste0(path_default, "plot_data_multiplier.rds"))
p1 <- data_mult %>% 
  ggplot(aes(x = HCRmult, y = f_rel, group = paper, 
             linetype = as.factor(cluster), colour = as.factor(cluster))) +
  theme_paper + 
  geom_line(size = 0.2, show.legend = FALSE) +
  ylim(0, 1.5) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "", y = expression(italic(F/F[MSY]))) +
  theme(strip.text = element_blank())
p2 <- data_mult %>% 
  ggplot(aes(x = HCRmult, y = ssb_rel, group = paper, 
             linetype = as.factor(cluster), colour = as.factor(cluster))) +
  theme_paper + 
  geom_line(size = 0.2, show.legend = FALSE) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "", y = expression(italic(SSB/B[MSY]))) +
  theme(strip.text = element_blank())
p3 <- data_mult %>% 
  ggplot(aes(x = HCRmult, y = risk_collapse, group = paper, 
             linetype = as.factor(cluster), colour = as.factor(cluster))) +
  theme_paper + 
  geom_line(size = 0.2, show.legend = FALSE) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "", y = expression(collapse~risk))
p4 <- data_mult %>% 
  ggplot(aes(x = HCRmult, y = risk_blim, group = paper, 
             linetype = as.factor(cluster), colour = as.factor(cluster))) +
  theme_paper + 
  geom_line(size = 0.2, show.legend = FALSE) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "", y = expression(italic(B[lim])~risk)) +
  theme(strip.text = element_blank())
p5 <- data_mult %>% 
  ggplot(aes(x = HCRmult, y = yield_rel, group = paper, 
             linetype = as.factor(cluster), colour = as.factor(cluster))) +
  theme_paper + 
  geom_line(size = 0.2, show.legend = FALSE) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "catch rule multiplier", 
       y = expression(catch/MSY)) +
  theme(strip.text = element_blank())
p6 <- data_mult %>% 
  ggplot(aes(x = HCRmult, y = iav, group = paper, 
             linetype = as.factor(cluster), colour = as.factor(cluster))) +
  theme_paper + 
  geom_line(size = 0.2) +
  scale_linetype("cluster") + scale_colour_discrete("cluster") +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "", y = "ICV")

### extract legend
legend <- get_legend(p6)
### combine plots
plot_grid(plot_grid(p1, p2, p3, p4, p5, p6 + theme(legend.position = "none"), 
                    nrow = 2, align = "hv"),
          legend, ncol = 2, rel_widths = c(3, 0.2))

ggsave(filename = "output/plots/paper_revision/multiplier.png",
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper_revision/multiplier.jpeg", quality = 100,
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper_revision/multiplier.pdf",
       width = 17, height = 13, units = "cm", dpi = 600)



### ------------------------------------------------------------------------ ###
### catch constraints ####
### ------------------------------------------------------------------------ ###


stats_df <- readRDS(paste0(path_default, "plot_data_constraints.rds"))

### upper constraints
p1 <- stats_df %>% filter(upper < Inf &
                            lower == 0) %>%
  ggplot(aes(x = upper, y = ssb_rel, group = stock, 
             linetype = as.factor(cluster), colour = as.factor(cluster))) +
  theme_paper + 
  geom_line(size = 0.2, show.legend = FALSE) +
  geom_line(data = stats_df %>% filter(upper %in% c(1.5, Inf) &
                                         lower == 0),
            size = 0.05, show.legend = FALSE) +
  geom_point(data = stats_df %>% filter(upper %in% c(Inf) &
                                          lower == 0),
            size = 0.8, show.legend = FALSE) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "", y = expression(italic(SSB/B[MSY]))) + ylim(0, NA) +
  theme(strip.text = element_blank())
p2 <- stats_df %>% filter(upper < Inf &
                            lower == 0) %>%
  ggplot(aes(x = upper, y = risk_blim, group = stock, 
             linetype = as.factor(cluster), colour = as.factor(cluster))) +
  theme_paper + 
  geom_line(size = 0.2, show.legend = FALSE) +
  geom_line(data = stats_df %>% filter(upper %in% c(1.5, Inf) &
                                         lower == 0),
            size = 0.05, show.legend = FALSE) +
  geom_point(data = stats_df %>% filter(upper %in% c(Inf) &
                                          lower == 0),
             size = 0.8, show.legend = FALSE) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "upper catch constraint", y = expression(italic(B[lim])~risk)) +
  theme(strip.text = element_blank())
p3 <- stats_df %>% filter(upper < Inf &
                            lower == 0) %>%
  ggplot(aes(x = upper, y = yield_rel, group = stock, 
             linetype = as.factor(cluster), colour = as.factor(cluster))) +
  theme_paper + 
  geom_line(size = 0.2) +
  geom_line(data = stats_df %>% filter(upper %in% c(1.5, Inf) &
                                         lower == 0),
            size = 0.05, show.legend = FALSE) +
  geom_point(data = stats_df %>% filter(upper %in% c(Inf) &
                                          lower == 0),
             size = 0.8, show.legend = FALSE) +
  scale_colour_discrete("cluster") + scale_linetype("cluster") + 
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "", y = expression(catch/MSY)) +
  theme(legend.position = "right")

### lower constraints
p4 <- stats_df %>% filter(upper == 1.2) %>%
  ggplot(aes(x = lower, y = ssb_rel, group = stock, 
             linetype = as.factor(cluster), colour = as.factor(cluster))) +
  theme_paper + 
  geom_line(size = 0.2, show.legend = FALSE) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "", y = expression(italic(SSB/B[MSY]))) + ylim(0, NA) +
  theme(strip.text = element_blank()) + xlim(0, 1)
p5 <- stats_df %>% filter(upper == 1.2) %>%
  ggplot(aes(x = lower, y = risk_blim, group = stock, 
             linetype = as.factor(cluster), colour = as.factor(cluster))) +
  theme_paper + 
  geom_line(size = 0.2, show.legend = FALSE) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "lower catch constraint", y = expression(italic(B[lim])~risk)) +
  theme(strip.text = element_blank()) + xlim(0, 1)
p6 <- stats_df %>% filter(upper == 1.2) %>%
  ggplot(aes(x = lower, y = yield_rel, group = stock, 
             linetype = as.factor(cluster), colour = as.factor(cluster))) +
  theme_paper + 
  geom_line(size = 0.2) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "", y = expression(catch/MSY)) +
  theme(legend.position = "none") + xlim(0, 1)

### extract legend
legend <- get_legend(p3)
### combine plots
p_upper <- plot_grid(p1, p2, p3 + theme(legend.position = "none"), 
                     nrow = 1, align = "hv")
p_lower <- plot_grid(p4, p5, p6 + theme(legend.position = "none"),
                    nrow = 1, align = "hv")
plot_grid(plot_grid(p_upper, p_lower, ncol = 1, labels = c("A", "B"), 
                    label_fontfamily = "serif"),
          legend, rel_widths = c(3, 0.2))


ggsave(filename = "output/plots/paper_revision/constraints.png",
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper_revision/constraints.jpeg", quality = 100,
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper_revision/constraints.pdf",
       width = 17, height = 13, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### contribution of catch rule components ####
### ------------------------------------------------------------------------ ###


df_comps_perc <- readRDS(paste0(path_default, "plot_data_components.rds"))

### plot
p <- df_comps_perc %>%
  filter(component %in% c("r", "f", "b", "prod")) %>%
  ggplot() +
  geom_ribbon(aes(x = year, ymin = `5%`, ymax = `95%`, fill = "90%")) +
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`, fill = "50%")) +
  scale_fill_manual("confidence\ninterval", 
                    values = c("50%" = "grey50", "90%" = "grey80")) +
  geom_line(aes(x = year, y = `50%`, linetype = "median"), size = 0.2) +
  geom_hline(yintercept = 1, linetype = "dotted", size = 0.2) +
  scale_linetype("") +
  facet_grid(label ~ stk_label2, 
             scales = "fixed", labeller = label_parsed) +
  theme_paper +
  labs(x = "year", y = "component value") +
  ylim(0, NA)# +
  #geom_line(data = df_comps_perc %>% filter(percentile == "SSB/Bmsy"),
  #          aes(x = year, y = `50%`), size = 0.2,
  #          colour = "red", linetype = "dashed")
p

ggsave(filename = "output/plots/paper_revision/components.png",
       width = 8.5, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper_revision/components.jpeg", quality = 100,
       width = 8.5, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper_revision/components.pdf",
       width = 8.5, height = 10, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### data timing ####
### ------------------------------------------------------------------------ ###

res <- readRDS(paste0(path_default, "plot_data_timing.rds"))
### plot SSB 
res %>% filter(stock %in% c("ang", "pol", "ane", "her")) %>%
  ggplot(aes(x = year, y = data, linetype = timing, colour = as.factor(TAC))) +
  geom_line(size = 0.2) +
  theme_paper +
  facet_wrap(~ label2, labeller = label_parsed) +
  scale_colour_manual("TAC period", values = c("grey", "black")) +
  scale_linetype_discrete("relative timing\nof catch and\nindex") +
  labs(y = expression(italic(SSB/B[MSY])))

ggsave(filename = "output/plots/paper_revision/timing.png",
       width = 8.5, height = 6, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper_revision/timing.jpeg", quality = 100,
       width = 8.5, height = 6, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper_revision/timing.pdf",
       width = 8.5, height = 6, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### perfect information/knowledge ####
### ------------------------------------------------------------------------ ###

res <- readRDS(paste0(path_default, "plot_data_perfect.rds"))
### plot  
res %>% 
  filter(qname %in% c("rjc2", "meg", "lin", "ang", "had", "syc2",
                      "whg", "her"),
         year >= 100) %>% 
  ggplot(aes(x = year - 100, y = data, linetype = scenario)) +
  geom_hline(yintercept = 1, linetype = "solid", size = 0.2, colour = "grey") +
  geom_line(size = 0.2) +
  theme_paper +
  facet_wrap(~ label2, labeller = label_parsed, ncol = 2) +
  labs(x = "year", y = expression(SSB/B[MSY]), 
       colour = "scenario")

ggsave(filename = "output/plots/paper_revision/perfect.png",
       width = 8.5, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper_revision/perfect.jpeg", quality = 100,
       width = 8.5, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper_revision/perfect.pdf",
       width = 8.5, height = 10, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### some replicates ####
### ------------------------------------------------------------------------ ###

tmp <- readRDS(paste0(path_default, "plot_data_replicates.rds"))

tmp %>%
  ggplot(aes(x = year - 100)) +
  geom_col(aes(y = catch, fill = "catch")) +
  scale_fill_manual("", values = c("grey60", "green")) +
  geom_line(aes(y = ssb, colour = "SSB"), size = 0.2) + 
  scale_colour_manual("", values = c("black", "green")) +
  facet_grid(iter ~ label2, labeller = label_parsed) +
  #scale_colour_manual(values = c("#00BFC4", "#F8766D")) +
  theme_bw() +
  theme_paper +
  coord_cartesian(xlim = c(0, 20)) +
  labs(x = "year",
       y = expression(paste("relative value  ",
                            "(catch/MSY, SSB/",
                            ~B[MSY], ")")))# +
  #theme(strip.text.y = element_blank())

ggsave(filename = "output/plots/paper_revision/replicates.png",
       width = 8.5, height = 12, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper_revision/replicates.jpeg", quality = 100,
       width = 8.5, height = 12, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper_revision/replicates.pdf",
       width = 8.5, height = 12, units = "cm", dpi = 600)

