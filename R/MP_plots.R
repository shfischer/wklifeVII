library(FLCore)
library(ggplotFL)
library(tidyverse)
library(foreach)
library(cowplot)
library(Cairo)
library(doParallel)

cl <- makeCluster(10)
registerDoParallel(cl)

### path to final scenario
path_default <- "output/observation_error/new_baseline/"

### scenario defintions
scn_df <- read.csv("MP_scenarios.csv")

### load short/paper names
names_short <- read.csv("input/names_short.csv")
### add names and order by scenario
scn_df <- merge(scn_df, names_short)
scn_df <- scn_df[order(scn_df$scenario), ]

### load lhist
lhist <- readRDS("input/lhist_extended.rds")
### reference points
refpts <- readRDS("input/refpts.rds")

### theme for plotting
theme_paper <- theme_bw(base_size = 8, base_family = "serif")

### ------------------------------------------------------------------------ ###
### default catch rule 3.2.1  & fishing history ####
### ------------------------------------------------------------------------ ###

### find scenarios
### catch rule 3.2.1 default
scns <- scn_df %>% filter(scenario %in% 6515:6572)

### load default results
qnts <- foreach(scenario = scns$scenario, stock = as.character(scns$paper),
                fhist = scns$fhist,
                refpts_i = refpts[as.character(scns$stock)]) %do% {
  #browser()
  ### load
  qts_tmp <- readRDS(paste0("output/perfect_knowledge/combined/corrected/",
                            scenario, ".rds"))
  ### keep SSB, catch, fbar
  qts_tmp <- qts_tmp[c("ssb", "catch", "fbar")]
  ### calculate relative values
  qts_tmp$catch_rel <- qts_tmp$catch / c(refpts_i["msy", "yield"])
  qts_tmp$fbar_rel <- qts_tmp$fbar / c(refpts_i["msy", "harvest"])
  qts_tmp$ssb_rel <- qts_tmp$ssb / c(refpts_i["msy", "ssb"])
  ### calculate median
  qts_tmp <- lapply(qts_tmp, iterMedians)
  ### coerce into data frame
  qts_tmp <- as.data.frame(FLQuants(qts_tmp))
  qts_tmp$stock = stock
  qts_tmp$fhist = fhist
  return(qts_tmp)
  
}
qnts <- do.call(rbind, qnts)

### format for plotting
qnts_df <- qnts %>% spread(qname, data) %>%
  mutate(year = year - 100,
         group = ifelse(year < 1, "history", "projection"))
qnts_df <- qnts_df %>% bind_rows(qnts_df %>% filter(year == 0) %>%
                             mutate(group = "projection"))
qnts_df <- qnts_df %>% 
  mutate(group = factor(group, levels = c("history", "projection")))

# qnts_df %>%
#   ggplot(aes(x = year, y = ssb, group = stock)) +
#   geom_line() +
#   facet_wrap(~ fhist) +
#   theme_bw() +
#   geom_vline(xintercept = 0.5, linetype = "dashed")

### plot elements
plot_ssb <- qnts_df %>%
  ggplot(aes(x = year, y = ssb, group = stock)) +
  geom_line(size = 0.15) +
  facet_grid(fhist ~ group, scales = "free_x", space = "free_x") +
  theme_paper +
  scale_x_continuous(breaks = c(-25, 0, 25, 50, 75, 100), expand = c(0, 0)) +
  labs(y = "SSB") +
  theme(strip.text.y = element_blank(), 
        panel.spacing.x = unit(0, units = "cm"),
        plot.margin = unit(x = c(4, 6, 4, 4), units = "pt"))
plot_ssb_rel <- qnts_df %>%
  ggplot(aes(x = year, y = ssb_rel, group = stock)) +
  geom_line(size = 0.15) +
  facet_grid(fhist ~ group, scales = "free_x", space = "free_x") +
  theme_paper +
  scale_x_continuous(breaks = c(-25, 0, 25, 50, 75, 100), expand = c(0, 0)) +
  labs(y = expression(italic(SSB/B[MSY]))) +
  theme(strip.text.y = element_blank(), panel.spacing.x = unit(0, units = "cm"),
        plot.margin = unit(x = c(4, 6, 4, 4), units = "pt"))
plot_fbar <- qnts_df %>%
  ggplot(aes(x = year, y = fbar, group = stock)) +
  geom_line(size = 0.15) +
  facet_grid(fhist ~ group, scales = "free_x", space = "free_x") +
  theme_paper +
  scale_x_continuous(breaks = c(-25, 0, 25, 50, 75, 100), expand = c(0, 0)) +
  labs(y = "fishing mortality") +
  theme(strip.text.y = element_blank(), panel.spacing.x = unit(0, units = "cm"),
        plot.margin = unit(x = c(4, 6, 4, 4), units = "pt"))
plot_fbar_rel <- qnts_df %>%
  ggplot(aes(x = year, y = fbar_rel, group = stock)) +
  geom_line(size = 0.15) +
  facet_grid(fhist ~ group, scales = "free_x", space = "free_x") +
  theme_paper +
  scale_x_continuous(breaks = c(-25, 0, 25, 50, 75, 100), expand = c(0, 0)) +
  labs(y = expression(italic(F/F[MSY]))) +
  theme(strip.text.y = element_blank(), panel.spacing.x = unit(0, units = "cm"),
        plot.margin = unit(x = c(4, 6, 4, 4), units = "pt"))
plot_catch <- qnts_df %>%
  ggplot(aes(x = year, y = catch, group = stock)) +
  geom_line(size = 0.15) +
  facet_grid(fhist ~ group, scales = "free_x", space = "free_x") +
  theme_paper +
  scale_x_continuous(breaks = c(-25, 0, 25, 50, 75, 100), expand = c(0, 0)) +
  labs(y = "catch") +
  theme(panel.spacing.x = unit(0, units = "cm"))
plot_catch_rel <- qnts_df %>%
  ggplot(aes(x = year, y = catch_rel, group = stock)) +
  geom_line(size = 0.15) +
  facet_grid(fhist ~ group, scales = "free_x", space = "free_x") +
  theme_paper +
  scale_x_continuous(breaks = c(-25, 0, 25, 50, 75, 100), expand = c(0, 0)) +
  labs(y = "catch/MSY") +
  theme(panel.spacing.x = unit(0, units = "cm"))

### combine the three plots
plot_grid(plot_ssb, plot_fbar, plot_catch, rel_widths = c(1, 1, 1.1),
          nrow = 1)
### and save
ggsave(filename = "output/perfect_knowledge/plots/paper/trajectories.png",
       width = 17, height = 7, units = "cm", dpi = 600, type = "cairo")

ggsave(filename = "output/perfect_knowledge/plots/paper/trajectories.pdf",
       width = 17, height = 7, units = "cm", dpi = 600)

### same for relative values
### combine the three plots
plot_grid(plot_ssb_rel, plot_fbar_rel, plot_catch_rel, 
          rel_widths = c(1, 1, 1.1), nrow = 1)
### and save
ggsave(filename = "output/perfect_knowledge/plots/paper/trajectories_rel.png",
       width = 17, height = 7, units = "cm", dpi = 600, type = "cairo")

ggsave(filename = "output/perfect_knowledge/plots/paper/trajectories_rel.pdf",
       width = 17, height = 7, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### penalized regression ####
### ------------------------------------------------------------------------ ###
library(glmnet)

### load data
glmnet_res <- readRDS(paste0(path_default, "corrected/one-way/glmnet.rds"))

### plots
p1 <- ggplot(data = glmnet_res$data,
             aes(x = k, y = f_rel_new)) +
  geom_point(size = 0.3) + theme_paper + 
  lims(x = c(0, NA), y = c(0, NA)) +
  geom_abline(data = data.frame(
    intercept = c(glmnet_res$glmnet$f_rel_new["(Intercept)",],
                  lms[lms$var == "f_rel_new", "Intercept"]),
    slope = c(glmnet_res$glmnet$f_rel_new["k",],
              lms[lms$var == "f_rel_new", "k"]),
    model = c("lasso regression", "linear regression")),
    aes(intercept = intercept, slope = slope, linetype = model), size = 0.3, 
    show.legend = FALSE) +
  labs(y = expression(italic(F/F[MSY])), x = "") +
  theme(axis.title.x = element_blank())
p2 <- ggplot(data = glmnet_res$data,
             aes(x = k, y = ssb_rel_new)) +
  geom_point(size = 0.3) + theme_paper + 
  lims(x = c(0, NA), y = c(0, NA)) +
  geom_abline(data = data.frame(
    intercept = c(glmnet_res$glmnet$ssb_rel_new["(Intercept)",],
                  lms[lms$var == "ssb_rel_new", "Intercept"]),
    slope = c(glmnet_res$glmnet$ssb_rel_new["k",],
              lms[lms$var == "ssb_rel_new", "k"]),
    model = c("lasso regression", "linear regression")),
    aes(intercept = intercept, slope = slope, linetype = model), size = 0.3, 
    show.legend = FALSE) +
  labs(y = expression(italic(SSB/B[MSY])), x = "") +
  theme(axis.title.x = element_blank())
p3 <- ggplot(data = glmnet_res$data,
             aes(x = k, y = collapse_total_new)) +
  geom_point(size = 0.3) + theme_paper + 
  lims(x = c(0, NA), y = c(0, NA)) +
  geom_abline(data = data.frame(
    intercept = c(glmnet_res$glmnet$collapse_total_new["(Intercept)",],
                  lms[lms$var == "collapse_total_new", "Intercept"]),
    slope = c(glmnet_res$glmnet$collapse_total_new["k",],
              lms[lms$var == "collapse_total_new", "k"]),
    model = c("lasso regression", "linear regression")),
    aes(intercept = intercept, slope = slope, linetype = model), size = 0.3, 
    show.legend = FALSE) +
  labs(y = expression(collapse~risk), x = "") +
  theme(axis.title.x = element_blank())
p4 <- ggplot(data = glmnet_res$data,
             aes(x = k, y = ssb_below_blim_new)) +
  geom_point(size = 0.3) + theme_paper + 
  lims(x = c(0, NA), y = c(0, NA)) +
  geom_abline(data = data.frame(
    intercept = c(glmnet_res$glmnet$ssb_below_blim_new["(Intercept)",],
                  lms[lms$var == "ssb_below_blim_new", "Intercept"]),
    slope = c(glmnet_res$glmnet$ssb_below_blim_new["k",],
              lms[lms$var == "ssb_below_blim_new", "k"]),
    model = c("lasso regression", "linear regression")),
    aes(intercept = intercept, slope = slope, linetype = model), size = 0.3, 
    show.legend = FALSE) +
  labs(y = expression(italic(B[lim])~risk), x = "") +
  theme(axis.title.x = element_blank())
p5 <- ggplot(data = glmnet_res$data[, ],
             aes(x = k, y = yield_rel_MSY_new)) +
  geom_point(size = 0.3) + theme_paper + 
  lims(x = c(0, NA), y = c(0, NA)) +
  geom_abline(data = data.frame(
    intercept = c(glmnet_res$glmnet$yield_rel_MSY_new["(Intercept)",],
                  lms[lms$var == "yield_rel_MSY_new", "Intercept"]),
    slope = c(glmnet_res$glmnet$yield_rel_MSY_new["k",],
              lms[lms$var == "yield_rel_MSY_new", "k"]),
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
                                         lms[lms$var == "iav", "Intercept"]),
                           slope = c(glmnet_res$glmnet$iav["k",],
                                          lms[lms$var == "iav", "k"]),
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
ggsave(filename = "output/plots/paper_revision/glmnet.pdf",
       width = 8.5, height = 12, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### clustering ####
### ------------------------------------------------------------------------ ###
library(dtwclust)
library(ggdendro)

cluster <- readRDS(paste0(path, "corrected/one-way/cluster_plot_data.rds"))

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
        #axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black", angle = 90, hjust = 0.5),
        axis.text.x = element_text(size = 8),
        text = element_text(size = 8, family = "serif")) +
  labs(x = "stocks", y = "DTW distance")

plot_grid(p3, plot_grid(p1, p2, ncol = 2, rel_widths = c(1.7, 1), 
                        labels = c("B", "C"), label_fontfamily = "serif"),
          nrow = 2, rel_heights = c(1, 2), labels = c("A", ""), 
          label_fontfamily = "serif")

ggsave(filename = "output/plots/paper_revision/clustering.png",
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = "output/plots/paper_revision/clustering.pdf",
       width = 17, height = 13, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### multiplier ####
### ------------------------------------------------------------------------ ###

### load full stats
stats <- readRDS(file = "output/stats_scn_new.RDS")
### load cluster allocations for colour-coding
cl_alloc_ow <- readRDS("output/cluster_allocations_one_way.rds")
### merge
stats <- merge(stats, cl_alloc_ow[, c("4", "stock")])
### subset to scenarios for multiplier
### subset to 3.2.1 & one-way & multiplier
stats_df <- stats %>%
  filter(scenario %in% 4485:6804 & b_z == 1 &
           `4` != 0)

p1 <- stats_df %>% 
  ggplot(aes(x = HCRmult, y = f_rel, group = stock, 
             linetype = as.factor(`4`), colour = as.factor(`4`))) +
  theme_paper + 
  geom_line(size = 0.2, show.legend = FALSE) +
  #geom_line() +
  # scale_linetype_manual(values = c(`2` = "dotted", `3` = "dashed", 
  #                                  `4` = "solid", `1` = "solid")) +
  # scale_colour_manual(values = c(`2` = "black", `3` = "grey50", 
  #                                `4` = "black", `1` = "grey50")) +
  ylim(0, 1.5) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "", y = expression(italic(F/F[MSY]))) +
  theme(strip.text = element_blank())
p2 <- stats_df %>% 
  ggplot(aes(x = HCRmult, y = ssb_rel, group = stock, 
             linetype = as.factor(`4`), colour = as.factor(`4`))) +
  theme_paper + 
  geom_line(size = 0.2, show.legend = FALSE) +
  # scale_linetype_manual(values = c(`2` = "dotted", `3` = "dashed", 
  #                                  `4` = "solid", `1` = "solid")) +
  # scale_colour_manual(values = c(`2` = "black", `3` = "grey50", 
  #                                `4` = "black", `1` = "grey50")) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "", y = expression(italic(SSB/B[MSY]))) +
  theme(strip.text = element_blank())
p3 <- stats_df %>% 
  ggplot(aes(x = HCRmult, y = collapse_total, group = stock, 
             linetype = as.factor(`4`), colour = as.factor(`4`))) +
  theme_paper + 
  geom_line(size = 0.2, show.legend = FALSE) +
  # scale_linetype_manual(values = c(`2` = "dotted", `3` = "dashed", 
  #                                  `4` = "solid", `1` = "solid")) +
  # scale_colour_manual(values = c(`2` = "black", `3` = "grey50", 
  #                                `4` = "black", `1` = "grey50")) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "", y = expression(collapse~risk))
p4 <- stats_df %>% 
  ggplot(aes(x = HCRmult, y = ssb_below_blim_total, group = stock, 
             linetype = as.factor(`4`), colour = as.factor(`4`))) +
  theme_paper + 
  geom_line(size = 0.2, show.legend = FALSE) +
  # scale_linetype_manual(values = c(`2` = "dotted", `3` = "dashed", 
  #                                  `4` = "solid", `1` = "solid")) +
  # scale_colour_manual(values = c(`2` = "black", `3` = "grey50", 
  #                                `4` = "black", `1` = "grey50")) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "", y = expression(italic(B[lim])~risk)) +
  theme(strip.text = element_blank())
p5 <- stats_df %>% 
  ggplot(aes(x = HCRmult, y = catch_MSY_prop, group = stock, 
             linetype = as.factor(`4`), colour = as.factor(`4`))) +
  theme_paper + 
  geom_line(size = 0.2, show.legend = FALSE) +
  # scale_linetype_manual(values = c(`2` = "dotted", `3` = "dashed", 
  #                                  `4` = "solid", `1` = "solid")) +
  # scale_colour_manual(values = c(`2` = "black", `3` = "grey50", 
  #                                `4` = "black", `1` = "grey50")) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "catch rule multiplier", 
       y = expression(catch/MSY)) +
  theme(strip.text = element_blank())
p6 <- stats_df %>% 
  ggplot(aes(x = HCRmult, y = iav, group = stock, 
             linetype = as.factor(`4`), colour = as.factor(`4`))) +
  theme_paper + 
  geom_line(size = 0.2) +
  # scale_linetype_manual(values = c(`2` = "dotted", `3` = "dashed", 
  #                                  `4` = "solid", `1` = "solid"), "cluster") +
  # scale_colour_manual(values = c(`2` = "black", `3` = "grey50", 
  #                                `4` = "black", `1` = "grey50"), "cluster") +
  scale_linetype("cluster") + scale_colour_discrete("cluster") +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "", y = "ICV")

### extract legend
legend <- get_legend(p6)
### combine plots
plot_grid(plot_grid(p1, p2, p3, p4, p5, p6 + theme(legend.position = "none"), 
                    nrow = 2, align = "hv"),
          legend, ncol = 2, rel_widths = c(3, 0.2))

ggsave(filename = paste0("output/perfect_knowledge/plots/paper/",
                         "multiplier_colour.png"),
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = paste0("output/perfect_knowledge/plots/paper/",
                         "multiplier_colour.pdf"),
       width = 17, height = 13, units = "cm", dpi = 600)



### ------------------------------------------------------------------------ ###
### catch constraints ####
### ------------------------------------------------------------------------ ###

### subset to scenarios for catch constraints
### subset to 3.2.1 & one-way & multiplier
stats_df <- stats %>%
  filter(scenario %in% 1237:4484)

### upper constraints
p1 <- stats_df %>% filter(upper_constraint < Inf &
                            lower_constraint == 0) %>%
  ggplot(aes(x = upper_constraint, y = ssb_rel, group = stock, 
             linetype = as.factor(`4`), colour = as.factor(`4`))) +
  theme_paper + 
  geom_line(size = 0.2, show.legend = FALSE) +
  geom_line(data = stats_df %>% filter(upper_constraint %in% c(1.5, Inf) &
                                         lower_constraint == 0),
            size = 0.05, show.legend = FALSE) +
  geom_point(data = stats_df %>% filter(upper_constraint %in% c(Inf) &
                                          lower_constraint == 0),
            size = 0.8, show.legend = FALSE) +
  # scale_linetype_manual(values = c(`2` = "dotted", `3` = "dashed", 
  #                                  `4` = "solid", `1` = "solid")) +
  # scale_colour_manual(values = c(`2` = "black", `3` = "grey50", 
  #                                `4` = "black", `1` = "grey50")) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "", y = expression(italic(SSB/B[MSY]))) + ylim(0, NA) +
  theme(strip.text = element_blank())
p2 <- stats_df %>% filter(upper_constraint < Inf &
                            lower_constraint == 0) %>%
  ggplot(aes(x = upper_constraint, y = ssb_below_blim_total, group = stock, 
             linetype = as.factor(`4`), colour = as.factor(`4`))) +
  theme_paper + 
  geom_line(size = 0.2, show.legend = FALSE) +
  geom_line(data = stats_df %>% filter(upper_constraint %in% c(1.5, Inf) &
                                         lower_constraint == 0),
            size = 0.05, show.legend = FALSE) +
  geom_point(data = stats_df %>% filter(upper_constraint %in% c(Inf) &
                                          lower_constraint == 0),
             size = 0.8, show.legend = FALSE) +
  # scale_linetype_manual(values = c(`2` = "dotted", `3` = "dashed", 
  #                                  `4` = "solid", `1` = "solid")) +
  # scale_colour_manual(values = c(`2` = "black", `3` = "grey50", 
  #                                `4` = "black", `1` = "grey50")) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "upper catch constraint", y = expression(italic(B[lim])~risk)) +
  theme(strip.text = element_blank())
p3 <- stats_df %>% filter(upper_constraint < Inf &
                            lower_constraint == 0) %>%
  ggplot(aes(x = upper_constraint, y = catch_MSY_prop, group = stock, 
             linetype = as.factor(`4`), colour = as.factor(`4`))) +
  theme_paper + 
  geom_line(size = 0.2) +
  geom_line(data = stats_df %>% filter(upper_constraint %in% c(1.5, Inf) &
                                         lower_constraint == 0),
            size = 0.05, show.legend = FALSE) +
  geom_point(data = stats_df %>% filter(upper_constraint %in% c(Inf) &
                                          lower_constraint == 0),
             size = 0.8, show.legend = FALSE) +
  # scale_linetype_manual(values = c(`2` = "dotted", `3` = "dashed", 
  #                                  `4` = "solid", `1` = "solid"), "cluster") +
  # scale_colour_manual(values = c(`2` = "black", `3` = "grey50", 
  #                                `4` = "black", `1` = "grey50"), "cluster") +
  scale_colour_discrete("cluster") + scale_linetype("cluster") + 
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "", y = expression(catch/MSY)) +
  theme(legend.position = "right")

### lower constraints
p4 <- stats_df %>% filter(upper_constraint == 1.2) %>%
  ggplot(aes(x = lower_constraint, y = ssb_rel, group = stock, 
             linetype = as.factor(`4`), colour = as.factor(`4`))) +
  theme_paper + 
  geom_line(size = 0.2, show.legend = FALSE) +
  # scale_linetype_manual(values = c(`2` = "dotted", `3` = "dashed", 
  #                                  `4` = "solid", `1` = "solid")) +
  # scale_colour_manual(values = c(`2` = "black", `3` = "grey50", 
  #                                `4` = "black", `1` = "grey50")) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "", y = expression(italic(SSB/B[MSY]))) + ylim(0, NA) +
  theme(strip.text = element_blank()) + xlim(0, 1)
p5 <- stats_df %>% filter(upper_constraint == 1.2) %>%
  ggplot(aes(x = lower_constraint, y = ssb_below_blim_total, group = stock, 
             linetype = as.factor(`4`), colour = as.factor(`4`))) +
  theme_paper + 
  geom_line(size = 0.2, show.legend = FALSE) +
  # scale_linetype_manual(values = c(`2` = "dotted", `3` = "dashed", 
  #                                  `4` = "solid", `1` = "solid")) +
  # scale_colour_manual(values = c(`2` = "black", `3` = "grey50", 
  #                                `4` = "black", `1` = "grey50")) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "lower catch constraint", y = expression(italic(B[lim])~risk)) +
  theme(strip.text = element_blank()) + xlim(0, 1)
p6 <- stats_df %>% filter(upper_constraint == 1.2) %>%
  ggplot(aes(x = lower_constraint, y = catch_MSY_prop, group = stock, 
             linetype = as.factor(`4`), colour = as.factor(`4`))) +
  theme_paper + 
  geom_line(size = 0.2) +
  # scale_linetype_manual(values = c(`2` = "dotted", `3` = "dashed", 
  #                                  `4` = "solid", `1` = "solid"), "cluster") +
  # scale_colour_manual(values = c(`2` = "black", `3` = "grey50", 
  #                                `4` = "black", `1` = "grey50"), "cluster") +
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

ggsave(filename = paste0("output/perfect_knowledge/plots/paper/",
                         "constraints_colour.png"),
       width = 17, height = 13, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = paste0("output/perfect_knowledge/plots/paper/",
                         "constraints_colour.pdf"),
       width = 17, height = 13, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### contribution of catch rule components ####
### ------------------------------------------------------------------------ ###

### go through scenarios
df_comps <- foreach(scenario = scns$scenario[1:2], paper = scns$paper[1:2],
                    refpts_i = refpts[scns$stock[1:2]],
                    .combine = rbind, .packages = "FLCore") %do% {
  #browser()
  ### load stock
  stk_tmp <- readRDS(paste0("output/perfect_knowledge/combined/", scenario, 
                           ".rds"))
  
  ### load quants
  quants <- readRDS(paste0("output/perfect_knowledge/combined/corrected/",
                          scenario, ".rds"))
  
  ### get tracking object
  comps <- stk_tmp@tracking[c("HCR3.2.1r", "HCR3.2.1f", "HCR3.2.1b"),
                           ac(seq(from = 99, to = 200, by = 2))]
  # comps <- 1
  #   apply(comps, c(2:6), prod)
  
  ### find collapses and remove values
  for (year_i in dimnames(comps)$year) {
   for (iter_i in dimnames(comps)$iter) {
     if (isTRUE(quants$ssb[, year_i,,,, iter_i] == 0)) {
       comps[, year_i,,,, iter_i] <- NA
     }
   }
  }

  ### coerce into df for plotting
  df_comps <- as.data.frame(comps)
  df_comps$metric <- substr(x = df_comps$metric, start = 9, stop = 9)
  
  ### add relative SSB/Bmsy
  ssb_rel <- (quants$ssb / c(refpts_i["msy", "ssb"]))[, ac(99:198)]
  df_ssb <- as.data.frame(ssb_rel)
  names(df_ssb)[1] <- "metric"
  df_ssb$metric <- "SSB"
  ### median ssb
  ssb_med <- apply(ssb_rel, 2, median)
  df_ssb_med <- as.data.frame(ssb_med)
  names(df_ssb_med)[1] <- "metric"
  df_ssb_med$metric <- "SSB_median"
  
  df_res <- rbind(df_comps, df_ssb, df_ssb_med)
  ### stock names in paper format
  df_res$paper <- paper
  
  return(df_res)
  
}

### log transform
# df_comps$log_data <- log(df_comps$data)
### add product
df_comps <- df_comps %>% spread(metric, data)
df_comps$prod <- with(df_comps, r * f * b)
### gather again
df_comps <- df_comps %>% 
  gather(key = "component", value = "data", r, f, b, prod, SSB, SSB_median)

# ### add sum
# df_comps <- df_comps %>% spread(metric, data)
# df_comps$sum <- with(df_comps, r + f + b)
# ### gather again
# df_comps <- df_comps %>% 
#   gather(key = "component", value = "data", r, f, b, sum)

### change years to start from 1
df_comps <- df_comps %>%
  mutate(year = year - 98) %>%
  left_join(names_short %>% select(paper, stock)) %>%
  left_join(lhist %>% select(k, stock))
df_comps <- df_comps %>% filter(!is.na(data))
df_comps <- df_comps %>%
  mutate(percentile = ifelse(component != "SSB_median", "median", 
                             "SSB/Bmsy"))
df_comps_perc <- df_comps %>%
  group_by(paper, year, percentile, component) %>%
  summarise(`5%` = quantile(data, probs = 0.05, na.rm = TRUE),
            `25%` = quantile(data, probs = 0.25, na.rm = TRUE),
            `50%` = quantile(data, probs = 0.50, na.rm = TRUE),
            `75%` = quantile(data, probs = 0.75, na.rm = TRUE),
            `95%` = quantile(data, probs = 0.95, na.rm = TRUE))
df_comps_perc <- df_comps_perc %>%
  mutate(label = factor(component, levels = c("r", "f", "b", "prod", "SSB",
                                                  "SSB_median"),
                            labels = c("italic(r)",
                                       "italic(f)",
                                       "italic(b)",
                                       "italic(r~f~b)",
                                       "SSB_rel", 
                                       "italic(r~f~b)")))

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
  facet_grid(label ~ paper, 
             scales = "fixed", labeller = label_parsed) +
  theme_paper +
  labs(x = "year", y = "component value") +
  ylim(0, NA)# +
  #geom_line(data = df_comps_perc %>% filter(percentile == "SSB/Bmsy"),
  #          aes(x = year, y = `50%`), size = 0.2,
  #          colour = "red", linetype = "dashed")
p

ggsave(filename = paste0("output/perfect_knowledge/plots/paper/",
                         "components.png"),
       width = 8.5, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = paste0("output/perfect_knowledge/plots/paper/",
                         "components.pdf"),
       width = 8.5, height = 10, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### data timing ####
### ------------------------------------------------------------------------ ###

### subset to scenarios
scns <- scn_df %>% filter(scenario %in% c(6863:7210))

### load quants and calculate SSB/Bmsy
qnts <- foreach(scenario = scns$scenario, stock = scns$stock,
                  .packages = "FLCore", .export = "refpts") %dopar% {
  ### load SSB
  ssb_i <- readRDS(paste0("output/perfect_knowledge/combined/corrected/",
                          scenario, ".rds"))$ssb
  ### calculate relative SSB
  ssb_i <- ssb_i / c(refpts[[stock]]["msy", "ssb"])
  ### use median only
  ssb_i <- iterMedians(ssb_i)
  ### return
  return(ssb_i)
}
### format
names(qnts) <- scns$scenario
qnts <- FLQuants(qnts)
qnts_df <- as.data.frame(qnts)
qnts_df <- qnts_df %>%
  mutate(scenario = as.numeric(as.character(qname)),
         year = year - 100) %>%
  filter(year >= 1) %>%
  left_join(scns) %>%
  left_join(lhist %>% select(stock, k))
  

### plot SSB for all stocks
qnts_df %>% filter(fhist == "one-way" & 
                     paper %in% c("ang", "pol", "ane", "her")) %>%
  mutate(timing = paste(lst_catch, lst_idx),
         wrap = paste0("italic(k)==", k, "~", paper),
         TAC = ifelse(TAC == 1, "annual", "biennial"),
         timing = ifelse(timing == "0 1", "0 +1", timing)) %>%
  ggplot(aes(x = year, y = data, linetype = timing, colour = as.factor(TAC))) +
  geom_line(size = 0.2) +
  theme_paper +
  facet_wrap(~ wrap, labeller = label_parsed) +
  scale_colour_manual("TAC period", values = c("grey", "black")) +
  scale_linetype_discrete("relative timing\nof catch and\nindex") +
  labs(y = expression(italic(SSB/B[MSY])))
ggsave(filename = paste0("output/perfect_knowledge/plots/paper/",
                         "timing.png"),
       width = 8.5, height = 6, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = paste0("output/perfect_knowledge/plots/paper/",
                         "timing.pdf"),
       width = 8.5, height = 6, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### perfect information/knowledge ####
### ------------------------------------------------------------------------ ###

### 6515:6572 default
### 8719:8776 perfect information

### load scenarios
scns <- scn_df %>% filter(scenario %in% c(6515:6572, 8719:8776) &
                            fhist == "one-way")

### load quants and calculate SSB/Bmsy
qnts <- foreach(scenario = scns$scenario, stock = scns$stock,
                .packages = "FLCore", .export = "refpts") %dopar% {
  ### load SSB
  ssb_i <- readRDS(paste0("output/perfect_knowledge/combined/corrected/",
                          scenario, ".rds"))$ssb
  ### calculate relative SSB
  ssb_i <- ssb_i / c(refpts[[stock]]["msy", "ssb"])
  ### use median only
  ssb_i <- iterMedians(ssb_i)
  ### return
  return(ssb_i)
}
### format
names(qnts) <- scns$scenario
qnts <- FLQuants(qnts)
qnts_df <- as.data.frame(qnts)
qnts_df <- qnts_df %>%
  mutate(scenario = as.numeric(as.character(qname)),
         year = year - 100) %>%
  ### perfect info has only 50 years
  filter(year >= 1 & year <= 50) %>%
  left_join(scns) %>%
  left_join(lhist %>% select(stock, k)) %>%
  mutate(info = ifelse(is.na(OM_scn), "default", "perfect\ninformation"),
         wrap = paste0("italic(k)==", k, "~", paper)) %>%
  filter(paper %in% c("ang", "rjc", "meg", "had", "lin", "mut", "whg", "her"))


### plot SSB for all stocks
qnts_df %>%
  ggplot(aes(x = year, y = data, linetype = info)) +
  geom_line(size = 0.2) +
  geom_hline(yintercept = 1, linetype = "dotted", size = 0.2) +
  theme_paper +
  facet_wrap(~ wrap, labeller = label_parsed, ncol = 2) +
  scale_linetype("scenario") +
  labs(y = expression(SSB/B[MSY]))
ggsave(filename = paste0("output/perfect_knowledge/plots/paper/",
                         "perfect.png"),
       width = 8.5, height = 10, units = "cm", dpi = 600, type = "cairo")
ggsave(filename = paste0("output/perfect_knowledge/plots/paper/",
                         "perfect.pdf"),
       width = 8.5, height = 10, units = "cm", dpi = 600)
