library(FLCore)
library(ggplotFL)
library(Cairo)
library(tidyverse)

library(doParallel)
cl <- makeCluster(4)

### ------------------------------------------------------------------------ ###
### group ####
### ------------------------------------------------------------------------ ###

res_df <- readRDS("output/stats_scn_new.RDS") ### quants
refpts <- readRDS("input/refpts.rds")
pars <- readRDS("input/all_stocks_repfts_lhpar.rds")

### groups
pars$stock <- as.character(pars$stock)
g1 <- c("her-nis", "san-ns4", "sardina_pilchardus", "zeus_faber")
g2 <- c("ane-pore", "lem-nsea", "scophthalmus_rhombus", "whg-7e-k",
        "chelidonichtys_lucerna")
g3 <- c("anarchias_lupus", "arg-comb-ex5.", "had-iris", "mut-comb", "nep-2829",
        "ple-celt", "syc27.67", "syc27.8c", "tur-nsea",
        "spondyliosoma_cantharus")
g4 <- c("ang-78ab", "ang-ivvi", "ang-78ab_2", "lin-comb", "pol-nsea", 
        "rjc.27.afg", "rjc.27.347d", "sdv.27.nea")
g5 <- c("meg-4a6a")
g6 <- c("smn-con")

pars$group <- NA
pars$group <- unlist(lapply(pars$stock, function(x) {
  if (x %in% g1) return(1)
  else if (x %in% g2) return(2)
  else if (x %in% g3) return(3)
  else if (x %in% g4) return(4)
  else if (x %in% g5) return(5)
  else if (x %in% g6) return(6)
}))

### load SSBs
df_ssb <- foreach(scenario = 6515:6572, .packages = c("FLCore"), 
                  .export = c("res_df"), .errorhandling = "pass") %do% {
  ### load quants
  qts <- readRDS(paste0("output/perfect_knowledge/combined/3.2.1_quants/",
                        scenario, ".rds"))[["SSB"]]
  qts <- iterMedians(qts) ### median over all iterations
  ### coerce into data frame
  df_tmp <- as.data.frame(qts)
  ### add scenario definitions
  scn <- res_df[res_df$scenario == scenario, ]
  df_tmp <- cbind(df_tmp, fhist = scn$fhist, scenario = scenario,
                  stock = scn$stock, group = pars$group[pars$stock == scn$stock])
  return(df_tmp)
}
df_ssb <- do.call(rbind, df_ssb)
### plot
ggplot(data = df_ssb, aes(x = year, y = data, colour = as.factor(group),
                          group = as.factor(stock))) +
  geom_vline(xintercept = 100, colour = "grey") + geom_line() +
  theme_bw() +
  facet_wrap(~ fhist, scales = "free_y") +
  labs(y = "") +
  scale_colour_discrete("group")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/correlations",
                         "/default_ssb_all_stocks_groups.png"),
       width = 30, height = 18, units = "cm", dpi = 300, type = "cairo-png")
### plot each group separately
ggplot(data = df_ssb, aes(x = year, y = data, colour = as.factor(group),
                          group = as.factor(stock))) +
  geom_vline(xintercept = 100, colour = "grey") + geom_line() +
  theme_bw() +
  facet_grid(fhist ~ group) +
  labs(y = "SSB") +
  scale_colour_discrete("group")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/correlations",
                         "/default_ssb_all_stocks_groups2.png"),
       width = 20, height = 12, units = "cm", dpi = 100, type = "cairo-png")

### format pars for plotting
pars$stock <- factor(pars$stock, 
                     levels = pars$stock[order(pars$group,
                                               -(pars$K))])
pars$group <- as.factor(pars$group)
pars_df <- gather(pars, key = "parameter", value = "value",
                  c("L_inf", "K", "a50", "M_mat", "M", "MK", "t0", "a", "b", 
                    "max_age", "LFeFmsy", "LFeMK", "F_MSY", "L_c"))
pars_df$input <- "input"
pars_df$input[pars_df$parameter %in% c("max_age", "LFeFmsy", "LFeMK",
                                       "F_MSY", "L_c", "M", 
                                       "M_mat", "MK")] <- "derived"
pars_df$input <- factor(pars_df$input, levels = unique(pars_df$input))
pars_df$parameter <- factor(pars_df$parameter,
                            levels = unique(pars_df$parameter))
### plot
ggplot(pars_df, aes(x = stock, y = value, fill = group)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_wrap(~ input + parameter, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/correlations",
                         "/groups_pars.png"),
       width = 30, height = 18, units = "cm", dpi = 100, type = "cairo-png")



### ------------------------------------------------------------------------ ###
### stats vs. tuning: multiplier ####
### ------------------------------------------------------------------------ ###

### load full stats
stats <- readRDS(file = "output/stats_scn_new.RDS")
### load extended lhist
lhist <- readRDS("input/lhist_extended.rds")
### merge
stats_df <- merge(stats, lhist, all = TRUE)
### load cluster allocations
cl_alloc_ow <- readRDS("output/cluster_allocations_one_way.rds")
cl_alloc_rc <- readRDS("output/cluster_allocations_roller_coaster.rds")
### merge
stats_df <- merge(stats_df, cl_alloc_ow[, c("stock", "4")], all = TRUE)

### subset to 3.2.1 & one-way & multiplier
stats_df2 <- subset(stats_df, scenario %in% 4485:6804 &
                   b_z == 1 & fhist == "one-way")

ggplot(data = stats_df2 %>% 
         gather(key = "parameter", value = "value", 
                collapse_total, catch_MSY_prop), 
       aes(x = k, y = value, group = stock)) +
  theme_bw() + geom_line(position = position_dodge(width = 0.01)) + 
  geom_point(aes(colour = as.factor(HCRmult)),
             position = position_dodge(width = 0.01)) +
  facet_wrap(~ parameter, scales = "free_y")
### colour-code cluster
ggplot(data = stats_df2 %>% 
         gather(key = "parameter", value = "value", 
                collapse_total, catch_MSY_prop), 
       aes(x = HCRmult, y = value, group = stock, colour = as.factor(`4`))) +
  theme_bw() + geom_line() + geom_point() +
  facet_wrap(~ parameter, scales = "free_y")
### colour code k
ggplot(data = stats_df2 %>% 
         gather(key = "parameter", value = "value", 
                collapse_total, catch_MSY_prop), 
       aes(x = HCRmult, y = value, group = stock, colour = k)) +
  theme_bw() + geom_line() + geom_point() +
  facet_wrap(~ parameter, scales = "free_y")
### separate ks in separate plots
ggplot(data = stats_df2 %>% 
         gather(key = "parameter", value = "value", 
                ssb_rel, f_rel, collapse_iter_old, iav, catch_MSY_prop), 
       aes(x = HCRmult, y = value, group = stock)) +
  theme_bw() + geom_line() + geom_point() +
  facet_grid(parameter ~ k, scales = "free_y")
### separate by stock
ggplot(data = stats_df2 %>% 
         gather(key = "parameter", value = "value", 
                ssb_rel, f_rel, collapse_iter_old, iav, catch_MSY_prop), 
       aes(x = HCRmult, y = value, group = stock)) +
  theme_bw() + geom_line() + geom_point() +
  facet_grid(parameter ~ stock, scales = "free_y")

### plot for report
p1 <- ggplot(data = subset(stats_df, scenario %in% 4485:6804 & b_z == 1), 
             aes(x = HCRmult, y = f_rel, group = stock,
                 colour = as.factor(`4`))) +
  theme_bw() + geom_line(show.legend = FALSE) + ylim(0, 1.5) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = " ", y = "F/Fmsy") +
  theme(strip.text = element_blank())
p2 <- ggplot(data = subset(stats_df, scenario %in% 4485:6804 & b_z == 1), 
             aes(x = HCRmult, y = ssb_rel, group = stock,
                 colour = as.factor(`4`))) +
  theme_bw() + geom_line(show.legend = FALSE) + 
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = " ", y = "B/Bmsy") +
  theme(strip.text = element_blank())
p3 <- ggplot(data = subset(stats_df, scenario %in% 4485:6804 & b_z == 1), 
             aes(x = HCRmult, y = collapse_total, group = stock,
                 colour = as.factor(`4`))) +
  theme_bw() + geom_line(show.legend = FALSE) + 
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") + ylim(0, 1) +
  labs(x = " ", y = "collapse risk")
p4 <- ggplot(data = subset(stats_df, scenario %in% 4485:6804 & b_z == 1), 
             aes(x = HCRmult, y = ssb_below_blim_total, group = stock,
                 colour = as.factor(`4`))) +
  theme_bw() + geom_line(show.legend = FALSE) + 
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") + ylim(0, 1) +
  labs(x = " ", y = "Blim risk") +
  theme(strip.text = element_blank())
p5 <- ggplot(data = subset(stats_df, scenario %in% 4485:6804 & b_z == 1), 
             aes(x = HCRmult, y = catch_MSY_prop, group = stock,
                 colour = as.factor(`4`))) +
  theme_bw() + geom_line(show.legend = FALSE) + 
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "multiplier", y = "yield/MSY yield") +
  theme(strip.text = element_blank())
p6 <- ggplot(data = subset(stats_df, scenario %in% 4485:6804 & b_z == 1), 
             aes(x = HCRmult, y = iav, group = stock,
                 colour = as.factor(`4`))) +
  theme_bw() + geom_line() + 
  scale_colour_discrete("cluster") + 
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") + ylim(0, NA) + 
  labs(x = " ", y = "catch iav") +
  theme(legend.position = c(0.3, -0.09),
        legend.direction = "horizontal", legend.margin = margin(0,0,0,0))

library(cowplot)
# p_upper <- plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(1, 1, 1.1))
# p_lower <- plot_grid(p4, p5, p6, nrow = 1, rel_widths = c(1, 1, 1.5))

plot_grid(p1, p2, p3, p4, p5, p6,
          nrow = 2, align = "hv")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/", 
                         "multiplier_exponent/all_multiplier2x.png"),
       width = 25, height = 20, units = "cm", dpi = 300, type = "cairo")


### find multiplier where risk is below 5%
### for WKLIFE8 conclusion/recommendations

### reduce to scenarios where only multiplier is changed
df_mult <- subset(stats_df, scenario %in% 4485:6804 & b_z == 1)
### remove cluster 1 stocks (collapsed)
df_mult <- df_mult %>% filter(`4` != 1)


ggplot(data = df_mult, 
       aes(x = HCRmult, y = ssb_below_blim_total, group = stock,
           colour = as.factor(`4`))) +
  theme_bw() + geom_line() + 
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "multiplier", y = "Blim risk")

### median per cluster
df_mult %>% group_by(`4`, fhist, HCRmult) %>%
  summarise(risk = median(ssb_below_blim_total)) %>%
  ggplot(aes(x = HCRmult, y = risk,
             colour = as.factor(`4`))) +
  theme_bw() + geom_line() + 
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "multiplier", y = "Blim risk")

### find highest multiplier where risk < 0.05 per stock
df_mult %>% group_by(stock, fhist) %>%
  filter(ssb_below_blim_total < 0.05) %>% 
  filter(ssb_below_blim_total == max(ssb_below_blim_total)) %>% 
  ggplot(aes(x = k, y = HCRmult, colour = stock)) +
  theme_bw() + #geom_point() + 
  geom_jitter(width = 0, height = 0.01) + 
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "k", y = "multiplier")

### create "new" groups
df_mult %>% 
  select(stock, ssb_below_blim_total, fhist, HCRmult, k, `4`) %>%
  mutate(k19 = ifelse(k < 0.2, "low", "med"),
         k20 = ifelse(k <= 0.2, "low", "med"),
         korig = ifelse(`4` == 2, "low", "med")) %>%
  gather(key = "calculation", value = "group", c("k19", "k20", "korig")) %>%
  group_by(calculation, fhist, HCRmult, group) %>% 
  summarise(risk = median(ssb_below_blim_total)) %>%
  ggplot(aes(x = HCRmult, y = risk, colour = group)) +
  geom_line() +
  theme_bw() +
  facet_grid(fhist ~ calculation) +
  geom_hline(yintercept = 0.05)

### median stocks and fishing scenarios
df_mult %>% 
  select(stock, ssb_below_blim_total, HCRmult, k, `4`) %>%
  mutate(k19 = ifelse(k < 0.2, "low", "med"),
         k20 = ifelse(k <= 0.2, "low", "med"),
         korig = ifelse(`4` == 2, "low", "med")) %>%
  gather(key = "calculation", value = "group", c("k19", "k20", "korig")) %>%
  group_by(calculation, HCRmult, group) %>% 
  summarise(risk = median(ssb_below_blim_total)) %>%
  arrange(calculation, group, risk) %>%
  ggplot(aes(x = HCRmult, y = risk, colour = group)) +
  geom_line() +
  theme_bw() +
  facet_grid( ~ calculation) +
  geom_hline(yintercept = 0.05)
 
### plot range 
df_mult %>% 
  select(stock, ssb_below_blim_total, HCRmult, k, `4`, fhist) %>%
  mutate(k19 = ifelse(k < 0.2, "0.08-0.19", "0.20-0.32")#,
         #k20 = ifelse(k <= 0.2, "low", "med")#,
         #korig = ifelse(`4` == 2, "low", "med")
         ) %>%
  gather(key = "calculation", value = "k range", c("k19"#, "k20"#, "korig"
                                                 )) %>%
  group_by(calculation, HCRmult, `k range`) %>% 
  #filter(calculation == "k19" & HCRmult == 0.8) %>%
  ggplot() +
  #geom_line(aes(x = HCRmult, y = ssb_below_blim_total, colour = group,
  #              group = interaction(stock, fhist))) +
  geom_boxplot(aes(x = (HCRmult), y = ssb_below_blim_total, 
             colour = `k range`, 
             group = interaction(HCRmult, `k range`))) +
  stat_summary(aes(x = HCRmult, y = ssb_below_blim_total,
                   colour = `k range`),
                fun.y = median, geom = "line") +
  theme_bw() +
  #facet_grid( ~ calculation) +
  #facet_wrap( ~ calculation) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  labs(x = "multiplier", y = "Blim risk")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/", 
                         "wklife8/multiplier_risks.png"),
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo")


### ------------------------------------------------------------------------ ###
### stats vs. tuning: constraint ####
### ------------------------------------------------------------------------ ###


ggplot(data = subset(stats_df, scenario %in% 1237:4484 &
                       lower_constraint == 0 & fhist == "one-way") %>% 
         gather(key = "parameter", value = "value", 
                collapse_total, catch_MSY_prop), 
       aes(x = k, y = value, group = stock)) +
  theme_bw() + geom_line() + 
  geom_point(aes(colour = as.factor(upper_constraint)), position = "jitter") +
  facet_wrap(~ parameter, scales = "free_y")
ggplot(data = subset(stats_df, scenario %in% 1237:4484 &
                       lower_constraint == 0 & fhist == "one-way") %>% 
         gather(key = "parameter", value = "value", 
                ssb_rel, f_rel, collapse_total, iav, catch_MSY_prop) %>%
         mutate(upper_constraint = ifelse(is.finite(upper_constraint), 
                                          upper_constraint, 1.55)), 
       aes(x = upper_constraint, y = value, group = stock)) +
  theme_bw() + geom_line() + geom_point() +
  facet_grid(parameter ~ k, scales = "free_y")

df_plot1 <- subset(stats_df, scenario %in% 1237:4484 &
                     lower_constraint == 0 & upper_constraint < Inf)
df_plot2 <- subset(stats_df, scenario %in% 1237:4484 &
                   lower_constraint == 0 & upper_constraint %in% c(1.5, Inf)) %>%
  mutate(upper_constraint = ifelse(upper_constraint > 1.5, 
                                   1.52, upper_constraint))
df_plot3 <- subset(df_plot2, upper_constraint != 1.5)

### plot for report
p1 <- ggplot(data = df_plot1, 
             aes(x = upper_constraint, y = f_rel, group = stock,
                 colour = as.factor(`4`))) +
  theme_bw(base_size = 8) + geom_line() + ylim(0, 1.5) +
  geom_point(data = df_plot3, size = 0.2) +
  geom_line(data = df_plot2, linetype = "dotted", size = 0.1) +
  scale_colour_discrete(guide = FALSE) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = " ", y = "F/Fmsy") +
  theme(strip.text = element_blank())
p2 <- ggplot(data = df_plot1, 
             aes(x = upper_constraint, y = ssb_rel, group = stock,
                 colour = as.factor(`4`))) +
  theme_bw(base_size = 8) + geom_line() +
  geom_point(data = df_plot3, size = 0.2) +
  geom_line(data = df_plot2, linetype = "dotted", size = 0.1) +
  scale_colour_discrete(guide = FALSE) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = " ", y = "B/Bmsy") +
  theme(strip.text = element_blank())
p3 <- ggplot(data = df_plot1, 
             aes(x = upper_constraint, y = collapse_total, group = stock,
                 colour = as.factor(`4`))) +
  theme_bw(base_size = 8) + geom_line() + ylim(0, 1) +
  geom_point(data = df_plot3, size = 0.2) +
  geom_line(data = df_plot2, linetype = "dotted", size = 0.1) +
  scale_colour_discrete(guide = FALSE) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = " ", y = "collapse risk")
p4 <- ggplot(data = df_plot1, 
             aes(x = upper_constraint, y = ssb_below_blim_total, group = stock,
                 colour = as.factor(`4`))) +
  theme_bw(base_size = 8) + geom_line() + ylim(0, 1) +
  geom_point(data = df_plot3, size = 0.2) +
  geom_line(data = df_plot2, linetype = "dotted", size = 0.1) +
  scale_colour_discrete(guide = FALSE) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = " ", y = "Blim risk") +
  theme(strip.text = element_blank())
p5 <- ggplot(data = df_plot1, 
             aes(x = upper_constraint, y = catch_MSY_prop, group = stock,
                 colour = as.factor(`4`))) +
  theme_bw(base_size = 8) + geom_line() + 
  geom_point(data = df_plot3, size = 0.2) +
  geom_line(data = df_plot2, linetype = "dotted", size = 0.1) +
  scale_colour_discrete(guide = FALSE) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "upper constraint", y = "yield/MSY yield") +
  theme(strip.text = element_blank())
p6 <- ggplot(data = df_plot1, 
             aes(x = upper_constraint, y = iav, group = stock,
                 colour = as.factor(`4`))) +
  theme_bw(base_size = 8) + geom_line() +
  geom_point(data = df_plot3, size = 0.2) +
  geom_line(data = df_plot2, linetype = "dotted", size = 0.1) +
  scale_colour_discrete("cluster") +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "", y = "catch iav") +
  theme(legend.position = c(0.3, -0.09),
        legend.direction = "horizontal", legend.margin = margin(0,0,0,0))
### combine into single plot
plot_grid(p1, p2, p3, p4, p5, p6,
          nrow = 2, align = "hv")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/", 
                         "constraints/all_upper_limit2_coloured.png"),
       width = 20, height = 16, units = "cm", dpi = 300, type = "cairo")


### plot effect of lower constraint in combination with 1.2 upper limit
### plot for report
df_plot4 <- subset(stats_df, scenario %in% 1237:4484 &
                     upper_constraint == 1.2)
p1 <- ggplot(data = df_plot4, 
             aes(x = lower_constraint, y = f_rel, group = stock,
                 colour = as.factor(`4`))) +
  theme_bw(base_size = 8) + geom_line() + xlim(0, 1) +
  scale_colour_discrete(guide = FALSE) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = " ", y = "F/Fmsy") +
  theme(strip.text = element_blank())
p2 <- ggplot(data = df_plot4, 
             aes(x = lower_constraint, y = ssb_rel, group = stock,
                 colour = as.factor(`4`))) +
  theme_bw(base_size = 8) + geom_line() + xlim(0, 1) +
  scale_colour_discrete(guide = FALSE) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = " ", y = "B/Bmsy") +
  theme(strip.text = element_blank())
p3 <- ggplot(data = df_plot4, 
             aes(x = lower_constraint, y = collapse_total, group = stock,
                 colour = as.factor(`4`))) +
  theme_bw(base_size = 8) + geom_line() + xlim(0, 1) +
  scale_colour_discrete(guide = FALSE) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = " ", y = "collapse risk")
p4 <- ggplot(data = df_plot4, 
             aes(x = lower_constraint, y = ssb_below_blim_total, group = stock,
                 colour = as.factor(`4`))) +
  theme_bw(base_size = 8) + geom_line() + xlim(0, 1) +
  scale_colour_discrete(guide = FALSE) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = " ", y = "Blim risk") +
  theme(strip.text = element_blank())
p5 <- ggplot(data = df_plot4, 
             aes(x = lower_constraint, y = catch_MSY_prop, group = stock,
                 colour = as.factor(`4`))) +
  theme_bw(base_size = 8) + geom_line() + xlim(0, 1) +
  scale_colour_discrete(guide = FALSE) +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "lower constraint", y = "yield/MSY yield") +
  theme(strip.text = element_blank())
p6 <- ggplot(data = df_plot4, 
             aes(x = lower_constraint, y = iav, group = stock,
                 colour = as.factor(`4`))) +
  theme_bw(base_size = 8) + geom_line() + xlim(0, 1) +
  scale_colour_discrete("cluster") +
  facet_wrap(~ fhist, ncol = 1, strip.position = "right") +
  labs(x = "", y = "catch iav") +
  theme(legend.position = c(0.3, -0.09),
        legend.direction = "horizontal", legend.margin = margin(0,0,0,0))
### combine into single plot
plot_grid(p1, p2, p3, p4, p5, p6,
          nrow = 2, align = "hv")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/", 
                         "constraints/all_lower_limit2_coloured.png"),
       width = 20, height = 16, units = "cm", dpi = 300, type = "cairo")

### ------------------------------------------------------------------------ ###
### change in timing, explore performance ####
### ------------------------------------------------------------------------ ###

### load stats
stats <- readRDS("output/stats_corrected_scn.RDS")
### subset to scenarios
stats_t <- stats %>% filter(scenario %in% c(6863:7210))
### reference points
refpts <- readRDS("input/refpts.rds")
### load cluster allocations
cl_alloc_ow <- readRDS("output/cluster_allocations_one_way.rds")
cl_alloc_rc <- readRDS("output/cluster_allocations_roller_coaster.rds")
### stock abreviations
short_names <- read.csv(file = "input/names_short.csv")
### life history parameters
lhist <- readRDS("input/lhist_extended.rds")

### load quants and calculate B/Bmsy
qnts_t <- foreach(scenario = stats_t$scenario, stock = stats_t$stock) %do% {
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
names(qnts_t) <- stats_t$scenario
qnts_t <- FLQuants(qnts_t)
qnts_t <- as.data.frame(qnts_t)
names(qnts_t)[8] <- "scenario"
qnts_t <- merge(qnts_t, stats_t)

### add short names and vB k
qnts_t <- merge(qnts_t, short_names[, c("stock", "short")])
qnts_t <- merge(qnts_t, lhist)

### plot SSB for all stocks
qnts_t %>% filter(fhist == "one-way") %>%
  mutate(timing = paste(lst_catch, lst_idx),
         wrap = paste0("k=", k, " ", short)) %>%
  ggplot(aes(x = year, y = data, colour = timing, linetype = as.factor(TAC))) +
  geom_line() +
  theme_bw() +
  facet_wrap(~ wrap) +
  scale_linetype_discrete("TAC period\n[years]") +
  scale_colour_discrete("relative timing of\ncatch and index") +
  labs(y = "SSB/Bmsy")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/", 
                         "timing/all_SSB.png"),
       width = 25, height = 16, units = "cm", dpi = 300, type = "cairo")

### create factor levels for stock names based on k
name_levels <- merge(lhist[, c("stock", "k")], 
                     short_names[, c("short", "stock")])[, c("short", "k")]
name_levels <- as.character(name_levels$short[rev(order(name_levels$k))])
### make some plots
stats_t <- stats_t %>%
  full_join(cl_alloc_ow[, c("stock", "4")]) %>%
  mutate(key = as.factor(paste(lst_idx, lst_catch, TAC))) %>%
  mutate(key = factor(key, levels(key)[c(2, 1, 4, 3, 6, 5)])) %>% 
  gather(key = "stat", value = "value", f_rel, ssb_rel, risk_collapse, 
         risk_blim, yield_relFmsy, catch_iav) %>%
  mutate(stat = factor(stat, levels = c("f_rel", "ssb_rel", "risk_collapse", 
                                        "risk_blim", "yield_relFmsy", 
                                        "catch_iav"),
                       labels = c("F/Fmsy", "B/Bmsy", "collapse risk",
                                  "Blim risk", "yield/MSY yield", 
                                  "catch iav"))) %>%
  rename("cluster" = "4") %>%
  mutate(cluster = paste("cluster", cluster)) %>%
  left_join(short_names) %>%
  left_join(lhist)
stats_t$short <- factor(stats_t$short, levels = name_levels)
  
### plot performance stats for all stocks
stats_t %>% filter(fhist == "one-way") %>%
  ggplot(aes(x = short, y = value, fill = key)) +
  facet_grid(stat ~ cluster, scales = "free", space = "free_x") +
  # facet_wrap(~ cluster, scales = "free") +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "stock", y = "") +
  scale_fill_discrete("rel timing\nof index and\ncatch, and\nTAC period")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/", 
                         "timing/all_performance.png"),
       width = 25, height = 16, units = "cm", dpi = 300, type = "cairo")

### ------------------------------------------------------------------------ ###
### WKLIFE8 ####
### ------------------------------------------------------------------------ ###
### new OMs with different natural mortalities and recruitment uncertainty

### new scenarios
scns <- 7907:8718

### scenario definitions
stats <- readRDS("output/stats_corrected_scn_wklife8.RDS")
### lhist
names_short <- read.csv("input/names_short.csv")[, c("stock", "short")]
lhist <- readRDS("input/lhist_extended.rds")

stats <- stats[stats$scenario %in% scns, ]
stats <- stats %>% mutate(
  M_type = ifelse(grepl(pattern = "^gis*", x = OM_scn), "Gislason", "Lorenzen"),
  rec_sd = as.numeric(lapply(strsplit(OM_scn, "_"), "[[", 2)),
  M1 = as.numeric(lapply(strsplit(OM_scn, "_"), "[[", 3)),
  M2 = as.numeric(gsub(x = lapply(strsplit(OM_scn, "_"), "[[", 4),
                       pattern = "/", replacement = ""))
  ) %>%
  left_join(names_short) %>%
  left_join(lhist[, c("linf", "k", "t0", "a", "b", "a50", "l50", "stock")])

### load relative ssb values
ssbs <- foreach(stats_i = split(stats, seq(nrow(stats))),
                .packages = "FLCore", .combine = rbind) %dopar% {
                  
  ### load quants
  quants <- readRDS(paste0("output/perfect_knowledge/combined/corrected/",
                           stats_i$scenario, ".rds"))
  ### get median SSB relative to Bmsy
  SSB <- iterMedians(quants$ssb_rel)
  ### coerce into df
  SSB <- as.data.frame(SSB)
  ### add scenario defintion
  SSB <- cbind(SSB, stats_i)
  
  return(SSB)
  
}

### Gislason only
ssbs %>% filter(fhist == "one-way" & M_type == "Gislason") %>%
  mutate(wrap = paste0("k=", k, " ", short)) %>%
  ggplot(aes(x = year, y = data, group = OM_scn, colour = as.factor(rec_sd))) +
  geom_line() +
  facet_wrap(~ wrap) +
  theme_bw() +
  labs(y = "SSB/Bmsy", title = "Gislason mortality") + 
  geom_hline(yintercept = 1) +
  scale_color_discrete("rec sd")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/", 
                         "M/1_Gislason_rec_sd_ow.png"),
       width = 25, height = 16, units = "cm", dpi = 200, type = "cairo")
ssbs %>% filter(fhist == "roller-coaster" & M_type == "Gislason") %>%
  mutate(wrap = paste0("k=", k, " ", short)) %>%
  ggplot(aes(x = year, y = data, group = OM_scn, colour = as.factor(rec_sd))) +
  geom_line() +
  facet_wrap(~ wrap) +
  theme_bw() +
  labs(y = "SSB/Bmsy", title = "Gislason mortality") + 
  geom_hline(yintercept = 1) +
  scale_color_discrete("rec sd")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/", 
                         "M/2_Gislason_rec_sd_rc.png"),
       width = 25, height = 16, units = "cm", dpi = 200, type = "cairo")


### plot all OM scenarios
ssbs %>% filter(fhist == "one-way") %>%
  mutate(wrap = paste0("k=", k, " ", short)) %>%
  ggplot(aes(x = year, y = data, group = OM_scn)) +
  geom_line() +
  facet_wrap(~ wrap) +
  theme_bw() +
  labs(y = "SSB/Bmsy") + 
  geom_hline(yintercept = 1)
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/", 
                         "M/3_all_ow.png"),
       width = 25, height = 16, units = "cm", dpi = 200, type = "cairo")
ssbs %>% filter(fhist == "roller-coaster") %>%
  mutate(wrap = paste0("k=", k, " ", short)) %>%
  ggplot(aes(x = year, y = data, group = OM_scn)) +
  geom_line() +
  facet_wrap(~ wrap) +
  theme_bw() +
  labs(y = "SSB/Bmsy") + 
  geom_hline(yintercept = 1)
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/", 
                         "M/4_all_rc.png"),
       width = 25, height = 16, units = "cm", dpi = 200, type = "cairo")

### separate mortality functions
ssbs %>% filter(fhist == "one-way") %>%
  mutate(wrap = paste0("k=", k, " ", short)) %>%
  ggplot(aes(x = year, y = data, group = OM_scn, colour = M_type)) +
  geom_line() +
  facet_wrap(~ wrap) +
  theme_bw() +
  labs(y = "SSB/Bmsy") + 
  geom_hline(yintercept = 1) +
  scale_color_discrete("M function")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/", 
                         "M/5_all_M_type_ow.png"),
       width = 25, height = 16, units = "cm", dpi = 200, type = "cairo")
ssbs %>% filter(fhist == "roller-coaster") %>%
  mutate(wrap = paste0("k=", k, " ", short)) %>%
  ggplot(aes(x = year, y = data, group = OM_scn, colour = M_type)) +
  geom_line() +
  facet_wrap(~ wrap) +
  theme_bw() +
  labs(y = "SSB/Bmsy") + 
  geom_hline(yintercept = 1) +
  scale_color_discrete("M function")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/", 
                         "M/6_all_M_type_rc.png"),
       width = 25, height = 16, units = "cm", dpi = 200, type = "cairo")


### Lorenzen & rec sd = 0.3
ssbs %>% filter(fhist == "one-way" & M_type == "Lorenzen" &
                  rec_sd == 0.3) %>%
  mutate(wrap = paste0("k=", k, " ", short)) %>%
  ggplot(aes(x = year, y = data, group = OM_scn,
             colour = paste(M1, M2))) +
  geom_line() +
  facet_wrap(~ wrap) +
  theme_bw() +
  labs(y = "SSB/Bmsy", title = "Lorenzen mortality, rec sd = 0.3") + 
  geom_hline(yintercept = 1) +
  scale_color_discrete("M2 M1")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/", 
                         "M/7_Lorenzen_recsd03_ow.png"),
       width = 25, height = 16, units = "cm", dpi = 200, type = "cairo")
ssbs %>% filter(fhist == "roller-coaster" & M_type == "Lorenzen" &
                  rec_sd == 0.3) %>%
  mutate(wrap = paste0("k=", k, " ", short)) %>%
  ggplot(aes(x = year, y = data, group = OM_scn,
             colour = paste(M1, M2))) +
  geom_line() +
  facet_wrap(~ wrap) +
  theme_bw() +
  labs(y = "SSB/Bmsy", title = "Lorenzen mortality, rec sd = 0.3") + 
  geom_hline(yintercept = 1) +
  scale_color_discrete("M2 M1")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/", 
                         "M/8_Lorenzen_recsd03_rc.png"),
       width = 25, height = 16, units = "cm", dpi = 200, type = "cairo")



### absolute SSB
### load relative ssb values
ssbs_abs <- foreach(stats_i = split(stats, seq(nrow(stats))),
                .packages = "FLCore", .combine = rbind) %dopar% {
                  
  ### load quants
  quants <- readRDS(paste0("output/perfect_knowledge/combined/corrected/",
                           stats_i$scenario, ".rds"))
  ### get median SSB relative to Bmsy
  SSB <- iterMedians(quants$ssb)
  ### coerce into df
  SSB <- as.data.frame(SSB)
  ### add scenario defintion
  SSB <- cbind(SSB, stats_i)
  
  return(SSB)
  
}
### plot all OM scenarios
ssbs_abs %>% filter(fhist == "one-way") %>%
  mutate(wrap = paste0("k=", k, " ", short)) %>%
  ggplot(aes(x = year, y = data, group = OM_scn)) +
  geom_line() +
  facet_wrap(~ wrap) +
  theme_bw() +
  labs(y = "SSB") + 
  geom_hline(yintercept = 1)
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/", 
                         "M/3_all_ow.png"),
       width = 25, height = 16, units = "cm", dpi = 300, type = "cairo")