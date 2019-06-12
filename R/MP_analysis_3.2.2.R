library(FLCore)
library(ggplotFL)
library(Cairo)
library(tidyverse)

library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)

### ------------------------------------------------------------------------ ###
### compare Fproxy settings: real Fmsy and LF=M ####
### ------------------------------------------------------------------------ ###


### load stats
stats <- readRDS("output/stats_corrected_scn.RDS")
### subset to scenarios
stats_t <- stats %>% filter(scenario %in% c(7211:7906))
### reference points
refpts <- readRDS("input/refpts.rds")
### stock abreviations
short_names <- read.csv(file = "input/names_short.csv")
### life history parameters
lhist <- readRDS("input/lhist_extended.rds")

### load quants for plotting
qnts_t <- foreach(scenario = as.numeric(as.character(stats_t$scenario)), 
                  stock = as.character(stats_t$stock),
                  .combine = rbind, 
                  .errorhandling = "pass") %dopar% {
  #if (scenario == 7212) browser()
  #browser()
  ### load SSB
  stk <- readRDS(paste0("output/perfect_knowledge/combined/corrected/",
                          scenario, ".rds"))
  ### calculate relative SSB
  ssb_rel <- stk$ssb / c(refpts[[stock]]["msy", "ssb"])
  
  ### combine quants into FLQuants object
  qnts_tmp <- FLQuants(SSB = stk$ssb, SSB_rel = ssb_rel, Fbar = stk$fbar,
                       catch = stk$catch)
  ### use median only
  qnts_tmp <- lapply(qnts_tmp, iterMedians)
  ### coerce into data frame
  qnts_tmp <- as.data.frame(qnts_tmp)
  ### add stock name
  qnts_tmp$stock <- stock
  qnts_tmp$scenario <- scenario
  ### return
  return(qnts_tmp)
}
### add stats and lhist
qnts_t <- qnts_t %>%
  left_join(stats_t) %>%
  left_join(short_names) %>%
  left_join(lhist)

### plot time series for some quants with default parametrization
qnts_t %>% filter(TAC == 2 & lst_idx == -1 & lst_catch == -1 &
                    FproxyMSY_type == "FproxyMSY") %>%
  mutate(qname = as.character(qname)) %>%
  mutate(qname = as.character(ifelse(qname == "SSB_rel", "SSB/Bmsy", qname))) %>%
  ggplot(aes(x = year, y = data, group = short, 
             colour = k >= 0.35)) +
  geom_line() +
  theme_bw() +
  facet_grid(qname ~ fhist, scales = "free") +
  labs(y = "") +
  scale_colour_discrete("k >= 0.35")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/", 
                         "all_quants_default.png"),
       width = 25, height = 16, units = "cm", dpi = 300, type = "cairo")

### plot stats for default parametrization
stats_plot <- stats_t %>%
  gather(key = "stat", value = "value", f_rel, ssb_rel, risk_collapse, 
         risk_blim, yield_relFmsy, catch_iav) %>%
  mutate(stat = factor(stat, levels = c("f_rel", "ssb_rel", "risk_collapse", 
                                        "risk_blim", "yield_relFmsy", 
                                        "catch_iav"),
                       labels = c("F/Fmsy", "B/Bmsy", "collapse risk",
                                  "Blim risk", "yield/MSY yield", 
                                  "catch iav"))) %>%
  left_join(short_names) %>%
  left_join(lhist)
stats_plot$short <- factor(stats_plot$short, levels = name_levels)
stats_plot %>% filter(TAC == 2 & lst_idx == -1 & lst_catch == -1 &
                     FproxyMSY_type == "FproxyMSY") %>%
  ggplot(aes(x = short, y = value)) +
  facet_grid(stat ~ fhist, scales = "free", space = "free_x") +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "stock", y = "")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/", 
                         "timing/all_performance.png"),
       width = 25, height = 16, units = "cm", dpi = 300, type = "cairo")

### plot SSB for all stocks: compare proxies
qnts_t %>% filter(fhist == "one-way" & TAC == 2 & lst_idx == -1 &
                    lst_catch == -1 & qname == "SSB_rel") %>%
  mutate(timing = paste(lst_catch, lst_idx),
         wrap = paste0("k=", k, " ", short),
         FproxyMSY_type = ifelse(FproxyMSY_type == "FproxyMSY", "Fmsy",
                                 "LF=M")) %>%
  ggplot(aes(x = year, y = data, colour = FproxyMSY_type)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~ wrap) +
  scale_colour_discrete("Fproxy") +
  labs(y = "SSB/Bmsy") +
  geom_hline(yintercept = 1)
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/", 
                         "proxies/all_SSB.png"),
       width = 25, height = 16, units = "cm", dpi = 300, type = "cairo")
### roller-coaster
qnts_t %>% filter(fhist == "roller-coaster" & TAC == 2 & lst_idx == -1 &
                    lst_catch == -1 & qname == "SSB_rel") %>%
  mutate(timing = paste(lst_catch, lst_idx),
         wrap = paste0("k=", k, " ", short),
         FproxyMSY_type = ifelse(FproxyMSY_type == "FproxyMSY", "Fmsy",
                                 "LF=M")) %>%
  ggplot(aes(x = year, y = data, colour = FproxyMSY_type)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~ wrap) +
  scale_colour_discrete("Fproxy") +
  labs(y = "SSB/Bmsy") +
  geom_hline(yintercept = 1)
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/", 
                         "proxies/all_SSB_rc.png"),
       width = 25, height = 16, units = "cm", dpi = 300, type = "cairo")


### plot SSB: compare timing
qnts_t %>% filter(fhist == "one-way" & qname == "SSB_rel" &
                    FproxyMSY_type == "FproxyMSY") %>%
  mutate(timing = paste(lst_idx),
         wrap = paste0("k=", k, " ", short),
         FproxyMSY_type = ifelse(FproxyMSY_type == "FproxyMSY", "Fmsy",
                                 "LF=M")) %>%
  ggplot(aes(x = year, y = data, colour = timing, linetype = as.factor(TAC))) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_line() +
  theme_bw() +
  facet_wrap(~ wrap) +
  scale_linetype_discrete("TAC period\n[years]") +
  scale_colour_discrete("relative timing\nof index") +
  labs(y = "SSB/Bmsy")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/", 
                         "timing/all_SSB_ow.png"),
       width = 25, height = 16, units = "cm", dpi = 300, type = "cairo")
### same for roller-coaster
qnts_t %>% filter(fhist == "roller-coaster" & qname == "SSB_rel" &
                    FproxyMSY_type == "FproxyMSY") %>%
  mutate(timing = paste(lst_idx),
         wrap = paste0("k=", k, " ", short),
         FproxyMSY_type = ifelse(FproxyMSY_type == "FproxyMSY", "Fmsy",
                                 "LF=M")) %>%
  ggplot(aes(x = year, y = data, colour = timing, linetype = as.factor(TAC))) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_line() +
  theme_bw() +
  facet_wrap(~ wrap) +
  scale_linetype_discrete("TAC period\n[years]") +
  scale_colour_discrete("relative timing\nof index") +
  labs(y = "SSB/Bmsy")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/", 
                         "timing/all_SSB_rc.png"),
       width = 25, height = 16, units = "cm", dpi = 300, type = "cairo")

### plot SSB, one-way, Fmsy proxy: compare timing
qnts_t %>% filter(fhist == "one-way" & qname == "SSB_rel" &
                    FproxyMSY_type == "FproxyMSY") %>%
  mutate(timing = paste(lst_idx),
         wrap = paste0("k=", k, " ", short),
         FproxyMSY_type = ifelse(FproxyMSY_type == "FproxyMSY", "Fmsy",
                                 "LF=M")) %>%
  ggplot(aes(x = year, y = data, colour = timing, linetype = as.factor(TAC))) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_line() +
  theme_bw() +
  facet_wrap(~ wrap) +
  scale_linetype_discrete("TAC period\n[years]") +
  scale_colour_discrete("relative timing\nof index") +
  labs(y = "SSB/Bmsy")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/", 
                         "timing/all_SSB_ow.png"),
       width = 25, height = 16, units = "cm", dpi = 300, type = "cairo")
### same for roller-coaster
qnts_t %>% filter(fhist == "roller-coaster" & qname == "SSB_rel" &
                    FproxyMSY_type == "FproxyMSY") %>%
  mutate(timing = paste(lst_idx),
         wrap = paste0("k=", k, " ", short),
         FproxyMSY_type = ifelse(FproxyMSY_type == "FproxyMSY", "Fmsy",
                                 "LF=M")) %>%
  ggplot(aes(x = year, y = data, colour = timing, linetype = as.factor(TAC))) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_line() +
  theme_bw() +
  facet_wrap(~ wrap) +
  scale_linetype_discrete("TAC period\n[years]") +
  scale_colour_discrete("relative timing\nof index") +
  labs(y = "SSB/Bmsy")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/", 
                         "timing/all_SSB_rc.png"),
       width = 25, height = 16, units = "cm", dpi = 300, type = "cairo")

### plot SSB, one-way, LF=M proxy: compare timing
qnts_t %>% filter(fhist == "one-way" & qname == "SSB_rel" &
                    FproxyMSY_type == "FproxyMSY_MK") %>%
  mutate(timing = paste(lst_idx),
         wrap = paste0("k=", k, " ", short),
         FproxyMSY_type = ifelse(FproxyMSY_type == "FproxyMSY", "Fmsy",
                                 "LF=M")) %>%
  ggplot(aes(x = year, y = data, colour = timing, linetype = as.factor(TAC))) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_line() +
  theme_bw() +
  facet_wrap(~ wrap) +
  scale_linetype_discrete("TAC period\n[years]") +
  scale_colour_discrete("relative timing\nof index") +
  labs(y = "SSB/Bmsy")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/", 
                         "timing/all_SSB_LFeM_ow.png"),
       width = 25, height = 16, units = "cm", dpi = 300, type = "cairo")
### same for roller-coaster
qnts_t %>% filter(fhist == "roller-coaster" & qname == "SSB_rel" &
                    FproxyMSY_type == "FproxyMSY_MK") %>%
  mutate(timing = paste(lst_idx),
         wrap = paste0("k=", k, " ", short),
         FproxyMSY_type = ifelse(FproxyMSY_type == "FproxyMSY", "Fmsy",
                                 "LF=M")) %>%
  ggplot(aes(x = year, y = data, colour = timing, linetype = as.factor(TAC))) +
  geom_hline(yintercept = 1, colour = "grey") +
  geom_line() +
  theme_bw() +
  facet_wrap(~ wrap) +
  scale_linetype_discrete("TAC period\n[years]") +
  scale_colour_discrete("relative timing\nof index") +
  labs(y = "SSB/Bmsy")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/", 
                         "timing/all_SSB_LFeM_rc.png"),
       width = 25, height = 16, units = "cm", dpi = 300, type = "cairo")



### create factor levels for stock names based on k
name_levels <- merge(lhist[, c("stock", "k")], 
                     short_names[, c("short", "stock")])[, c("short", "k")]
name_levels <- as.character(name_levels$short[rev(order(name_levels$k))])
### make some plots
stats_plot <- stats_t %>%
  #full_join(cl_alloc_ow[, c("stock", "4")]) %>%
  mutate(key = as.factor(paste(lst_idx, TAC))) %>%
  mutate(key = factor(key, levels(key)[c(2, 1, 4, 3, 6, 5)])) %>% 
  gather(key = "stat", value = "value", f_rel, ssb_rel, risk_collapse, 
         risk_blim, yield_relFmsy, catch_iav) %>%
  mutate(stat = factor(stat, levels = c("f_rel", "ssb_rel", "risk_collapse", 
                                        "risk_blim", "yield_relFmsy", 
                                        "catch_iav"),
                       labels = c("F/Fmsy", "B/Bmsy", "collapse risk",
                                  "Blim risk", "yield/MSY yield", 
                                  "catch iav"))) %>%
  left_join(short_names) %>%
  left_join(lhist)
stats_plot$short <- factor(stats_plot$short, levels = name_levels)

### plot performance stats for all stocks
stats_plot %>% filter(fhist == "one-way" & FproxyMSY_type == "FproxyMSY") %>%
  ggplot(aes(x = short, y = value, fill = key)) +
  #facet_wrap(stat ~ , scales = "free", space = "free_x") +
  facet_grid(stat ~ ., scales = "free") +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "stock", y = "") +
  scale_fill_discrete("rel timing\nof index and\nTAC period")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/", 
                         "timing/all_performance_ow.png"),
       width = 25, height = 16, units = "cm", dpi = 300, type = "cairo")
### same for roller-coaster
stats_plot %>% filter(fhist == "roller-coaster" & 
                        FproxyMSY_type == "FproxyMSY") %>%
  ggplot(aes(x = short, y = value, fill = key)) +
  #facet_wrap(stat ~ , scales = "free", space = "free_x") +
  facet_grid(stat ~ ., scales = "free") +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "stock", y = "") +
  scale_fill_discrete("rel timing\nof index and\nTAC period")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/", 
                         "timing/all_performance_rc.png"),
       width = 25, height = 16, units = "cm", dpi = 300, type = "cairo")
