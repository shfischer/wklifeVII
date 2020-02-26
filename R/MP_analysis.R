### ------------------------------------------------------------------------ ###
### script for analysing MSE results ####
### ------------------------------------------------------------------------ ###
### preparation for plotting, penalized regression, clustering, etc.

library(mseDL)
library(tidyverse)
library(doParallel)

stocks <- read.csv("input/stock_list_full2.csv", stringsAsFactors = FALSE)
lhist <- read.csv("input/lhist_extended.csv", stringsAsFactors = FALSE)
OM_scns <- read.csv("input/OM_scns.csv", stringsAsFactors = FALSE)
### path to final scenario
path_default <- "output/observation_error/new_baseline/"
### stock labels
stk_labels <- stocks %>%
  select(stock_short, k) %>%
  mutate(stock = stock_short,
         label = paste0("k = ", k, ", ", stock),
         label2 = paste0("italic(k)==", k, "~~", stock))
### reference points
refpts <- readRDS("input/refpts_paper.rds")$new_baseline

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

### check terminal SSB relative to SSB0
depletion <- sapply(seq_along(names(files)), function(x) {
  median(res_corrected[[x]]$ssb[, ac(200)])/1000
})
names(depletion) <- names(files)

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



### ------------------------------------------------------------------------ ###
### clustering ####
### ------------------------------------------------------------------------ ###

library(dtwclust)
library(ggdendro)

### OM scenario
OM_scn <- 38

### get files
path <- paste0("output/observation_error/", OM_scns$id[[OM_scn]], "/")
files <- paste0(path_default, "/corrected/one-way/", stocks$stock, ".rds")
names(files) <- stocks$stock_short
res_corrected <- lapply(files, readRDS)

### extract SSB/Bmsy
SSBrel <- lapply(res_corrected, "[[", "ssb_rel")
SSBrel <- lapply(SSBrel, iterMedians)
SSBrel <- lapply(SSBrel, window, start = 101)
SSBmat <- lapply(SSBrel, function(x) as.data.frame(t(x[, drop = TRUE])))

df <- as.data.frame(FLQuants(SSBrel))
ggplot(data = df, aes(x = year, y = data)) +
  geom_line() +
  facet_wrap(~ qname)

### hierarchical clustering
cl_one_way <- tsclust(SSBmat, 
                      type = "hierarchical",
                      k = 10, distance = "dtw", seed = 1,
                      trace = TRUE, error.check = TRUE)
plot(cl_one_way)
### save
saveRDS(cl_one_way, file = paste0(path_default, "corrected/one-way/cluster.rds"))

clusters <- cl_one_way
### find all possible clusters and allocations of stocks to them
cl_alloc <- cutree(clusters, k = 1:29)

### stock SSB medians
stock_medians <- clusters@datalist
stock_medians <- lapply(stock_medians, function(x) {
  apply(x, 2, median)
})
stock_medians <- as.data.frame(do.call(rbind, stock_medians))
stock_medians$stock <- names(clusters@datalist)

### extract medians and create "centroids"
df_plot <- foreach(k = colnames(cl_alloc)) %do% {
  #browser()
  ### load SSB medians per stock
  res <- stock_medians
  ### allocate to cluster
  #identical(res$stock, rownames(cl_alloc))
  res$cluster <- cl_alloc[, k]
  ### cluster median
  ### go through each cluster
  res_add <- lapply(split(res$stock, res$cluster), function(stocks) {
    ### extract SSBs 
    tmp <- res[stocks, ]
    tmp$stock <- 0
    ### calculate median
    apply(tmp, MARGIN = 2, mean)
  })
  ### combine
  res_add <- as.data.frame(do.call(rbind, res_add))
  res_add$stock <- "cluster"
  res_add$cluster <- as.numeric(rownames(res_add))
  ### add 
  res <- rbind(res, res_add)
  res$cluster <- as.numeric(res$cluster)
  ### add number of clusters
  res$k = k
  return(res)
}
df_plot <- do.call(rbind, df_plot)
df_plot <- df_plot %>% 
  gather(key = "year", value = "value", "101":"200") %>%
  mutate(year = as.numeric(year) - 100,
         k = as.numeric(k))

### load lhist - k
lhist <- read.csv("input/stock_list_full2.csv", stringsAsFactors = FALSE)
lhist_k <- lhist[, c("stock_short", "k")]
names(lhist_k)[1] <- "stock"
lhist_k <- lhist_k[order(lhist_k$k), ]
lhist_k$k_pos <- seq(nrow(lhist_k))
### merge k and cluster allocations
cl_alloc_k <- as.data.frame(cl_alloc)
cl_alloc_k$stock <- rownames(cl_alloc_k)
cl_alloc_k <- merge(cl_alloc_k, lhist_k)
cl_alloc_k <- gather(cl_alloc_k, key = "n_cluster", value = "cluster",
                     `1`:`29`)
cl_alloc_k$column = "k"
cl_alloc_k$label <- " "

### data for dendrogram
dat_dend <- dendro_data(clusters)
dat_dend$labels$label <- as.character(dat_dend$labels$label)

### add cluster to stock labels
cl_alloc_dend <- as.data.frame(cl_alloc)
cl_alloc_dend$label <- row.names(cl_alloc_dend)
dat_dend$labels <- full_join(dat_dend$labels, cl_alloc_dend)


### save data for plotting
saveRDS(list(curves = df_plot, dend = dat_dend, bars = cl_alloc_k),
        file = paste0(path_default, "corrected/one-way/plot_data_cluster.rds"))
cluster <- readRDS(paste0(path_default, "corrected/one-way/plot_data_cluster.rds"))

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


### ------------------------------------------------------------------------ ###
### regression model ####
### ------------------------------------------------------------------------ ###
library(glmnet)

### get stats
stats <- read.csv(paste0(path_default, "stats.csv"), stringsAsFactors = FALSE)

### select data for glmnet
glm_data <- stats %>% 
  #filter(fhist == "one-way") %>%
  select(stock, paper, f_rel, ssb_rel, risk_collapse, 
         risk_blim, yield_rel, iav, fhist) %>%
  full_join(lhist)

### select predictor and response variables
var_pred_all <- c("a", "b", "linf", "k", "t0", "a50", "alpha", "beta", "spr0",
              "lopt", "r", "rc", "m", "mk", "fmsym", "bmsyb0")
var_pred_primary <- c("a", "b", "linf", "k", "t0", "a50")
var_resp <- c("f_rel", "ssb_rel", "risk_collapse", 
              "risk_blim", "yield_rel", "iav")

### lasso regression
ow_primary <- cv.glmnet(
  x = as.matrix(glm_data[glm_data$fhist == "one-way", var_pred_primary]),
  y = as.matrix(glm_data[glm_data$fhist == "one-way", var_resp]),
  foldid = 1:29, nlambda = 1000, family = "mgaussian", alpha = 1)
ow_full <- cv.glmnet(
  x = as.matrix(glm_data[glm_data$fhist == "one-way", var_pred_all]),
  y = as.matrix(glm_data[glm_data$fhist == "one-way", var_resp]),
  foldid = 1:29, nlambda = 1000, family = "mgaussian", alpha = 1)
rc_primary <- cv.glmnet(
  x = as.matrix(glm_data[glm_data$fhist == "roller-coaster", var_pred_primary]),
  y = as.matrix(glm_data[glm_data$fhist == "roller-coaster", var_resp]),
  foldid = 1:29, nlambda = 1000, family = "mgaussian", alpha = 1)

### coefficients
coef(ow_primary, s = "lambda.min")
coef(ow_full, s = "lambda.min")
coef(ow_primary, s = "lambda.1se")
coef(ow_full, s = "lambda.1se")
coef(rc_primary, s = "lambda.1se")
### mean squared error
min(ow_primary$cvm)
min(ow_full$cvm)

### fit lm separately to each of the response variables
lms_ow <- lapply(var_resp, function(x) {
  lm(glm_data[, x] ~ glm_data[, "k"])
})
names(lms_ow) <- var_resp
lapply(lms_ow, summary)
lms <- lapply(var_resp, function(x) {
  tmp <- lm(glm_data[, x] ~ glm_data[, "k"])
  tmp <- as.data.frame(t(data.frame(coef(tmp))))
  colnames(tmp) <- c("Intercept", "k")
  tmp$var <- x
  tmp
})
names(lms) <- var_resp
lms <- do.call(rbind, lms)
res_glmnet <- coef(ow_full, s = "lambda.1se")
### save glmnet results
saveRDS(object = list(data = glm_data, glmnet = res_glmnet, lms = lms, 
                      var_resp = var_resp, var_pred = var_pred_primary), 
        file = paste0(path_default, "corrected/one-way/glmnet.rds"))

glmnet_res <- readRDS(paste0(path_default, "corrected/one-way/glmnet.rds"))

### plots
p1 <- ggplot(data = glmnet_res$data,
             aes(x = k, y = f_rel)) +
  geom_point(size = 0.3) + theme_paper + 
  lims(x = c(0, NA), y = c(0, NA)) +
  geom_abline(data = data.frame(
    intercept = c(glmnet_res$glmnet$f_rel["(Intercept)",],
                  lms[lms$var == "f_rel", "Intercept"]),
    slope = c(glmnet_res$glmnet$f_rel_new["k",],
              lms[lms$var == "f_rel", "k"]),
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
                  lms[lms$var == "ssb_rel", "Intercept"]),
    slope = c(glmnet_res$glmnet$ssb_rel["k",],
              lms[lms$var == "ssb_rel", "k"]),
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
                  lms[lms$var == "risk_collapse", "Intercept"]),
    slope = c(glmnet_res$glmnet$risk_collapse["k",],
              lms[lms$var == "risk_collapse", "k"]),
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
                  lms[lms$var == "risk_blim", "Intercept"]),
    slope = c(glmnet_res$glmnet$risk_blim["k",],
              lms[lms$var == "risk_blim", "k"]),
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
                  lms[lms$var == "yield_rel", "Intercept"]),
    slope = c(glmnet_res$glmnet$yield_rel["k",],
              lms[lms$var == "yield_rel", "k"]),
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

### plot saved in MP_plots.R

### check elastic-net regularisation
### test alpha (lasso vs. ridge)
foldid <- 1:29
alphas <- seq(from = 0, to = 1, 0.01)
names(alphas) <- alphas
fits <- foreach(x = alphas, .errorhandling = "pass") %dopar% {
  cv.glmnet(
    foldid = 1:29, nlambda = 1000, alpha = x, family = "mgaussian",
    x = as.matrix(glm_data[glm_data$fhist == "one-way", var_pred_all]),
    y = as.matrix(glm_data[glm_data$fhist == "one-way", var_resp]))
}
names(fits) <- alphas
### extract mean-sqared error
mses <- lapply(alphas, function(x) {
  data.frame(lambda = fits[[as.character(x)]]$lambda, 
             mse = fits[[as.character(x)]]$cvm, alpha = x)
})
mses <- do.call(rbind, mses)
ggplot(data = mses[], 
       aes(x = log(lambda), y = mse, colour = alpha,
           group = alpha)) +
  geom_line() + theme_bw() +
  labs(x = "log(Lambda)", y = "Mean-Squared Error")
ggsave(filename = "output/perfect_knowledge/plots/3.2.1/glmnet/alphas_mgauss_revised.png",
       width = 20, height = 15, units = "cm", dpi = 300, type = "cairo")
### find alpha with lowest error
a <- mses$alpha[which.min(mses$mse)]
### coefficients
coef(fits[[as.character(a)]])
### mse
min(fits[[as.character(a)]]$cvm)
### before regularisation
min(ow_full$cvm)



### ------------------------------------------------------------------------ ###
### multiplier ####
### ------------------------------------------------------------------------ ###

### get cluster allocations 1-4
cl_alloc <- readRDS(paste0(path_default, "corrected/one-way/cluster.rds"))
cl_alloc <- cutree(cl_alloc, k = 4)
cl_alloc <- data.frame(paper = names(cl_alloc), cluster = unlist(cl_alloc))

### get stats
stats <- read.csv(paste0(path_default, "stats.csv"), stringsAsFactors = FALSE)
stats_tuning <- read.csv(paste0(path_default, "stats_tuning.csv"), 
                         stringsAsFactors = FALSE)

### data for plot
data_mult <- stats_tuning %>% 
  filter(grepl(x = dir, pattern = "HCRmult*")) %>%
  bind_rows(stats) %>%
  mutate(HCRmult = ifelse(is.na(HCRmult), 1, HCRmult)) %>%
  left_join(cl_alloc) %>%
  mutate(dir = ifelse(is.na(dir), "", dir))
saveRDS(object = data_mult, file = paste0(path_default, 
                                          "plot_data_multiplier.rds"))

### check when collapse occurs
collapse <- foreach(scn = split(data_mult, seq(nrow(data_mult)))) %do% {
  tmp <- readRDS(paste0(path_default, "corrected/", scn$fhist, "/",
                        scn$dir, "/", scn$stock, ".rds"))
  duration <- head(which(iterMedians(tmp$ssb)[, ac(101:200)] == 0), 1)
  duration <- ifelse(length(duration) > 0, duration, NA)
  return(duration)
}
data_mult$collapse <- unlist(collapse)
data_mult %>% 
  filter(paper == "ane" & fhist == "one-way")
data_mult %>% 
  filter(paper == "sar" & fhist == "roller-coaster")
data_mult %>% 
  select(paper, fhist, HCRmult, collapse) %>%
  arrange(paper, fhist, HCRmult) %>% View()

data_mult %>% 
  select(paper, fhist, HCRmult, collapse) %>%
  arrange(paper, fhist, HCRmult) %>%
  filter(!is.na(collapse))
### collapse at lowest HCRmult
data_mult %>% 
  select(paper, fhist, HCRmult, collapse) %>%
  arrange(paper, fhist, HCRmult) %>%
  filter(!is.na(collapse)) %>%
  arrange(HCRmult)
### find HCRmult with catch maximum
data_mult %>% 
  select(paper, fhist, HCRmult, yield_rel) %>%
  full_join(cl_alloc) %>%
  group_by(cluster, paper, fhist) %>%
  summarise(max_yield = max(yield_rel),
            HCRmult = HCRmult[which.max(yield_rel)]) %>%
  arrange(fhist, cluster, desc(HCRmult), paper) %>%
  print(n = Inf)


### ------------------------------------------------------------------------ ###
### catch constraints ####
### ------------------------------------------------------------------------ ###

stats <- read.csv(paste0(path_default, "stats.csv"), stringsAsFactors = FALSE)
stats_tuning <- read.csv(paste0(path_default, "stats_tuning.csv"), 
                         stringsAsFactors = FALSE)
data_constr <- stats_tuning %>% 
  filter(!is.na(upper) | !is.na(lower)) %>%
  bind_rows(stats) %>%
  mutate(upper = ifelse(is.na(upper), Inf, upper),
         lower = ifelse(is.na(lower), 0, lower)) %>%
  left_join(cl_alloc)
saveRDS(object = data_constr, 
        file = paste0(path_default, "plot_data_constraints.rds"))

### find upper with catch maximum
data_constr %>%
  filter(lower == 0) %>%
  select(paper, fhist, upper, yield_rel) %>%
  full_join(cl_alloc) %>%
  group_by(cluster, paper, fhist) %>%
  summarise(max_yield = max(yield_rel),
            upper = upper[which.max(yield_rel)]) %>%
  arrange(fhist, cluster, desc(upper), paper) %>%
  print(n = Inf)

### ------------------------------------------------------------------------ ###
### timing ####
### ------------------------------------------------------------------------ ###

timing_scns <- 
  data.frame(dir = c("", "interval-1", "interval-1_lstidx-0",
                     "interval-1_lstidx-1_lstcatch-0", "lstidx-0", 
                     "lstidx-1_lstcatch-0"),
             catch = c(-1, -1, -1, 0, -1, 0),
             idx = c(-1, -1, 0, 1, 0, 1),
             interval = c(2, 1, 1, 1, 2, 2))

### get quants
res <- lapply(seq(nrow(timing_scns)), function(x) {
  qnts <- lapply(seq_along(stocks$stock), function(y) {
    tmp <- readRDS(paste0("output/observation_error/new_baseline/corrected/",
                          "one-way/", timing_scns$dir[x], "/", stocks$stock[y], 
                          ".rds"))
    iterMedians(tmp$ssb_rel)
  })
  names(qnts) <- stocks$stock_short
  qnts <- as.data.frame(FLQuants(qnts))
  qnts$interval <- timing_scns$interval[x]
  qnts$idx <- timing_scns$idx[x]
  qnts$catch <- timing_scns$catch[x]
  return(qnts)
})
res <- do.call(rbind, res)
res <- res %>% mutate(stock = qname) %>%
  full_join(stk_labels) %>% 
  mutate(label = as.factor(label),
         label2 = as.factor(label2)) %>%
  filter(year >= 0)
res$label <- factor(res$label, levels = levels(res$label)[c(1:19, 21, 20, 22:29)])
res$label2 <- factor(res$label2, levels = levels(res$label2)[c(1:19, 21, 20, 22:29)])
res <- res %>%
  mutate(timing = paste(catch, idx),
         TAC = ifelse(interval == 1, "annual", "biennial"),
         timing = paste(catch, idx),
         timing = ifelse(timing == "0 1", "0 +1", timing))

saveRDS(object = res, 
        file = paste0(path_default, "plot_data_timing.rds"))


### plot all stocks
res %>% 
  filter(k >= 0.32) %>%
  ggplot(aes(x = year, y = data, linetype = timing, colour = as.factor(TAC))) +
  geom_line(size = 0.2) +
  theme_paper +
  facet_wrap(~ label2, labeller = label_parsed) +
  scale_colour_manual("TAC period", values = c("grey", "black")) +
  scale_linetype_discrete("relative timing\nof catch and\nindex") +
  labs(y = expression(italic(SSB/B[MSY])))

### ------------------------------------------------------------------------ ###
### contribution of catch rule components ####
### ------------------------------------------------------------------------ ###

### get results
df_comps <- lapply(c("her", "pol"), function(stock) {
  
  stk_long <- stocks$stock[stocks$stock_short == stock]
  res_mse <- readRDS(paste0(path_default, "one-way/", stk_long, ".rds"))
  res_corrected <- readRDS(paste0(path_default, "corrected/one-way/", stk_long, ".rds"))
  ### get tracking object
  comps <- res_mse@tracking[c("HCR3.2.1r", "HCR3.2.1f", "HCR3.2.1b"),
                           ac(seq(from = 99, to = 200, by = 2))]
  ### find collapses and remove values
  for (year_i in dimnames(comps)$year) {
    comps[, year_i,,,, which(res_corrected$ssb[, year_i ] == 0)] <- NA
  }
  
  ### coerce into df for plotting
  df_comps <- as.data.frame(comps)
  df_comps$metric <- substr(x = df_comps$metric, start = 9, stop = 9)
  
  ### add relative SSB/Bmsy
  ssb_rel <- res_corrected$ssb_rel[, ac(99:198)]
  df_ssb <- as.data.frame(ssb_rel)
  names(df_ssb)[1] <- "metric"
  df_ssb$metric <- "SSB"
  ### median ssb
  ssb_med <- apply(ssb_rel, 2, median)
  df_ssb_med <- as.data.frame(ssb_med)
  names(df_ssb_med)[1] <- "metric"
  df_ssb_med$metric <- "SSB_median"
  df_res <- rbind(df_comps, df_ssb, df_ssb_med)
  df_res$stock <- stock
  return(df_res)
})
df_comps <- do.call(rbind, df_comps)
### add product
df_comps <- df_comps %>% spread(metric, data)
df_comps$prod <- with(df_comps, r * f * b)
### gather again
df_comps <- df_comps %>% 
  gather(key = "component", value = "data", r, f, b, prod, SSB, SSB_median)
### change years to start from 1
df_comps <- df_comps %>%
  mutate(year = year - 98) %>%
  left_join(stk_labels) %>%
  mutate(stk_label = label,
         stk_label2 = label2,
         label = NULL, label2 = NULL)
df_comps <- df_comps %>% filter(!is.na(data))
df_comps <- df_comps %>%
  mutate(percentile = ifelse(component != "SSB_median", "median", 
                             "SSB/Bmsy"))
df_comps_perc <- df_comps %>%
  group_by(stock, stk_label2, year, percentile, component) %>%
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
### save
saveRDS(object = df_comps_perc, 
        file = paste0(path_default, "components_data_plot.rds"))

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
  ylim(0, NA)
p



### ------------------------------------------------------------------------ ###
### perfect information ####
### ------------------------------------------------------------------------ ###

scns <- c(38, 53)
values <- c("default", "perfect\ninformation")
res <- lapply(seq_along(scns), function(i_scn) {#browser()
  tmp <- lapply(seq_along(stocks$stock), function(i_stock) {
    qnt_i <- readRDS(paste0("output/",
                            ifelse(OM_scns$obs_error[scns[i_scn]], 
                                   "observation_error", "perfect_knowledge"), 
                            "/", OM_scns$id[scns[i_scn]], "/",
                            "corrected/one-way/", stocks$stock[i_stock], ".rds"))
    iterMedians(qnt_i$ssb_rel)
  })
  names(tmp) <- stocks$stock_short
  tmp <- as.data.frame(FLQuants(tmp))
  tmp$scenario <- values[i_scn]
  return(tmp)
})
res <- do.call(rbind, res)
res <- res %>% mutate(scenario = as.factor(scenario),
               stock = qname) %>%
  full_join(stk_labels) %>% 
  mutate(label = as.factor(label),
         label2 = as.factor(label2))
res$label <- factor(res$label, levels = levels(res$label)[c(1:19, 21, 20, 22:29)])
res$label2 <- factor(res$label2, levels = levels(res$label2)[c(1:19, 21, 20, 22:29)])


### save
saveRDS(object = res, 
        file = paste0(path_default, "plot_data_perfect.rds"))

### plot all stocks
res %>% 
  ggplot(aes(x = year - 100, y = data, linetype = scenario)) +
  geom_hline(yintercept = 1, linetype = "dotted", size = 0.2) +
  geom_line() +
  theme_paper +
  facet_wrap(~ label2, labeller = label_parsed) +
  labs(x = "year", y = expression(SSB/B[MSY]), 
       colour = "scenario") +
  geom_vline(xintercept = 0.5, colour = "darkgrey")


### ------------------------------------------------------------------------ ###
### individual replicates ####
### ------------------------------------------------------------------------ ###
### to show poor performance for higher k stocks

### load MSE results
files <- paste0(path_default, "corrected/", 
                rep(c("one-way", "roller-coaster"), each = 29), 
                "/", stocks$stock, ".rds")
names(files) <- rep(stocks$stock_short, 2)
res_corrected <- lapply(files, readRDS)


tmp_ssb <- as.data.frame(FLQuants(her = res_corrected$her$ssb_rel,
                                  pol = res_corrected$pol$ssb_rel))
tmp_catch <- as.data.frame(FLQuants(her = res_corrected$her$catch_rel,
                                    pol = res_corrected$pol$catch_rel))
tmp <- rbind(cbind(tmp_ssb, quant = "ssb"),
             cbind(tmp_catch, quant = "catch"))
tmp2 <- tmp %>% group_by(year, qname, quant) %>%
  summarise(data = median(data)) %>%
  mutate(iter = "median")
tmp <- bind_rows(tmp, tmp2)
tmp <- tmp %>%
  mutate(stock = qname) %>%
  left_join(stk_labels) %>%
  filter(year >= 100 & iter %in% c("median", 3, 4, 6, 13, 22)) %>%
  spread(quant, data)
tmp <- tmp %>%
  mutate(iter = replace(iter, iter == 3, "replicate~1"),
         iter = replace(iter, iter == 4, "replicate~2"),
         iter = replace(iter, iter == 6, "replicate~3"),
         iter = replace(iter, iter == 13, "replicate~4"),
         iter = replace(iter, iter == 22, "replicate~5"),)

saveRDS(object = tmp, file = paste0(path_default, "plot_data_replicates.rds"))



### ------------------------------------------------------------------------ ###
### supplementary: sensitivity - SSB ####
### ------------------------------------------------------------------------ ###

scns <- 38:40
label <- "steepness"
values <- c(0.75, 0.9, 0.6)

### function for comparing SSB trajectories from different OM sets
SSB_comparison <- function(scns, label, values, stks_subset, 
                           nrow = NULL, ncol = NULL) {
  
  ### stock labels
  stk_labels <- stocks %>%
    select(stock_short, k) %>%
    mutate(stock = stock_short,
           label = paste0("k = ", k, ", ", stock),
           label2 = paste0("italic(k)==", k, "~~", stock))
  
  res <- lapply(seq_along(scns), function(i_scn) {#browser()
    tmp <- lapply(seq_along(stocks$stock), function(i_stock) {
      qnt_i <- readRDS(paste0("output/observation_error/", 
                              OM_scns$id[scns[i_scn]], "/",
                              "corrected/one-way/", stocks$stock[i_stock], ".rds"))
      iterMedians(qnt_i$ssb_rel)
    })
    names(tmp) <- stocks$stock_short
    tmp <- as.data.frame(FLQuants(tmp))
    tmp$scenario <- values[i_scn]
    return(tmp)
  })
  res <- do.call(rbind, res)
  res <- res %>% mutate(scenario = as.factor(scenario),
                 stock = qname) %>%
    full_join(stk_labels) %>% 
    mutate(label = as.factor(label),
           label2 = as.factor(label2))
  res$label <- factor(res$label, levels = levels(res$label)[c(1:19, 21, 20, 22:29)])
  res$label2 <- factor(res$label2, levels = levels(res$label2)[c(1:19, 21, 20, 22:29)])
  
  ### subset stocks
  if (!missing(stks_subset)) {
    res <- res %>% filter(stock %in% stks_subset)
  }
  
  ### plot  
  p <- res %>% 
    ggplot(aes(x = year - 100, y = data, colour = scenario)) +
    geom_line() +
    theme_bw() +
    facet_wrap(~ label2, labeller = label_parsed, nrow = nrow, ncol = ncol) +
    labs(x = "year", y = expression("median SSB/B"["MSY"]), 
         colour = label) +
    geom_vline(xintercept = 0.5, colour = "darkgrey")
  return(p)
}

### h levels
SSB_comparison(scns = 38:40, label = "steepness", values = c(0.75, 0.9, 0.6))
ggsave(filename = paste0("output/plots/paper_revision/supplementary_material/",
                         "SSB_h_levels.png"),
       width = 20, height = 15, units = "cm", dpi = 300, type = "cairo-png")

### h functional relationships
SSB_comparison(scns = c(38, 41, 42), label = "steepness", 
               values = c("h=0.75", "h~k", "h~L50/Linf"))
ggsave(filename = paste0("output/plots/paper_revision/supplementary_material/",
                         "SSB_h_functional.png"),
       width = 20, height = 15, units = "cm", dpi = 300, type = "cairo-png")

### h from Myers
SSB_comparison(scns = c(38, 43), label = "steepness", nrow = 3,
               values = c(0.75, "Myers et\nal. (1999)"),
               stks_subset = c("ang3", "smn", "ang", "ang2", "pol", "had",
                               "sbb", "ple", "whg", "lem", "ane", "sar", "her"))
ggsave(filename = paste0("output/plots/paper_revision/supplementary_material/",
                         "SSB_h_Myers.png"),
       width = 20, height = 10, units = "cm", dpi = 300, type = "cairo-png")

### recruitment sd
SSB_comparison(scns = c(38, 29, 44), label = "recruitment\nvariability\n(sd)", 
               values = c("0.6", "0.3", "0.9"))
ggsave(filename = paste0("output/plots/paper_revision/supplementary_material/",
                         "SSB_rec_sd.png"),
       width = 20, height = 15, units = "cm", dpi = 300, type = "cairo-png")

### index sd
SSB_comparison(scns = c(38, 45, 46), label = "index\nuncertainty\n(sd)", 
               values = c(0.2, 0.4, 0.6))
ggsave(filename = paste0("output/plots/paper_revision/supplementary_material/",
                         "SSB_idx_sd.png"),
       width = 20, height = 15, units = "cm", dpi = 300, type = "cairo-png")

### length frequency sd
SSB_comparison(scns = c(38, 47, 48), 
               label = "length\nfrequency\nuncertainty\n(sd)", 
               values = c(0.2, 0.4, 0.6))
ggsave(filename = paste0("output/plots/paper_revision/supplementary_material/",
                         "SSB_lngth_sd.png"),
       width = 20, height = 15, units = "cm", dpi = 300, type = "cairo-png")

### length frequency and index sd
SSB_comparison(scns = c(38, 49), 
               label = "index and\nlength\nfrequency\nuncertainty\n(sd)", 
               values = c(0.2, 0.4))
ggsave(filename = paste0("output/plots/paper_revision/supplementary_material/",
                         "SSB_idx_lngth_sd.png"),
       width = 20, height = 15, units = "cm", dpi = 300, type = "cairo-png")

### selectivity
SSB_comparison(scns = c(38, 50, 52), 
               label = "selectivity", 
               values = c("default", "before\nmaturity",
                          "after\nmaturity"))
ggsave(filename = paste0("output/plots/paper_revision/supplementary_material/",
                         "SSB_selectivity.png"),
       width = 20, height = 15, units = "cm", dpi = 300, type = "cairo-png")

### selectivity, the 2nd
SSB_comparison(scns = c(38, 54:57), sort = TRUE,
               label = "selectivity", 
               values = c("default", "age +1", "age +2", "age +3", "age +4"))
ggsave(filename = paste0("output/plots/paper_revision/supplementary_material/",
                         "SSB_selectivity_shift.png"),
       width = 20, height = 15, units = "cm", dpi = 300, type = "cairo-png")



### ------------------------------------------------------------------------ ###
### supplementary: sensitivity - stats ####
### ------------------------------------------------------------------------ ###

### load stats from scenarios
stats_all <- lapply(c(29, c(38:52, 54:57)), function(x) {
  OM_scn <- OM_scns[OM_scns$idSEQ == x, ]
  tmp <- read.csv(paste0("output/",
                         ifelse(OM_scn$obs_error, "observation_error",
                                "perfect_knowledge"), "/",
                         OM_scn$id, "/",
                         "stats.csv"))
  tmp$id <- OM_scn$id
  tmp$idSEQ <- OM_scn$idSEQ
  if (x == 43) tmp <- tmp %>% filter(paper %in% c("her", "pol", "smn", "lem",
                                                  "ple", "whg", "had", "ang",
                                                  "ang2", "sar", "sbb", "ane",
                                                  "ang3"))
  return(tmp)
})
stats_all <- do.call(rbind, stats_all)
stats_all <- stats_all %>%
  select(f_rel, ssb_rel, risk_collapse, risk_blim, yield_rel, iav,
         id, paper, fhist) %>%
  rename("F/F[MSY]" = f_rel, 
         "SSB/B[MSY]" = ssb_rel, 
         "collapse~risk" = risk_collapse, 
         "B[lim]~risk" = risk_blim, 
         "catch/MSY" = yield_rel, 
         "ICV" = iav) %>%
  gather(key = "stat", value = "value", 1:6)
stats_all$stat <- as.factor(stats_all$stat)
stats_all$stat <- factor(stats_all$stat,
                         levels = levels(stats_all$stat)[c(4, 6, 3, 1, 2, 5)])
#levels(stats_all$stat)
stats_all <- stats_all %>% spread(id, value)

p_base <- stats_all %>%
  filter(fhist == "one-way") %>%
  ggplot() +
  stat_function(fun = function(x){x}, colour = "darkgrey") +
  facet_wrap(~ stat, scales = "free", nrow = 1, 
             labeller = label_parsed) +
  theme_bw() +
  xlim(0, NA) + ylim(0, NA)
ph0.9 <- p_base + 
  geom_point(aes(x = new_baseline, y = rev_h0.9), size = 0.5) +
  labs(x = "default,  steepness h = 0.75", y = "\nh = 0.9")
ph0.6 <- p_base + 
  geom_point(aes(x = new_baseline, y = rev_h0.6), size = 0.5) +
  labs(x = "default,  steepness h = 0.75", y = "\nh = 0.6")
ph_k <- p_base + 
  geom_point(aes(x = new_baseline, y = rev_h_k), size = 0.5) +
  labs(x = "default,  steepness h = 0.75", y = "\nh ~ k")
ph_l50linf <- p_base + 
  geom_point(aes(x = new_baseline, y = rev_h_l50linf), size = 0.5) +
  labs(x = "default,  steepness h = 0.75", y = "\nh ~ L50/Linf")
ph_Myers <- p_base + 
  geom_point(aes(x = new_baseline, y = rev_h_Myers), size = 0.5) +
  labs(x = "default,  steepness h = 0.75", y = "\nh ~ Myers et al. (1999)")

p_recSD0.3 <- p_base + 
  geom_point(aes(x = new_baseline, y = test_idx0.2), size = 0.5) +
  labs(x = "default, recruitment sd = 0.6", y = "\nrecruitment sd = 0.3")
p_recSD0.9 <- p_base + 
  geom_point(aes(x = new_baseline, y = rev_sR0.9), size = 0.5) +
  labs(x = "default, recruitment sd = 0.6", y = "\nrecruitment sd = 0.9")

p_idxSD0.4 <- p_base + 
  geom_point(aes(x = new_baseline, y = rev_idxSD0.4), size = 0.5) +
  labs(x = "default, index sd = 0.2", y = "\nindex sd = 0.4")
p_idxSD0.6 <- p_base + 
  geom_point(aes(x = new_baseline, y = rev_idxSD0.6), size = 0.5) +
  labs(x = "default, index sd = 0.2", y = "\nindex sd = 0.6")
p_lengthSD0.4 <- p_base + 
  geom_point(aes(x = new_baseline, y = rev_lengthSD0.4), size = 0.5) +
  labs(x = "default, length sd = 0.2", y = "\nlength sd = 0.4")
p_lengthSD0.6 <- p_base + 
  geom_point(aes(x = new_baseline, y = rev_lengthSD0.6), size = 0.5) +
  labs(x = "default, length sd = 0.2", y = "\nlength sd = 0.6")
p_idxSD0.4_lengthSD0.4 <- p_base + 
  geom_point(aes(x = new_baseline, y = rev_idxSD0.4_lengthSD0.4), size = 0.5) +
  labs(x = "default, length sd = 0.2, index sd = 0.2", 
       y = "length sd = 0.4,\nindex sd = 0.4")

p_sel_before <- p_base + 
  geom_point(aes(x = new_baseline, y = rev_sel_before), size = 0.5) +
  labs(x = "default selectivity", y = "selectivity\nbefore maturity")
p_sel_after <- p_base + 
  geom_point(aes(x = new_baseline, y = rev_sel_after), size = 0.5) +
  labs(x = "default selectivity", y = "selectivity\nafter maturity")

p_sel_1 <- p_base + 
  geom_point(aes(x = new_baseline, y = rev_sel_1), size = 0.5) +
  labs(x = "default selectivity", y = "selectivity\nage +1")
p_sel_2 <- p_base + 
  geom_point(aes(x = new_baseline, y = rev_sel_2), size = 0.5) +
  labs(x = "default selectivity", y = "selectivity\nage +2")
p_sel_3 <- p_base + 
  geom_point(aes(x = new_baseline, y = rev_sel_3), size = 0.5) +
  labs(x = "default selectivity", y = "selectivity\nage +3")
p_sel_4 <- p_base + 
  geom_point(aes(x = new_baseline, y = rev_sel_4), size = 0.5) +
  labs(x = "default selectivity", y = "selectivity\nage +4")

### recruitment sensitivity
plot_grid(ph0.9, ph0.6, ph_l50linf, ph_k, ph_Myers, p_recSD0.3, p_recSD0.9,
          align = "vh", ncol = 1, labels = "auto")
ggsave(filename = paste0("output/plots/paper_revision/supplementary_material/",
                         "stats_steepness.png"),
       width = 25, height = 35, units = "cm", dpi = 300, type = "cairo-png")
### survey uncertainty
plot_grid(p_idxSD0.4, p_idxSD0.6, p_lengthSD0.4, p_lengthSD0.6, 
          p_idxSD0.4_lengthSD0.4,
          align = "vh", ncol = 1, labels = "auto")
ggsave(filename = paste0("output/plots/paper_revision/supplementary_material/",
                         "stats_uncertainty.png"),
       width = 25, height = 25, units = "cm", dpi = 300, type = "cairo-png")
### selectivity
plot_grid(p_sel_before, p_sel_after,
          align = "vh", ncol = 1, labels = "auto")
ggsave(filename = paste0("output/plots/paper_revision/supplementary_material/",
                         "stats_selectivity.png"),
       width = 25, height = 10, units = "cm", dpi = 300, type = "cairo-png")
### selectivity, he 2nd
plot_grid(p_sel_1, p_sel_2, p_sel_3, p_sel_4,
          align = "vh", ncol = 1, labels = "auto")
ggsave(filename = paste0("output/plots/paper_revision/supplementary_material/",
                         "stats_selectivity_shift.png"),
       width = 25, height = 20, units = "cm", dpi = 300, type = "cairo-png")

