### ------------------------------------------------------------------------ ###
### cluster analysis of 3.2.1 results ####
### ------------------------------------------------------------------------ ###

### load packages
library(FLCore)
library(ggplotFL)
library(Cairo)
library(tidyr)
library(dtwclust)
library(ggdendro)


### set up parallel environment
library(parallel)
workers <- makeCluster(detectCores())
invisible(clusterEvalQ(workers, {
  library(dtwclust)
  RcppParallel::setThreadOptions(1L)
}))
library(doParallel)
registerDoParallel(workers)

# registerDoSEQ()


### ------------------------------------------------------------------------ ###
### load data ####
### ------------------------------------------------------------------------ ###

### scenario definitions
res_df <- readRDS("output/stats_scn_new.RDS")
### reference points
refpts <- readRDS("input/refpts.rds")

### find & sort scenarios
scns <- res_df[res_df$scenario %in% 6515:6572, ]
scns <- scns[order(scns$fhist), ]

# ### load/extract quants
# qts <- foreach(scenario = scns$scenario, .packages = c("FLCore"),
#                .export = c("scns", "refpts"), .errorhandling = "pass") %dopar% {
#   #browser()
#   ### load "corrected" quants
#   qts_tmp <- readRDS(paste0("output/perfect_knowledge/combined/3.2.1_quants/",
#                            scenario, ".rds"))
#   # qts_tmp <- readRDS(paste0("/gpfs/afmcefas/simonf/output/combined/",
#   #                           "3.2.1_quants/", scenario, ".rds"))
#   
#   ### load stock
#   stk_tmp <- readRDS(paste0("output/perfect_knowledge/combined/",
#                            scenario, ".rds"))
#   # stk_tmp <- readRDS(paste0("/gpfs/afmcefas/simonf/output/combined/",
#   #                           scenario, ".rds"))
#   
#   ### refpts
#   refpts_tmp <- refpts[[as.character(scns$stock[scns$scenario == scenario])]]
#   
#   ### calucalate values relative to MSY and combine all quants
#   list(### "corrected" values
#    SSB = qts_tmp$SSB, fbar = qts_tmp$fbar,
#    Catch = qts_tmp$Catch, Rec = qts_tmp$Rec,
#    ### original values
#    SSB_raw = ssb(stk_tmp), fbar_raw = fbar(stk_tmp),
#    Catch_raw = catch(stk_tmp), Rec_raw = rec(stk_tmp),
#    ### "corrected" values, relative to MSY
#    SSB_rel = qts_tmp$SSB / c(refpts_tmp["msy", "ssb"]),
#    fbar_rel = qts_tmp$fbar / c(refpts_tmp["msy", "harvest"]),
#    Catch_rel = qts_tmp$Catch / c(refpts_tmp["msy", "yield"]),
#    Rec_rel = qts_tmp$Rec / c(refpts_tmp["msy", "rec"]),
#    ### original values, relative to MSY
#    SSB_raw_rel = ssb(stk_tmp) / c(refpts_tmp["msy", "ssb"]),
#    fbar_raw_rel = fbar(stk_tmp) / c(refpts_tmp["msy", "harvest"]),
#    Catch_raw_rel = catch(stk_tmp) / c(refpts_tmp["msy", "yield"]),
#    Rec_raw_rel = rec(stk_tmp) / c(refpts_tmp["msy", "rec"])
#   )
#   
#   }
# names(qts) <- scns$stock
# ### save
# saveRDS(qts, file = "output/perfect_knowledge/clustering/def_scns.rds")
### load
qts <- readRDS("output/perfect_knowledge/clustering/def_scns.rds")

### extract SSBs
# cl_SSB <- lapply(qts, function(x) {
#   as.data.frame(t(x[["SSB"]][, drop = TRUE]))
# })
# cl_SSB_raw <- lapply(qts, function(x) {
#   as.data.frame(t(x[["SSB_raw"]][, drop = TRUE]))
# })
# cl_SSB_rel <- lapply(qts, function(x) {
#   as.data.frame(t(x[["SSB_rel"]][, drop = TRUE]))
# })
# cl_SSB_raw_rel <- lapply(qts, function(x) {
#   as.data.frame(t(x[["SSB_raw_rel"]][, drop = TRUE]))
# })
### extract medians
# cl_SSB_rel_median <- lapply(qts, function(x) {
#   as.data.frame(t(iterMedians(x[["SSB_rel"]])[, drop = TRUE]))
# })



### create list with input objects for cluster analysis
### extract required SSB quants
cl_input <- c("SSB", "SSB_raw", "SSB_rel", "SSB_raw_rel")
names(cl_input) <- unlist(cl_input)
cl_input <- lapply(cl_input, function(x) {
  tmp <- lapply(qts, "[[", x)
  names(tmp) <- names(qts)
  tmp
})
### split into fishing scenarios
cl_input <- list("one-way" = lapply(cl_input, "[", 1:29),
                  "roller-coaster" = lapply(cl_input, "[", 30:58))
### add iter medians
cl_input <- list("iter" = cl_input,
                  "median" = lapply(cl_input, function(x) {
                    lapply(x, function(y) {
                      lapply(y, iterMedians)
                    })
                  }))
### convert FLQuants into matrices
cl_input <- lapply(cl_input, function(x) {
  lapply(x, function(y) {
    lapply(y, function(z) {
      lapply(z, function(w) {
        as.data.frame(t(w[, drop = TRUE]))
      })
    })
  })
})

### ------------------------------------------------------------------------ ###
### hierarchical cluster analysis for all list elements ####
### ------------------------------------------------------------------------ ###

### "loop" through elements
res <- foreach(stat = names(cl_input), 
               .packages = c("FLCore", "dtwclust", "foreach", "tidyr", "dplyr",
                             "ggplot2"),
               .export = "cl_input") %:%
  foreach(fhist = names(cl_input$iter)) %:% 
  foreach(quant = names(cl_input$iter$`one-way`), 
          .errorhandling = "stop") %dopar% {
    #browser()
    
    ### -------------------------------------------------------------------- ###
    ### run cluster analysis ####
    ### -------------------------------------------------------------------- ###
    cluster <- tsclust(cl_input[[stat]][[fhist]][[quant]], 
                       type = "hierarchical",
                       k = 10, distance = "dtw", seed = 1,
                       trace = TRUE, error.check = TRUE)
    
    ### -------------------------------------------------------------------- ###
    ### plot dendrogram ####
    ### -------------------------------------------------------------------- ###
    
    png(filename = paste0("output/perfect_knowledge/plots/3.2.1/cluster/",
                          "hierarchical/", stat, "_", fhist, "_", quant,
                          "_dendrogram.png"),
        width = 25, height = 20, units = "cm", res = 300)
    plot(cluster, main = paste(stat, fhist, quant, sep = " - "))
    dev.off()
    
    ### -------------------------------------------------------------------- ###
    ### plot cluster SSB trends ####
    ### -------------------------------------------------------------------- ###
    
    ### find all possible clusters and allocations of stocks to them
    cl_alloc <- cutree(cluster, k = 1:29)
    
    ### stock SSB medians
    stock_medians <- cluster@datalist
    stock_medians <- lapply(stock_medians, function(x) {
      apply(x, 2, median)
    })
    stock_medians <- as.data.frame(do.call(rbind, stock_medians))
    stock_medians$stock <- names(cluster@datalist)
    
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
        ### calculate median
        apply(tmp [, ac(75:200)], MARGIN = 2, median)
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
    df_plot <- gather(df_plot, key = "year", value = "value", "75":"200")
    df_plot$year <- as.numeric(df_plot$year)
    df_plot$k <- as.numeric(df_plot$k)
    ### plot (in junks)
    for (junks in split(1:29, ceiling(1:29/5))) {
      ggplot(df_plot[df_plot$k %in% junks, ], 
             aes(x = year, y = value, group = stock, 
                 linetype = stock == "cluster", colour = stock == "cluster")) +
        geom_line() + 
        facet_grid(k ~ cluster) + theme_bw() +
        scale_linetype_manual("data", values = c("dotted", "solid"), 
                              labels = c("stock", "cluster")) +
        scale_colour_manual("data", values = c("grey", "black"), 
                            labels = c("stock", "cluster")) +
        theme(panel.grid = element_blank()) +
        labs(y = "SSB/SSBmsy", title = paste(stat, fhist, quant, sep = " - "))
      ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/cluster/",
                               "hierarchical/", stat, "_", fhist, "_", quant,
                               "_centroids_", head(junks, 1), "-", 
                               tail(junks, 1), ".png"),
             width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
    }
    
    return(cluster)
    
}

### ------------------------------------------------------------------------ ###
### same again, but now only from year 101 onwards ####
### ------------------------------------------------------------------------ ###

### remove early years
### convert FLQuants into matrices
cl_input2 <- lapply(cl_input, function(x) {
  lapply(x, function(y) {
    lapply(y, function(z) {
      lapply(z, function(w) {
        w[!names(w) %in% ac(75:100)]
      })
    })
  })
})


### "loop" through elements
res2 <- foreach(stat = names(cl_input2), 
               .packages = c("FLCore", "dtwclust", "foreach", "tidyr", "dplyr",
                             "ggplot2"),
               .export = "cl_input2") %:%
  foreach(fhist = names(cl_input2$iter)) %:% 
  foreach(quant = names(cl_input2$iter$`one-way`), 
          .errorhandling = "stop") %dopar% {
  #browser()
  
  ### -------------------------------------------------------------------- ###
  ### run cluster analysis ####
  ### -------------------------------------------------------------------- ###
  cluster <- tsclust(cl_input2[[stat]][[fhist]][[quant]], 
                     type = "hierarchical",
                     k = 10, distance = "dtw", seed = 1,
                     trace = TRUE, error.check = TRUE)
  
  ### -------------------------------------------------------------------- ###
  ### plot dendrogram ####
  ### -------------------------------------------------------------------- ###
  
  png(filename = paste0("output/perfect_knowledge/plots/3.2.1/cluster/",
                        "hierarchical/excl_history/", stat, "_", fhist, "_",
                        quant, "_dendrogram.png"),
      width = 25, height = 20, units = "cm", res = 300)
  plot(cluster, main = paste(stat, fhist, quant, sep = " - "))
  dev.off()
  
  ### -------------------------------------------------------------------- ###
  ### plot cluster SSB trends ####
  ### -------------------------------------------------------------------- ###
  
  ### find all possible clusters and allocations of stocks to them
  cl_alloc <- cutree(cluster, k = 1:29)
  
  ### stock SSB medians
  stock_medians <- cluster@datalist
  stock_medians <- lapply(stock_medians, function(x) {
    apply(x, 2, median)
  })
  stock_medians <- as.data.frame(do.call(rbind, stock_medians))
  stock_medians$stock <- names(cluster@datalist)
  
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
      ### calculate median
      apply(tmp[, ac(101:200)], MARGIN = 2, median)
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
  df_plot <- gather(df_plot, key = "year", value = "value", "101":"200")
  df_plot$year <- as.numeric(df_plot$year)
  df_plot$k <- as.numeric(df_plot$k)
  ### plot (in junks)
  for (junks in split(1:29, ceiling(1:29/5))) {
    ggplot(df_plot[df_plot$k %in% junks, ], 
           aes(x = year, y = value, group = stock, 
               linetype = stock == "cluster", colour = stock == "cluster")) +
      geom_line() + 
      facet_grid(k ~ cluster) + theme_bw() +
      scale_linetype_manual("data", values = c("dotted", "solid"), 
                            labels = c("stock", "cluster")) +
      scale_colour_manual("data", values = c("grey", "black"), 
                          labels = c("stock", "cluster")) +
      theme(panel.grid = element_blank()) +
      labs(y = "SSB/SSBmsy", title = paste(stat, fhist, quant, sep = " - "))
    ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/cluster/",
                             "hierarchical/excl_history/", stat, "_", fhist, 
                             "_", quant, "_centroids_", head(junks, 1), "-", 
                             tail(junks, 1), ".png"),
           width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")
  }
  
  return(cluster)
  
}

### ------------------------------------------------------------------------ ###
### final clustering ####
### ------------------------------------------------------------------------ ###
### hierarchical
### use median SSB relative to Bmsy
### one-way & roller-coaster
cl_one_way <- tsclust(cl_input[["median"]][["one-way"]][["SSB_rel"]], 
                      type = "hierarchical",
                      k = 10, distance = "dtw", seed = 1,
                      trace = TRUE, error.check = TRUE)
cl_roller_coaster <- tsclust(cl_input[["median"]][["roller-coaster"]][["SSB_rel"]], 
                      type = "hierarchical",
                      k = 10, distance = "dtw", seed = 1,
                      trace = TRUE, error.check = TRUE)

### select 4 clusters: get cluster allocation for stocks
cl_alloc_ow <- as.data.frame(cutree(cl_one_way, k = 1:29))
cl_alloc_rc <- as.data.frame(cutree(cl_roller_coaster, k = 1:29))
cl_alloc_ow$stock <- rownames(cl_alloc_ow)
cl_alloc_rc$stock <- rownames(cl_alloc_rc)

### save
saveRDS(cl_one_way, file = "output/clusters_final_ow.rds")
saveRDS(cl_roller_coaster, file = "output/clusters_final_rc.rds")

### save cluster allocations
saveRDS(object = cl_alloc_ow, file = "output/cluster_allocations_one_way.rds")
saveRDS(object = cl_alloc_rc, file = "output/cluster_allocations_roller_coaster.rds")

### groups 3 & 4 replaced in one-way compared to roller-coaster
### -> replace
cl_alloc_ow4 <- cl_alloc_ow[, c("stock", "4")]
cl_alloc_rc4 <- cl_alloc_rc[, c("stock", "4")]
cl_alloc_rc4$`4` <- sapply(cl_alloc_rc4$`4`, function(x) {
  switch(x, "1" = 1, "2" = 2, "3" = 4, "4" = 3)
})
res <- data.frame(stock = cl_alloc_ow$stock,
           ow = cl_alloc_ow$`4`,
           rc = sapply(cl_alloc_rc$`4`, function(x) {
             switch(x, "1" = 1, "2" = 2, "3" = 4, "4" = 3)
           }))
res$is <- res$ow == res$rc
sum(res$is)
### 21 of 29 have same allocation

### ------------------------------------------------------------------------ ###
### final clustering: compare with lhist and check odd ones ####
### ------------------------------------------------------------------------ ###


### check clusters, compare with lhist
lhist <- readRDS("input/lhist_extended.rds")
allocs <- cl_alloc_ow[, 3, drop = FALSE]
allocs$stock <- rownames(allocs)
rownames(allocs) <- NULL
allocs <- merge(allocs, lhist[, c("stock", "k")])
allocs <- allocs[order(allocs$`3`, allocs$k), ]

### get scenario definitions
source("MP_scenarios.R")
scn_df <- scn_df[scn_df$scenario %in% 6515:6572 & scn_df$fhist == "one-way", ]

### merge scenario into allocations
allocs <- merge(allocs, scn_df[, c("stock", "scenario")])

### cluster 3 only (cluster 2/middle in plots)
allocs2 <- allocs[allocs$`3` == "3", ]

### load quants
allocs_qnts <- foreach(scenario = allocs2$scenario) %do% {
  qnts <- readRDS("output/perfect_knowledge/combined/3.2.1_quants/6522.rds")
}

### ------------------------------------------------------------------------ ###
### final clustering - plot clusters and k ####
### ------------------------------------------------------------------------ ###
library(cowplot)

### find all possible clusters and allocations of stocks to them
cl_alloc <- cutree(cl_one_way, k = 1:29)

### stock SSB medians
stock_medians <- cl_one_way@datalist
stock_medians <- lapply(stock_medians, function(x) {
  apply(x, 2, median)
})
stock_medians <- as.data.frame(do.call(rbind, stock_medians))
stock_medians$stock <- names(cl_one_way@datalist)

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
df_plot <- gather(df_plot, key = "year", value = "value", "75":"200")
df_plot$year <- as.numeric(df_plot$year)
df_plot$k <- as.numeric(df_plot$k)

### load lhist - k
lhist <- readRDS("input/lhist_extended.rds")
lhist_k <- lhist[, c("stock", "k")]
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
#cl_alloc_k$cluster <- paste("cluster ", cl_alloc_k$cluster)
#cl_alloc_k$n_cluster <- paste(cl_alloc_k$n_cluster, "clusters")

### data for dendrogram
dat_dend <- dendro_data(cl_one_way)
### replace some names
dat_dend$labels$label <- as.character(dat_dend$labels$label)
dat_dend$labels$label[c(2, 7, 8, 11, 12, 13)] <- c("bll", "sar", "jnd", "wlf",
                                                   "gut", "sbb")

### plot
p1 <- ggplot(df_plot[df_plot$k %in% 1:4, ], 
       aes(x = year, y = value, group = stock, 
           linetype = stock == "cluster", alpha = stock == "cluster",
           size = stock == "cluster", colour = as.factor(cluster))) +
  geom_line() + 
  facet_grid(ifelse(k > 1, 
                    paste(k, "clusters"), 
                    paste(k, "cluster")) ~ paste("cluster:", cluster)) +
  theme_bw() +
  scale_linetype_manual("", values = c("dotted", "solid"), 
                        labels = c("stock", "cluster")) +
  scale_alpha_manual("", values = c(0.7, 1), 
                     labels = c("stock", "cluster")) +
  scale_size_manual("", values = c(0.25, 1), 
                    labels = c("stock", "cluster")) +
  scale_colour_discrete(guide = FALSE) +
  theme(panel.grid = element_blank()) +
  labs(y = "SSB/Bmsy") +
  theme(legend.position = c(0.87, 0.87),
        legend.background = element_blank(),
        legend.key = element_blank(),
        strip.text.y = element_blank())

p2 <- ggplot(data = cl_alloc_k[cl_alloc_k$n_cluster %in% 1:4, ],
       aes(x = k_pos, y = k, fill = as.factor(cluster))) +
  geom_bar(stat = "identity", colour = "black", size = 0.1) +
  # scale_fill_manual("cluster", values = c("grey1", "white", 
  #                                         "grey33", "grey67")) +
  scale_fill_discrete("cluster") +
  #facet_wrap(~ n_cluster, ncol = 1, strip.position = "right") +
  facet_grid(ifelse(n_cluster > 1, 
                    paste(n_cluster, "clusters"), 
                    paste(n_cluster, "cluster")) ~ label) + 
  theme_bw() + labs(x = "stocks") +
  theme(axis.text.x = element_text(colour = "white"),
        axis.ticks.x = element_blank())

p3 <- ggdendrogram(dat_dend) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text.y = element_blank(), axis.title.x = element_blank(),
        #axis.title.x = element_text(colour = "black"),
        axis.title.y = element_text(colour = "black", angle = 90, hjust = -2),
        axis.text.x = element_text(size = 9),
        text = element_text(size = 11)) +
  labs(x = "stocks", y = "DTW distance")

plot_grid(p3, plot_grid(p1, p2, ncol = 2, rel_widths = c(1.7, 1), 
                        labels = c("B", "C")),
          nrow = 2, rel_heights = c(1, 2), labels = c("A", ""))

ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/cluster/",
                         "hierarchical/summary_one-way_coloured.png"),
       width = 20, height = 16, units = "cm", dpi = 300, type = "cairo")


### same for roller-coaster
### find all possible clusters and allocations of stocks to them
cl_alloc <- cutree(cl_roller_coaster, k = 1:29)

### stock SSB medians
stock_medians <- cl_roller_coaster@datalist
stock_medians <- lapply(stock_medians, function(x) {
  apply(x, 2, median)
})
stock_medians <- as.data.frame(do.call(rbind, stock_medians))
stock_medians$stock <- names(cl_one_way@datalist)

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
df_plot <- gather(df_plot, key = "year", value = "value", "75":"200")
df_plot$year <- as.numeric(df_plot$year)
df_plot$k <- as.numeric(df_plot$k)

### merge k and cluster allocations
cl_alloc_k <- as.data.frame(cl_alloc)
cl_alloc_k$stock <- rownames(cl_alloc_k)
cl_alloc_k <- merge(cl_alloc_k, lhist_k)
cl_alloc_k <- gather(cl_alloc_k, key = "n_cluster", value = "cluster",
                     `1`:`29`)
cl_alloc_k$column = "k"
#cl_alloc_k$cluster <- paste("cluster ", cl_alloc_k$cluster)
#cl_alloc_k$n_cluster <- paste(cl_alloc_k$n_cluster, "clusters")
### plot
p1 <- ggplot(df_plot[df_plot$k %in% 1:4, ], 
             aes(x = year, y = value, group = stock, 
                 linetype = stock == "cluster", colour = stock == "cluster",
                 size = stock == "cluster")) +
  geom_line() + 
  facet_grid(ifelse(k > 1, 
                    paste(k, "clusters"), 
                    paste(k, "cluster")) ~ paste("cluster:", cluster)) +
  theme_bw() +
  scale_linetype_manual("", values = c("dotted", "solid"), 
                        labels = c("stock", "cluster")) +
  scale_size_manual("", values = c(0.25, 1), 
                    labels = c("stock", "cluster")) +
  scale_colour_manual("", values = c("grey40", "black"), 
                      labels = c("stock", "cluster")) +
  theme(panel.grid = element_blank()) +
  labs(y = "SSB/Bmsy") +
  theme(legend.position = c(0.87, 0.87),
        legend.background = element_blank(),
        legend.key = element_blank(),
        strip.text.y = element_blank())

p2 <- ggplot(data = cl_alloc_k[cl_alloc_k$n_cluster %in% 1:4, ],
             aes(x = k_pos, y = k, fill = as.factor(cluster))) +
  geom_bar(stat = "identity", colour = "black", size = 0.1) +
  scale_fill_manual("cluster", values = c("grey1", "white", 
                                          "grey33", "grey67")) +
  #facet_wrap(~ n_cluster, ncol = 1, strip.position = "right") +
  facet_grid(ifelse(n_cluster > 1, 
                    paste(n_cluster, "clusters"), 
                    paste(n_cluster, "cluster")) ~ column) + 
  theme_bw() + labs(x = "stocks") +
  theme(axis.text.x = element_text(colour = "white"),
        axis.ticks.x = element_blank())

plot_grid(p1, p2, ncol = 2, rel_widths = c(1.7, 1))
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/cluster/",
                         "hierarchical/summary_roller-coaster.png"),
       width = 20, height = 14, units = "cm", dpi = 300, type = "cairo")

### ------------------------------------------------------------------------ ###
### dendrogram ####
### ------------------------------------------------------------------------ ###


library(latticeExtra)
### distance matrix
dat_dist <- cl_one_way@distmat
### dendrogramm
dat_dend <- as.dendrogram(cl_one_way)

### get order from dendrogram
dat_order <- order.dendrogram(dat_dend)

# png(filename = "output/perfect_knowledge/clustering/dist_dend.png", 
#     width = 18, height = 18, units = "cm", res = 300, type = "cairo")
levelplot(
  dat_dist[(dat_order), (dat_order)],
  aspect = "fill",
  scales = list(x = list(rot = 90)),
  colorkey = list(space = "left"),
  #col.regions = colorRamps::green2red(29),
  #col.regions = gray(seq(from = 1, to = 0, length.out = 29)),
  col.regions = (grDevices::heat.colors(29)),
  xlab = NULL, ylab = NULL, 
  legend =
    list(
      top =
        list(fun = dendrogramGrob,
             args =
               list(x = dat_dend,
                    side = "top",
                    size = 5)))
)


### ------------------------------------------------------------------------ ###
### hierarchical cluster: rel SSB with iterations ####
### ------------------------------------------------------------------------ ###
c2_iter <- tsclust(cl_SSB_rel[1:29], type = "hierarchical",
                   k = 10, distance = "dtw", seed = 1,
                   trace = TRUE, error.check = TRUE)
### raw not corrected SSBs (i.e. stocks recover after collapse)
# c2_iter_raw <- tsclust(cl_SSB_raw_rel[1:29], type = "hierarchical",
#                    k = 10, distance = "dtw", seed = 1,
#                    trace = TRUE, error.check = TRUE)
plot(c2_iter)
ggdendrogram(c2_iter) + #coord_cartesian(ylim = c(20000, 100000)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
png(filename = "output/perfect_knowledge/clustering/hierarchical_iter.png",
    width = 20, height = 15, units = "cm", res = 300, type = "cairo")
plot(c2_iter)
dev.off()


dist_mat <- as.data.frame(c2_iter@distmat)
dist_mat$stock <- rownames(dist_mat)
dist_mat %>% gather(key = "stock2", value = "value", 1:29) %>%
  ggplot(aes(x = stock, y = stock2, fill = value)) + theme_bw() +
  geom_raster() + 
  labs(x = "", y = "") + scale_fill_continuous("distance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25))
ggsave(filename = "output/perfect_knowledge/clustering/hierarchical_dist.png",
       width = 20, height = 15, units = "cm", dpi = 600, type = "cairo")



### find all possible clusters and allocations of stocks to them
cl_alloc <- cutree(c2_iter, k = 1:29)

### stock SSB medians
stock_medians <- do.call(rbind, cl_SSB_rel_median)
stock_medians$stock <- rownames(stock_medians)
stock_medians <- stock_medians[1:29, ]

### extract medians and create "centroids"
df_plot <- foreach(k = colnames(cl_alloc)) %do% {
  ### load SSB medians per stock
  res <- stock_medians
  ### allocate to cluster
  #identical(res$stock, rownames(cl_alloc))
  res$cluster <- cl_alloc[, k]
  ### cluster median
  ### go through each cluster
  res_add <- lapply(split(res$stock, res$cluster), function(stocks) {
    ### extract SSBs 
    tmp <- lapply(qts[1:29][stocks], function(x) {
      as.data.frame(t(x[["SSB_rel"]][, drop = TRUE]))
    })
    ### combine into single table
    tmp <- do.call(rbind, tmp)
    ### calculate median
    apply(tmp, MARGIN = 2, median)
  })
  ### combine
  res_add_cluster <- names(res_add)
  res_add <- as.data.frame(do.call(rbind, res_add))
  res_add$stock <- "cluster"
  res_add$cluster <- res_add_cluster
  ### add 
  res <- rbind(res, res_add)
  res$cluster <- as.numeric(res$cluster)
  ### add number of clusters
  res$k = k
  return(res)
}
df_plot <- do.call(rbind, df_plot)
df_plot <- gather(df_plot, key = "year", value = "value", "75":"200")
df_plot$year <- as.numeric(df_plot$year)
df_plot$k <- as.numeric(df_plot$k)
### plot
ggplot(df_plot[df_plot$k %in% 25:29, ], 
       aes(x = year, y = value, group = stock, linetype = stock == "cluster",
           colour = stock == "cluster")) +
  geom_line() + 
  facet_grid(k ~ cluster) + theme_bw() +
  scale_linetype_manual("data", values = c("dotted", "solid"), 
                        labels = c("stock", "cluster")) +
  scale_colour_manual("data", values = c("grey", "black"), 
                      labels = c("stock", "cluster")) +
  theme(panel.grid = element_blank()) +
  labs(y = "SSB/SSBmsy")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/cluster/iter/", 
                         "hierarchical_cluster_25-29.png"),
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")


### ------------------------------------------------------------------------ ###
### plot lhist grouped by cluster ####
### ------------------------------------------------------------------------ ###
### load params
pars <- readRDS("input/all_stocks_repfts_lhpar.rds")
### format for plotting
pars$stock <- as.character(pars$stock)
pars$stock <- factor(pars$stock, 
                     levels = pars$stock[order(-(pars$K))])
#pars$group <- as.factor(pars$group)
# pars_df <- gather(pars, key = "parameter", value = "value",
#                   c("L_inf", "K", "a50", "M_mat", "M", "MK", "t0", "a", "b", 
#                     "max_age", "LFeFmsy", "LFeMK", "F_MSY", "L_c"))
# pars_df$input <- "input"
# pars_df$input[pars_df$parameter %in% c("max_age", "LFeFmsy", "LFeMK",
#                                        "F_MSY", "L_c", "M", 
#                                        "M_mat", "MK")] <- "derived"
# pars_df$input <- factor(pars_df$input, levels = unique(pars_df$input))
# pars_df$parameter <- factor(pars_df$parameter,
#                             levels = unique(pars_df$parameter))

### merge clusters and params
cl_alloc <- cbind(as.data.frame(cl_alloc), stock = rownames(cl_alloc))
df_plot <- merge(pars, cl_alloc)
### replicate content for plotting
df_plot <- gather(df_plot, key = "n_cluster", value = "cluster", "1":"29")
df_plot$cluster <- as.factor(df_plot$cluster)
df_plot$n_cluster <- as.factor(as.numeric(df_plot$n_cluster))

### plot
ggplot(df_plot[df_plot$n_cluster %in% as.character(1:10), ], 
       aes(x = stock, y = K, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_wrap(~ n_cluster, nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/cluster/iter/",
                         "/hierarchical_K_cl17-29.png"),
       width = 30, height = 18, units = "cm", dpi = 300, type = "cairo-png")


### ------------------------------------------------------------------------ ###
### hierarchical cluster: rel SSB - medians only ####
### ------------------------------------------------------------------------ ###
c2_median <- tsclust(cl_SSB_rel_median[1:29], type = "hierarchical",
                     k = 10, distance = "dtw", seed = 1,
                     trace = TRUE, error.check = TRUE)

plot(c2_median)
png(filename = paste0("output/perfect_knowledge/plots/3.2.1/cluster/median/",
                      "hierarchical.png"),
    width = 20, height = 15, units = "cm", res = 300, type = "cairo")
plot(c2_median)
dev.off()


dist_mat <- as.data.frame(c2_median@distmat)
dist_mat$stock <- rownames(dist_mat)
dist_mat %>% gather(key = "stock2", value = "value", 1:29) %>%
  ggplot(aes(x = stock, y = stock2, fill = value)) + theme_bw() +
  geom_raster() + 
  labs(x = "", y = "") + scale_fill_continuous("distance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25))
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/cluster/median/",
                         "hierarchical_dist.png"),
       width = 20, height = 15, units = "cm", dpi = 600, type = "cairo")

### find all possible clusters and allocations of stocks to them
cl_alloc <- cutree(c2_median, k = 1:29)

### extract medians and create "centroids"
df_plot <- foreach(k = colnames(cl_alloc)) %do% {
  ### load SSB medians per stock
  res <- stock_medians
  ### allocate to cluster
  #identical(res$stock, rownames(cl_alloc))
  res$cluster <- cl_alloc[, k]
  ### cluster median
  ### go through each cluster
  res_add <- lapply(split(res$stock, res$cluster), function(stocks) {
    ### extract SSBs 
    tmp <- lapply(qts[1:29][stocks], function(x) {
      as.data.frame(t(x[["SSB_rel"]][, drop = TRUE]))
    })
    ### combine into single table
    tmp <- do.call(rbind, tmp)
    ### calculate median
    apply(tmp, MARGIN = 2, median)
  })
  ### combine
  res_add_cluster <- names(res_add)
  res_add <- as.data.frame(do.call(rbind, res_add))
  res_add$stock <- "cluster"
  res_add$cluster <- res_add_cluster
  ### add 
  res <- rbind(res, res_add)
  res$cluster <- as.numeric(res$cluster)
  ### add number of clusters
  res$k = k
  return(res)
}
df_plot <- do.call(rbind, df_plot)
df_plot <- gather(df_plot, key = "year", value = "value", "75":"200")
df_plot$year <- as.numeric(df_plot$year)
df_plot$k <- as.numeric(df_plot$k)
### plot
ggplot(df_plot[df_plot$k %in% 25:29, ], 
       aes(x = year, y = value, group = stock, linetype = stock == "cluster",
           colour = stock == "cluster")) +
  geom_line() + 
  facet_grid(k ~ cluster) + theme_bw() +
  scale_linetype_manual("data", values = c("dotted", "solid"), 
                        labels = c("stock", "cluster")) +
  scale_colour_manual("data", values = c("grey", "black"), 
                      labels = c("stock", "cluster")) +
  theme(panel.grid = element_blank()) +
  labs(y = "SSB/SSBmsy")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/cluster/median/", 
                         "hierarchical_cluster_25-29.png"),
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo")

### plot lhist grouped by cluster ####
### merge clusters and params
cl_alloc <- cbind(as.data.frame(cl_alloc), stock = rownames(cl_alloc))
df_plot <- merge(pars, cl_alloc)
### replicate content for plotting
df_plot <- gather(df_plot, key = "n_cluster", value = "cluster", "1":"29")
df_plot$cluster <- as.factor(df_plot$cluster)
df_plot$n_cluster <- as.factor(as.numeric(df_plot$n_cluster))

### plot
ggplot(df_plot[df_plot$n_cluster %in% as.character(1:10), ], 
       aes(x = stock, y = K, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_wrap(~ n_cluster, nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/cluster/median/",
                         "/hierarchical_K_cl17-29.png"),
       width = 30, height = 18, units = "cm", dpi = 300, type = "cairo-png")




### ------------------------------------------------------------------------ ###
### new plot, idea from Laurie
### plot distance matrix and add dendrogram
library(latticeExtra)
### distance matrix
dat_dist <- c2_median@distmat
### dendrogramm
dat_dend <- as.dendrogram(c2_median)

### get order from dendrogram
dat_order <- order.dendrogram(dat_dend)

png(filename = "output/perfect_knowledge/clustering/dist_dend.png", 
    width = 18, height = 18, units = "cm", res = 300, type = "cairo")
levelplot(
  dat_dist[(dat_order), (dat_order)],
  aspect = "fill",
  scales = list(x = list(rot = 90)),
  colorkey = list(space = "left"),
  #col.regions = colorRamps::green2red(29),
  #col.regions = gray(seq(from = 1, to = 0, length.out = 29)),
  col.regions = (grDevices::heat.colors(29)),
  xlab = NULL, ylab = NULL, 
  legend =
    list(
      top =
        list(fun = dendrogramGrob,
             args =
               list(x = dat_dend,
                    side = "top",
                    size = 5)))
)
dev.off()





### ------------------------------------------------------------------------ ###
### dead ends from here ####
### ------------------------------------------------------------------------ ###


# 
# install.packages("dtwclust")
# library(dtwclust)
# data("uciCT")
# 
# require("TSclust")
# proxy::pr_DB$set_entry(FUN = diss.ACF, names = c("ACFD"),
#                        loop = TRUE, type = "metric", distance = TRUE,
#                        description = "Autocorrelation-based distance")
# 
# proxy::dist(CharTraj[3:8], method = "ACFD", upper = TRUE)
# 
# library(TSdist)
# # Register the Fourier distance
# proxy::pr_DB$set_entry(FUN = FourierDistance, names = c("fourier"),
#                        loop = TRUE, type = "metric", distance = TRUE,
#                        description = "Distance with Fourier coefficients")
# # Fourier distance requires equal length
# data <- reinterpolate(CharTraj, new.length = max(lengths(CharTraj)))
# # Partitional clustering
# tsclust(data[1L:10L], k = 2L, distance = "fourier", seed = 838)
# 
# 
# 
# # Reinterpolate to same length
# data <- reinterpolate(CharTraj, new.length = max(lengths(CharTraj)))
# # Calculate the DTW distances between all elements
# system.time(D1 <- proxy::dist(data[1L:5L], data[6L:100L],
#                               method = "dtw_basic",
#                               window.size = 20L))
# 
# # Nearest neighbors
# NN1 <- apply(D1, 1L, which.min)
# # Calculate the distance matrix with dtw_lb
# system.time(D2 <- dtw_lb(data[1L:5L], data[6L:100L],
#                          window.size = 20L))
# 
# # Nearest neighbors
# NN2 <- apply(D2, 1L, which.min)
# # Same results?
# all(NN1 == NN2)
# 
# # Exclude a series as an example
# database <- data[-100L]
# classify_series <- function(query) {
#   # Nearest neighbor
#   nn <- which.min(dtw_lb(database, query,
#                          window.size = 18L,
#                          nn.margin = 2L))
#   # Return a label
#   CharTrajLabels[nn]
# }
# # 100-th series is a Z character
# classify_series(data[100L])
# 
# hc_sbd <- tsclust(CharTraj, type = "h", k = 20L,
#                   preproc = zscore, seed = 899,
#                   distance = "sbd", centroid = shape_extraction,
#                   control = hierarchical_control(method = "average"))
# # By default, the dendrogram is plotted in hierarchical clustering
# plot(hc_sbd)
# # The series and the obtained prototypes can be plotted too
# plot(hc_sbd, type = "sc")
# 
# # Focusing on the first cluster
# plot(hc_sbd, type = "series", clus = 1L)
# plot(hc_sbd, type = "centroids", clus = 1L)
# 
# 
# 
# 
# # Calculate autocorrelation up to 50th lag
# acf_fun <- function(dat, ...) {
#   lapply(dat, function(x) {
#     as.numeric(acf(x, lag.max = 50, plot = FALSE)$acf)
#   })
# }
# # Fuzzy c-means
# fc <- tsclust(CharTraj[1:20], type = "f", k = 4L,
#               preproc = acf_fun, distance = "L2",
#               seed = 42)
# # Fuzzy membership matrix
# fc@fcluster
# plot(fc, series = CharTraj[1:20], type = "series")
# 
# 


library(FLCore)
library(ggplotFL)
library(Cairo)
library(tidyr)
library(dtwclust)

library(parallel)
workers <- makeCluster(detectCores())
invisible(clusterEvalQ(workers, {
  library(dtwclust)
  RcppParallel::setThreadOptions(1L)
}))
library(doParallel)
registerDoParallel(workers)

# registerDoSEQ()

### ------------------------------------------------------------------------ ###
### group ####
### ------------------------------------------------------------------------ ###

res_df <- readRDS("output/stats_scn_new.RDS") ### quants
refpts <- readRDS("input/refpts.rds")


### find scenarios
scns <- res_df[res_df$scenario %in% 6515:6572 &
               res_df$fhist == "one-way", ]
# ### load SSBs
# ssbs <- foreach(scenario = scns$scenario, .packages = c("FLCore"), 
#                   .export = c("res_df"), .errorhandling = "pass") %do% {
#   ### load quants
#   qts <- readRDS(paste0("output/perfect_knowledge/combined/3.2.1_quants/",
#                         scenario, ".rds"))[["SSB"]]
#   iterMedians(qts)
# }
# names(ssbs) <- scns$stock
# ### ------------------------------------------------------------------------ ###
# ### ssb relative to Bmsy ####
# ### load refpts
# refpts <- readRDS("input/refpts.rds")
### extract Bmsy
Bmsy <- lapply(refpts, function(x) {
  c(x["msy", "ssb"])
})
# ### calculate SSB/Bmsy
# SSBmsy <- lapply(names(ssbs), function(x) {
#   c(ssbs[[x]]) / c(Bmsy[[x]])
# })
# names(SSBmsy) <- names(ssbs)
# 
# 
# cluster <- tsclust(SSBmsy, 
#                        type = "fuzzy", k = 5, distance = "dtw",
#                        seed = 42)
# plot(cluster)
# cluster
# 
# ssb_cluster <- tsclust(lapply(ssbs, c), 
#                        type = "fuzzy", k = 5, distance = "dtw",
#                        seed = 42)
# ssb_cluster
# plot(ssb_cluster)
# 
# 
# ### find numbers of clusters...
# clusters <- tsclust(SSBmsy, 
#                    type = "fuzzy", k = 2:10, distance = "dtw",
#                    seed = 42)
# names(clusters) <- paste0("k_", 2:10)
# sapply(clusters, cvi, type = "internal")
# 


### try iteration-wise
ssbs <- foreach(scenario = scns$scenario, .packages = c("FLCore"), 
                .export = c("res_df"), .errorhandling = "pass") %do% {
    ### load quants
    # qts <- readRDS(paste0("output/perfect_knowledge/combined/3.2.1_quants/",
    #                       scenario, ".rds"))[["SSB"]]
    qts <- readRDS(paste0("/gpfs/afmcefas/simonf/output/combined/3.2.1_quants/",
                          scenario, ".rds"))[["SSB"]]
    as.data.frame(t(qts[, drop = TRUE]))
}
names(ssbs) <- scns$stock[seq(length(ssbs))]
SSBmsy <- lapply(names(ssbs), function(x) {
  ssbs[[x]] / c(Bmsy[[x]])
})
names(SSBmsy) <- names(ssbs)
#SSBmsy <- do.call(rbind, SSBmsy)

### ------------------------------------------------------------------------ ###
### see if there are different clusters within stocks
### ------------------------------------------------------------------------ ###
### number of repetitions
nrep = 100
### number of clusters
k = 2:5


### perform cluster analysis for all stock
system.time({
cls_stks <- foreach(ssb = SSBmsy, .errorhandling = "pass") %do% {
  tsclust(ssb, 
          type = "partitional", k = k, distance = "dtw",
          seed = 1, control = partitional_control(nrep = nrep))
}
})
### check if all runs succeded
unlist(lapply(cls_stks, is))
names(cls_stks) <- names(SSBmsy)
### save results
saveRDS(cls_stks, 
        file = "output/perfect_knowledge/clustering/dtw_stocks_single.rds")
#cls_stks <- readRDS("output/perfect_knowledge/clustering/dtw_stocks_single.rds")


### return to serial execution
stopCluster(workers)
registerDoSEQ()
### calculate statistics: cluster validity indices
cl_stats <- foreach(stock_i = cls_stks) %:%
  foreach(part_i = seq(length(k) * nrep), .combine = rbind) %do% {
    
    ### calculate
    res <- cvi(stock_i[[part_i]])
    ### format and add description
    res <- data.frame(t(res),
                      k = stock_i[[part_i]]@k)

}
### add stock names
cl_stats <- lapply(seq_along(cl_stats), function(x) {
  cbind(cl_stats[[x]], names(ssbs)[x])
})
### combine all stocks
cl_stats <- do.call(rbind, cl_stats)
### save
saveRDS(cl_stats, 
       file = "output/perfect_knowledge/clustering/dtw_stocks_single_stats.rds")
# cl_stats <- readRDS(paste0("output/perfect_knowledge/clustering/",
#                            "dtw_stocks_single_stats.rds"))

names(cl_stats)[9] <- "stock"
cl_stats$run <- rep(seq(nrep), length(k))

### format for plotting
df_plot <- gather(cl_stats, key = "index", value = "value", -k, -stock, -run)
#df_plot$run <- seq(nrow(cl_stats))

. <- foreach(stock = unique(df_plot$stock), .packages = "ggplot2") %do% {
  ### plot
  ggplot(data = df_plot[df_plot$stock == stock, ]) +
    geom_line(aes(x = k, y = value, group = run),
              colour = "grey", size = 0.2) +
    stat_summary(aes(x = k, y = value),
                 fun.y = "median", geom = "line", colour = "black",
                 size = 1) +
    facet_wrap(~ index, scales = "free_y") +
    theme_bw() +
    labs(x = "number of clusters", y = "index", title = stock)
  ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/cluster/",
                           "single/", stock, ".png"),
         width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
  
}


### ------------------------------------------------------------------------ ###
### cluster stocks - one matrix for each stock ####
### ------------------------------------------------------------------------ ###

### cluster
cl_stks <- tsclust(SSBmsy, type = "partitional",
                   k = 2:15, distance = "dtw", seed = 1)
sapply(cl_stks, cvi)
saveRDS(cl_stks, file = "output/perfect_knowledge/clustering/stks_mt_cl.rds")
cl_stks <- readRDS("output/perfect_knowledge/clustering/tmp.rds")


### cluster with absolute SSB values
cl_stks2 <- tsclust(ssbs, type = "partitional",
                   k = 2:15, distance = "dtw", seed = 1)
sapply(cl_stks2, cvi)


table(cl_stks[[2]]@cluster)
ct <- lapply(cl_stks[[2]]@centroids, function(x) {
  apply(x, 2, median)
})
ct <- t(do.call(rbind, ct))
ct <- gather(data.frame(ct))
ct$year <- 75:200
ggplot(ct, aes(x = year, y = value, colour = key)) +
  geom_line()


blubb <- dist(SSBmsy[[1]], method = "dtw")
hc <- hclust(blubb, method="average")

### ------------------------------------------------------------------------ ###
### cluster stocks - median SSBs ####
### ------------------------------------------------------------------------ ###

### calculate median
SSBmedian <- lapply(SSBmsy, function(x) {
  apply(x, 2, median)
})

### cluster
cl_median <- tsclust(SSBmedian, type = "partitional",
                     k = 2:15, distance = "dtw", seed = 1)
sapply(cl_median, cvi)


### ------------------------------------------------------------------------ ###
### cluster stocks - each iteration independently ####
### ------------------------------------------------------------------------ ###


clusters <- 2
system.time({
  ssb_cluster <- tsclust(do.call(rbind, SSBmsy), 
                         type = "partitional", k = clusters, distance = "dtw",
                         seed = 1)
})

stats2 <- sapply(ssb_cluster2, cvi)
stats2 <- as.data.frame(stats2)
names(stats2) <- clusters
stats2$stat <- row.names(stats2)

stats2 <- gather(stats2, key = "n_cluster", value = "value", -"stat")
stats2$n_cluster <- as.numeric(stats2$n_cluster)
ggplot(stats2, aes(x = n_cluster, y = value)) +
  geom_line() +
  facet_wrap(~ stat, scales = "free_y")
sapply(ssb_cluster2, cvi, b = rep(1:3, each = 500))


ssb_cluster2 <- tsclust(SSBmsy, 
                       type = "partitional", k = 3, distance = "dtw",
                       seed = 42)
ssb_cluster
plot(ssb_cluster)
res <- as.data.frame(ssb_cluster@fcluster)
plot(res$cluster_1)


### ------------------------------------------------------------------------ ###
### load data ####
### ------------------------------------------------------------------------ ###

### scenario definitions
res_df <- readRDS("output/stats_scn_new.RDS")
### reference points
refpts <- readRDS("input/refpts.rds")

### find & sort scenarios
scns <- res_df[res_df$scenario %in% 6515:6572, ]
scns <- scns[order(scns$fhist), ]

### load/extract quants
qts <- foreach(scenario = scns$scenario, .packages = c("FLCore"), 
               .export = c("scns", "refpts"), .errorhandling = "pass") %dopar% {
  #browser()
  ### load "corrected" quants
  qts_tmp <- readRDS(paste0("output/perfect_knowledge/combined/3.2.1_quants/",
                            scenario, ".rds"))
  # qts_tmp <- readRDS(paste0("/gpfs/afmcefas/simonf/output/combined/",
  #                           "3.2.1_quants/", scenario, ".rds"))
  
  ### load stock
  stk_tmp <- readRDS(paste0("output/perfect_knowledge/combined/",
                            scenario, ".rds"))
  # stk_tmp <- readRDS(paste0("/gpfs/afmcefas/simonf/output/combined/",
  #                           scenario, ".rds"))
  
  ### refpts
  refpts_tmp <- refpts[[as.character(scns$stock[scns$scenario == scenario])]]
  
  ### calucalate values relative to MSY and combine all quants
  list(### "corrected" values 
       SSB = qts_tmp$SSB, fbar = qts_tmp$fbar, 
       Catch = qts_tmp$Catch, Rec = qts_tmp$Rec,
       ### original values 
       SSB_raw = ssb(stk_tmp), fbar_raw = fbar(stk_tmp), 
       Catch_raw = catch(stk_tmp), Rec_raw = rec(stk_tmp),
       ### "corrected" values, relative to MSY
       SSB_rel = qts_tmp$SSB / c(refpts_tmp["msy", "ssb"]), 
       fbar_rel = qts_tmp$fbar / c(refpts_tmp["msy", "harvest"]), 
       Catch_rel = qts_tmp$Catch / c(refpts_tmp["msy", "yield"]), 
       Rec_rel = qts_tmp$Rec / c(refpts_tmp["msy", "rec"]),
       ### original values, relative to MSY
       SSB_raw_rel = ssb(stk_tmp) / c(refpts_tmp["msy", "ssb"]), 
       fbar_raw_rel = fbar(stk_tmp) / c(refpts_tmp["msy", "harvest"]), 
       Catch_raw_rel = catch(stk_tmp) / c(refpts_tmp["msy", "yield"]), 
       Rec_raw_rel = rec(stk_tmp) / c(refpts_tmp["msy", "rec"])
       )
  
}
names(qts) <- scns$stock
### save
saveRDS(qts, file = "output/perfect_knowledge/clustering/def_scns.rds")
qts <- readRDS("output/perfect_knowledge/clustering/def_scns.rds")

### extract SSBs
cl_SSB <- lapply(qts, function(x) {
  as.data.frame(t(x[["SSB"]][, drop = TRUE]))
})
cl_SSB_raw <- lapply(qts, function(x) {
  as.data.frame(t(x[["SSB_raw"]][, drop = TRUE]))
})
cl_SSB_rel <- lapply(qts, function(x) {
  as.data.frame(t(x[["SSB_rel"]][, drop = TRUE]))
})
cl_SSB_raw_rel <- lapply(qts, function(x) {
  as.data.frame(t(x[["SSB_raw_rel"]][, drop = TRUE]))
})

### ------------------------------------------------------------------------ ###
### cluster - median SSBs ####
### ------------------------------------------------------------------------ ###
### extract medians
cl_SSB_rel_median <- lapply(qts, function(x) {
  as.data.frame(t(iterMedians(x[["SSB_rel"]])[, drop = TRUE]))
})

### cluster
### loop through ks, otherwise the computation will fail, bug in dtwclust...
cluster <- lapply(2:10, function(x) {
  tsclust(cl_SSB_rel_median[1:29], type = "hierarchical",
                   k = x, distance = "dtw", seed = 1,
                   control = partitional_control(nrep = 1000),
                   trace = TRUE, error.check = TRUE)
})
#cluster <- unlist(cluster)
plot(cluster[[3]])
clusterX <- lapply(2:10, function(x) {
  tsclust(cl_SSB_rel_median[1:29], type = "hierarchical",
          k = x, distance = "dtw", seed = 1,
          control = partitional_control(nrep = 1000),
          trace = TRUE, error.check = TRUE)
})
plot(c2)
c2 <- tsclust(cl_SSB_rel_median[1:29], type = "hierarchical",
              k = 2, distance = "dtw", seed = 1,
              trace = TRUE, error.check = TRUE)
c2_iter <- tsclust(cl_SSB_rel[1:29], type = "hierarchical",
                   k = 10, distance = "dtw", seed = 1,
                   trace = TRUE, error.check = TRUE)
c2_iter2 <- tsclust(cl_SSB_rel[30:58], type = "hierarchical",
                   k = 10, distance = "dtw", seed = 1,
                   trace = TRUE, error.check = TRUE)

# c2_iter2 <- tsclust(cl_SSB_rel[1:29], type = "hierarchical",
#                    k = 10, distance = "dtw", seed = 2,
#                    trace = TRUE, error.check = TRUE)
### check centroids
identical(c2_iter@centroids$zeus_faber, c2_iter@datalist$zeus_faber)
### cluster allocations
cutree(tree = c2_iter, k = 5)
### plot distance matrix
dist_mat <- as.data.frame(c2_iter@distmat)
dist_mat$stock <- rownames(dist_mat)
dist_mat %>% gather(key = "stock2", value = "value", 1:29) %>%
  ggplot(aes(x = stock, y = stock2, fill = value)) + theme_bw() +
  geom_raster() + 
  labs(x = "", y = "") + scale_fill_continuous("distance") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25))

c3 <- tsclust(cl_SSB_rel_median[1:29], type = "hierarchical",
              k = 3, distance = "dtw", seed = 1,
              trace = TRUE, error.check = TRUE)
cX <- tsclust(cl_SSB_rel_median[1:29], type = "hierarchical",
              k = 2:20, distance = "dtw", seed = 1,
              trace = TRUE, error.check = TRUE)
plot(c(cX[[5]]@centroids[[1]]), type = "l")
lines(c(cX[[5]]@centroids[[2]]))
lines(c(cX[[5]]@centroids[[3]]))
lines(c(cX[[5]]@centroids[[4]]))
lines(c(cX[[5]]@centroids[[5]]))
lines(c(cX[[5]]@centroids[[6]]))

sapply(cX, cvi)

### calculate statistics: cluster validity indices
cl_stats <- foreach(k = cluster) %do% {
  res <- foreach(i = k, .combine = rbind,
                    .packages = "dtwclust") %dopar% {
    ### calculate
    res <- cvi(i)
    ### format and add description
    res <- data.frame(t(res), k = i@k)
  }
  names(res) <- c("Sil [max]", "SF [max]", "CH [max]", "DB [min]",
                       "DBstar [min]", "D [max]", "COP [min]", "k")
  return(res)
}
cl_stats2 <- do.call(rbind, cl_stats)

### check number of index results for all indices
lapply(cl_stats, function(x) {
  apply(x, 2, function(y) {
    length(unique(y))
  })
})
hist(cl_stats[[9]]$`Sil [max]`)
which.max(cl_stats2$`Sil [max]`)
plot()

cl_stats2 %>% gather(key = "index", value = "value", 1:7) %>%
  ggplot(aes(x = k, y = value, group = k)) +
  geom_boxplot() +
  facet_wrap(~ index, scale = "free_y")

### only 5 possible results for k=2
### 8 for k = 3
### plot
cl_stats[cl_stats$k == 3, ] %>% 
  gather(key = "index", value = "value", 1:7) %>%
  ggplot(aes(x = value)) + theme_bw() +
  geom_bar() + facet_wrap(~ index, scales = "free")
### clear result:
### 5 of 7 indices give same result, best index in ~ 44.9% of runs
best <- lapply(names(cl_stats)[-8], function(x) {
  if (grepl(x = x, pattern = "\\[max\\]")) 
    return(which(cl_stats[, x] == max(unique(cl_stats[, x]))))
  else if (grepl(x = x, pattern = "\\[min\\]")) 
    return((which(cl_stats[, x] == min(unique(cl_stats[, x])))))
})
names(best) <- names(cl_stats)[-8]
lapply(best, length)
identical(best[[1]], best[[3]])
identical(best[[1]], best[[4]])
identical(best[[1]], best[[5]])
identical(best[[1]], best[[6]])

### extract centroids
centroids <- lapply(cluster[best[[1]]], function(x) {
  tmp <- x@centroids #### extract centroids
  ### sort by first element, because order of centroids can be different
  tmp[order(unlist(lapply(tmp, "[[", 1)))]
})
### check if they are all the same
### by compariing the first element to all other ones
all(lapply(centroids, function(x) {
  all.equal(x, centroids[[1]])
}))
### all identical
head(best$`Sil [max]`)

# 
# cl_SSB_res <- tsclust(cl_SSB[1:29], type = "partitional",
#                    k = 2:10, distance = "dtw", seed = 1)
# cl_SSB_rel_res <- tsclust(cl_SSB_rel[1:29], type = "partitional",
#                       k = 2:10, distance = "dtw", seed = 1)
# cl_SSB_raw_res <- tsclust(cl_SSB_raw[1:29], type = "partitional",
#                       k = 2:10, distance = "dtw", seed = 1)
# cl_SSB_raw_rel_res <- tsclust(cl_SSB_raw_rel[1:29], type = "partitional",
#                       k = 2:10, distance = "dtw", seed = 1)
# 
# res <- list(sapply(cl_SSB_res, cvi),
#             sapply(cl_SSB_rel_res, cvi),
#             sapply(cl_SSB_raw_res, cvi),
#             sapply(cl_SSB_raw_rel_res, cvi))
# 
# ### same for roller-coaster
# cl_SSB_res_ <- tsclust(cl_SSB[30:58], type = "partitional",
#                       k = 2:10, distance = "dtw", seed = 1)
# cl_SSB_rel_res_ <- tsclust(cl_SSB_rel[30:58], type = "partitional",
#                           k = 2:10, distance = "dtw", seed = 1)
# cl_SSB_raw_res_ <- tsclust(cl_SSB_raw[30:58], type = "partitional",
#                           k = 2:10, distance = "dtw", seed = 1)
# cl_SSB_raw_rel_res_ <- tsclust(cl_SSB_raw_rel[30:58], type = "partitional",
#                               k = 2:10, distance = "dtw", seed = 1)
# 
# res_ <- list(sapply(cl_SSB_res_, cvi),
#             sapply(cl_SSB_rel_res_, cvi),
#             sapply(cl_SSB_raw_res_, cvi),
#             sapply(cl_SSB_raw_rel_res_, cvi))

### ------------------------------------------------------------------------ ###
### cluster - use relative SSB ####
### ------------------------------------------------------------------------ ###

### number of clusters
clusters <- 2:26
### number of repetitions
#nrep <- 50

### cluster with 100 repetitions
### split into 3 part due to bug in tsclust...
### not more than 41 repetitions possible, otherwise error message
### seems to be related with the creation of the random seed length in tsclust()
cl_large1 <- tsclust(cl_SSB_rel[1:29], type = "partitional",
                     k = clusters, distance = "dtw", seed = 1,
                     control = partitional_control(nrep = 20),
                     trace = TRUE, error.check = TRUE)
cl_large2 <- tsclust(cl_SSB_rel[1:29], type = "partitional",
                     k = clusters, distance = "dtw", seed = 2,
                     control = partitional_control(nrep = 40),
                     trace = TRUE, error.check = TRUE)
cl_large3 <- tsclust(cl_SSB_rel[1:29], type = "partitional",
                     k = clusters, distance = "dtw", seed = 3,
                     control = partitional_control(nrep = 40),
                     trace = TRUE, error.check = TRUE)

### combine into single object
cl_large <- unlist(c(cl_large1, cl_large2, cl_large3), recursive = FALSE)
rm(cl_large1, cl_large2, cl_large3)

### save
saveRDS(cl_large, file = paste0("/gpfs/afmcefas/simonf/output/clustering/",
                                "clusters_matrix.rds"))
# cl_large <- readRDS("/gpfs/afmcefas/simonf/output/clustering/clusters_matrix.rds")

### repeat with many repetitions, plot cvis and centroids!

### calculate statistics: cluster validity indices
cl_stats <- foreach(i = cl_large, .combine = rbind) %dopar% {
    
  ### calculate
  res <- cvi(i)
  ### format and add description
  res <- data.frame(t(res), k = i@k)
    
}
### add run numbers
cl_stats$run <- c(rep(1:20, 25), rep(21:60, 25), rep(61:100, 25))

saveRDS(cl_stats, file = paste0("/gpfs/afmcefas/simonf/output/clustering/",
                                "clusters_matrix_stats.rds"))
cl_stats <- readRDS(paste0("output/perfect_knowledge/clustering/",
                           "clusters_matrix_stats.rds"))

### create dummy
# cl_stats <- rbind(cl_stats, cl_stats, cl_stats)
# cl_stats$run <- rep(1:3, each = 9)
# cl_stats$run <- rep(seq(nrep), length(k))

### format for plotting
df_plot <- gather(cl_stats, key = "index", value = "value", -k, -run)
### dummy
df_plot$value <- df_plot$value * rlnorm(nrow(df_plot))

### index levels
# df_plot$index <- factor(df_plot$index, 
#                         levels = c("Sil", "D", "CH", "SF", "COP","DB", 
#                                    "DBstar"),
#                         labels = c("Sil [max]", "D [max]", "CH [max]" , 
#                                    "SF [max]",
#                                    "COP [min]", "DB [min]", "DBstar [min]"))

#. <- foreach(stock = unique(df_plot$stock), .packages = "ggplot2") %do% {
  ### plot
  ggplot(data = df_plot) +
    geom_line(aes(x = k, y = value, group = run),
              colour = "grey", size = 0.2) +
    stat_summary(aes(x = k, y = value),
                 fun.y = "median", geom = "line", colour = "black",
                 size = 1) +
    facet_wrap(~ index, scales = "free_y") +
    theme_bw() +
    labs(x = "number of clusters", y = "index")
  ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/cluster/",
                           "matrix_stats.png"),
         width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
  
#}

### standardize indices
library(dplyr)
blubb <- group_by(df_plot, run, index) %>% 
  summarise(value_max = max(value))
blubb2 <- full_join(df_plot, blubb) %>% mutate(value_rel = value / value_max)

ggplot(data = blubb2) +
  geom_line(aes(x = k, y = value_rel, group = run),
            colour = "grey", size = 0.2) +
  stat_summary(aes(x = k, y = value_rel),
               fun.y = "mean", geom = "line", colour = "black",
               size = 1) +
  facet_wrap(~ index, scales = "free_y") +
  theme_bw() +
  labs(x = "number of clusters", y = "index") +
  ylim(0, NA)
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/cluster/",
                         "matrix_stats_rel.png"),
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")

### ------------------------------------------------------------------------ ###
### centroids

### remove input data to save memory
cl_large_red <- foreach(cluster = cl_large) %dopar% {
  tmp <- cluster
  tmp@datalist <- list()
  tmp
}

cl_large_red <- foreach(cluster = cl_large_red) %dopar% {
  tmp <- cluster
  tmp@centroids <- lapply(tmp@centroids, function(x) {
    apply(x, 2, median)
  })
  tmp
}

saveRDS(cl_large_red, file = paste0("/gpfs/afmcefas/simonf/output/clustering/",
                                "clusters_matrix_reduced.rds"))
cl_large_red <- readRDS(paste0("output/perfect_knowledge/clustering/",
                           "clusters_matrix_reduced.rds"))

### extract centroids (medians)
centroids <- foreach(cluster = cl_large_red, run = seq_along(cl_large_red),
                     .combine = rbind) %dopar% {
  tmp <- data.frame(do.call(rbind, cluster@centroids))
  tmp$cluster <- seq(nrow(tmp)) ### cluster
  tmp$k <- cluster@k ### number of clusters
  tmp$run <- run ### run ID
  tmp
}
### run per k
centroids <- group_by(centroids, k) %>% 
  mutate(run2 = rep(1:100, each = mean(k)))

### format
centroids <- gather(centroids, key = "year", value = "value", -cluster, -k, 
                    -run, -run2)
centroids$year <- as.numeric(gsub(pattern = "X", replacement = "",
                                  x = centroids$year))


### plot centroids
### individually for each tested number of clusters
for (k in unique(centroids$k)) {
  p <- ggplot(data = centroids[centroids$k %in% k, ],
              aes(x = year, y = value, colour = as.factor(cluster))) +
    geom_line() +
    facet_wrap(~ run2) +
    theme_bw() +
    scale_color_discrete("cluster") +
    labs(x = "year", y = "SSB / Bmsy")
  ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/cluster/",
                           "clusters_", k, ".png"), plot = p,
         width = 30, height = 20, units = "cm", dpi = 300, type = "cairo-png")
  
}
### ------------------------------------------------------------------------ ###
### find best clustering for specific number of clusters
res <- foreach(k = unique(cl_stats$k), .combine = rbind,
               .packages = c("tidyr", "ggplot2", "dplyr")) %dopar% {
  tmp <- cl_stats[cl_stats$k == k, ]
  ### find best run (max)
  best_max <- apply(tmp[, c("Sil", "D", "CH", "SF")], 2, which.max)
  ### find best run (min)
  best_min <- apply(tmp[, c("COP", "DB", "DBstar")], 2, which.min)
  ### combine
  res <- c(best_max, best_min)#data.frame(t())
  ### identify run IDs
  res_df <- data.frame(t(tmp$run[res]))
  names(res_df) <- c("Sil", "D", "CH", "SF", "COP", "DB", "DBstar")
  res_df$k <- k
  
  ### plot distribution
  tmp2 <- gather(tmp, key = "key", value = "value", -k, -run)
  p <- ggplot(tmp2, aes(x = value)) +
    geom_histogram(bins = 100) +
    facet_wrap(~ key, scales = "free") +
    theme_bw()
  ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/cluster/",
                           "distribution/clusters_", k, ".png"), plot = p,
         width = 30, height = 20, units = "cm", dpi = 300, type = "cairo-png")
  
  return(res_df)
  
}


### format for plotting
tmp2 <- gather(tmp, key = "key", value = "value", -k, -run)
ggplot(tmp2, aes(x = value)) +
  geom_histogram(bins = 100) +
  facet_wrap(~ key, scales = "free") +
  theme_bw()

tmp3 <- gather(cl_stats, key = "key", value = "value", -k, -run)
ggplot(tmp3, aes(x = value)) +
  geom_histogram(bins = 100) +
  facet_grid(k ~ key, scales = "free") +
  theme_bw()

### ------------------------------------------------------------------------ ###

# 
# 
# library(kml)
# data("epipageShort", package = "kml")
# head(epipageShort)
# #epipageShort[, 3:6] <- imputation(as.matrix(epipageShort[, 3:6]))
# 
# cldSDQ <- cld(epipageShort, timeInData = 3:6)
# 
# #kml(cldSDQ, nbRedraw = 2, toPlot = "both")
# 
# kml(cldSDQ)
# epipageShort$clusters <- getClusters(cldSDQ, 4)
# 
# 
# tmp_inp <- cld(SSBmsy)
# kml(tmp_inp, nbClusters = 2:10)
# try(choice(tmp_inp))
# plot(as.numeric(getClusters(tmp_inp, 3)))
# plotAllCriterion(tmp_inp)