library(FLCore)
library(ggplotFL)
library(doParallel)
library(dplyr)
library(tidyr)


### parallelization
workers <- makeCluster(detectCores())
registerDoParallel(workers)


### ------------------------------------------------------------------------ ###
### correlation f <-> stock ####
### ------------------------------------------------------------------------ ###

res_df <- readRDS("output/stats_scn_new.RDS") ### quants
refpts <- readRDS("input/refpts.rds")


### find & sort scenarios
scns <- res_df[res_df$scenario %in% 6515:6572, ]
scns <- scns[order(scns$fhist), ]


### load/extract quants
qts <- foreach(scenario = scns$scenario, 
               .packages = c("FLCore", "foreach"),
               .export = c("scns"), .errorhandling = "pass") %dopar% {
  #browser()
  ### load "corrected" quants
  qts_tmp <- readRDS(paste0("output/perfect_knowledge/combined/3.2.1_quants/",
                           scenario, ".rds"))
 
  ### "loop" through iterations
  tmp <- foreach(i = 1:500, .packages = "FLCore", .combine = rbind) %do% {
    
    ### subset to iteration
    ### and year range
    qts_i <- window(FLCore::iter(qts_tmp, i), end = 100)
    
    ### correlation test
    ### first for inverted F, then SSB
    cor_f <- cor.test(c(qts_i$HCR3.2.1f), 1/c(qts_i$fbar))
    cor_SSB <- cor.test(c(qts_i$HCR3.2.1f), c(qts_i$SSB))
    
    ### data frame with results
    data.frame(cor_f = cor_f$estimate,
               p_f = cor_f$p.value,
               cor_SSB = cor_SSB$estimate,
               p_SSB = cor_SSB$p.value,
               iter = i)
    
  }
  
  cbind(tmp, scenario = scenario, 
        stock = scns$stock[scns$scenario == scenario],
        fhist = scns$fhist[scns$scenario == scenario])
  
}
names(qts) <- scns$stock

lapply(lapply(qts, "[[", "cor_f"), summary)

qts_df <- do.call(rbind, qts)
qts_df <- gather(qts_df, key, value, cor_f:p_SSB) %>%
  mutate(stat = ifelse(grepl(x = key, pattern = "^cor_"), 
                        "correlation", "p"),
         quant = ifelse(grepl(x = key, pattern = "*_SSB$"), 
                        "f ~ SSB", "f ~ 1/fbar"))


### plot
ggplot(data = qts_df,
       aes(x = stock, y = value)) +
  facet_grid(stat ~ quant, scales = "free_y") +
  geom_boxplot(size = 0.3, outlier.size = 0.5) +
  labs(y = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/correlations/",
                         "component_f_vs_stock_quant_correlation.png"),
       width = 15, height = 15, units = "cm", dpi = 300, type = "cairo-png")


### ------------------------------------------------------------------------ ###
### age structure - one-way vs. roller-coaster ####
### ------------------------------------------------------------------------ ###

### load/extract quants
qts <- foreach(scenario = scns$scenario, 
               .packages = c("FLCore"),
               .errorhandling = "pass") %dopar% {

  ### load stock
  tmp <- readRDS(paste0("output/perfect_knowledge/combined/",
                scenario, ".rds"))
  iterMedians(stock.n(tmp))
}
names(qts) <- paste(scns$stock, scns$fhist, sep = "__")
### convert into data frame
qts_df <- as.data.frame(FLQuants(qts))
### add descriptors
qts_df$stock <- sapply(strsplit(x = as.character(qts_df$qname), split = "__"), 
                       "[[", 1)
qts_df$fhist <- sapply(strsplit(x = as.character(qts_df$qname), split = "__"), 
                       "[[", 2)

### plot age structure at beginning of historical period
qts_df %>% 
  filter(year == 100) %>%
  ggplot(aes(x = age, y = data, fill = fhist, colour = fhist)) +
  #geom_bar(stat = "identity", position = "dodge") + 
  geom_line(size = 0.5) + geom_point(size = 0.5) +
  facet_wrap(~ stock, scales = "free") +
  theme_bw(base_size = 8) +
  labs(y = "stock numbers at age") +
  xlim(0, NA) + ylim(0, NA)
ggsave(filename = paste0("input/OMs_age_structure_year100.png"),
       width = 20, height = 15, units = "cm", dpi = 300, type = "cairo-png")
