library(FLCore)
library(ggplotFL)
library(Cairo)
library(reshape2)

### set up parallel computing environment
library(doParallel)
cl <- makeCluster(28)
registerDoParallel(cl)
getDoParWorkers()
getDoParName()

### ------------------------------------------------------------------------ ###
### list files ####
### ------------------------------------------------------------------------ ###

path_res <- "output/perfect_knowledge/"
files <- list.files(path_res)

### keep only MSR results
files <- files[grepl(pattern = "^[0-9]{1,}_[0-9]{1,}\\.rds$", x = files)]

### split files name into parts (scenario, part, extension)
files_desc <- strsplit(x = files, split = "_|[.]")

### create data frame with file info
df <- data.frame(scenario = unlist(lapply(files_desc, function(x){x[[1]]})),
                 part = unlist(lapply(files_desc, function(x){x[[2]]})),
                 file_name = files)
### coerce into numeric
df$scenario <- as.numeric(as.character(df$scenario))
df$part <- as.numeric(as.character(df$part))

### sort
df <- df[order(df$scenario, df$part), ]

### show number of parts for each scenario
df2 <- data.frame(table(df$scenario))
names(df2) <- c("scenario", "count")
df2$scenario <- as.numeric(as.character(df2$scenario))

### list with all scenarios
df3 <- data.frame(scenario = 1:180)

### merge
df4 <- merge(df2, df3, all = TRUE)

### show missing scenarios and scenarios with parts<10
df4[is.na(df4$count) | (!is.na(df4$count) & df4$count != 10), ]

# wrong <- c(5,11,17,23,29,35,41,47,53,59,65,71,77,83,89,95,101,107,113,119,125,131,137,143,149,155,161,167,173,179)
# seq(from = 5, to = 180, by = 6)

### ------------------------------------------------------------------------ ###
### combine parts ####
### ------------------------------------------------------------------------ ###


### function that combines list of object over iter dimension
ibind <- function(object){
  
  if (!is.list(object)) stop("object has to be a list")
  
  ### get iterations per list element
  iters_list <- unlist(lapply(object, function(x){dims(x)$iter}))
  
  ### create list with iterations for each part
  iter_end <- cumsum(iters_list)
  iter_start <- iter_end - iters_list + 1
  iters <- lapply(seq_along(iter_start), function(x){
    seq(iter_start[x], iter_end[x])
  })
  
  ### create template with correct iter dimension
  template <- propagate(FLCore::iter(object[[1]], 1), sum(iters_list))
  
  ### fill with values from list
  for (part in seq_along(iters)) {
    
    FLCore::iter(template, iters[[part]]) <- object[[part]]
    
  }
  
  return(template)
  
}

### load parts and combine them 
scenarios <- df2$scenario#df2$scenario[df2$count == 10]
### only 3.2.2 scenarios
#scenarios <- 181:210


res <- foreach(scenario = scenarios, 
        .packages = c("FLCore"), 
        .export = ls()) %dopar% {
  
  ### get parts
  parts_i <- nrow(df[df$scenario == scenario, ])
  stk_list <- lapply(1:parts_i,
    function(part){
      readRDS(paste0(path_res, df$file_name[df$scenario == scenario &
                                            df$part == part]))
    })
  
  ### fill stock
  stk_temp <- ibind(stk_list)
  
  ### do the same for the attributes
  attr(stk_temp, "lhpar") <- ibind(lapply(stk_list, function(x){
    attr(x, "lhpar")
  }))
  attr(stk_temp, "refpts") <- ibind(lapply(stk_list, function(x){
    attr(x, "refpts")
  }))
  attr(stk_temp, "tracking") <- ibind(lapply(stk_list, function(x){
    attr(x, "tracking")
  }))
  attr(stk_temp, "catch_len") <- ibind(lapply(stk_list, function(x){
    attr(x, "catch_len")
  }))
  
  ### set name
  name(stk_temp) <- as.character(scenario)
  
  ### save
  saveRDS(stk_temp, file = paste0(path_res, "combined/", scenario, ".rds"))
  
  ### return stk
  return(stk_temp)

}

### set names to scenario number
names(res) <- scenarios

### save stock list
#res2 <- readRDS(file = paste0("output/perfect_knowledge/0_all.rds"))
#res <- c(res2, res)
#res <- res[order(as.numeric(names(res)))]
#saveRDS(res, file = paste0("output/perfect_knowledge/0_all.rds"))
#res <- readRDS(file = paste0("output/perfect_knowledge/0_all.rds"))



### ------------------------------------------------------------------------ ###
### load results from all available (combined) scenarios ####
### ------------------------------------------------------------------------ ###
### get available files
files <- list.files("output/perfect_knowledge/combined/")
### get scenario numbers
scenarios <- sort(unlist(lapply(files, function(x){
  as.numeric(unlist(strsplit(x, split = "\\."))[1])
})))
### read results
res <- foreach(scenario = scenarios) %dopar% {
  readRDS(paste0("output/perfect_knowledge/combined/", scenario, ".rds"))
}
names(res) <- scenarios

### ------------------------------------------------------------------------ ###
### add results for 3.2.1 f:c, only first part available
### ------------------------------------------------------------------------ ###
# scenarios <- seq(5, 180, by = 6)
# 
# ### load files
# res_add <- foreach(scenario = scenarios, .packages = "FLCore") %dopar% {
#   stk_tmp <- readRDS(paste0("output/perfect_knowledge/", scenario, "_1.rds"))
#   name(stk_tmp) <- as.character(scenario)
#   stk_tmp
# }
# names(res_add) <- scenarios
# 
# ### add
# res <- c(res, res_add)


### ------------------------------------------------------------------------ ###
### plot results for all scenarios ####
### ------------------------------------------------------------------------ ###

#scenarios <- c(181:210, 331:390)
. <- lapply(scenarios, function(x){
  stk_tmp <- readRDS(paste0("output/perfect_knowledge/combined/", x, ".rds"))
  p <- plot(stk_tmp)
  ggsave(filename = paste0("output/perfect_knowledge/plots/", x,
                           ".png"),
         width = 15, height = 15, units = "cm", dpi = 100, type = "cairo-png",
         plot = p)
})


### ------------------------------------------------------------------------ ###
### calculate Blim ####
### ------------------------------------------------------------------------ ###
### definition:
### SSB where recruitment is impaired by 30%

BevHolt <- function(a, b, ssb) {
  return(a * ssb / (b + ssb))
}

### calculate Blim
Blim <- lapply(1:30, function(x){
  (function(){
    ### load input objects
    load(paste0("input/stocks/", x, ".RData"))
    
    #params(sr.om)
    
    ### ssb from 0 to B0
    ssbs <- seq(0, 1000, 0.01)
    ### recruitment with Beverton & Holt
    rec <- BevHolt(a = c(params(sr.om)["a"]), b = c(params(sr.om)["b"]),
                   ssb = ssbs)

    ssbs[tail(which(rec <= c(max(rec) * 0.7)), 1)]
      
  })()
})

### 162.79 for all stocks

### ------------------------------------------------------------------------ ###
### statistics ####
### ------------------------------------------------------------------------ ###
B_lim <- 162.79

# risks_blim <- lapply(seq_along(res), function(i){
#   n_iter <- dim(res[[i]])[6]
#   ### risks of falling below B_lim
#   apply(ssb(res[[i]])[, ac(101:200)] < B_lim, c(2), function(x){
#     length(which(x)) / n_iter
#   })
# })
stats <- foreach(scenario = scenarios, .packages = "FLCore", .export = "B_lim") %dopar% {
#stats <- foreach(stk_tmp = res, .packages = "FLCore", .export = "B_lim") %dopar% {
  #lapply(seq_along(res), function(i){
  ### load stk
  stk_tmp <- readRDS(paste0("output/perfect_knowledge/combined/", scenario, ".rds"))
  res_temp <- list()
  n_iter <- dim(stk_tmp)[6]
  
  ### proportion where SSB was below B_lim
  res_temp$ssb_below_blim_total <- 
    sum(ssb(stk_tmp)[, ac(101:200)] < B_lim) / (n_iter*100)
  ### proportion of iterations where SSB dropped below B_lim
  res_temp$ssb_below_blim_iter <- 
    sum(yearSums(ssb(stk_tmp)[, ac(101:200)] < B_lim) > 0) / n_iter
  
  ### stock collapse = ssb < 1
  res_temp$collapse_total <- 
    sum(ssb(stk_tmp)[, ac(100:200)] < 1) / (n_iter*100)
  ### proportion of iterations with collapse
  res_temp$collapse_iter <- 
    sum(yearSums(ssb(stk_tmp)[, ac(101:200)] < 1) > 0) / n_iter
  
  ### how frequently is max F reached?
  res_temp$fmaxed_total <- sum(fbar(stk_tmp)[, ac(101:200)] == 5) / (n_iter*100)
  ### in how many iterations did this happen?
  res_temp$fmaxed_iter <- sum(yearSums(fbar(stk_tmp)[, ac(101:200)] == 5) > 0)/
                          n_iter
  
  ### yield
  res_temp$yield <- mean(catch(stk_tmp[,ac(101:200)]))
  res_temp$rel_yield <- mean(catch(stk_tmp[,ac(101:200)])) / mean(catch(stk_tmp[,ac(75:100)]))
  ### scenario definition
  #res_temp$scenario <- names(res)[i]
  
  return(as.data.frame(res_temp))
  
}#)
stats <- do.call(rbind, stats)
stats$scenario <- as.numeric(as.character(names(res)))

### combine with results from earlier scenarios
### load existing results
stats_add <- readRDS("output/stats.RDS")
### combine and skip old results, if newer ones exist
stats <- rbind(stats_add[!stats_add$scenario %in% stats$scenario, ],
               stats)

### save
saveRDS(stats, file = "output/stats.RDS")
write.csv(stats, file = "output/stats.csv")

### ------------------------------------------------------------------------ ###
### combine results with scenario definitions ####
### ------------------------------------------------------------------------ ###
stats <- readRDS("output/stats.RDS")
source("MP_scenarios.R")

### merge
res_df <- merge(stats, scn_df, all = TRUE)
### sort
res_df <- res_df[order(res_df$scenario), ]

### save
saveRDS(res_df, file = "output/stats_scn.RDS")
write.csv(res_df, file = "output/stats_scn.csv")
#

res_df <- readRDS("output/stats_scn.RDS")

### duplicate results df for easier plotting
res_df2 <- res_df
res_df2$HCRmult[is.na(res_df2$HCRmult)] <- 1
res_df2$w[is.na(res_df2$w)] <- 1.4
res_df2$b_w[is.na(res_df2$b_w)] <- 1.4

### ------------------------------------------------------------------------ ###
### plot groups of stocks #####
### ------------------------------------------------------------------------ ###

plot_lst <- function(select, name = ""){
  stk_lst <- lapply(select$scenario, function(x){
    readRDS(paste0("output/perfect_knowledge/combined/", x, ".rds"))
  })
  # stk_lst <- res[as.character(select$scenario)]
  names(stk_lst) <- select$stock
  ### median over all iterations
  stk_lst <- lapply(stk_lst, function(x){
    qapply(x, iterMedians)
  })
  stk_lst <- FLStocks(stk_lst)
  p <- plot(stk_lst) + theme_bw() +
    xlab("year")
  #return(p)
  ggsave(filename = paste0("output/perfect_knowledge/plots/", name, ".png"),
         width = 18, height = 15, units = "cm", dpi = 300, type = "cairo-png",
         plot = p)
}

### plot some ...
### 3.2.1
### one-way, single components
plot_lst(res_df[res_df$options == "option_r:a" & res_df$fhist == "one-way", ],
         name = "3.2.1/factor_r.a_one-way")
plot_lst(res_df[res_df$options == "option_r:b" & res_df$fhist == "one-way", ],
         name = "3.2.1/factor_r.b_one-way")
plot_lst(res_df[res_df$options == "option_f:a" & res_df$fhist == "one-way", ],
         name = "3.2.1/factor_f.a_one-way")
plot_lst(res_df[res_df$options == "option_f:b" & res_df$fhist == "one-way", ],
         name = "3.2.1/factor_f.b_one-way")
plot_lst(res_df[res_df$options == "option_f:c" & res_df$fhist == "one-way", ],
         name = "3.2.1/factor_f.c_one-way")
plot_lst(res_df[res_df$options == "option_b:a" & res_df$fhist == "one-way", ],
         name = "3.2.1/factor_b.a_one-way")
### roller-coaster
plot_lst(res_df[res_df$options == "option_r:a" & res_df$fhist == "roller-coaster", ],
         name = "3.2.1/factor_r.a_roller-coaster")
plot_lst(res_df[res_df$options == "option_r:b" & res_df$fhist == "roller-coaster", ],
         name = "3.2.1/factor_r.b_roller-coaster")
plot_lst(res_df[res_df$options == "option_f:a" & res_df$fhist == "roller-coaster", ],
         name = "3.2.1/factor_f.a_roller-coaster")
plot_lst(res_df[res_df$options == "option_f:b" & res_df$fhist == "roller-coaster", ],
         name = "3.2.1/factor_f.b_roller-coaster")
plot_lst(res_df[res_df$options == "option_f:c" & res_df$fhist == "roller-coaster", ],
         name = "3.2.1/factor_f.c_roller-coaster")
plot_lst(res_df[res_df$options == "option_b:a" & res_df$fhist == "roller-coaster", ],
         name = "3.2.1/factor_b.a_roller-coaster")

### with actual values from OM for factor f
plot_lst(res_df[res_df$options == "option_f:a perfect_knowledge:TRUE" & res_df$fhist == "one-way", ],
         name = "3.2.1/perfect_factor_f.a_one-way")
plot_lst(res_df[res_df$options == "option_f:b perfect_knowledge:TRUE" & res_df$fhist == "one-way", ],
         name = "3.2.1/perfect_factor_f.b_one-way")
plot_lst(res_df[res_df$options == "option_f:c perfect_knowledge:TRUE" & res_df$fhist == "one-way", ],
         name = "3.2.1/perfect_factor_f.c_one-way")
plot_lst(res_df[res_df$options == "option_f:a perfect_knowledge:TRUE" & res_df$fhist == "roller-coaster", ],
         name = "3.2.1/perfect_factor_f.a_roller-coaster")
plot_lst(res_df[res_df$options == "option_f:b perfect_knowledge:TRUE" & res_df$fhist == "roller-coaster", ],
         name = "3.2.1/perfect_factor_f.b_roller-coaster")
plot_lst(res_df[res_df$options == "option_f:c perfect_knowledge:TRUE" & res_df$fhist == "roller-coaster", ],
         name = "3.2.1/perfect_factor_f.c_roller-coaster")

### M/K = 1.5
### with actual values from OM for factor f
plot_lst(res_df[res_df$options == "option_f:a MK:1.5" & res_df$fhist == "one-way", ],
         name = "3.2.1/perfect_factor_factor_f.a_MK1.5_one-way")
plot_lst(res_df[res_df$options == "option_f:a MK:1.5" & res_df$fhist == "roller-coaster", ],
         name = "3.2.1/perfect_factor_factor_f.a_MK1.5_roller-coaster")

### 3.2.2. perfect knowledge
plot_lst(res_df[res_df$catch_rule == "3.2.2" & res_df$fhist == "one-way" &
                res_df$TAC == 2 & res_df$uncertainty == "perfect_knowledge", ],
         name = "3.2.2/perfect_knowledge_one-way")
plot_lst(res_df[res_df$catch_rule == "3.2.2" & res_df$fhist == "roller-coaster" &
                  res_df$TAC == 2 & res_df$uncertainty == "perfect_knowledge", ],
         name = "3.2.2/perfect_knowledge_roller-coaster")
### 3.2.2 perfect knowledge annual TAC
plot_lst(res_df[res_df$catch_rule == "3.2.2" & res_df$fhist == "one-way" &
                res_df$TAC == 1 & res_df$uncertainty == "perfect_knowledge", ],
         name = "3.2.2/perfect_knowledge_annual_one-way")
plot_lst(res_df[res_df$catch_rule == "3.2.2" & res_df$fhist == "roller-coaster" &
                res_df$TAC == 1 & res_df$uncertainty == "perfect_knowledge", ],
         name = "3.2.2/perfect_knowledge_annual_roller-coaster")
### 3.2.2 with observation error
plot_lst(res_df[res_df$catch_rule == "3.2.2" & res_df$fhist == "one-way" &
                  res_df$TAC == 2 & res_df$uncertainty == "observation_error" &
                  is.na(res_df$HCRmult) & is.na(res_df$w), ],
         name = "3.2.2/observation_error_one-way")
plot_lst(res_df[res_df$catch_rule == "3.2.2" & res_df$fhist == "roller-coaster" &
                  res_df$TAC == 2 & res_df$uncertainty == "observation_error" &
                  is.na(res_df$HCRmult) & is.na(res_df$w), ],
         name = "3.2.2/observation_error_roller-coaster")
### 3.2.2 with observation error and HCR multiplier
plot_lst(res_df[res_df$catch_rule == "3.2.2" & res_df$fhist == "one-way" &
                res_df$TAC == 2 & res_df$uncertainty == "observation_error" &
                res_df$HCRmult == 0.5 & !is.na(res_df$HCRmult), ],
         name = "3.2.2/observation_error_HCRmult0.5_one-way")
plot_lst(res_df[res_df$catch_rule == "3.2.2" & res_df$fhist == "roller-coaster" &
                  res_df$TAC == 2 & res_df$uncertainty == "observation_error" &
                  res_df$HCRmult == 0.5 & !is.na(res_df$HCRmult), ],
         name = "3.2.2/observation_error_HCRmult0.5_roller-coaster")

### 3.2.1. combinations
plot_lst(res_df[res_df$options == "option_f:a perfect_knowledge:TRUE option_r:a option_b:a" & 
                  res_df$fhist == "one-way", ],
         name = "3.2.1/perfect_comb_f.a_r.a_one-way")
plot_lst(res_df[res_df$options == "option_f:a perfect_knowledge:TRUE option_r:b option_b:a" & 
                  res_df$fhist == "one-way", ],
         name = "3.2.1/perfect_comb_f.a_r.b_one-way")
plot_lst(res_df[res_df$options == "option_f:b perfect_knowledge:TRUE option_r:a option_b:a" & 
                  res_df$fhist == "one-way", ],
         name = "3.2.1/perfect_comb_f.b_r.a_one-way")
plot_lst(res_df[res_df$options == "option_f:b perfect_knowledge:TRUE option_r:b option_b:a" & 
                  res_df$fhist == "one-way", ],
         name = "3.2.1/perfect_comb_f.b_r.b_one-way")
plot_lst(res_df[res_df$options == "option_f:a perfect_knowledge:TRUE option_r:a option_b:a" & 
                  res_df$fhist == "roller-coaster", ],
         name = "3.2.1/perfect_comb_f.a_r.a_roller-coaster")
plot_lst(res_df[res_df$options == "option_f:a perfect_knowledge:TRUE option_r:b option_b:a" & 
                  res_df$fhist == "roller-coaster", ],
         name = "3.2.1/perfect_comb_f.a_r.b_roller-coaster")
plot_lst(res_df[res_df$options == "option_f:b perfect_knowledge:TRUE option_r:a option_b:a" & 
                  res_df$fhist == "roller-coaster", ],
         name = "3.2.1/perfect_comb_f.b_r.a_roller-coaster")
plot_lst(res_df[res_df$options == "option_f:b perfect_knowledge:TRUE option_r:b option_b:a" & 
                  res_df$fhist == "roller-coaster", ],
         name = "3.2.1/perfect_comb_f.b_r.b_roller-coaster")

### 3.2.1. combinations & noise
plot_lst(res_df[res_df$options == "option_f:a option_r:a option_b:a MK:1.5" & 
                  res_df$fhist == "one-way" & 
                  is.na(res_df$HCRmult) & is.na(res_df$upper_constraint), ],
         name = "3.2.1/noise_comb_f.a_r.a_one-way")
plot_lst(res_df[res_df$options == "option_f:a option_r:b option_b:a MK:1.5" & 
                  res_df$fhist == "one-way" & 
                  is.na(res_df$HCRmult) & is.na(res_df$upper_constraint), ],
         name = "3.2.1/noise_comb_f.a_r.b_one-way")
plot_lst(res_df[res_df$options == "option_f:a option_r:a option_b:a MK:1.5" & 
                  res_df$fhist == "roller-coaster" & 
                  is.na(res_df$HCRmult) & is.na(res_df$upper_constraint), ],
         name = "3.2.1/noise_comb_f.a_r.a_roller-coaster")
plot_lst(res_df[res_df$options == "option_f:a option_r:b option_b:a MK:1.5" & 
                  res_df$fhist == "roller-coaster" & 
                  is.na(res_df$HCRmult) & is.na(res_df$upper_constraint), ],
         name = "3.2.1/noise_comb_f.a_r.b_roller-coaster")


### ------------------------------------------------------------------------ ###
#### plotting ####
### ------------------------------------------------------------------------ ###


### plot all scenarios
scenarios <- res_df$scenario[grepl(pattern = "^option_r\\:*", x = res_df$options)]
. <- lapply(scenarios, function(x){
  plot_factor(scenarios = x, names = x, tracking_factors = "HCR3.2.1r",
              res_path = "output/perfect_knowledge/plots/factors/r/",
              save = TRUE, file_name = x,
              refpts = c(HCR3.2.1r = 1))
})

plot_factor <- function(scenarios, names, res_path = NULL, save = FALSE,
                        tracking_factors = NULL, file_name = "",
                        quantiles = c(0.05, 0.25, 0.50, 0.75, 0.95),
                        quants = c("ssb", "catch", "fbar"),
                        refpts = NULL, ncol = 2, median = FALSE){
  
  ### load stock(s)
  stk_lst <- lapply(scenarios, function(x){
    readRDS(paste0("output/perfect_knowledge/combined/", x, ".rds"))
  })
  names(stk_lst) <- names
  
  ### extract additional tracking factors, if requested
  ### and get quantiles
  if (!is.null(tracking_factors)) {
    ### go through stock list
    factors_lst <- lapply(stk_lst, function(stk_i){
      ### go through factors
      res <- lapply(tracking_factors, function(fact_i){
        ### extract factors
        res_temp <- attr(stk_i, "tracking")[fact_i]
        ### get quantiles
        res_temp <- quantile(res_temp, quantiles, na.rm = TRUE)
      })
      names(res) <- tracking_factors
      return(res)
    })
  }
  
  ### extract quants from stock list
  stk_quants <- lapply(stk_lst, function(stk_i){
    ### extract each requested quant
    res <- lapply(quants, function(quant){
      ### get slot
      res <- get(quant)(stk_i)
      ### quantiles
      res <- quantile(res, quantiles, na.rm = TRUE)
      ### modify dimnames for coercion
      names(dimnames(res))[1] <- "metric"
      return(res)
    })
    names(res) <- quants
    return(res)
  })
  
  ### combine the two list
  quant_lst <- lapply(seq_along(scenarios), function(scenario){
    if (exists("factors_lst")) {
      FLQuants(c(factors_lst[[scenario]], stk_quants[[scenario]]))
    } else {
      FLQuants(stk_quants[[scenario]])
    }
  })
  names(quant_lst) <- names
  
  ### coerce into data frames
  ### and reshape
  quant_lst <- lapply(seq_along(quant_lst), function(quant_i){
    res <- as.data.frame(quant_lst[[quant_i]])
    res <- res[!is.na(res$data), ]
    ### reshape
    res <- dcast(res, qname + year ~ iter, value.var = "data")
    ### add stock name
    res$stock <- names(quant_lst)[quant_i]
    return(res)
  })
  
  ### combine into single df
  data <- do.call(rbind, quant_lst)
  
  ### create df for refpts
  if (!is.null(refpts)) {
    df_ref <- data.frame(qname = names(refpts),
                         data = refpts)
  }
  
  ### more than 1 stock?
  mult_stk <- ifelse(length(unique(data$stock)) > 1, TRUE, FALSE)
  
  ### plot
  p <- ggplot(data, aes(x = year, y = `50%`, group = stock))
  
  # ### reference values
  if (!is.null(refpts)) {
    p <- p + 
      geom_hline(data = df_ref, aes(yintercept = data), colour = "black",
                 linetype = "dotted")
  }
  
  if (mult_stk){
    p <- p + geom_line(aes(colour = stock))
    if (!isTRUE(median)) {
      p <- p +
        geom_ribbon(aes(x = year, ymin = `5%`, ymax = `95%`, fill = stock), 
                    alpha = 0.15) +
        geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`, fill = stock), 
                    alpha = 0.25)
    }
  } else {
    p <- p + geom_line()
    if (!isTRUE(median)) {
      p <- p +
        geom_ribbon(aes(x = year, ymin = `5%`, ymax = `95%`), 
                    alpha = 0.15) +
        geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), 
                    alpha = 0.25)
    }
  }
  #geom_hline(data = data_refpts, aes(yintercept = y), colour = "red") +
  p <- p +
    geom_vline(xintercept = 100) +
    facet_wrap(~ qname, scales = "free_y", ncol = ncol) +
    ylim(0, NA) +
    theme_bw() +
    ylab("")

  
  # data_refpts <- rbind(data.frame(qname = "L_mean",
  #                                 y = mean(attr(stk, "refpts")["LFeFmsy"])),
  #                      data.frame(qname = "factor_f",
  #                                 y = 1))
  
  ### save
  if (isTRUE(save)){
    ggsave(filename = paste0(res_path, file_name, ".png"), plot = p,
           width = 30, height = 20, units = "cm", dpi = 300, type = "cairo-png")
  } else {
    p
  }
  
}


### ------------------------------------------------------------------------ ###
### extract refpts for all stocks ####
### ------------------------------------------------------------------------ ###
### load list with initial OM objects
# OM_list <- (function(){
#   load("om_list.RData")
#   return(OM_list)
# })()
# library(FLBRP)
# ### extract reference points from brp
# refpts_lst <- lapply(OM_list[1:15], function(x){
#   refpts(x$brp)
# })
# names(refpts_lst) <- unlist(lapply(strsplit(names(refpts_lst), split = "_"), 
#                                    function(x){
#   x[[1]]
# }))
# saveRDS(refpts_lst, file = "output/refpts.rds")
#refpts_lst <- readRDS("output/refpts.rds")
refpts_lst <- readRDS("input/refpts.rds")

### ------------------------------------------------------------------------ ###
### 3.2.1 perfect knowledge r ####
### ------------------------------------------------------------------------ ###

res_df[res_df$catch_rule == "3.2.1" & res_df$stock == "pol-nsea" & 
       res_df$fhist %in% c("one-way") &
       res_df$options %in% c("option_r:a", "option_r:b"), ]

p <- plot_factor(scenarios = c(7, 8), 
                 names = c("r:a", "r:b"), quants = c("rec", "ssb", "catch", "fbar"),
                 tracking_factors = c("HCR3.2.1r"),
                 res_path = "",
                 save = FALSE, file_name = "", ncol = 2)
#stk7 <- attr(readRDS("output/perfect_knowledge/combined/7.rds"), "refpts")
ref_MSY <- data.frame(qname = c("ssb", "catch", "fbar", "HCR3.2.1r"), 
                      data = c(refpts_lst$`pol-nsea`["msy", "ssb"],
                               refpts_lst$`pol-nsea`["msy", "yield"],
                               refpts_lst$`pol-nsea`["msy", "harvest"], 1))
ref_lim <- data.frame(qname = "ssb", data = B_lim)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted")
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "perfect_knowledge_r.a_pol-nsea_one_way.png"), plot = p,
       width = 18, height = 15, units = "cm", dpi = 300, type = "cairo-png")

### risk plot
### "risk" plots
df_plot <- res_df[res_df$catch_rule == "3.2.1" & 
                  res_df$uncertainty == "perfect_knowledge" &
                  res_df$options %in% c("option_r:a", "option_r:b"), ]
### reshape
df_plot <- melt(data = df_plot, 
                id.vars = c("scenario", "stock", "options", "fhist"),
                measure.vars = c("ssb_below_blim_total", "collapse_iter",
                                 "rel_yield"))
levels(df_plot$variable) <- c("p(SSB<B[lim])", 
                              ("iter~collapse"), "rel.~yield")
df_plot$options <- as.factor(df_plot$options)
#levels(df_plot$options) <- c()

p <- ggplot(data = df_plot, 
            aes(x = stock, y = value,  fill = stock)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(value, 2), y = value + 0.05), position = "dodge", size = 2) +
  facet_grid(variable + fhist ~ options, labeller = "label_parsed") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ylab("")
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "perfect_knowledge_r_risks.png"), plot = p,
       width = 18, height = 20, units = "cm", dpi = 300, type = "cairo-png")



### ------------------------------------------------------------------------ ###
### 3.2.1 perfect knowledge f ####
### ------------------------------------------------------------------------ ###

res_df[res_df$catch_rule == "3.2.1" & res_df$stock == "pol-nsea" & 
         res_df$fhist %in% c("one-way", "roller-coaster") &
         res_df$options %in% c("option_f:a perfect_knowledge:TRUE"), ]

p <- plot_factor(scenarios = c(214), 
                 names = c(""), quants = c("rec", "ssb", "catch", "fbar"),
                 tracking_factors = c("L_mean", "HCR3.2.1f"),
                 res_path = "",
                 save = FALSE, file_name = "", ncol = 2)
stk7 <- readRDS("output/perfect_knowledge/combined/7.rds")
ref_MSY <- data.frame(qname = c("ssb", "catch", "fbar", "HCR3.2.1f", "L_mean"), 
                      data = c(refpts_lst$`pol-nsea`["msy", "ssb"],
                               refpts_lst$`pol-nsea`["msy", "yield"],
                               refpts_lst$`pol-nsea`["msy", "harvest"], 
                               1, attr(stk7, "refpts")["LFeFmsy",1]))
ref_lim <- data.frame(qname = "ssb", data = B_lim)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "perfect_knowledge_f.a_LFeFmsy_pol-nsea_one_way.png"), plot = p,
       width = 18, height = 15, units = "cm", dpi = 300, type = "cairo-png")
### roller-coaster
p <- plot_factor(scenarios = c(259), 
                 names = c(""), quants = c("rec", "ssb", "catch", "fbar"),
                 tracking_factors = c("L_mean", "HCR3.2.1f"),
                 res_path = "",
                 save = FALSE, file_name = "", ncol = 2)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "perfect_knowledge_f.a_LFeFmsy_pol-nsea_roller-coaster.png"), plot = p,
       width = 18, height = 15, units = "cm", dpi = 300, type = "cairo-png")


### "risk" plots
df_plot <- res_df[res_df$catch_rule == "3.2.1" & 
                  res_df$uncertainty == "perfect_knowledge" &
                  res_df$options %in% c("option_f:a", "option_f:a perfect_knowledge:TRUE",
                                        "option_f:a MK:1.5"), ]
### reshape
df_plot <- melt(data = df_plot, 
                id.vars = c("scenario", "stock", "options", "fhist"),
                measure.vars = c("collapse_iter"))
df_plot$options <- as.factor(df_plot$options)
levels(df_plot$options) <- c("M&K from OM", "M/K=1.5", "LF=FMSY")

p <- ggplot(data = df_plot, 
            aes(x = stock, y = value,  fill = stock)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(value, 2), y = value + 0.05), position = "dodge", size = 2) +
  facet_grid(options ~ fhist) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ylab("proportion of collapsed iterations")
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "perfect_knowledge_f_a_risks_collapse.png"), plot = p,
       width = 20, height = 12, units = "cm", dpi = 300, type = "cairo-png")


### ------------------------------------------------------------------------ ###
### examples for combinations ####
### ------------------------------------------------------------------------ ###

res_df[res_df$catch_rule == "3.2.1" & res_df$stock == "pol-nsea" & 
       res_df$fhist == "one-way" &
       res_df$options == "option_f:a perfect_knowledge:TRUE option_r:a option_b:a", ]

p <- plot_factor(scenarios = 392, 
            names = c(""), quants = c("rec", "ssb", "catch", "fbar"),
            tracking_factors = c("L_mean", "HCR3.2.1r", "HCR3.2.1f", "HCR3.2.1b"),
            res_path = "",
            save = FALSE, file_name = "", ncol = 2)
stk392refs <- attr(readRDS("output/perfect_knowledge/combined/392.rds"), "refpts")
ref_MSY <- data.frame(qname = c("ssb", "catch", "fbar", "L_mean", "HCR3.2.1r", 
                                "HCR3.2.1f", "HCR3.2.1b"), 
                     data = c(refpts_lst$`pol-nsea`["msy", "ssb"],
                              refpts_lst$`pol-nsea`["msy", "yield"],
                              refpts_lst$`pol-nsea`["msy", "harvest"],
                              stk392refs["LFeFmsy", 1],
                              1, 1, 1))
ref_lim <- data.frame(qname = "ssb", data = B_lim)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "combs_pol-nsea_one_way_r.a_f.a_b.a.png"), plot = p,
       width = 15, height = 15, units = "cm", dpi = 300, type = "cairo-png")

### ------------------------------------------------------------------------ ###
### combinations risk plots ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### risk plots for default parametrization
df_plot <- res_df[res_df$catch_rule == "3.2.1" &
                    res_df$options %in% 
                    c("option_f:a perfect_knowledge:TRUE option_r:a option_b:a",
                      "option_f:a perfect_knowledge:TRUE option_r:b option_b:a",
                      "option_f:b perfect_knowledge:TRUE option_r:a option_b:a",
                      "option_f:b perfect_knowledge:TRUE option_r:b option_b:a"), ]
### reshape
df_plot <- melt(data = df_plot, 
                id.vars = c("scenario", "stock", "fhist", "options"),
                measure.vars = c("ssb_below_blim_total", "ssb_below_blim_iter",
                                 "collapse_total", "collapse_iter"))
df_plot <- melt(data = df_plot, 
                id.vars = c("scenario", "stock", "fhist", "options"),
                measure.vars = c("collapse_iter"))
### set factor levels for options
df_plot$options <- as.factor(df_plot$options)
levels(df_plot$options) <- c("f:a & r:a", "f:a & r:b", "f:b & r:a",
                             "f:b & r:b")

### number of collapse iterations
p <- ggplot(data = df_plot[df_plot$variable == "collapse_iter", ], 
            aes(x = stock, y = value,  fill = stock)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(value, 2), y = value + 0.05), position = "dodge", size = 2) + 
  facet_grid(options ~ fhist) +
  theme_bw() +
  ylim(0, NA) + 
  theme(axis.text.x = element_blank()) +
  ylab("proportion of collapsed iterations")
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "combs_risks_iter.png"), plot = p,
       width = 15, height = 15, units = "cm", dpi = 300, type = "cairo-png")
### risk Blim
p <- ggplot(data = df_plot[df_plot$variable == "ssb_below_blim_total", ], 
            aes(x = stock, y = value,  fill = stock)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(value, 2), y = value + 0.05), position = "dodge", size = 2) + 
  facet_grid(options ~ fhist) +
  theme_bw() +
  ylim(0, NA) + 
  theme(axis.text.x = element_blank()) +
  ylab("probability SSB < Blim")
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "combs_risks_blim.png"), plot = p,
       width = 15, height = 15, units = "cm", dpi = 300, type = "cairo-png")


### ------------------------------------------------------------------------ ###
### examples for combinations with noise ####
### ------------------------------------------------------------------------ ###

res_df[res_df$catch_rule == "3.2.1" & res_df$stock == "pol-nsea" & 
         res_df$fhist %in% c("one-way", "roller-coaster") &
         res_df$options %in% c("option_f:a option_r:a option_b:a MK:1.5",
                               "option_f:a perfect_knowledge:TRUE option_r:a option_b:a") &
         is.na(res_df$HCRmult) & is.na(res_df$upper_constraint), ]
### pol-nsea roller-coaster
p <- plot_factor(scenarios = c(407, 643), 
                 names = c("pol-nsea", "pol-nsea\n+noise"), quants = c("rec", "ssb", "catch", "fbar"),
                 tracking_factors = c("L_mean", "HCR3.2.1r", "HCR3.2.1f", "HCR3.2.1b"),
                 res_path = "",
                 save = FALSE, file_name = "", ncol = 2)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "noise_combs_pol-nsea_roller-coaster_r.a_f.a_b.a.png"), plot = p,
       width = 15, height = 15, units = "cm", dpi = 300, type = "cairo-png")
### zoom into f
p <- plot_factor(scenarios = c(407, 643), 
                 names = c("pol-nsea", "pol-nsea\n+noise"), quants = NULL,
                 tracking_factors = c("HCR3.2.1f"),
                 res_path = "",
                 save = FALSE, file_name = "", ncol = 2)
p <- p + geom_hline(yintercept = 1, linetype = "dashed")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "noise_combs_pol-nsea_roller-coaster_r.a_f.a_b.a_zoomf.png"),
       plot = p, width = 15, height = 15, units = "cm", dpi = 300,
       type = "cairo-png")

### one-way
p <- plot_factor(scenarios = c(392, 628), 
                 names = c("pol-nsea", "pol-nsea\n+noise"), quants = c("rec", "ssb", "catch", "fbar"),
                 tracking_factors = c("L_mean", "HCR3.2.1r", "HCR3.2.1f", "HCR3.2.1b"),
                 res_path = "",
                 save = FALSE, file_name = "", ncol = 2)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "noise_combs_pol-nsea_one_way_r.a_f.a_b.a.png"), plot = p,
       width = 15, height = 15, units = "cm", dpi = 300, type = "cairo-png")

### ------------------------------------------------------------------------ ###
### combinations risk plots with noise ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### risk plots for default parametrization
df_plot <- res_df[res_df$catch_rule == "3.2.1" &
                    res_df$options %in% 
                    c("option_f:a option_r:a option_b:a MK:1.5",
                      "option_f:a option_r:b option_b:a MK:1.5") &
                    is.na(res_df$HCRmult) & is.na(res_df$upper_constraint), ]

### reshape
df_plot <- melt(data = df_plot, 
                id.vars = c("scenario", "stock", "fhist", "options"),
                measure.vars = c("ssb_below_blim_total", "collapse_iter", 
                                 "rel_yield"))
levels(df_plot$variable) <- c("p(SSB<Blim)", 
                              ("iter collapse"), "rel. yield")

### set factor levels for options
df_plot$options <- as.factor(df_plot$options)
levels(df_plot$options) <- c("f:a & r:a & b:a", "f:a & r:b & b:a")

### 
p <- ggplot(data = df_plot, 
            aes(x = stock, y = value,  fill = stock)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(value, 2), y = value + 0.05), position = "dodge", size = 2) + 
  facet_grid(variable+options ~ fhist) +
  theme_bw() +
  ylim(0, NA) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "noise_combs_risks.png"), plot = p,
       width = 20, height = 22, units = "cm", dpi = 300, type = "cairo-png")

### risk Blim
p <- ggplot(data = df_plot[df_plot$variable == "ssb_below_blim_total", ], 
            aes(x = stock, y = value,  fill = stock)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(value, 2), y = value + 0.05), position = "dodge", size = 2) + 
  facet_grid(options ~ fhist) +
  theme_bw() +
  ylim(0, NA) + 
  theme(axis.text.x = element_blank()) +
  ylab("probability SSB < Blim")
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "combs_risks_blim.png"), plot = p,
       width = 15, height = 15, units = "cm", dpi = 300, type = "cairo-png")


### ------------------------------------------------------------------------ ###
### 3.2.1 combinations & noise & b_w ####
### ------------------------------------------------------------------------ ###

### pol-nsea only
### one-way
### f:a r:a
df <- res_df[res_df$catch_rule == "3.2.1" &
              res_df$stock == "pol-nsea" & res_df$fhist == "one-way" & 
              res_df$uncertainty == "observation_error" &
              !is.na(res_df$b_w) &
             grepl(x = res_df$options, pattern = "^option_f:a option_r:a option_b:a MK:1.5*"), ]
p <- plot_factor(scenarios = df$scenario, quants = c("rec", "ssb", "catch", "fbar"),
                 names = df$b_w, 
                 save = FALSE, ncol = 1, median = TRUE)
ref_MSY <- data.frame(qname = c("ssb", "catch", "fbar"), 
                      data = c(refpts_lst[["pol-nsea"]]["msy", "ssb"],
                               refpts_lst[["pol-nsea"]]["msy", "yield"],
                               refpts_lst[["pol-nsea"]]["msy", "harvest"]))
ref_lim <- data.frame(qname = "ssb", data = B_lim)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted") + 
  scale_colour_discrete("w")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "noise_combs_pol-nsea_b_w_one-way_f.a_r.a.png"), plot = p,
       width = 15, height = 12, units = "cm", dpi = 300, type = "cairo-png")
### f:a r:b
df <- res_df[res_df$catch_rule == "3.2.1" &
               res_df$stock == "pol-nsea" & res_df$fhist == "one-way" & 
               res_df$uncertainty == "observation_error" &
               !is.na(res_df$b_w) &
               grepl(x = res_df$options, pattern = "^option_f:a option_r:b option_b:a MK:1.5*"), ]
p <- plot_factor(scenarios = df$scenario, quants = c("rec", "ssb", "catch", "fbar"),
                 names = df$b_w, 
                 save = FALSE, ncol = 1, median = TRUE)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted") + 
  scale_colour_discrete("w")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "noise_combs_pol-nsea_b_w_one-way_f.a_r.b.png"), plot = p,
       width = 15, height = 12, units = "cm", dpi = 300, type = "cairo-png")
### roller-coaster
### f:a r:a
df <- res_df[res_df$catch_rule == "3.2.1" &
               res_df$stock == "pol-nsea" & res_df$fhist == "roller-coaster" & 
               res_df$uncertainty == "observation_error" &
               !is.na(res_df$b_w) &
               grepl(x = res_df$options, pattern = "^option_f:a option_r:a option_b:a MK:1.5*"), ]
p <- plot_factor(scenarios = df$scenario, quants = c("rec", "ssb", "catch", "fbar"),
                 names = df$b_w, 
                 save = FALSE, ncol = 1, median = TRUE)
ref_MSY <- data.frame(qname = c("ssb", "catch", "fbar"), 
                      data = c(refpts_lst[["pol-nsea"]]["msy", "ssb"],
                               refpts_lst[["pol-nsea"]]["msy", "yield"],
                               refpts_lst[["pol-nsea"]]["msy", "harvest"]))
ref_lim <- data.frame(qname = "ssb", data = B_lim)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted") + 
  scale_colour_discrete("w")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "noise_combs_pol-nsea_b_w_roller-coaster_f.a_r.a.png"), plot = p,
       width = 15, height = 12, units = "cm", dpi = 300, type = "cairo-png")
### f:a r:b
df <- res_df[res_df$catch_rule == "3.2.1" &
               res_df$stock == "pol-nsea" & res_df$fhist == "roller-coaster" & 
               res_df$uncertainty == "observation_error" &
               !is.na(res_df$b_w) &
               grepl(x = res_df$options, pattern = "^option_f:a option_r:b option_b:a MK:1.5*"), ]
p <- plot_factor(scenarios = df$scenario, quants = c("rec", "ssb", "catch", "fbar"),
                 names = df$b_w, 
                 save = FALSE, ncol = 1, median = TRUE)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted") + 
  scale_colour_discrete("w")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "noise_combs_pol-nsea_b_w_roller-coaster_f.a_r.b.png"), plot = p,
       width = 15, height = 12, units = "cm", dpi = 300, type = "cairo-png")

### risk plot
df_plot <- res_df[res_df$catch_rule == "3.2.1" &
               res_df$stock == "pol-nsea" & 
               res_df$uncertainty == "observation_error" &
               !is.na(res_df$b_w) &
                is.na(res_df$b_z) & is.na(res_df$upper_constraint), ]
df_plot$options <- substr(x = df_plot$options, start = 1, stop = 39)
df_plot$options <- as.factor(df_plot$options)
levels(df_plot$options) <- c("r:a & f:a & b:a", "r:b & f:a & b:a")
df_plot <- melt(data = df_plot, 
                id.vars = c("scenario", "stock", "fhist", "b_w", "options"),
                measure.vars = c("ssb_below_blim_total", "collapse_iter", 
                                 "rel_yield"))
levels(df_plot$variable) <- c("p(SSB<Blim)", 
                              ("iter collapse"), "rel. yield")

p <- ggplot(data = df_plot, 
            aes(x = b_w, y = value, colour = fhist)) +
  geom_line() + geom_point() +
  facet_grid(variable ~ options) +
  theme_bw() +
  ylim(0, NA) +
  xlab("w for Itrigger") + ylab("")
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "noise_combs_pol-nsea_risks_b_w.png"), plot = p,
       width = 15, height = 12.5, units = "cm", dpi = 300, type = "cairo-png")



### ------------------------------------------------------------------------ ###
### 3.2.1 combinations & noise & multiplier ####
### ------------------------------------------------------------------------ ###

### pol-nsea only
### one-way
### f:a r:a
df <- res_df[res_df$catch_rule == "3.2.1" &
               res_df$stock == "pol-nsea" & res_df$fhist == "one-way" & 
               res_df$uncertainty == "observation_error" &
               is.na(res_df$b_w) &
               grepl(x = res_df$options, 
                     pattern = "^option_f:a option_r:a option_b:a MK:1.5*") &
               is.na(res_df$b_z) & is.na(res_df$upper_constraint), ]
df$HCRmult[is.na(df$HCRmult)] <- 1
p <- plot_factor(scenarios = df$scenario, quants = c("rec", "ssb", "catch", "fbar"),
                 names = df$HCRmult, 
                 save = FALSE, ncol = 1, median = TRUE)
ref_MSY <- data.frame(qname = c("ssb", "catch", "fbar"), 
                      data = c(refpts_lst[["pol-nsea"]]["msy", "ssb"],
                               refpts_lst[["pol-nsea"]]["msy", "yield"],
                               refpts_lst[["pol-nsea"]]["msy", "harvest"]))
ref_lim <- data.frame(qname = "ssb", data = B_lim)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted") + 
  scale_colour_discrete("advice\nmultiplier")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "noise_combs_pol-nsea_multiplier_one-way_f.a_r.a.png"), plot = p,
       width = 15, height = 12, units = "cm", dpi = 300, type = "cairo-png")
### f:a r:b
df <- res_df[res_df$catch_rule == "3.2.1" &
               res_df$stock == "pol-nsea" & res_df$fhist == "one-way" & 
               res_df$uncertainty == "observation_error" &
               is.na(res_df$b_w) &
               grepl(x = res_df$options, 
                     pattern = "^option_f:a option_r:b option_b:a MK:1.5*") &
               is.na(res_df$b_z) & is.na(res_df$upper_constraint), ]
df$HCRmult[is.na(df$HCRmult)] <- 1
p <- plot_factor(scenarios = df$scenario, quants = c("rec", "ssb", "catch", "fbar"),
                 names = df$HCRmult, 
                 save = FALSE, ncol = 1, median = TRUE)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted") + 
  scale_colour_discrete("advice\nmultiplier")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "noise_combs_pol-nsea_multiplier_one-way_f.a_r.b.png"), plot = p,
       width = 15, height = 12, units = "cm", dpi = 300, type = "cairo-png")
### roller-coaster
### f:a r:a
df <- res_df[res_df$catch_rule == "3.2.1" &
               res_df$stock == "pol-nsea" & res_df$fhist == "roller-coaster" & 
               res_df$uncertainty == "observation_error" &
               is.na(res_df$b_w) &
               grepl(x = res_df$options, 
                     pattern = "^option_f:a option_r:a option_b:a MK:1.5*") &
               is.na(res_df$b_z) & is.na(res_df$upper_constraint), ]
df$HCRmult[is.na(df$HCRmult)] <- 1

p <- plot_factor(scenarios = df$scenario, quants = c("rec", "ssb", "catch", "fbar"),
                 names = df$HCRmult, 
                 save = FALSE, ncol = 1, median = TRUE)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted") + 
  scale_colour_discrete("advice\nmultiplier")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "noise_combs_pol-nsea_multiplier_roller-coaster_f.a_r.a.png"), plot = p,
       width = 15, height = 12, units = "cm", dpi = 300, type = "cairo-png")
### f:a r:b
df <- res_df[res_df$catch_rule == "3.2.1" &
               res_df$stock == "pol-nsea" & res_df$fhist == "roller-coaster" & 
               res_df$uncertainty == "observation_error" &
               is.na(res_df$b_w) &
               grepl(x = res_df$options, 
                     pattern = "^option_f:a option_r:b option_b:a MK:1.5*") &
               is.na(res_df$b_z) & is.na(res_df$upper_constraint), ]
df$HCRmult[is.na(df$HCRmult)] <- 1
p <- plot_factor(scenarios = df$scenario, quants = c("rec", "ssb", "catch", "fbar"),
                 names = df$HCRmult, 
                 save = FALSE, ncol = 1, median = TRUE)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted") + 
  scale_colour_discrete("advice\nmultiplier")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "noise_combs_pol-nsea_multiplier_roller-coaster_f.a_r.b.png"), plot = p,
       width = 15, height = 12, units = "cm", dpi = 300, type = "cairo-png")


### risk plot
df_plot <- res_df[res_df$catch_rule == "3.2.1" &
               res_df$stock == "pol-nsea" & 
               res_df$uncertainty == "observation_error" &
               is.na(res_df$b_w) &
               (grepl(x = res_df$options, 
                      pattern = "^option_f:a option_r:a option_b:a MK:1.5*") |
               grepl(x = res_df$options, 
                     pattern = "^option_f:a option_r:b option_b:a MK:1.5*")) &
              is.na(res_df$b_z) & is.na(res_df$upper_constraint), ]
df_plot$HCRmult[is.na(df_plot$HCRmult)] <- 1
df_plot$options <- as.factor(df_plot$options)
levels(df_plot$options) <- c("r:a & f:a & b:a", "r:b & f:a & b:a")
df_plot <- melt(data = df_plot, 
                id.vars = c("scenario", "stock", "fhist", "HCRmult", "options"),
                measure.vars = c("ssb_below_blim_total", "collapse_iter", 
                                 "rel_yield"))
levels(df_plot$variable) <- c("p(SSB<Blim)", 
                              ("iter collapse"), "rel. yield")

p <- ggplot(data = df_plot, 
            aes(x = HCRmult, y = value, colour = fhist)) +
  geom_line() + geom_point() +
  facet_grid(variable ~ options) +
  theme_bw() +
  ylim(0, NA) +
  xlab("advice multiplier")
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "noise_combs_pol-nsea_risks_iter_multiplier.png"), plot = p,
       width = 15, height = 12.5, units = "cm", dpi = 300, type = "cairo-png")


### ------------------------------------------------------------------------ ###
### individual plots for 3.2.2 ####
### ------------------------------------------------------------------------ ###

### get scenarios
scns <- res_df$scenario[res_df$catch_rule == "3.2.2"]

### loop through them
. <- lapply(scns, function(x){
  scn_p <- res_df[x, ]
  p <- plot_factor(scenarios = scn_p$scenario, quants = c("rec", "ssb", "catch", "fbar"),
                   names = c(""), save = FALSE, ncol = 1)
  ref_MSY <- data.frame(qname = c("ssb", "catch", "fbar"), 
                        data = c(refpts_lst[[ac(scn_p$stock)]]["msy", "ssb"],
                                 refpts_lst[[ac(scn_p$stock)]]["msy", "yield"],
                                 refpts_lst[[ac(scn_p$stock)]]["msy", "harvest"]))
  ref_lim <- data.frame(qname = "ssb", data = 162.79)
  p <- p +
    geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
    geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted")
  ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/", 
                           scn_p$scenario,".png"), plot = p,
         width = 15, height = 15, units = "cm", dpi = 300, type = "cairo-png")
})



### pol-nsea
### add uncertainty
res_df[res_df$catch_rule == "3.2.2" &
                res_df$TAC == 2 &
                res_df$stock == "pol-nsea" & res_df$fhist == "one-way", ]
p <- plot_factor(scenarios = c(182, 362), quants = c("rec", "ssb", "catch", "fbar"),
                 names = c("perfect\nknowledge", "observation\nerror"), 
                 save = FALSE, ncol = 1)
ref_MSY <- data.frame(qname = c("ssb", "catch", "fbar"), 
                      data = c(refpts_lst[["pol-nsea"]]["msy", "ssb"],
                               refpts_lst[["pol-nsea"]]["msy", "yield"],
                               refpts_lst[["pol-nsea"]]["msy", "harvest"]))
ref_lim <- data.frame(qname = "ssb", data = B_lim)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/",
                         "0_pol-nsea_uncertainty.png"), plot = p,
       width = 15, height = 15, units = "cm", dpi = 300, type = "cairo-png")


### her-nis
### HCR multiplier

df <- res_df2[res_df$catch_rule == "3.2.2" &
         res_df$stock == "her-nis" & res_df$fhist == "one-way" & 
         res_df$uncertainty == "observation_error", ]
p <- plot_factor(scenarios = df$scenario, quants = c("rec", "ssb", "catch", "fbar"),
                 names = df$HCRmult, 
                 save = FALSE, ncol = 1, median = TRUE)
ref_MSY <- data.frame(qname = c("ssb", "catch", "fbar"), 
                      data = c(refpts_lst[["her-nis"]]["msy", "ssb"],
                               refpts_lst[["her-nis"]]["msy", "yield"],
                               refpts_lst[["her-nis"]]["msy", "harvest"]))
ref_lim <- data.frame(qname = "ssb", data = 250)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted") + 
  scale_colour_discrete("HCR\nmultiplier")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/",
                         "0_her-nis_multiplier_one-way.png"), plot = p,
       width = 15, height = 15, units = "cm", dpi = 300, type = "cairo-png")

### same for roller-coaster
df <- res_df2[res_df$catch_rule == "3.2.2" &
                res_df$stock == "her-nis" & res_df$fhist == "roller-coaster" & 
                res_df$uncertainty == "observation_error", ]
p <- plot_factor(scenarios = df$scenario, quants = c("rec", "ssb", "catch", "fbar"),
                 names = df$HCRmult, 
                 save = FALSE, ncol = 1, median = TRUE)
ref_MSY <- data.frame(qname = c("ssb", "catch", "fbar"), 
                      data = c(refpts_lst[["her-nis"]]["msy", "ssb"],
                               refpts_lst[["her-nis"]]["msy", "yield"],
                               refpts_lst[["her-nis"]]["msy", "harvest"]))
ref_lim <- data.frame(qname = "ssb", data = 250)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted") + 
  scale_colour_discrete("HCR\nmultiplier")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/",
                         "0_her-nis_multiplier_roller-coaster.png"), plot = p,
       width = 15, height = 15, units = "cm", dpi = 300, type = "cairo-png")

### san-ns4
### HCR multiplier
df <- res_df2[res_df2$catch_rule == "3.2.2" &
                res_df2$stock == "san-ns4" & res_df2$fhist == "one-way" & 
                res_df2$uncertainty == "observation_error", ]
p <- plot_factor(scenarios = df$scenario, quants = c("rec", "ssb", "catch", "fbar"),
                 names = df$HCRmult, 
                 save = FALSE, ncol = 1, median = TRUE)
ref_MSY <- data.frame(qname = c("ssb", "catch", "fbar"), 
                      data = c(refpts_lst[["san-ns4"]]["msy", "ssb"],
                               refpts_lst[["san-ns4"]]["msy", "yield"],
                               refpts_lst[["san-ns4"]]["msy", "harvest"]))
ref_lim <- data.frame(qname = "ssb", data = 250)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted") + 
  scale_colour_discrete("HCR\nmultiplier")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/",
                         "0_san-ns4_multiplier_one-way.png"), plot = p,
       width = 15, height = 15, units = "cm", dpi = 300, type = "cairo-png")
### same for roller-coaster
df <- res_df2[res_df$catch_rule == "3.2.2" &
                res_df$stock == "san-ns4" & res_df$fhist == "roller-coaster" & 
                res_df$uncertainty == "observation_error", ]
p <- plot_factor(scenarios = df$scenario, quants = c("rec", "ssb", "catch", "fbar"),
                 names = df$HCRmult, 
                 save = FALSE, ncol = 1, median = TRUE)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted") + 
  scale_colour_discrete("HCR\nmultiplier")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/",
                         "0_san-ns4_multiplier_roller-coaster.png"), plot = p,
       width = 15, height = 15, units = "cm", dpi = 300, type = "cairo-png")


### "risk" plots
df_plot <- res_df2[res_df2$catch_rule == "3.2.2" & 
                   res_df2$uncertainty == "observation_error" & 
                   res_df2$stock %in% c("her-nis", "san-ns4") &
                   res_df2$w == 1.4, ]
### reshape
df_plot <- melt(data = df_plot, 
                id.vars = c("scenario", "stock", "fhist", "HCRmult"),
                measure.vars = c("ssb_below_blim_total", "collapse_total", 
                                 "rel_yield"))
levels(df_plot$variable) <- c("p(SSB<B[lim])", 
                              ("iter~collapse"), "rel.~yield")

p <- ggplot(data = df_plot, 
            aes(x = HCRmult, y = value,  colour = stock, 
                linetype = fhist, shape = fhist)) +
  geom_line() + geom_point() +
  facet_wrap(~ variable, labeller = "label_parsed") +
  theme_bw() +
  ylim(0, NA) +
  scale_shape_discrete("fishing\nhistory") +
  scale_linetype_discrete("fishing\nhistory") +
  xlab("catch advice multiplier") + ylab("")
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/",
                         "pelagics_risk_HCRmult.png"), plot = p,
       width = 15, height = 8, units = "cm", dpi = 300, type = "cairo-png")

### ------------------------------------------------------------------------ ###
### w
### san-ns4
### w
df <- res_df2[res_df2$catch_rule == "3.2.2" &
                res_df2$stock == "san-ns4" & res_df2$fhist == "one-way" & 
                res_df2$uncertainty == "observation_error" &
                res_df2$HCRmult == 1, ]
p <- plot_factor(scenarios = df$scenario, quants = c("rec", "ssb", "catch", "fbar"),
                 names = df$w, 
                 save = FALSE, ncol = 1, median = TRUE)
ref_MSY <- data.frame(qname = c("ssb", "catch", "fbar"), 
                      data = c(refpts_lst[["san-ns4"]]["msy", "ssb"],
                               refpts_lst[["san-ns4"]]["msy", "yield"],
                               refpts_lst[["san-ns4"]]["msy", "harvest"]))
ref_lim <- data.frame(qname = "ssb", data = B_lim)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted") + 
  scale_colour_discrete("w")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/",
                         "0_san-ns4_w_one-way.png"), plot = p,
       width = 15, height = 15, units = "cm", dpi = 300, type = "cairo-png")
### same for roller-coaster
df <- res_df2[res_df2$catch_rule == "3.2.2" &
                res_df2$stock == "san-ns4" & res_df2$fhist == "roller-coaster" & 
                res_df2$uncertainty == "observation_error" &
                res_df2$HCRmult == 1, ]
p <- plot_factor(scenarios = df$scenario, quants = c("rec", "ssb", "catch", "fbar"),
                 names = df$w, 
                 save = FALSE, ncol = 1, median = TRUE)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted") + 
  scale_colour_discrete("w")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/",
                         "0_san-ns4_w_roller-coaster.png"), plot = p,
       width = 15, height = 15, units = "cm", dpi = 300, type = "cairo-png")
### her-nis
df <- res_df2[res_df2$catch_rule == "3.2.2" &
                res_df2$stock == "her-nis" & res_df2$fhist == "one-way" & 
                res_df2$uncertainty == "observation_error" &
                res_df2$HCRmult == 1, ]
p <- plot_factor(scenarios = df$scenario, quants = c("rec", "ssb", "catch", "fbar"),
                 names = df$w, 
                 save = FALSE, ncol = 1, median = TRUE)
ref_MSY <- data.frame(qname = c("ssb", "catch", "fbar"), 
                      data = c(refpts_lst[["san-ns4"]]["msy", "ssb"],
                               refpts_lst[["san-ns4"]]["msy", "yield"],
                               refpts_lst[["san-ns4"]]["msy", "harvest"]))
ref_lim <- data.frame(qname = "ssb", data = B_lim)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted") + 
  scale_colour_discrete("w")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/",
                         "0_her-nis_w_one-way.png"), plot = p,
       width = 15, height = 15, units = "cm", dpi = 300, type = "cairo-png")
### same for roller-coaster
df <- res_df2[res_df2$catch_rule == "3.2.2" &
                res_df2$stock == "her-nis" & res_df2$fhist == "roller-coaster" & 
                res_df2$uncertainty == "observation_error" &
                res_df2$HCRmult == 1, ]
p <- plot_factor(scenarios = df$scenario, quants = c("rec", "ssb", "catch", "fbar"),
                 names = df$w, 
                 save = FALSE, ncol = 1, median = TRUE)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted") + 
  scale_colour_discrete("w")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/",
                         "0_her-nis_w_roller-coaster.png"), plot = p,
       width = 15, height = 15, units = "cm", dpi = 300, type = "cairo-png")

### "risk" plots
df_plot <- res_df2[res_df2$catch_rule == "3.2.2" & 
                   res_df2$uncertainty == "observation_error" & 
                   res_df2$stock %in% c("her-nis", "san-ns4") &
                   res_df2$HCRmult == 1, ]
### reshape
df_plot <- melt(data = df_plot, 
                id.vars = c("scenario", "stock", "fhist", "w"),
                measure.vars = c("ssb_below_blim_total", "collapse_total", 
                                 "rel_yield"))
levels(df_plot$variable) <- c("p(SSB<B[lim])", 
                              ("iter~collapse"), "rel.~yield")

p <- ggplot(data = df_plot, 
            aes(x = w, y = value,  colour = stock, 
                linetype = fhist, shape = fhist)) +
  geom_line() + geom_point() +
  facet_wrap(~ variable, labeller = "label_parsed") +
  theme_bw() +
  ylim(0, NA) +
  scale_shape_discrete("fishing\nhistory") +
  scale_linetype_discrete("fishing\nhistory") +
  xlab("catch advice multiplier") + ylab("")
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/",
                         "pelagics_risk_w.png"), plot = p,
       width = 15, height = 8, units = "cm", dpi = 300, type = "cairo-png")

### ------------------------------------------------------------------------ ###
### risk plots for default parametrization
df_plot <- subset(res_df, catch_rule == "3.2.2" & 
                    uncertainty == "observation_error" &
                    is.na(HCRmult) & is.na(w))
### reshape
df_plot <- melt(data = df_plot, 
                id.vars = c("scenario", "stock", "fhist"),
                measure.vars = c("ssb_below_blim_total", "collapse_total",
                                 "rel_yield"))
levels(df_plot$variable) <- c("p(SSB<Blim)", "iter collapse", "rel. yield")

p <- ggplot(data = df_plot, 
            aes(x = stock, y = value,  fill = stock)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(value, 3), y = value + 0.05), position = "dodge", 
            size = 2) + 
  facet_grid(variable ~ fhist) +
  theme_bw() +
  ylim(0, NA) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  #scale_shape_discrete("fishing\nhistory") +
  #scale_linetype_discrete("fishing\nhistory")
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.2/",
                         "all_risks.png"), plot = p,
       width = 20, height = 13, units = "cm", dpi = 300, type = "cairo-png")



### ------------------------------------------------------------------------ ###
### 3.1 SPiCT, individual plots ####
### ------------------------------------------------------------------------ ###

### example pol-nsea
res_df[res_df$catch_rule == "3.1" & res_df$stock == "pol-nsea", ]
stk613 <- readRDS("output/perfect_knowledge/combined/613.rds")
p <- plot_factor(scenarios = 613, 
                 names = c(""), quants = c("rec", "ssb", "catch", "fbar"),
                 tracking_factors = c("spict_f", "spict_fmsy", "spict_b", "spict_bmsy"),
                 res_path = "",
                 save = FALSE, file_name = "", ncol = 2)
ref_MSY <- data.frame(qname = c("ssb", "catch", "fbar"), 
                      data = c(refpts_lst$`pol-nsea`["msy", "ssb"],
                               refpts_lst$`pol-nsea`["msy", "yield"],
                               refpts_lst$`pol-nsea`["msy", "harvest"]
                               ))
ref_lim <- data.frame(qname = "ssb", data = B_lim)
p <- p +
  geom_hline(data = ref_MSY,  aes(yintercept = data), linetype = "dashed") +
  geom_hline(data = ref_lim, aes(yintercept = data), linetype = "dotted")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.1/",
                         "pol-nsea_roller-coaster.png"), plot = p,
       width = 15, height = 15, units = "cm", dpi = 300, type = "cairo-png")


### ------------------------------------------------------------------------ ###
#### fishing history plots ####
### ------------------------------------------------------------------------ ###

fhist_lst <- lapply(1:30, function(x){
  ### load stock
  tmp <- (function(){
    load(paste0("input/stocks/observation_error/", x, ".RData"))
    stk <- qapply(stk, iterMedians)
    stk[, ac(75:100)]
  })()
})
#data(wklife)
fhist_lst <- as(fhist_lst, "FLStocks")
### set names
data(wklife, package = "FLife")
stk_names <- ac(rep(wklife$stock, 2))
for(i in seq_along(fhist_lst)){
  name(fhist_lst[[i]]) <- stk_names[i]
}

### extract quants for plotting for each stock
lst_quants <- lapply(fhist_lst, function(x){
  ### extract quants
  res <- FLQuants(Rec = rec(x), SSB = ssb(x), Catch = catch(x), Harvest = fbar(x))
  ### coerce into data frame
  res <- as.data.frame(res)
  ### add stock name
  res$stock <- name(x)
  ### fishing history (dummy)
  res$fhist <- "one-way"
  return(res)
})
### add fishing history
lst_quants[16:30] <- lapply(lst_quants[16:30], function(x){
  x$fhist <- "roller-coaster"
  return(x)
})
lst_quants <- do.call(rbind, lst_quants)
### plot all stocks
p <- ggplot(data = lst_quants, aes(x = year, y = data, colour = stock)) + geom_line() +
  facet_grid(qname ~ fhist, scales = "free") +
  xlab("year") + ylab("") +
  theme_bw() +
  ylim(0, NA)
p
ggsave(filename = paste0("input/OMs_fhist.png"),
       width = 15, height = 12, units = "cm", dpi = 300, type = "cairo-png")


# plot(as(fhist_lst, "FLStocks")[1:15]) + theme_bw()
# ggsave(filename = paste0("input/OMs_one-way.png"),
#        width = 10, height = 15, units = "cm", dpi = 300, type = "cairo-png")
# plot(as(fhist_lst, "FLStocks")[16:30]) + theme_bw()
# ggsave(filename = paste0("input/OMs_roller-coaster.png"),
#        width = 10, height = 15, units = "cm", dpi = 300, type = "cairo-png")

### plot one example history with uncertainty
# tmp <- (function(){
#   load(paste0("input/stocks/observation_error/", 1, ".RData"))
#   stk
# })()
### extract
OM_list <- readRDS("../../wklifeFL/exec/oms.rds")
stk_ex <- OM_list[[17]]$stk_full
### extract quants
stk_ex <- FLQuants(Rec = rec(stk_ex), SSB = ssb(stk_ex), Catch = catch(stk_ex), 
                   Harvest = fbar(stk_ex))
### quantiles
stk_ex <- lapply(stk_ex, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                na.rm = TRUE)
### data frame
stk_ex <- as.data.frame(stk_ex)
### reshape
stk_ex <- dcast(stk_ex, qname + year ~ iter, value.var = "data")
### plot

p <- ggplot(data = stk_ex, 
            aes(x = year, y = `50%`)) + 
  geom_ribbon(aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.15) +
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`),  alpha = 0.25) +
  geom_line() +
  facet_wrap(~ qname, scales = "free", ncol = 1) +
  xlab("year") + ylab("") +
  theme_bw() +
  ylim(0, NA)
p
ggsave(filename = paste0("input/OMs_fhist_pol-nsea.png"),
       width = 15, height = 12, units = "cm", dpi = 300, type = "cairo-png")



# plot(OM_list[[17]]$stk_full, probs = c(0.05,0.25,0.5,0.75,0.95)) +
#   theme_bw() + xlab("year")
# ggsave(filename = paste0("input/OMs_example_pol-nsea_roller-coaster.png"),
#        width = 17.1, height = 15, units = "cm", dpi = 300, type = "cairo-png")

### ------------------------------------------------------------------------ ###
### SRR plot
### ------------------------------------------------------------------------ ###
df <- data.frame(SSB = seq(0, 1000, length.out = 10000))
df$rec <- (12.990 * df$SSB / (90.909 + df$SSB)) / 
  ((12.990 * 1000 / (90.909 + 1000)))
ggplot(df, aes(x = SSB, y = rec)) + geom_line() + theme_bw() +
  geom_hline(yintercept = 0.7, colour = "red") +
  geom_vline(xintercept = B_lim, colour = "red") +
  ylab("relative recruitment")
ggsave(filename = paste0("input/OMs_srr.png"),
       width = 10, height = 6, units = "cm", dpi = 300, type = "cairo-png")

### ------------------------------------------------------------------------ ###
### comparison of length reference points
### ------------------------------------------------------------------------ ###

### go through stocks and extract life history reference points
lh_list <- lapply(1:15, function(x){
  load(paste0("input/stocks/perfect_knowledge/", x, ".RData"))
  return(attr(stk, "refpts")[,,,,, 1])
})
lhpar_list <- lapply(1:15, function(x){
  load(paste0("input/stocks/perfect_knowledge/", x, ".RData"))
  return(attr(stk, "lhpar")[,,,,, 1])
})
### calc Lc
source("functions/fFun.R")
lc_list <- unlist(lapply(1:15, function(x){
  load(paste0("input/stocks/perfect_knowledge/", x, ".RData"))
  mean(calc_Lc(attr(stk, "catch_len")[,ac(100)]))
}))
### extract LFeFmsy
df_refs <- data.frame(stock_pos = 1:15)
df_refs$LFeFmsy <- unlist(lapply(lh_list, function(x){
  x["LFeFmsy", 1]
}))

### calculate length with M & K
df_refs$Lcalc <- unlist(lapply(1:15, function(x){
  lhpar <- lhpar_list[[x]][, 1]
  c((lhpar["L_inf"] + 2 * lhpar["M"]/lhpar["K"] * lc_list[x]) / 
    (1 + 2 * lhpar["M"]/lhpar["K"]))
}))
### calculate with M/K = 1.5
df_refs$LMK15 <- unlist(lapply(1:15, function(x){
  lhpar <- lhpar_list[[x]][, 1]
  c((lhpar["L_inf"] + 2 * 1.5 * lc_list[x]) / 
      (1 + 2 * 1.5))
}))
### add stock names
stk_names <- unique(res_df[, c("stk_pos2", "stock")])
df_refs <- merge(df_refs, stk_names, by.x = "stock_pos", by.y = "stk_pos2")
df_refs$stock_pos <- NULL
### reshape
library(reshape2)
library(ggplot2)
df_refs2 <- melt(data = df_refs, id.vars = c("stock"))
### plot
p <- ggplot(data = df_refs2, aes(x = stock, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() +
  ylab("reference length") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(filename = paste0("input/reference_lengths.png"),
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")



### ------------------------------------------------------------------------ ###
### find lh groups
### ------------------------------------------------------------------------ ###

### 3.2.1 combinations & perfect knowledge
select = res_df[res_df$options == "option_f:a perfect_knowledge:TRUE option_r:a option_b:a" & 
                  res_df$fhist == "one-way", ]

stk_lst <- lapply(select$scenario, function(x){
  readRDS(paste0("output/perfect_knowledge/combined/", x, ".rds"))
})
# stk_lst <- res[as.character(select$scenario)]
names(stk_lst) <- select$stock
### median over all iterations
stk_lst2 <- lapply(stk_lst, function(x){
  qapply(x, iterMedians)
})
stk_lst2 <- FLStocks(stk_lst2)
p <- plot(stk_lst2) + theme_bw()

### remove pelagics
plot(stk_lst2[!names(stk_lst2) %in% c("her-nis", "san-ns4")])

### find highest one
plot(stk_lst2["smn-con"])
### smn-con is the one with the heigh SSB hill

### plot remaining ones
plot(stk_lst2[!names(stk_lst2) %in% c("her-nis", "san-ns4", "smn-con")])

### group of 2
plot(stk_lst2[names(stk_lst2) %in% c("lin-comb", "meg-4a6a")])

### group
plot(stk_lst2[names(stk_lst2) %in% c("lem-nsea", "tur-nsea")])

### remaining ones
plot(stk_lst2[!names(stk_lst2) %in% c("her-nis", "san-ns4", "smn-con",
                                      "lin-comb", "meg-4a6a", "lem-nsea",
                                      "tur-nsea")])

### group
plot(stk_lst2[names(stk_lst2) %in% c("ple-celt", "mut-comb")])

### remaining ones
plot(stk_lst2[!names(stk_lst2) %in% c("her-nis", "san-ns4", "smn-con",
                                      "lin-comb", "meg-4a6a", "lem-nsea",
                                      "tur-nsea", "mut-comb", "ple-celt")])

### whg-7e-k, own group, between pelagics and remaining non-pelagics
plot(stk_lst2["whg-7e-k"])

### group
plot(stk_lst2[names(stk_lst2) %in% c("had-iris", "nep-2829")])

### group
plot(stk_lst2[names(stk_lst2) %in% c("ang-ivvi", "pol-nsea")])

### last one
### group
plot(stk_lst2[names(stk_lst2) %in% c("ang-78ab")])

### groups:
g1 <- c("her-nis", "san-ns4")
g2 <- c("smn-con")
g3 <- c("lin-comb", "meg-4a6a")
g4 <- c("lem-nsea", "tur-nsea")
g5 <- c("ple-celt", "mut-comb")
g6 <- c("had-iris", "nep-2829")
g7 <- c("ang-ivvi", "pol-nsea")
g8 <- c("ang-78ab")
g9 <- c("whg-7e-k")

ssb_list <- lapply(stk_lst2, ssb)
ssb_list <- as.data.frame(ssb_list)
### add group
ssb_list$group <- NA
ssb_list$group[ssb_list$qname %in% g1] <- paste(g1, collapse = "_")
ssb_list$group[ssb_list$qname %in% g2] <- paste(g2, collapse = "_")
ssb_list$group[ssb_list$qname %in% g3] <- paste(g3, collapse = "_")
ssb_list$group[ssb_list$qname %in% g4] <- paste(g4, collapse = "_")
ssb_list$group[ssb_list$qname %in% g5] <- paste(g5, collapse = "_")
ssb_list$group[ssb_list$qname %in% g6] <- paste(g6, collapse = "_")
ssb_list$group[ssb_list$qname %in% g7] <- paste(g7, collapse = "_")
ssb_list$group[ssb_list$qname %in% g8] <- paste(g8, collapse = "_")
ssb_list$group[ssb_list$qname %in% g9] <- paste(g9, collapse = "_")

ggplot(ssb_list[!ssb_list$qname %in% c("her-nis", "san-ns4"), ],
       aes(x = year, y = data, colour = group, group = qname)) +
  geom_line() +
  theme_bw()


### with noise
select <- res_df[res_df$options == "option_f:a option_r:a option_b:a MK:1.5" & 
                   res_df$fhist == "one-way" & is.na(res_df$HCRmult) &
                   is.na(res_df$upper_constraint) &
                   is.na(res_df$b_z), ]
stk_lst <- lapply(select$scenario, function(x){
  readRDS(paste0("output/perfect_knowledge/combined/", x, ".rds"))
})
# stk_lst <- res[as.character(select$scenario)]
names(stk_lst) <- select$stock
### median over all iterations
stk_lst2 <- lapply(stk_lst, function(x){
  qapply(x, iterMedians)
})
stk_lst2 <- FLStocks(stk_lst2)
plot(stk_lst2) + theme_bw()

### remove pelagics
plot(stk_lst2[!names(stk_lst2) %in% c("her-nis", "san-ns4")])

### group
plot(stk_lst2[names(stk_lst2) %in% c("had-iris", "nep-2829", "ple-celt")])
plot(stk_lst2[names(stk_lst2) %in% c("smn-con")])
plot(stk_lst2[names(stk_lst2) %in% c("meg-4a6a")])
plot(stk_lst2[names(stk_lst2) %in% c("ang-ivvi", "pol-nsea", "ang-78ab",
                                     "lin-comb")])
plot(stk_lst2[names(stk_lst2) %in% c("whg-7e-k", "lem-nsea")])
plot(stk_lst2[names(stk_lst2) %in% c("tur-nsea", "mut-comb")])
### groups
g1 <-c("had-iris", "nep-2829", "ple-celt")
g2 <- c("smn-con")
g3 <- c("meg-4a6a")
g4 <- c("ang-ivvi", "pol-nsea", "ang-78ab",
        "lin-comb")
g5 <- c("whg-7e-k", "lem-nsea")
g6 <- c("tur-nsea", "mut-comb")

ssb_list <- lapply(stk_lst2, ssb)
ssb_list <- as.data.frame(ssb_list)
### add group
ssb_list$group <- NA
ssb_list$group[ssb_list$qname %in% g1] <- paste(g1, collapse = "_")
ssb_list$group[ssb_list$qname %in% g2] <- paste(g2, collapse = "_")
ssb_list$group[ssb_list$qname %in% g3] <- paste(g3, collapse = "_")
ssb_list$group[ssb_list$qname %in% g4] <- paste(g4, collapse = "_")
ssb_list$group[ssb_list$qname %in% g5] <- paste(g5, collapse = "_")
ssb_list$group[ssb_list$qname %in% g6] <- paste(g6, collapse = "_")

p <- ggplot(ssb_list[!ssb_list$qname %in% c("her-nis", "san-ns4"), ],
       aes(x = year, y = data, colour = group, group = qname)) +
  geom_line(size = 2) +
  theme_bw()

### plot for one-way trip
### same for roller-coaster
### with noise
select <- res_df[res_df$options == "option_f:a option_r:a option_b:a MK:1.5" & 
                      res_df$fhist == "one-way" & is.na(res_df$HCRmult) &
                      res_df$uncertainty == "observation_error" &
                      is.na(res_df$upper_constraint), ]
stk_lst <- lapply(select$scenario, function(x){
  readRDS(paste0("output/perfect_knowledge/combined/", x, ".rds"))
})
# stk_lst <- res[as.character(select$scenario)]
names(stk_lst) <- select$stock
### median over all iterations
stk_lst2 <- lapply(stk_lst, function(x){
  qapply(x, iterMedians)
})
stk_lst2 <- FLStocks(stk_lst2)
plot(stk_lst2) + theme_bw()
### apply same grouping as for one-way
ssb_list <- lapply(stk_lst2, ssb)
ssb_list <- as.data.frame(ssb_list)
### add group
ssb_list$group <- NA
ssb_list$group[ssb_list$qname %in% g1] <- paste(g1, collapse = " &\n")
ssb_list$group[ssb_list$qname %in% g2] <- paste(g2, collapse = " &\n")
ssb_list$group[ssb_list$qname %in% g3] <- paste(g3, collapse = " &\n")
ssb_list$group[ssb_list$qname %in% g4] <- paste(g4, collapse = " &\n")
ssb_list$group[ssb_list$qname %in% g5] <- paste(g5, collapse = " &\n")
ssb_list$group[ssb_list$qname %in% g6] <- paste(g6, collapse = " &\n")
ssb_list$group[!is.na(ssb_list$group)] <- paste0(ssb_list$group[!is.na(ssb_list$group)], "\n")
### add identifier whether stock has been used for further exploration
ssb_list$used <- FALSE
ssb_list$used[ssb_list$qname %in% c("pol-nsea", "nep-2829", "tur-nsea",
                                    "whg-7e-k")] <- TRUE

### add catch
ssb_list$catch <- as.data.frame(lapply(stk_lst2, catch))$data
ssb_list$fbar <- as.data.frame(lapply(stk_lst2, fbar))$data
ssb_list$rec <- as.data.frame(lapply(stk_lst2, rec))$data
ssb_list$ssb <- as.data.frame(lapply(stk_lst2, ssb))$data

### reshape
ssb_list_res <- melt(data = ssb_list, 
                        id.vars = c("year", "qname", "group", "used"),
                        measure.vars = c("ssb", "catch", "rec", "fbar"))
### plot ssb
p <- ggplot(ssb_list_res[!ssb_list$qname %in% c("her-nis", "san-ns4") &
                              ssb_list_res$variable == "ssb", ],
            aes(x = year, y = value, colour = group, group = qname)) +
  geom_line(size = 1) +
  theme_bw() +
  facet_wrap(~variable, scales = "free_y") +
  geom_vline(xintercept = 100) +
  #scale_linetype_manual(values = c("TRUE" = 6, "FALSE" = 1)) + 
  labs(y = "")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "groups_ssb_one-way.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")

### same for roller-coaster
### with noise
select_rc <- res_df[res_df$options == "option_f:a option_r:a option_b:a MK:1.5" & 
                    res_df$fhist == "roller-coaster" & is.na(res_df$HCRmult) &
                    res_df$uncertainty == "observation_error" &
                    is.na(res_df$upper_constraint), ]
stk_lst_rc <- lapply(select_rc$scenario, function(x){
  readRDS(paste0("output/perfect_knowledge/combined/", x, ".rds"))
})
# stk_lst <- res[as.character(select$scenario)]
names(stk_lst_rc) <- select_rc$stock
### median over all iterations
stk_lst2_rc <- lapply(stk_lst_rc, function(x){
  qapply(x, iterMedians)
})
stk_lst2_rc <- FLStocks(stk_lst2_rc)
plot(stk_lst2_rc) + theme_bw()
### apply same grouping as for one-way
ssb_list_rc <- lapply(stk_lst2_rc, ssb)
ssb_list_rc <- as.data.frame(ssb_list_rc)
### add group
ssb_list_rc$group <- NA
ssb_list_rc$group[ssb_list_rc$qname %in% g1] <- paste(g1, collapse = " &\n")
ssb_list_rc$group[ssb_list_rc$qname %in% g2] <- paste(g2, collapse = " &\n")
ssb_list_rc$group[ssb_list_rc$qname %in% g3] <- paste(g3, collapse = " &\n")
ssb_list_rc$group[ssb_list_rc$qname %in% g4] <- paste(g4, collapse = " &\n")
ssb_list_rc$group[ssb_list_rc$qname %in% g5] <- paste(g5, collapse = " &\n")
ssb_list_rc$group[ssb_list_rc$qname %in% g6] <- paste(g6, collapse = " &\n")
ssb_list_rc$group[!is.na(ssb_list_rc$group)] <- paste0(ssb_list_rc$group[!is.na(ssb_list_rc$group)], "\n")

### add catch
ssb_list_rc$catch <- as.data.frame(lapply(stk_lst2_rc, catch))$data
ssb_list_rc$fbar <- as.data.frame(lapply(stk_lst2_rc, fbar))$data
ssb_list_rc$rec <- as.data.frame(lapply(stk_lst2_rc, rec))$data
ssb_list_rc$ssb <- as.data.frame(lapply(stk_lst2_rc, ssb))$data

### add identifier whether stock has been used for further exploration
ssb_list_rc$used <- FALSE
ssb_list_rc$used[ssb_list_rc$qname %in% c("pol-nsea", "nep-2829", "tur-nsea",
                                    "whg-7e-k")] <- TRUE

### reshape
ssb_list_rc_res <- melt(data = ssb_list_rc, 
                  id.vars = c("year", "qname", "group", "used"),
                  measure.vars = c("ssb", "catch", "rec", "fbar"))
### plot ssb
p <- ggplot(ssb_list_rc_res[!ssb_list_rc$qname %in% c("her-nis", "san-ns4") &
                               ssb_list_rc_res$variable == "ssb", ],
             aes(x = year, y = value, colour = group, group = qname)) +
    geom_line(size = 1) +
    theme_bw() +
    facet_wrap(~variable, scales = "free_y") +
    geom_vline(xintercept = 100) +
    #scale_linetype_manual(values = c("TRUE" = 6, "FALSE" = 1)) + 
    labs(y = "")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "groups_ssb_roller-coaster.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")

### plot all quants
(p <- ggplot(ssb_list_rc_res[!ssb_list_rc$qname %in% c("her-nis", "san-ns4"), ],
             aes(x = year, y = value, colour = group, group = qname)) +
  geom_line(size = 2) +
  theme_bw() +
  facet_wrap(~variable, scales = "free_y") )
### plot fbar
(pfbar <- ggplot(ssb_list_rc_res[!ssb_list_rc$qname %in% c("her-nis", "san-ns4") &
                         ssb_list_rc_res$variable == "fbar", ],
       aes(x = year, y = value, colour = group, group = qname)) +
  geom_line(size = 2) +
  theme_bw() +
  facet_wrap(~variable, scales = "free_y"))
pfbar + ylim(0, 0.5)
 

p2 <- ggplot(ssb_list_rc[!ssb_list_rc$qname %in% c("her-nis", "san-ns4"), ],
            aes(x = year, y = ssb, colour = group, group = qname)) +
  geom_line(size = 2) +
  theme_bw()
p2

### load life-history parameter input
library(FLife)
data(wklife)
wklife$sex <- as.character(wklife$sex)
wklife$sex[15] <- "M"; wklife$a[15] <- 0.00028; wklife$b[15] <- 3.229; 
wklife$linf[15] <- 70; wklife$l50[15] <- 28.4; wklife$k[15] <- 0.2
wklife$group <- NA
wklife$group[wklife$stock %in% g1] <- paste(g1, collapse = "_")
wklife$group[wklife$stock %in% g2] <- paste(g2, collapse = "_")
wklife$group[wklife$stock %in% g3] <- paste(g3, collapse = "_")
wklife$group[wklife$stock %in% g4] <- paste(g4, collapse = "_")
wklife$group[wklife$stock %in% g5] <- paste(g5, collapse = "_")
wklife$group[wklife$stock %in% g6] <- paste(g6, collapse = "_")

### get more lhist params
### M
wklife$M <- unlist(lapply(stk_lst, function(x){
  median(c(attr(x, "lhpar")["M"]))
}))
### max age
wklife$amax <- unlist(lapply(stk_lst, function(x){
  median(c(attr(x, "lhpar")["max_age"]))
}))
### l50
wklife$l50[12] <- 42.3
### a50
wklife$a50 <- unlist(lapply(stk_lst, function(x){
  median(c(attr(x, "lhpar")["a50"]))
}))
wklife$l50linf <- with(wklife, l50/linf)
wklife$MK <- with(wklife, M/k)

### calculate mature M
wklife$M_mat <-  unlist(lapply(1:15, function(x){
  weighted.mean(x = m(stk_lst[[x]])[,1,,,, 1], mat(stk_lst[[x]])[,1,,,, 1])
}))
### M/K
wklife$MK_mat <- with(wklife, M_mat / k)

### modify and save
wklife2 <- wklife
wklife2$group <- as.factor(wklife2$group)
levels(wklife2$group) <- 1:length(levels(wklife2$group))
wklife2 <- wklife2[order(wklife2$group), c("name", "common", "area", "stock",
                                           "a", "b", "linf", "l50", "a50",
                                           "t0", "k", "amax", "M_mat", 
                                           "MK_mat", "l50linf", "group")]
write.csv(x = wklife2, file = paste0("output/perfect_knowledge/plots/3.2.1",
                                     "/subgroup/lh.csv"), row.names = FALSE)


#stk_list
wklife_df <- melt(data = wklife, 
                id.vars = c("stock", "group"),
                measure.vars = c("amax", "linf", "l50linf", "l50", 
                                 "k", "MK_mat"))
ggplot(data = wklife_df,
       aes(x = value, y = 1, fill = group, colour = group, shape = group)) +
  #geom_bar(stat = "identity") +
  geom_point(size = 10, alpha = 0.5) +
  #geom_jitter(size = 5, alpha = 0.5, width = 0, height = 0.01) +
  facet_wrap(~ variable, scale = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
# ggplot(data = wklife_df[wklife_df$variable %in% c("MK", "MK_mat"), ],
#             aes(x = stock, y = value, fill = group)) +
#   geom_bar(stat = "identity") +
#   facet_wrap(~ variable, scale = "free") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# ggplot(data = wklife_df[wklife_df$variable %in% c("MK", "MK_mat"), ],
#             aes(x = value, y = 1, fill = group, colour = group)) +
#   #geom_bar(stat = "identity") +
#   geom_point(size = 10, alpha = 0.5) +
#   #geom_jitter(size = 5, alpha = 0.5, width = 0, height = 0.01) +
#   facet_wrap(~ variable, scale = "free") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

### grouping visible in values?
wklife[order(wklife$l50linf), c("stock", "MK_mat", "l50linf", "group")]


plot(stk_lst2[c("ang-ivvi", "pol-nsea", "ang-78ab", "lin-comb")])
plot(stk_lst2[c("had-iris", "ple-celt", "nep-2829")])

### stock for futher analysis
stocks <- c(2,6,11,15)

### ------------------------------------------------------------------------ ###
### risk plots 3.2.1 combs for reduced stocks
### ------------------------------------------------------------------------ ###

### multiplier
df_plot <- res_df[res_df$catch_rule == "3.2.1" &
                    res_df$stk_pos2 %in% c(2,6,11,15) & 
                    res_df$uncertainty == "observation_error" &
                    res_df$options == "option_f:a option_r:a option_b:a MK:1.5" &
                    is.na(res_df$upper_constraint) &
                    is.na(res_df$b_w), ]
df_plot$HCRmult[is.na(df_plot$HCRmult)] <- 1
### reshape
df_plot <- melt(data = df_plot, 
                id.vars = c("scenario", "stock", "fhist", "HCRmult"),
                measure.vars = c("ssb_below_blim_total", "collapse_iter", 
                                 "rel_yield"))
levels(df_plot$variable) <- c("p(SSB < Blim)", "iteration collapse risk", "relative yield")

p <- ggplot(data = df_plot, 
            aes(x = HCRmult, y = value, colour = stock)) +
  geom_line() + geom_point() +
  theme_bw() +
  ylim(0, NA) +
  facet_grid(variable ~ fhist, scale = "free") +
  xlab("advice multiplier") + ylab("")
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_multiplier.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")


### z
df_plot <- res_df[res_df$catch_rule == "3.2.1" &
                    res_df$stk_pos2 %in% c(2,6,11,15) & 
                    res_df$uncertainty == "observation_error" &
                    grepl(x = res_df$options, 
                          pattern = "b_z") &
                    is.na(res_df$upper_constraint) &
                    is.na(res_df$HCRmult), ]
### reshape
df_plot <- melt(data = df_plot, 
                id.vars = c("scenario", "stock", "fhist", "b_z"),
                measure.vars = c("ssb_below_blim_total", "collapse_iter", 
                                 "rel_yield"))
levels(df_plot$variable) <- c("p(SSB < Blim)", "iteration collapse risk", "relative yield")

p <- ggplot(data = df_plot, 
            aes(x = b_z, y = value, colour = stock)) +
  geom_line() + geom_point() +
  theme_bw() +
  ylim(0, NA) +
  facet_grid(variable ~ fhist, scale = "free") +
  xlab("exponent of component b") + ylab("")
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_exponent_z.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")

### add relative values
df_plot2 <- aggregate(value ~ stock + fhist + variable, data = df_plot,
                      FUN = max)
names(df_plot2)[length(df_plot2)] <- "value_sum"
df_plot <- merge(df_plot, df_plot2, all = TRUE)
df_plot$value_rel <- with(df_plot, value/value_sum)
### plot relative values
p <- ggplot(data = df_plot, 
            aes(x = b_z, y = value_rel, colour = stock)) +
  geom_line() + geom_point() +
  theme_bw() +
  facet_grid(variable ~ fhist, scale = "free") +
  xlab("exponent of component b") + ylab("relative values")
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_exponent_z_rel.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")







### TAC constraint (upper)
df_plot <- res_df[res_df$catch_rule == "3.2.1" &
                    res_df$stk_pos2 %in% c(2,6,11,15) & 
                    res_df$uncertainty == "observation_error" &
                    res_df$options == "option_f:a option_r:a option_b:a MK:1.5" &
                    is.na(res_df$b_z) &
                    is.na(res_df$HCRmult), ]
df_plot$upper_constraint[is.na(df_plot$upper_constraint)] <- 3.05
### add rows without constraint
df_add <- df_plot[df_plot$upper_constraint == 3.05, ]
df_add$ssb_below_blim_total <- df_add$collapse_iter <- df_add$rel_yield <- NA
df_add$upper_constraint <- 3.025
df_plot <- rbind(df_plot, df_add)
 
### reshape
df_plot <- melt(data = df_plot, 
                id.vars = c("scenario", "stock", "fhist", "upper_constraint"),
                measure.vars = c("ssb_below_blim_total", "collapse_iter",
                                 "rel_yield"))
levels(df_plot$variable) <- c("p(SSB < Blim)", "iteration collapse risk", "relative yield")

p <- ggplot(data = df_plot, 
            aes(x = upper_constraint, y = value, colour = stock,
                shape = upper_constraint > 3)) +
  geom_line() + geom_point() +
  theme_bw() +
  ylim(0, NA) +
  facet_grid(variable ~ fhist, scale = "free") +
  xlab("upper constraint") + ylab("risk") + 
  xlim(1, NA) +
  scale_shape_discrete("", labels = c("constraint", "without\nconstraint"))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_upper_constraint.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")


### ------------------------------------------------------------------------ ###
### risk plots 3.2.1 reduced stocks, combinations of modifications:constraints
### ------------------------------------------------------------------------ ###

### load scenario results
df_plot <- res_df[c(913:1032), ]

### reshape
df_plot2 <- melt(data = df_plot, 
                id.vars = c("scenario", "stock", "lower_constraint", 
                            "upper_constraint", "fhist"),
                measure.vars = c("ssb_below_blim_total", "collapse_iter", 
                                 "rel_yield"))
levels(df_plot2$variable) <- c("p(SSB < Blim)", "iteration collapse", "relative yield")

p <- ggplot(data = df_plot2, 
            aes(x = lower_constraint, y = upper_constraint)) +
  scale_fill_gradient("value", limits = c(0,1),
                      low = "green", high = "red") +
  labs(x = "lower limit", y = "upper limit") +
  facet_grid(fhist+variable ~ stock) +
  geom_raster(aes(fill = value)) +
  geom_text(aes(label = round(value, 2))) +
  theme_bw()
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_limits_combs_grid.png"), plot = p,
       width = 25, height = 15, units = "cm", dpi = 300, type = "cairo-png")
### line plots
p <- ggplot(data = df_plot2, 
             aes(x = lower_constraint, y = value, linetype = fhist, 
                 colour = as.factor(upper_constraint),
                 shape = fhist)) +
  geom_line() + geom_point() + 
  labs(x = "lower limit") +
  facet_grid(variable ~ stock, scales = "free") +
  theme_bw() +
  scale_color_discrete("upper limit") +
  scale_shape_discrete("fishing\nhistory") + 
  scale_linetype_discrete("fishing\nhistory") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_limits_combs.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
### one-way only
p <- ggplot(data = df_plot2[df_plot2$fhist == "one-way", ], 
            aes(x = lower_constraint, y = value, 
                colour = as.factor(upper_constraint))) +
  geom_line() + geom_point() + 
  labs(x = "lower limit") +
  facet_grid(variable ~ stock, scales = "free") +
  theme_bw() +
  scale_color_discrete("upper limit\none-way") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_limits_combs_one-way.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
### roller-coaster only
p <- ggplot(data = df_plot2[df_plot2$fhist == "roller-coaster", ], 
            aes(x = lower_constraint, y = value, 
                colour = as.factor(upper_constraint))) +
  geom_line() + geom_point() + 
  labs(x = "lower limit") +
  facet_grid(variable ~ stock, scales = "free") +
  theme_bw() +
  scale_color_discrete("upper limit\nroller-coaster") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_limits_combs_roller-coaster.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
### same, but inverted
p <- ggplot(data = df_plot2, 
            aes(x = upper_constraint, y = value, linetype = fhist, 
                colour = as.factor(lower_constraint),
                shape = fhist)) +
  geom_line() + geom_point() + 
  labs(x = "upper limit") +
  facet_grid(variable ~ stock, scales = "free") +
  theme_bw() +
  scale_color_discrete("lower limit") +
  scale_shape_discrete("fishing\nhistory") + 
  scale_linetype_discrete("fishing\nhistory") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_limits_combs2.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
### one-way only
p <- ggplot(data = df_plot2[df_plot2$fhist == "one-way", ], 
            aes(x = upper_constraint, y = value,
                colour = as.factor(lower_constraint))) +
  geom_line() + geom_point() + 
  labs(x = "upper limit") +
  facet_grid(variable ~ stock, scales = "free") +
  theme_bw() +
  scale_color_discrete("lower limit\none-way") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_limits_combs2_one_way.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
### roller-coaster only
p <- ggplot(data = df_plot2[df_plot2$fhist == "roller-coaster", ], 
            aes(x = upper_constraint, y = value,
                colour = as.factor(lower_constraint))) +
  geom_line() + geom_point() + 
  labs(x = "upper limit") +
  facet_grid(variable ~ stock, scales = "free") +
  theme_bw() +
  scale_color_discrete("lower limit\nroller-coaster") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_limits_combs2_roller-coaster.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")


### ------------------------------------------------------------------------ ###
### risk plots 3.2.1 reduced stocks, combinations of modifications: z and x
### ------------------------------------------------------------------------ ###

### load scenario results
df_plot <- res_df[res_df$scenario %in% c(1033:1152), ]

### reshape
df_plot2 <- melt(data = df_plot, 
                 id.vars = c("scenario", "stock", "HCRmult", 
                             "b_z", "fhist"),
                 measure.vars = c("ssb_below_blim_total", "collapse_iter", 
                                  "rel_yield"))
levels(df_plot2$variable) <- c("p(SSB < Blim)", "iteration collapse", "relative yield")

p <- ggplot(data = df_plot2, 
            aes(x = HCRmult, y = b_z)) +
  scale_fill_gradient("value", limits = c(0,1),
                      low = "green", high = "red") +
  labs(x = "advice multiplier", y = "b exponent") +
  facet_grid(fhist+variable ~ stock) +
  geom_raster(aes(fill = value)) +
  geom_text(aes(label = round(value, 2))) +
  theme_bw()
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_mult_b_combs_grid.png"), plot = p,
       width = 25, height = 15, units = "cm", dpi = 300, type = "cairo-png")
### line plot
p <- ggplot(data = df_plot2, 
            aes(x = HCRmult, y = value, linetype = fhist, 
                colour = as.factor(b_z),
                shape = fhist)) +
  geom_line() + geom_point() + 
  labs(x = "advice multiplier") +
  facet_grid(variable ~ stock, scales = "free") +
  theme_bw() +
  scale_color_discrete("b exponent") +
  scale_shape_discrete("fishing\nhistory") + 
  scale_linetype_discrete("fishing\nhistory") +
  scale_x_continuous(breaks = c(.8, 0.9, 1))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_mult_z_combs.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
### for one-way only
p <- ggplot(data = df_plot2[df_plot2$fhist == "one-way", ], 
            aes(x = HCRmult, y = value, colour = as.factor(b_z))) +
  geom_line() + geom_point() + 
  labs(x = "advice multiplier") +
  facet_grid(variable ~ stock, scales = "free") +
  theme_bw() +
  scale_color_discrete("b exponent\none-way") +
  scale_x_continuous(breaks = c(.8, 0.9, 1))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_mult_z_combs_one-way.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
### for roller-coaster only
p <- ggplot(data = df_plot2[df_plot2$fhist == "roller-coaster", ], 
            aes(x = HCRmult, y = value, colour = as.factor(b_z))) +
  geom_line() + geom_point() + 
  labs(x = "advice multiplier") +
  facet_grid(variable ~ stock, scales = "free") +
  theme_bw() +
  scale_color_discrete("b exponent\nroller-coaster") +
  scale_x_continuous(breaks = c(.8, 0.9, 1))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_mult_z_combs_roller-coaster.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
### same, but inverted
p <- ggplot(data = df_plot2, 
            aes(x = b_z, y = value, linetype = fhist, 
                colour = as.factor(HCRmult),
                shape = fhist)) +
  geom_line() + geom_point() + 
  labs(x = "b exponent") +
  facet_grid(variable ~ stock, scales = "free") +
  theme_bw() +
  scale_color_discrete("advice multiplier") +
  scale_shape_discrete("fishing\nhistory") + 
  scale_linetype_discrete("fishing\nhistory") +
  scale_x_continuous(breaks = 1:3)
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_mult_z_combs2.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
### for one-way only
p <- ggplot(data = df_plot2[df_plot2$fhist == "one-way", ], 
            aes(x = b_z, y = value, colour = as.factor(HCRmult))) +
  geom_line() + geom_point() + 
  labs(x = "b exponent") +
  facet_grid(variable ~ stock, scales = "free") +
  theme_bw() +
  scale_color_discrete("advice multiplier\none-way") +
  scale_x_continuous(breaks = 1:3)
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_mult_z_combs2_one-way.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
### for roller-coaster only
p <- ggplot(data = df_plot2[df_plot2$fhist == "roller-coaster", ], 
            aes(x = b_z, y = value, colour = as.factor(HCRmult))) +
  geom_line() + geom_point() + 
  labs(x = "b exponent") +
  facet_grid(variable ~ stock, scales = "free") +
  theme_bw() +
  scale_color_discrete("advice multiplier\nroller-coaster") +
  scale_x_continuous(breaks = 1:3)
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_mult_z_combs2_roller_coaster.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")

### same plots, but with standardised values
### add relative values
### multiplier
df_plot_x <- aggregate(value ~ stock + fhist + variable + HCRmult, 
                       data = df_plot2, FUN = max)
names(df_plot_x)[length(df_plot_x)] <- "value_max_HCRmult"
df_plot <- merge(df_plot_x, df_plot2, all = TRUE)
df_plot$value_rel_HCRmult <- with(df_plot, value/value_max_HCRmult)
### exponent z
df_plot_z <- aggregate(value ~ stock + fhist + variable + b_z, 
                       data = df_plot2, FUN = max)
names(df_plot_z)[length(df_plot_z)] <- "value_max_b_z"
df_plot <- merge(df_plot_z, df_plot, all = TRUE)
df_plot$value_rel_b_z <- with(df_plot, value/value_max_b_z)
### correct NaNs
df_plot$value_rel_HCRmult[is.na(df_plot$value_rel_HCRmult)] <- 1

### multiplier on x-axis
p <- ggplot(data = df_plot, 
            aes(x = HCRmult, y = value_rel_b_z, linetype = fhist, 
                colour = as.factor(b_z),
                shape = fhist)) +
  geom_line() + geom_point() + 
  labs(x = "advice multiplier") +
  facet_grid(variable ~ stock, scales = "free") +
  theme_bw() +
  scale_color_discrete("b exponent") +
  scale_shape_discrete("fishing\nhistory") + 
  scale_linetype_discrete("fishing\nhistory") +
  scale_x_continuous(breaks = c(.8, 0.9, 1)) +
  ylab("standardised values")
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_mult_z_combs_rel.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
### exponent on x-axis
p <- ggplot(data = df_plot, 
            aes(x = b_z, y = value_rel_HCRmult, linetype = fhist, 
                colour = as.factor(HCRmult),
                shape = fhist)) +
  geom_line() + geom_point() + 
  labs(x = "b exponent") +
  facet_grid(variable ~ stock, scales = "free") +
  theme_bw() +
  scale_color_discrete("advice multiplier") +
  scale_shape_discrete("fishing\nhistory") + 
  scale_linetype_discrete("fishing\nhistory") +
  scale_x_continuous(breaks = 1:3) +
  ylab("standardised values")
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_mult_z_combs2_rel.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")

### ------------------------------------------------------------------------ ###
### plot default ra fa ba & noise ssb and catch
### ------------------------------------------------------------------------ ###
res_plot <- res_df[res_df$options == "option_f:a option_r:a option_b:a MK:1.5" & 
                     res_df$fhist == "one-way" &
                     is.na(res_df$HCRmult) & is.na(res_df$upper_constraint), ]

### load stocks
res_stks <- lapply(res_plot$scenario, function(x){
  stk <- readRDS(paste0("output/perfect_knowledge/combined/", x, ".rds"))
})
res_stks <- lapply(seq_along(res_stks), function(x){
  res <- FLQuants(SSB = apply(ssb(res_stks[[x]]), 1:2, median),
                  catch = apply(catch(res_stks[[x]]), 1:2, median))
  res <- as.data.frame(res)
  res$year <- res$year - 100
  res$stock <- ac(res_plot$stock[x])
  return(res)
})
res_stks <- do.call(rbind, res_stks)
res_stks$qname <- as.factor(res_stks$qname)
levels(res_stks$qname) <- c("stock size", "catch")
### plot
p <- ggplot(data = res_stks[res_stks$year >= -1, ],
            aes(x = year, y = data, colour = stock)) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~ qname, scales = "free", nrow = 1) +
  ylab("") + 
  theme_bw()
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/default.png"), 
       plot = p,
       width = 10, height = 4, units = "cm", dpi = 300, type = "cairo-png")

### ------------------------------------------------------------------------ ###
### more stocks: plot list of stocks
### not used for WKLIFE
### ------------------------------------------------------------------------ ###
### 3.2.1 f:a r:a b:a perfect knowledge
plot_lst(res_df[1153:1166, ],
         name = "more/3.2.1_perfect_one-way")
plot_lst(res_df[1167:1180, ],
         name = "more/3.2.1_perfect_roller_coaster")
### 3.2.1 f:a r:a b:a observation error & assumptions
plot_lst(res_df[1181:1194, ],
         name = "more/3.2.1_error_one-way")
plot_lst(res_df[1195:1208, ],
         name = "more/3.2.1_error_roller_coaster")
### 3.2.2 observation error
plot_lst(res_df[1209:1222, ],
         name = "more/3.2.2_error_one-way")
plot_lst(res_df[1223:1236, ],
         name = "more/3.2.2_error_roller_coaster")

### risk plots: 3.2.1 obs error all stocks ####
### ------------------------------------------------------------------------ ###
### risk plots for default parametrization
df_plot <- res_df[res_df$catch_rule == "3.2.1" &
                    res_df$options == "option_f:a option_r:a option_b:a MK:1.5" &
                    is.na(res_df$HCRmult) & is.na(res_df$upper_constraint), ]

### reshape
df_plot <- melt(data = df_plot, 
                id.vars = c("scenario", "stock", "fhist", "options", "stk_pos"),
                measure.vars = c("ssb_below_blim_total", "collapse_iter", 
                                 "rel_yield"))
levels(df_plot$variable) <- c("p(SSB<Blim)", 
                              ("iter collapse"), "rel. yield")
### set stock group
df_plot$group[df_plot$stk_pos <= 30] <- "wklife"
df_plot$group[is.na(df_plot$group)] <- "new"

### plot
p <- ggplot(data = df_plot, 
            aes(x = stock, y = value,  fill = stock)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
  geom_text(aes(label = round(value, 2), y = value + 0.05), position = "dodge", size = 2) + 
  facet_grid(variable + fhist ~ group, scale = "free_x") +
  theme_bw() +
  ylim(0, NA) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/more/",
                         "risks_3.2.1_error.png"), plot = p,
       width = 20, height = 22, units = "cm", dpi = 300, type = "cairo-png")

### risk plots: 3.2.2 obs error all stocks ####
### ------------------------------------------------------------------------ ###
### 
df_plot <- subset(res_df, catch_rule == "3.2.2" & 
                    uncertainty == "observation_error" &
                    is.na(HCRmult) & is.na(w))

### reshape
df_plot <- melt(data = df_plot, 
                id.vars = c("scenario", "stock", "fhist", "options", "stk_pos"),
                measure.vars = c("ssb_below_blim_total", "collapse_iter", 
                                 "rel_yield"))
levels(df_plot$variable) <- c("p(SSB<Blim)", 
                              ("iter collapse"), "rel. yield")
### set stock group
df_plot$group[df_plot$stk_pos <= 30] <- "wklife"
df_plot$group[is.na(df_plot$group)] <- "new"

### plot
p <- ggplot(data = df_plot, 
            aes(x = stock, y = value,  fill = stock)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
  geom_text(aes(label = round(value, 2), y = value + 0.05), position = "dodge", size = 2) + 
  facet_grid(variable + fhist ~ group, scale = "free_x") +
  theme_bw() +
  ylim(0, NA) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/more/",
                         "risks_3.2.2_error.png"), plot = p,
       width = 20, height = 22, units = "cm", dpi = 300, type = "cairo-png")


### ------------------------------------------------------------------------ ###
### older attempts
### ------------------------------------------------------------------------ ###

# ### show all files in trial_runs folder
# files <- list.files("trial_runs/")
# 
# ### keep only RData files with number in front
# files <- files[grepl(pattern = "[0-9]{1,}.RData", x = files)]
# 
# ### sort according to number
# files <- files[order(an(gsub(pattern = "[^0-9]", replacement = "", x = files)))]
# 
# 
# ### load files
# stks <- lapply(files, function(x){
#   (function(){
#     load(paste0("trial_runs/", x))
#     return(stk)
#   })()
# })
# 
# ### coerce into FLStocks
# stks <- FLStocks(stks)
# 
# ### correct names
# names(stks)[1:12] <- paste0("lemon_sole_oneway_", names(stks)[1:12])
# names(stks)[13:length(stks)] <- paste0("lemon_sole_rollercoaster_", 
#                                        names(stks)[13:length(stks)])
# 
# 
# plot(stks) + theme(legend.position = "bottom")
# plot(stks[c(6,12+6)]) + theme(legend.position = "bottom")
# 
# ### check if catch matches targeted catch
# which(catch(stks[[1]])[,ac(101:200)] / (attr(stks[[1]], "tracking")["IEM",ac(100:199)]) > 1.001)
# ### catch in year 107, iterations 1:6
# 
# ### example
# catch(stks[[1]])[,ac(101:109),,,,1]
# attr(stks[[1]], "tracking")["IEM",ac(100:108),,,,1]
# plot(stks[[1]][,ac(100:110),,,,1:6]) + geom_vline(xintercept = 107, colour = "red")
# 


#stk210 <- readRDS(paste0("output/perfect_knowledge/combined/", 204, ".rds"))
### ------------------------------------------------------------------------ ###
### older attempts
### ------------------------------------------------------------------------ ###




### ------------------------------------------------------------------------ ###
### 3.2.1 all stocks - constraints ####
### ------------------------------------------------------------------------ ###

### load scenario results
df_plot <- res_df[res_df$scenario %in% 1237:2164, ]

### reshape
df_plot2 <- melt(data = df_plot, 
                 id.vars = c("scenario", "stock", "lower_constraint", 
                             "upper_constraint", "fhist"),
                 measure.vars = c("ssb_below_blim_total", "collapse_iter", 
                                  "rel_yield"))
levels(df_plot2$variable) <- c("p(SSB < Blim)", "iteration collapse", "relative yield")

p <- ggplot(data = df_plot2[df_plot2$stock == unique(df_plot2$stock)[1],], 
            aes(x = lower_constraint, y = upper_constraint)) +
  scale_fill_gradient("value", limits = c(0,1),
                      low = "green", high = "red") +
  labs(x = "lower limit", y = "upper limit") +
  facet_grid(fhist+variable ~ stock) +
  geom_raster(aes(fill = value)) +
  geom_text(aes(label = round(value, 2))) +
  theme_bw()
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_limits_combs_grid.png"), plot = p,
       width = 25, height = 15, units = "cm", dpi = 300, type = "cairo-png")
### line plots
p <- ggplot(data = df_plot2[df_plot2$stock == unique(df_plot2$stock)[1],], 
            aes(x = lower_constraint, y = value, linetype = fhist, 
                colour = as.factor(upper_constraint),
                shape = fhist)) +
  geom_line() + geom_point() + 
  labs(x = "lower limit") +
  facet_grid(variable ~ stock, scales = "free") +
  theme_bw() +
  scale_color_discrete("upper limit") +
  scale_shape_discrete("fishing\nhistory") + 
  scale_linetype_discrete("fishing\nhistory") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_limits_combs.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
### one-way only
p <- ggplot(data = df_plot2[df_plot2$fhist == "one-way", ], 
            aes(x = lower_constraint, y = value, 
                colour = as.factor(upper_constraint))) +
  geom_line() + geom_point() + 
  labs(x = "lower limit") +
  facet_grid(variable ~ stock, scales = "free") +
  theme_bw() +
  scale_color_discrete("upper limit\none-way") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_limits_combs_one-way.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
### roller-coaster only
p <- ggplot(data = df_plot2[df_plot2$fhist == "roller-coaster", ], 
            aes(x = lower_constraint, y = value, 
                colour = as.factor(upper_constraint))) +
  geom_line() + geom_point() + 
  labs(x = "lower limit") +
  facet_grid(variable ~ stock, scales = "free") +
  theme_bw() +
  scale_color_discrete("upper limit\nroller-coaster") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_limits_combs_roller-coaster.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
### same, but inverted
p <- ggplot(data = df_plot2, 
            aes(x = upper_constraint, y = value, linetype = fhist, 
                colour = as.factor(lower_constraint),
                shape = fhist)) +
  geom_line() + geom_point() + 
  labs(x = "upper limit") +
  facet_grid(variable ~ stock, scales = "free") +
  theme_bw() +
  scale_color_discrete("lower limit") +
  scale_shape_discrete("fishing\nhistory") + 
  scale_linetype_discrete("fishing\nhistory") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_limits_combs2.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
### one-way only
p <- ggplot(data = df_plot2[df_plot2$fhist == "one-way", ], 
            aes(x = upper_constraint, y = value,
                colour = as.factor(lower_constraint))) +
  geom_line() + geom_point() + 
  labs(x = "upper limit") +
  facet_grid(variable ~ stock, scales = "free") +
  theme_bw() +
  scale_color_discrete("lower limit\none-way") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_limits_combs2_one_way.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
### roller-coaster only
p <- ggplot(data = df_plot2[df_plot2$fhist == "roller-coaster", ], 
            aes(x = upper_constraint, y = value,
                colour = as.factor(lower_constraint))) +
  geom_line() + geom_point() + 
  labs(x = "upper limit") +
  facet_grid(variable ~ stock, scales = "free") +
  theme_bw() +
  scale_color_discrete("lower limit\nroller-coaster") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/subgroup/",
                         "risks_limits_combs2_roller-coaster.png"), plot = p,
       width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
