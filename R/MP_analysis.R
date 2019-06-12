library(FLCore)
library(ggplotFL)
library(Cairo)
library(reshape2)
library(dplyr)

### set up parallel computing environment
library(doParallel)
cl <- makeCluster(parallel::detectCores())
registerDoParallel(cl)
getDoParWorkers()
getDoParName()

### ------------------------------------------------------------------------ ###
### functions ####
### ------------------------------------------------------------------------ ###

### select R scripts from functions folder
load_files <- list.files("functions/")
load_files <- load_files[grepl(pattern = "*.R$", x = load_files)]
### source the scripts
invisible(lapply(paste0("functions/", load_files), source))

source("MP_scenarios.R")

### ------------------------------------------------------------------------ ###
### find result files ####
### ------------------------------------------------------------------------ ###

### get file names
#path_res <- "output/perfect_knowledge/"
path_res <- "/gpfs/afmcefas/simonf/output/"
files <- list.files(path = path_res, pattern = "^[0-9]{1,}_[0-9]{1,}\\.rds$")

### split files name into parts (scenario, part, extension)
files_desc <- strsplit(x = files, split = "_|[.]")

### create data frame with file info
df <- data.frame(scenario = as.numeric(sapply(files_desc, "[[", 1)),
                 part = as.numeric(sapply(files_desc, "[[", 2)),
                 file_name = files)

### sort
df <- df[order(df$scenario, df$part), ]

### show number of parts for each scenario
count(df, scenario) %>%
  filter(n == 10)

### keep only scenarios that have fininshed, i.e. 10 parts
df_process <- df %>% 
  full_join(count(df, scenario)) %>%
  filter(n == 10)

### ------------------------------------------------------------------------ ###
### combine parts ####
### ------------------------------------------------------------------------ ###

. <- foreach(scenario = unique(df_process$scenario), 
               .packages = c("FLCore"), 
               .export = ls()) %dopar% {
  
  ### get parts results
  stk_list <- lapply(df$file_name[df$scenario == scenario], function(file_i){
    readRDS(paste0(path_res, file_i))
  })
  
  ### combine parts
  stk_tmp <- do.call(stock_ibind, stk_list)
  
  ### set name
  name(stk_tmp) <- as.character(scenario)
  
  ### save
  saveRDS(stk_tmp, file = paste0(path_res, "combined/", scenario, ".rds"))
  
}

### remove original files
file.remove(paste0(path_res, df_process$file_name))

### ------------------------------------------------------------------------ ###
### "correct" SSB and catch ####
### ------------------------------------------------------------------------ ###
### once a stock collapsed (SSB < 1), it stays collapsed

### select scenarios for full analysis
scns <- list.files(path = paste0(path_res, "combined"), 
                   pattern = "[0-9]+\\.rds")
scns <- gsub(x = scns, pattern = "\\.rds", replacement = "")
scns <- sort(as.numeric(scns))

### subset to new ones...
scns <- scns[scns > 7906]

dir.create(paste0(path_res, "combined/corrected"))

### load refpts
refpts <- readRDS("input/refpts_wklife8.rds")

### loop through all scenarios
res <- foreach(scenario = scns,
               .packages = "FLCore",
               .export = c("path_res", "refpts", "scn_df")) %dopar% {
  #browser()
  ### scenario definitions
  OM_scn_i <- gsub(x = scn_df$OM_scn[scn_df$scenario == scenario], 
                   pattern = "/", replacement = "")
  stock <- scn_df$stock[scn_df$scenario == scenario]
                 
  ### load results
  stk_tmp <- readRDS(paste0(path_res, "combined/", scenario, ".rds"))
  
  ### get reference points
  refpts_i <- refpts[[OM_scn_i]][[stock]]
  
  ### find first F=5 for each iteration
  ### return numeric year or NA if F not maxed out
  fmax <- apply(fbar(stk_tmp), 6, FUN = function(x){
    #browser()
    as.numeric(dimnames(x)$year[min(which(x >= 4.999999999))])
  })
  
  ### extract SSB and catch and Fbar
  SSB_old <- SSB <- ssb(stk_tmp)
  catch_old <- catch <- catch(stk_tmp)
  fbar_old <- fbar <- fbar(stk_tmp)
  ### assume collapse/extinction if SSB drops below 1 (0.1% of B0)
  for (i in which(is.finite(fmax))) {
    ### years to check
    years <- NA
    years <- seq(from = min(fmax[,,,,, i] + 1, range(stk_tmp)[["maxyear"]]),
                 to = dims(SSB)$maxyear)
    ### find first year with SSB < 1
    min_year <- years[min(which(SSB[, ac(years),,,, i] < 1))]
    ### years to correct
    if (length(min_year) > 0) {
      if (!is.na(min_year)) {
        years_chg <- years[years >= min_year]
       
        ### change SSB
        SSB[, ac(years_chg),,,, i] <- 0
        ### and catch
        catch[, ac(years_chg),,,, i] <- 0
        ### and Fbar
        fbar[, ac(years_chg),,,, i] <- 0
      }
    }
  }
  
  ### calculate values relative to MSY
  ssb_rel <- SSB / refpts_i["msy", "ssb"]
  fbar_rel <- fbar / refpts_i["msy", "harvest"]
  
  ### save SSB and catch, both versions
  saveRDS(object = list(ssb = SSB, catch = catch, ssb_old = SSB_old,
                        catch_old = catch_old, fbar = fbar, 
                        fbar_old = fbar_old, fbar_rel = fbar_rel,
                        ssb_rel = ssb_rel),
          file = paste0(path_res, "combined/corrected/", scenario, ".rds"))
  
}

### ------------------------------------------------------------------------ ###
### fish at Fmsy ####
### ------------------------------------------------------------------------ ###
### as reference scenario, for all OMs

### load OM scenario definitions
OM_scns <- readRDS("input/OM_scns.rds")

### get stock names corresponding to stk_pos
stk_name_pos <- unique(scn_df[, c("stk_pos", "stk_pos2", "fhist", "stock")])

### fish at Fmsy
res <- foreach(OM_scn = OM_scns$id,
               .final = function(x) {
                 names(x) <- OM_scns$id
                 return(x)}) %:% 
  foreach(stk_pos = stk_name_pos$stk_pos, 
          stk_pos2 = stk_name_pos$stk_pos2,
          stock = as.character(stk_name_pos$stock),
          .export = c("refpts", "scn_df"),
          .packages = c("FLash"), .errorhandling = "pass",
          .final = function(x) {
            names(x) <- as.character(stk_name_pos$stock)
            return(x)
          }) %dopar% {
                 # if (stk_pos > 2) stop()
                 #browser()
  ### load OM
  load(paste0("input/stocks/perfect_knowledge/", OM_scn, "/", stk_pos, 
              ".RData"))
  ### get refpts
  refpts_i <- refpts[[OM_scn]][[stock]]
  
  ### fish at FMSY
  ctrl <- fwdControl(data.frame(year = 101:dims(stk)$maxyear, 
                               quantity = "f",
                               val = c(refpts_i["msy", "harvest"])))
  stk_fwd <- fwd(stk, ctrl = ctrl, sr = sr.om, sr.residuals = sr.om.res,
                 sr.residuals.mult = TRUE, maxF = 5)
  
  return(stk_fwd)
  
}
### save
saveRDS(object = res, file = "input/stocks/perfect_knowledge/fished_msy.rds")

### get catch and ssb over projection period
catch_MSY <- lapply(res, function(x) {
  lapply(x, function(y) {
    apply(window(catch(y), start = 101), 6, sum)
  })
})
ssb_MSY <- lapply(res, function(x) {
  lapply(x, function(y) {
    apply(window(ssb(y), start = 101), 6, sum)
  })
})

### save them
saveRDS(object = catch_MSY, file = "input/all_stocks_catch_at_Fmsy_wklife8.rds")
saveRDS(object = ssb_MSY, file = "input/all_stocks_ssb_at_Fmsy_wklife8.rds")

### ------------------------------------------------------------------------ ###
### stats ####
### ------------------------------------------------------------------------ ###

### load values achieved when stocks fished at MSY
catch_msy <- readRDS("input/all_stocks_catch_at_Fmsy_wklife8.rds")
ssb_msy <- readRDS("input/all_stocks_ssb_at_Fmsy_wklife8.rds")
### reference points
refpts <- readRDS("input/refpts_wklife8.rds")
B_lim <- 162.79

### go through requested scenarios
stats <- foreach(scenario = scns,
                 stk_pos = scn_df$stk_pos[scn_df$scenario %in% scns],
                 stk_pos2 = scn_df$stk_pos2[scn_df$scenario %in% scns],
                 stock = as.character(scn_df$stock[scn_df$scenario %in% scns]),
                 OM_scn = gsub(x = scn_df$OM_scn[scn_df$scenario %in% scns],
                               pattern = "/", replacement = ""),
                 .packages = "FLCore",
                 .export = c("B_lim", "catch_msy", "ssb_msy", 
                             "refpts")) %dopar% {
                               
  ### load quants
  quants <- readRDS(paste0(path_res, "combined/corrected/", scenario, ".rds"))
 
  ### subset to simulation period
  quants <- window(FLQuants(quants), start = 101)
  
  ### get catch and SSB when fished at MSY
  catch_msy_i <- catch_msy[[OM_scn]][[stk_pos]]
  ssb_msy_i <- ssb_msy[[OM_scn]][[stk_pos]]
  ### reference points
  refpts_i <- refpts[[OM_scn]][[stock]]
  
  ### list for storing results
  res <- list()
  ### number of iterations for current scenario
  n_iter <- dim(quants[[1]])[6]
  n_years <- dim(quants[[1]])[2]
  
  ### proportion where SSB was below B_lim
  res$risk_blim <- 
    sum(quants$ssb < B_lim) / (n_iter*n_years)
  ### proportion of iterations where SSB dropped below B_lim
  res$risk_blim_iter <-
    sum(yearSums(quants$ssb < B_lim) > 0) / n_iter
  
  ### stock collapse = ssb < 1
  res$risk_collapse <- sum(quants$ssb < 1) / (n_iter*n_years)
  ### proportion of iterations with collapse
  res$risk_collapse_iter <-
    sum(yearSums(quants$ssb < 1) > 0) / n_iter
  
  ### risk SSB < half Bmsy
  res$risk_0.5Bmsy <- 
    sum(quants$ssb < c(refpts_i["msy", "ssb"]/2)) / (n_iter*n_years)
  
  ### risk SSB < Bmsy
  res$risk_Bmsy <- 
    sum(quants$ssb < c(refpts_i["msy", "ssb"])) / (n_iter*n_years)
  
  ### yield and SSB relative to when fished at MSY
  res$yield_relFmsy <- median(yearSums(quants$catch) / catch_msy_i)
  res$ssb_relFmsy <- median(yearSums(quants$ssb) / ssb_msy_i)
  
  ### Fbar and SSB relative to MSY reference points
  res$ssb_rel <- median(yearMeans(quants$ssb / refpts_i["msy", "ssb"]))
  res$f_rel <- median(yearMeans(quants$fbar / refpts_i["msy", "harvest"]))
 
  ### catch iav
  res$catch_iav <- c(iav(object = quants$catch, period = 2, start = 101,
                         summary_per_iter = mean, summary = median))
 
  return(res)
 
}

stats <- do.call(rbind, lapply(stats, data.frame))
stats$scenario <- scns

### save
saveRDS(stats, file = "output/stats_corrected_wklife8.RDS")

### combine results with scenario definitions
source("MP_scenarios.R")

### merge
stats <- merge(stats, scn_df, all = TRUE)
### sort
stats <- stats[order(stats$scenario), ]

### save
saveRDS(stats, file = "output/stats_corrected_scn_wklife8.RDS")
stats <- readRDS("output/stats_corrected_scn_wklife8.RDS")



### ------------------------------------------------------------------------ ###
### extract refpts from OMs ####
### ------------------------------------------------------------------------ ###

# ### FLBRP refpts already exist
# refpts <- readRDS("input/refpts.rds")
# 
# ### get stock names
# stk_names <- names(refpts)
# ### stock positions
# stk_pos <- seq_along(stk_names)
# stk_pos[stk_pos > 15] <- stk_pos[stk_pos > 15] + 15
# 
# ### extract & save refpts used in simulation
# refpts_OM_error <- lapply(stk_pos, function(x){
#   load(paste0("input/stocks/observation_error/", x, ".RData"))
#   return(attr(stk, "refpts")[,,,,, 1])
# })
# names(refpts_OM_error) <- stk_names
# saveRDS(refpts_OM_error, file = "input/refpts_OM_error.rds")
# ### same without error model
# refpts_OM_det <- lapply(stk_pos, function(x){
#   load(paste0("input/stocks/perfect_knowledge//", x, ".RData"))
#   return(attr(stk, "refpts")[,,,,, 1])
# })
# names(refpts_OM_det) <- stk_names
# saveRDS(refpts_OM_det, file = "input/refpts_OM_det.rds")

### ------------------------------------------------------------------------ ###
### load results from all available (combined) scenarios ####
### ------------------------------------------------------------------------ ###
### get available files
# files <- list.files("output/perfect_knowledge/combined/")
# ### get scenario numbers
# scenarios <- sort(unlist(lapply(files, function(x){
#   as.numeric(unlist(strsplit(x, split = "\\."))[1])
# })))
# ### read results
# res <- foreach(scenario = scenarios) %dopar% {
#   readRDS(paste0("output/perfect_knowledge/combined/", scenario, ".rds"))
# }
# names(res) <- scenarios

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
# . <- lapply(scenarios, function(x){
#   stk_tmp <- readRDS(paste0("output/perfect_knowledge/combined/", x, ".rds"))
#   p <- plot(stk_tmp)
#   ggsave(filename = paste0("output/perfect_knowledge/plots/", x,
#                            ".png"),
#          width = 15, height = 15, units = "cm", dpi = 100, type = "cairo-png",
#          plot = p)
# })

### plot all scenarios from external HDD
files <- list.files("D:/WKLIFEVII/github/wklifeVII/R/output/perfect_knowledge/combined/")
scenarios <- sort(unlist(lapply(files, function(x){
  as.numeric(unlist(strsplit(x, split = "\\."))[1])
})))

### loop through all scenarios and plot
. <- foreach(scenario = scenarios, .packages = "ggplotFL") %dopar% {
  stk_tmp <- readRDS(paste0("D:/WKLIFEVII/github/wklifeVII/R/output/",
                            "perfect_knowledge/combined/", scenario, ".rds"))
  p <- plot(stk_tmp, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  ggsave(filename = paste0("output/perfect_knowledge/plots/all/", scenario,
                           ".png"),
         width = 15, height = 15, units = "cm", dpi = 100, type = "cairo-png",
         plot = p)
}

### ------------------------------------------------------------------------ ###
### calculate Blim ####
### ------------------------------------------------------------------------ ###
### definition:
### SSB where recruitment is impaired by 30%

# BevHolt <- function(a, b, ssb) {
#   return(a * ssb / (b + ssb))
# }
# 
# ### calculate Blim
# Blim <- lapply(1:30, function(x){
#   (function(){
#     ### load input objects
#     load(paste0("input/stocks/", x, ".RData"))
#     
#     ### ssb from 0 to B0
#     ssbs <- seq(0, 1000, 0.01)
#     ### recruitment with Beverton & Holt
#     rec <- BevHolt(a = c(params(sr.om)["a"]), b = c(params(sr.om)["b"]),
#                    ssb = ssbs)
#     
#     ssbs[tail(which(rec <= c(max(rec) * 0.7)), 1)]
#     
#   })()
# })

### 162.79 for all stocks

### ------------------------------------------------------------------------ ###
### statistics ####
### ------------------------------------------------------------------------ ###
B_lim <- 162.79

stats <- foreach(scenario = scenarios, .packages = "FLCore",
                 .export = "B_lim") %dopar% {
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
  
  if (mult_stk) {
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
  if (isTRUE(save)) {
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
### "correct" SSB and catch ####
### ------------------------------------------------------------------------ ###

library(FLCore)
### set up parallel computing environment
library(doParallel)
cl <- makeCluster(28)
registerDoParallel(cl)

### find available scenarios
files <- list.files("/gpfs/afmcefas/simonf/output/combined/")
### sort by scenarios
scenarios <- lapply(files, function(x) {
  tmp <- strsplit(x = x, split = "\\.")
  return(as.numeric(tmp[[1]][1]))
})
scenarios <- sort(unlist(scenarios))

### once SSB down, stay down ...
### load scenario 6803
# stk <- readRDS("D:/WKLIFEVII/github/wklifeVII/R/output/perfect_knowledge/combined/6803.rds")

### loop through all scenarios
res <- foreach(scenario = scenarios, .packages = "FLCore") %dopar% {
  
  ### load results
  stk_tmp <- readRDS(paste0("/gpfs/afmcefas/simonf/output/combined/",
                            scenario, ".rds"))
  
  ### find first F=5 for each iteration
  ### return numeric year or NA if F not maxed
  fmax <- apply(fbar(stk_tmp), 6, FUN = function(x){
    #browser()
    as.numeric(dimnames(x)$year[min(which(x >= 4.999999999))])
  })
  
  ### extract SSB and catch
  SSB_old <- SSB <- ssb(stk_tmp)
  catch_old <- catch <- catch(stk_tmp)
  ### assume collapse/extinction if SSB drops below 1 (0.1% of B0)
  for(i in which(is.finite(fmax))) {
    ### years to check
    years <- NA
    years <- seq(from = min(fmax[,,,,, i]+1, 200),
                 to = dims(SSB)$maxyear)
    ### find first year with SSB < 1
    min_year <- years[min(which(SSB[, ac(years),,,, i] < 1))]
    ### years to correct
    if (length(min_year) > 0) {
      if (!is.na(min_year)) {
        years_chg <- years[years >= min_year]
        
        ### change SSB
        SSB[, ac(years_chg),,,, i] <- 0
        ### and catch
        catch[, ac(years_chg),,,, i] <- 0
      }
    }
  }
  
  ### return SSB and catch
  ### both versions
  return(list(ssb = SSB, catch = catch, ssb_old = SSB_old,
              catch_old = catch_old))
  
}
### name scenarios
names(res) <- res_df[!is.na(res_df$ssb_below_blim_total), "scenario"]

### save
saveRDS(object = res, "output/quants_1_6804.rds")

stopCluster(cl)

### ------------------------------------------------------------------------ ###
### calculate statistics with new quants ####
### ------------------------------------------------------------------------ ###

cl <- makeCluster(28)
registerDoParallel(cl)

### define Blim
B_lim <- 162.79

### load quants
res <- readRDS("output/quants_1_6804.rds")

#res <- res[5001:length(res)]

stats <- foreach(scenario = res, .packages = "FLCore",
                 .export = "B_lim") %dopar% {
                   
  ### list for storing results
  res_temp <- list()
  ### number of iterations for current scenario
  n_iter <- dim(scenario[[1]])[6]
  
  ### proportion where SSB was below B_lim
  res_temp$ssb_below_blim_total_old <-
   sum(scenario$ssb_old[, ac(101:200)] < B_lim) / (n_iter*100)
  ### proportion of iterations where SSB dropped below B_lim
  res_temp$ssb_below_blim_iter_old <-
   sum(yearSums(scenario$ssb_old[, ac(101:200)] < B_lim) > 0) / n_iter
  
  ### stock collapse = ssb < 1
  res_temp$collapse_total_old <-
   sum(scenario$ssb_old[, ac(100:200)] < 1) / (n_iter*100)
  ### proportion of iterations with collapse
  res_temp$collapse_iter_old <-
   sum(yearSums(scenario$ssb_old[, ac(101:200)] < 1) > 0) / n_iter
  
  # ### how frequently is max F reached?
  # res_temp$fmaxed_total <- sum(fbar(stk_tmp)[, ac(101:200)] == 5) / (n_iter*100)
  # ### in how many iterations did this happen?
  # res_temp$fmaxed_iter <- sum(yearSums(fbar(stk_tmp)[, ac(101:200)] == 5) > 0)/
  #   n_iter
  
  ### yield
  res_temp$yield_old <- mean(scenario$catch_old[,ac(101:200)])
  res_temp$rel_yield_old <- mean(scenario$catch_old[,ac(101:200)]) /
   mean(scenario$catch_old[,ac(75:100)])
  
  ### same but for "corrected" quants
  ### proportion where SSB was below B_lim
  res_temp$ssb_below_blim_total <-
   sum(scenario$ssb[, ac(101:200)] < B_lim) / (n_iter*100)
  ### proportion of iterations where SSB dropped below B_lim
  res_temp$ssb_below_blim_iter <-
   sum(yearSums(scenario$ssb[, ac(101:200)] < B_lim) > 0) / n_iter
  
  ### stock collapse = ssb < 1
  res_temp$collapse_total <-
   sum(scenario$ssb[, ac(100:200)] < 1) / (n_iter*100)
  ### proportion of iterations with collapse
  res_temp$collapse_iter <-
   sum(yearSums(scenario$ssb[, ac(101:200)] < 1) > 0) / n_iter
  
  ### yield
  res_temp$yield <- mean(scenario$catch[,ac(101:200)])
  res_temp$rel_yield <- median(apply(scenario$catch[, ac(101:200)], 6, mean) /
                                apply(scenario$catch[, ac(75:100)], 6, mean))
  
  ### scenario definition
  #res_temp$scenario <- names(res)[i]
  
  return(as.data.frame(res_temp))
  
}


stats <- do.call(rbind, stats)
stats$scenario <- scenarios

### save
saveRDS(stats, file = "output/stats_new.RDS")

### combine results with scenario definitions
source("MP_scenarios.R")

### merge
res_df <- merge(stats, scn_df, all = TRUE)
### sort
res_df <- res_df[order(res_df$scenario), ]

### save
saveRDS(res_df, file = "output/stats_scn_new.RDS")
res_df <- readRDS("output/stats_scn_new.RDS")

stopCluster(cl)

### ------------------------------------------------------------------------ ###
### 3.2.1 all stocks - constraints ####
### ------------------------------------------------------------------------ ###

library(reshape)

### load scenario results
df_plot <- res_df[res_df$scenario %in% 1237:4484, ]

### reshape
df_plot2 <- melt(data = df_plot,
                 id.vars = c("scenario", "stock", "lower_constraint",
                             "upper_constraint", "fhist"),
                 measure.vars = c("ssb_below_blim_total", "collapse_iter",
                                  "rel_yield"))
# levels(df_plot2$variable) <- c("p(SSB < Blim)", "iteration collapse",
#                                "relative yield")
levels(df_plot2$variable) <- c("p(SSB<B[lim])",
                               ("p(iter~collapse)"), "rel.~yield")

### modify for line plots
df_lower <- df_plot2
df_tmp <- df_lower[df_lower$lower_constraint == 0, ]
df_tmp$lower_constraint <- 0.49
df_tmp$value <- NA
df_lower$lower_constraint[df_lower$lower_constraint == 0] <- 0.48
df_lower <- rbind(df_lower, df_tmp)
df_lower$x_lower <- "limit"
df_lower$x_lower[df_lower$lower_constraint == 0.48] <- "no limit"

df_upper <- df_plot2
df_tmp <- df_upper[df_upper$upper_constraint == Inf, ]
df_tmp$upper_constraint <- 1.51
df_tmp$value <- NA
df_upper$upper_constraint[df_upper$upper_constraint == Inf] <- 1.52
df_upper <- rbind(df_upper, df_tmp)
df_upper$x_upper <- "limit"
df_upper$x_upper[df_upper$upper_constraint == 1.52] <- "no limit"

### loop through stocks
. <- foreach(stock = unique(df_plot2$stock), .packages = "ggplot2",
             .export = c("df_lower", "df_upper")) %dopar% {
               
               ### line plots - lower limit
               p <- ggplot(data = df_lower[df_lower$stock == stock, ],
                           aes(x = lower_constraint, y = value,
                               colour = as.factor(upper_constraint),
                               shape = x_lower)) +
                 geom_line() + geom_point() +
                 labs(x = "lower limit") +
                 facet_grid(variable ~ fhist, scales = "free", labeller = "label_parsed") +
                 theme_bw() +
                 scale_color_discrete("upper limit") +
                 scale_shape_discrete("lower limit") +
                 theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
               p
               ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/constraints/",
                                        stock, "_lower_limit.png"), plot = p,
                      width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
               
               ### line plots - upper limit
               p <- ggplot(data = df_upper[df_upper$stock == stock, ],
                           aes(x = upper_constraint, y = value,
                               colour = as.factor(lower_constraint),
                               shape = x_upper)) +
                 geom_line() + geom_point() +
                 labs(x = "upper limit") +
                 facet_grid(variable ~ fhist, scales = "free", labeller = "label_parsed") +
                 theme_bw() +
                 scale_color_discrete("lower limit", guide = guide_legend(order = 1)) +
                 scale_shape_discrete("upper limit", guide = guide_legend(order = 2)) +
                 theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
               p
               ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/constraints/",
                                        stock, "_upper_limit.png"), plot = p,
                      width = 15, height = 10, units = "cm", dpi = 300, type = "cairo-png")
               
             }

### line plots - all stocks, lower limit only
p <- ggplot(data = df_lower[df_lower$upper_constraint == Inf, ],
            aes(x = lower_constraint, y = value,
                colour = as.factor(stock),
                shape = x_lower)) +
  geom_line() + geom_point() +
  labs(x = "lower limit") +
  facet_grid(variable ~ fhist, scales = "free", labeller = "label_parsed") +
  theme_bw() +
  scale_color_discrete("stock") +
  scale_shape_discrete("lower limit") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/constraints/",
                         "all_lower_limit.png"), plot = p,
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo-png")

### line plots - all stocks, upper limit only
p <- ggplot(data = df_upper[df_upper$lower_constraint == 0, ],
            aes(x = upper_constraint, y = value,
                colour = as.factor(stock),
                shape = x_upper)) +
  geom_line() + geom_point() +
  labs(x = "upper limit") +
  facet_grid(variable ~ fhist, scales = "free", labeller = "label_parsed") +
  theme_bw() +
  scale_color_discrete("stock") +
  scale_shape_discrete("upper limit") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/constraints/",
                         "all_upper_limit.png"), plot = p,
       width = 30, height = 20, units = "cm", dpi = 300, type = "cairo-png")


### look at values...
upper <- spread(data = df_plot2[df_plot2$lower_constraint == 0, -1],
                upper_constraint, value, drop = TRUE)
upper$check <- with(upper, `Inf` / `1.1`)
upper[, c("stock", "fhist", "variable", "check")]

# blubb1 <- df_plot2[df_plot2$upper_constraint %in% c(Inf, 1.5) &
#            df_plot2$lower_constraint == 0 & df_plot2$stock == "her-nis", ]
# spread(blubb1, upper_constraint, value, -scenario)
### ------------------------------------------------------------------------ ###
### 3.2.1 all stocks - multiplier & exponent ####
### ------------------------------------------------------------------------ ###

library(reshape)

### load scenario results
df_plot <- res_df[res_df$scenario %in% 4485:6804, ]

### reshape
df_plot2 <- gather(data = df_plot, key = "variable", value = "value",
                   ssb_below_blim_total, collapse_iter, rel_yield,
                   factor_key = TRUE)
levels(df_plot2$variable) <- c("p(SSB<B[lim])",
                               ("p(iter~collapse)"), "rel.~yield")

### values relative to value at HCRmult 1
df_merge <- df_plot2[df_plot2$HCRmult == 1,
                     c("stock", "fhist", "variable", "value", "b_z")]
names(df_merge)[names(df_merge) == "value"] <- "value_1"
### merge
df_merged <- merge(df_plot2, df_merge)
### calculate relative values
df_merged$value_rel <- with(df_merged, value / value_1)
### Infs
df_merged[df_merged$value_rel == Inf & df_merged$variable == "p(iter~collapse)" &
            !is.na(df_merged$value_rel), ]
### rel values at multiplier 0.5
risk0.5 <- subset(df_merged, variable == "p(iter~collapse)" & HCRmult == 0.5 & b_z == 1, c("stock", "fhist", "value_rel"))
summary(risk0.5)
summary(risk0.5[risk0.5$fhist == "one-way", ])
summary(risk0.5[risk0.5$fhist == "roller-coaster", ])

### calculate values relative to max
df_plot3 <- df_plot2 %>% group_by(stock, fhist, variable, b_z) %>%
  summarise(value_max = max(value, na.rm = TRUE)) %>%
  left_join(df_plot2) %>%
  mutate(value_rel = value / value_max)



### loop through stocks
. <- foreach(stock = unique(df_plot2$stock), .packages = "ggplot2") %dopar% {
  browser()
  ### line plots - multiplier
  p <- ggplot(data = df_plot2[df_plot2$stock == stock, ],
              aes(x = HCRmult, y = value,
                  colour = as.factor(b_z))) +
    geom_line() + geom_point() +
    labs(x = "multiplier") +
    facet_grid(variable ~ fhist, scales = "free", labeller = "label_parsed") +
    theme_bw() +
    scale_color_discrete("exponent b") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  p
  ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                           "multiplier_exponent/", stock, "_multiplier.png"),
         plot = p, width = 15, height = 10, units = "cm", dpi = 300,
         type = "cairo-png")
  
  ### line plots - exponent
  p <- ggplot(data = df_plot2[df_plot2$stock == stock, ],
              aes(x = b_z, y = value,
                  colour = as.factor(HCRmult))) +
    geom_line() + geom_point() +
    labs(x = "exponent") +
    facet_grid(variable ~ fhist, scales = "free", labeller = "label_parsed") +
    theme_bw() +
    scale_color_discrete("multiplier") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  p
  ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                           "multiplier_exponent/", stock, "_exponent.png"),
         plot = p, width = 15, height = 10, units = "cm", dpi = 300,
         type = "cairo-png")
  
  ### grid
  p <- ggplot(data = df_plot2[df_plot2$stock == stock, ],
              aes(y = HCRmult, x = b_z)) +
    scale_fill_gradient("value", limits = c(0, 1),
                        low = "green", high = "red") +
    labs(y = "advice multiplier", x = "exponent b") +
    facet_grid(variable ~ fhist, labeller = "label_parsed") +
    geom_raster(aes(fill = value)) +
    geom_text(aes(label = round(value, 2)), size = 2) +
    theme_bw()
  p
  ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                           "multiplier_exponent/", stock, "_grid.png"),
         plot = p, width = 15, height = 10, units = "cm", dpi = 300,
         type = "cairo-png")
  
}

### line plots - all stocks - multiplier
p <- ggplot(data = df_plot2[df_plot2$b_z == 1, ],
            aes(x = HCRmult, y = value,
                colour = as.factor(stock))) +
  geom_line() + geom_point() +
  labs(x = "multiplier") +
  facet_grid(variable ~ fhist, scales = "free", labeller = "label_parsed") +
  theme_bw() +
  scale_color_discrete("stock") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "multiplier_exponent/all_multiplier.png"),
       plot = p, width = 30, height = 20, units = "cm", dpi = 300,
       type = "cairo-png")

### line plots - all stocks - exponent
p <- ggplot(data = df_plot2[df_plot2$HCRmult == 1, ],
            aes(x = b_z, y = value,
                colour = as.factor(stock))) +
  geom_line() + geom_point() +
  labs(x = "exponent z") +
  facet_grid(variable ~ fhist, scales = "free", labeller = "label_parsed") +
  theme_bw() +
  scale_color_discrete("stock") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "multiplier_exponent/all_exponent.png"),
       plot = p, width = 30, height = 20, units = "cm", dpi = 300,
       type = "cairo-png")

### line plots - all stocks - multiplier - relative
p <- ggplot(data = df_plot3[df_plot3$b_z == 1, ],
            aes(x = HCRmult, y = value_rel,
                colour = as.factor(stock))) +
  geom_line() + geom_point() +
  labs(x = "multiplier") +
  facet_grid(variable ~ fhist, scales = "free", labeller = "label_parsed") +
  theme_bw() +
  scale_color_discrete("stock") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(y = "relative value")
p
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "multiplier_exponent/all_multiplier_rel.png"),
       plot = p, width = 30, height = 20, units = "cm", dpi = 300,
       type = "cairo-png")
p2 <- ggplot(data = blubb[blubb$b_z == 1, ],
             aes(x = HCRmult, y = value_rel,
                 colour = as.factor(stock))) +
  geom_line() + geom_point() +
  labs(x = "multiplier") +
  facet_grid(variable ~ fhist, scales = "free", labeller = "label_parsed") +
  theme_bw() +
  scale_color_discrete("stock") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(y = "relative value")
p2




### look at values: multiplier ...
HCRmult <- spread(data = df_plot2[df_plot2$b_z == 1, -1:-10],
                  HCRmult, value, drop = TRUE)
upper$check <- with(upper, `Inf` / `1.1`)
upper[, c("stock", "fhist", "variable", "check")]

blubb1 <- df_plot2[df_plot2$upper_constraint %in% c(Inf, 1.5) &
                     df_plot2$lower_constraint == 0 & df_plot2$stock == "her-nis", ]
spread(blubb1, upper_constraint, value, -scenario)



### plot 3.2.1 default
4485:6804
res[1:4484] <- NA

res_df2 <- res_df[4485:6804, ]
res_df2 <- res_df2[res_df2$HCRmult == 1 & res_df2$b_z == 1, ]

quant_lst <- res[res_df2$scenario]
quant_lst <- lapply(quant_lst, function(x){
  x[[1]]
})
ssbs <- FLQuants(quant_lst)
names(ssbs) <- res_df2$stock
plot(ssbs[c(1:22, 24:29)])



### ------------------------------------------------------------------------ ###
### 3.2.1 correlate with life-history ####
### ------------------------------------------------------------------------ ###
library(tidyverse)

### load stats
res_df <- readRDS("output/stats_scn_new.RDS")

### subset to default settings
res_df2 <- res_df[res_df$scenario %in% 4485:6804, ]
res_df2 <- res_df2[res_df2$HCRmult == 1 & res_df2$b_z == 1, ]
### 58 lines, 29 stocks & 2 fishing histories

### load lhist
lhist <- read.csv("input/stock_list_full2.csv", header = TRUE)
# lhist <- lhist[, c("stock", "a", "b", "lmax", "linf", "l50", "a50", "t0",
#                    "k", "M", "M_mat", "MK")]
#lhist <- read.csv("output/perfect_knowledge/plots/3.2.1/subgroup/lh.csv")

### load more life-history parameters
### find stock positions
stocks <- unique(res_df[res_df$fhist == "one-way", ][, c("stock", "stk_pos")])
### loop through all stocks
lhist_lst <- lapply(stocks$stk_pos, function(stock) {
  ### load data
  load(paste0("input/stocks/perfect_knowledge/", stock, ".RData"))
  ### extract lhist
  lhist_tmp <- attr(stk, "lhpar")[, 1]
  ### create data frame with values
  tmp <- c(lhist_tmp)
  names(tmp) <- dimnames(lhist_tmp)$params
  if (stock == 12) {tmp <- c(tmp[1:16], l50=NA, tmp[17:18])}
  return(tmp)
})
lhist_lst <- as.data.frame(do.call(rbind, lhist_lst))
lhist_lst$stock <- stocks$stock

lhist_full <- merge(lhist[, c("stock", "M_mat", "MK")],
                    lhist_lst[, c("L_inf", "K", "t0", "a", "b", "a50", "l50",
                                  "M", "max_age", "stock")], by.x = "stock",
                    by.y = "stock")

### merge with stats
res_df2 <- merge(res_df2, lhist_full)

### reshape
res_df3 <- gather(data = res_df2, key = variable, value = value,
                  c("ssb_below_blim_total", "collapse_iter", "rel_yield"))
### sort variables
res_df3$variable <- as.factor(res_df3$variable)
res_df3$variable <- factor(res_df3$variable, levels(res_df3$variable)[c(3,1,2)])

### plot risks
ggplot(data = res_df3, aes(x = stock, y = value, fill = stock)) +
  theme_bw() +
  facet_grid(variable ~ fhist, scales = "free_y") +
  geom_bar(stat = "identity", show.legend = FALSE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25))

### plot correlations
ggplot(data = res_df3, aes(x = max_age, y = value, colour = stock)) +
  theme_bw() +
  facet_grid(variable ~ fhist, scales = "free_y") +
  #geom_bar(stat = "identity", show.legend = FALSE) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25))

### plot risk ~ risk
ggplot(data = res_df2, aes(x = ssb_below_blim_total, y = collapse_iter,
                           colour = stock)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ fhist) +
  stat_smooth(method = "lm", colour = "black", se = FALSE, fill = "grey90") +
  geom_abline(intercept = 0, slope = 1)

### plot risk ~ yield
ggplot(data = res_df2, aes(x = ssb_below_blim_total, y = rel_yield,
                           colour = stock)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ fhist) +
  stat_smooth(method = "lm", colour = "black", se = TRUE, fill = "grey90")


### ------------------------------------------------------------------------ ###
### 3.2.1 try grouping ####
### ------------------------------------------------------------------------ ###

### load quants
# res <- readRDS("output/quants_1_6804.rds")

### stats
res_df <- readRDS("output/stats_scn_new.RDS")

### find default scenarios
scenarios <- res_df[res_df$scenario %in% 4485:6804 &
                      res_df$HCRmult == 1 & res_df$b_z == 1, "scenario"]

### extract and save default 3.2.1 runs for all stocks/fhist
# quants <- res[ac(scenarios)]
# saveRDS(quants, file = "output/quants_3.2.1_default.rds")
quants <- readRDS("output/quants_3.2.1_default.rds")

### coerce into data frames, extract median for each scenario
quants_df <- lapply(scenarios, function(x) {
  df_tmp <- lapply(names(quants[[ac(x)]]), function(y){
    cbind(as.data.frame(apply(quants[[ac(x)]][[y]], 2, median)), quant = y)
  })
  df_tmp <- do.call(rbind, df_tmp)
  df_tmp <- cbind(df_tmp, res_df[res_df$scenario == x, ])
  df_tmp
})
quants_df <- do.call(rbind, quants_df)

ggplot(data = quants_df, aes(x = year, y = data, colour = stock)) +
  facet_grid(quant ~ fhist) +
  geom_line() +
  theme_bw()

### plot for all stocks
for (i in 1:29) {
  p <- ggplot(data = quants_df[quants_df$stk_pos2 == i, ],
              aes(x = year, y = data, colour = stock)) +
    facet_grid(quant ~ fhist) +
    geom_line() +
    theme_bw()
  print(p)
  a <- readline()
}

### ------------------------------------------------------------------------ ###
### compare MSY length reference points
### ------------------------------------------------------------------------ ###
cl <- makeCluster(4)
registerDoParallel(cl)

### unique stocks
stocks <- unique(res_df[, c("stk_pos", "stock")])

### load functions for calculating Lc
source("functions/fFun.R")

### loop through all stocks
pars <- foreach(stock = stocks$stk_pos, .packages = "FLCore") %dopar% {
  
  ### load results
  load(paste0("input/stocks/perfect_knowledge/", stock, ".RData"))
  ### calculate Lc
  L_c <- median(calc_Lc(attr(stk, "catch_len")[,ac(100)]))
  ### life history data
  lhpar <- apply(stk@lhpar, 1, median)
  ### add l50, if missing
  if (!"l50" %in% dimnames(lhpar)$params) {
    lhpar <- rbind2(lhpar[1:16], FLPar(l50 = NA), lhpar[17:18])
  }
  ### reference points
  refpts <- apply(stk@refpts, 1, median)
  ### calculate mature M
  M_mat <- weighted.mean(x = c(m(stk)[,1,,,, 1]), c(mat(stk)[,1,,,, 1]))
  ### M/K
  MK <- M_mat / c(lhpar["K"])
  ### bind everything together
  tmp <- rbind2(lhpar, refpts, FLPar(L_c = L_c, M_mat = M_mat, MK = MK))
  ### return as data frame
  tmp <- t(data.frame(tmp))
  row.names(tmp) <- NULL
  return(tmp)
  
}
pars <- as.data.frame(do.call(rbind, pars))
pars$stk_pos <- 1:nrow(pars)

### merge with stock details
pars <- merge(pars,
              unique(res_df[, c("stk_pos", "stk_pos2", "stock", "fhist")]))
### keep only one row per stock
pars <- pars[match(unique(pars$stock), pars$stock), ]


### calculate LFeM with M/K = 1.5
pars$LFeMK <- with(pars, (L_inf + 2 * 1.5 * L_c) / (1 + 2 * 1.5))

### save
saveRDS(object = pars, "input/all_stocks_repfts_lhpar.rds")

stopCluster(cl)

### plotting
plot(pars$LFeMK ~ pars$LFeFmsy)

pars_df <- gather(data = pars, key = "method", value = "length",
                  c("LFeFmsy", "LFeMK", "LFeM"))
ggplot(data = pars_df, aes(x = stock, y = length, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "correlations/length_ref_comparison.png"),
       width = 20, height = 15, units = "cm", dpi = 300,
       type = "cairo-png")

### correlate with risks
### subset to default settings
res_df2 <- res_df[res_df$scenario %in% 4485:6804 &
                    res_df$HCRmult == 1 & res_df$b_z == 1, ]
### merge with reference values
res_df3 <- merge(res_df2, pars, by = "stock", all = TRUE)
### absolute deviation
ggplot(data = res_df3,
       aes(x = ssb_below_blim_total,
           y = LFeFmsy - LFeMK)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ fhist.x) +
  theme_bw() +
  labs(x = "p(SSB<Blim)")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "correlations/risk_dev_length_ref_absolute.png"),
       width = 20, height = 15, units = "cm", dpi = 300,
       type = "cairo-png")
### relative deviation
ggplot(data = res_df3,
       aes(x = ssb_below_blim_total,
           y = (LFeFmsy - LFeMK) / LFeFmsy)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ fhist.x) +
  theme_bw() +
  labs(x = "p(SSB<Blim)")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "correlations/risk_dev_length_ref_rel.png"),
       width = 20, height = 15, units = "cm", dpi = 300,
       type = "cairo-png")
### absolute values of relative deviation
ggplot(data = res_df3,
       aes(x = ssb_below_blim_total,
           y = abs((LFeFmsy - LFeMK) / LFeFmsy))) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ fhist.x) +
  theme_bw() +
  labs(x = "p(SSB<Blim)")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "correlations/risk_dev_length_ref_absrel.png"),
       width = 20, height = 15, units = "cm", dpi = 300,
       type = "cairo-png")

### ------------------------------------------------------------------------ ###
### 3.2.1 correlate stats with lhpars
### ------------------------------------------------------------------------ ###

### select required columns
df_cor <- res_df3[, c(1, 9, 12, 14, 17, 30:32, 46:48, 51:54, 57:59, 62)]
### reshape
df_cor <- gather(df_cor, key = "par", value = "par_value",
                 L_inf:LFeMK)
df_cor <- gather(df_cor, key = "stat", value = "stat_value",
                 ssb_below_blim_total:rel_yield)
### plot correlations
ggplot(data = df_cor, aes(x = par_value, y = stat_value)) +
  geom_point() +
  facet_grid(stat + fhist.x ~ par, scales = "free") +
  theme_bw() +
  geom_smooth(method = "lm") +
  ylim(0, 1) +
  labs(x = "parameter value", y = "stats value")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "correlations/correlations.png"),
       width = 29.7, height = 21, units = "cm", dpi = 200,
       type = "cairo-png")


### ------------------------------------------------------------------------ ###
### 3.2.1 default - try to understand behaviour ####
### ------------------------------------------------------------------------ ###

### load stocks and stats
res_df <- readRDS("output/stats_scn_new.RDS") ### quants
quants <- readRDS("output/quants_3.2.1_default.rds") ### stats

### stocks
stocks <- foreach(stock = 6515:6572) %dopar% {
  # readRDS(paste0("D:/WKLIFEVII/github/wklifeVII/R/output/perfect_knowledge",
  #                "/combined/", stock, ".rds"))
  readRDS(paste0("C:/Users/SF02/OneDrive - CEFAS/WKLIFEVII/github/wklifeVII/",
                 "R/output/perfect_knowledge/combined/", stock, ".rds"))
}
names(stocks) <- 6515:6572

stocks_desc <- res_df[6515:6572, c("scenario", "stock", "fhist")]
### Zeus faber 6551

plot(stocks[["6554"]])
plot(stocks[["6551"]])
plot(quants[["6551"]]$ssb)

### correlate SSB at beginning of simulation vs. risk
SSB_status <- foreach(stock = stocks, .packages = "FLCore") %dopar% {
  c(iterMedians(ssb(stock)[, ac(100)]))
}
res_df$SSB_status <- NA
res_df$SSB_status[res_df$scenario %in% 6515:6572] <- unlist(SSB_status)

ggplot(data = res_df, aes(x = SSB_status, y = collapse_iter)) +
  geom_point() +
  facet_wrap(~ fhist)
plot(head(unlist(SSB_status), 15) ~
       head(res_df$ssb_below_blim_total[res_df$scenario %in% 6515:6572], 15))
### no correlation whatsoever...


stock <- stocks[["6551"]]
plot(stocks[["6551"]]@tracking[c("HCR3.2.1r", "HCR3.2.1f", "HCR3.2.1b"),,,,, 1]) +
  geom_hline(yintercept = 1)# + ylim(0, 2)
plot(stocks[["6516"]]@tracking[c("HCR3.2.1r", "HCR3.2.1f", "HCR3.2.1b"),,,,, 1]) + geom_hline(yintercept = 1)

HCR3.2.1r <- foreach(stock = stocks, .packages = "FLCore") %dopar% {
  c(iterMedians(stock@tracking["HCR3.2.1r", ac(99)]))
}
HCR3.2.1f <- foreach(stock = stocks, .packages = "FLCore") %dopar% {
  c(iterMedians(stock@tracking["HCR3.2.1f", ac(99)]))
}
HCR3.2.1b <- foreach(stock = stocks, .packages = "FLCore") %dopar% {
  c(iterMedians(stock@tracking["HCR3.2.1b", ac(99)]))
}
HCR3.2.1fvar <- foreach(stock = stocks, .packages = "FLCore") %dopar% {
  #iterMedians(apply(stock@tracking["HCR3.2.1f"], 6, var, na.rm = TRUE))
  iterMedians(apply(abs(stock@tracking["HCR3.2.1f"] - 1), 6, median, na.rm = TRUE))
}

### create df with stats ...
df_default <- res_df[res_df$scenario %in% 6515:6572, ]
df_default$SSB_stats <- unlist(SSB_status)
df_default$HCR3.2.1r <- unlist(HCR3.2.1r)
df_default$HCR3.2.1f <- unlist(HCR3.2.1f)
df_default$HCR3.2.1b <- unlist(HCR3.2.1b)
df_default$HCR3.2.1fvar <- unlist(HCR3.2.1fvar)

df_def_plot <- gather(df_default, key = "stats", value = "value",
                      SSB_stats:HCR3.2.1fvar)
ggplot(df_def_plot, aes(x = value, y = ssb_below_blim_total)) +
  geom_point() +
  theme_bw() +
  facet_grid(fhist ~ stats, scale = "free")

stock <- 6551
plot_trial <- function(stock, yrs, iters){
  track <- FLQuants(HCR_r = stock@tracking[c("HCR3.2.1r")],
                    HCR_f = stock@tracking[c("HCR3.2.1f")],
                    HCR_b = stock@tracking[c("HCR3.2.1b")])
  track <- lapply(track, function(x) {
    names(dimnames(x))[1] <- "age"
    dimnames(x)[1] <- "all"
    x
  })
  qts <- FLQuants(Rec = rec(stock),
                  SSB = ssb(stock),
                  catch = catch(stock),
                  fbar = fbar(stock),
                  HCR_r = track$HCR_r,
                  HCR_f = track$HCR_f,
                  HCR_b = track$HCR_b)
  qts_plot <- FLCore::iter(window(qts, start = head(yrs, 1),
                                  end = tail(yrs, 1)),
                           iters)
  plot(qts_plot) + theme_bw() +
    geom_hline(data = data.frame(yintercept = 1,
                                 qname = c("HCR_r", "HCR_f", "HCR_b")),
               aes(yintercept = yintercept),
               linetype = "dashed", colour = "grey")
}
plot_trial(stocks[["6551"]], yrs = 99:200, iters = 1)

### loop through stocks & quants
cor_stats <- foreach(stock = stocks, quant = quants,
                     .packages = "FLCore", .combine = rbind) %dopar% {
  
  ### extract values of r, f & b from HCR
  vals <- lapply(c(r = "r", f = "f", b = "b"), function(x) {
   stock@tracking[paste0("HCR3.2.1", x)]
  })
  # fr <- stock@tracking[c("HCR3.2.1r"), ,,,, ]
  # ff <- stock@tracking[c("HCR3.2.1f"), ,,,, ]
  # fb <- stock@tracking[c("HCR3.2.1b"), ,,,, ]
  
  ### remove years after collapse
  ### (defined as SSB<1 after hitting F=5)
  ### replace with NAs
  vals <- lapply(vals, function(x){
    x@.Data[quant$ssb@.Data == 0] <- NA
    return(x)
  })
  
  ### contribution to total deviation from 1
  vals_1 <- lapply(vals, function(x) abs(x - 1)) ### difference from 1
  vals_sum <- vals_1$r + vals_1$f + vals_1$b ### sum of differences
  contr <- lapply(vals_1, function(x) x/vals_sum) ### contribution
  # cr <- abs(fr-1)/(abs(fr-1) + abs(ff-1) + abs(fb-1))
  # cf <- abs(ff-1)/(abs(fr-1) + abs(ff-1) + abs(fb-1))
  # cb <- abs(fb-1)/(abs(fr-1) + abs(ff-1) + abs(fb-1))
  
  ### remove years after collapse
  ### (defined as SSB<1 after hitting F=5)
  ### replace with NAs
  # contr <- lapply(contr, function(x){
  #  x@.Data[quant$ssb@.Data == 0] <- NA
  #  return(x)
  # })
  # cr@.Data[quant$ssb@.Data == 0] <- NA
  # cf@.Data[quant$ssb@.Data == 0] <- NA
  # cb@.Data[quant$ssb@.Data == 0] <- NA
  
  ### mean contribution per iter over all years
  contr_mean <- lapply(contr, yearMeans)
  ### median over iterations
  contr_median <- lapply(contr_mean, iterMedians)
  
  ### factor = r * f * b
  ### remove one, calculate factor, calculate sum of squares
  factor <- vals$r * vals$f * vals$b
  ### models
  factor_rf <- vals$r * vals$f
  factor_rb <- vals$r * vals$b
  factor_fb <- vals$f * vals$b
  ### sum of squares, median of that
  SS_rf <- iterMedians(apply((factor - factor_rf)^2, 6, sum, na.rm = TRUE))
  SS_rb <- iterMedians(apply((factor - factor_rb)^2, 6, sum, na.rm = TRUE))
  SS_fb <- iterMedians(apply((factor - factor_fb)^2, 6, sum, na.rm = TRUE))

  return(c(unlist(contr_median), SS_b = SS_rf, SS_f = SS_rb, SS_r = SS_fb,
           SS_b_rel = SS_rf/(SS_rf + SS_rb + SS_fb),
           SS_f_rel = SS_rb/(SS_rf + SS_rb + SS_fb),
           SS_r_rel = SS_fb/(SS_rf + SS_rb + SS_fb)))
  
}

### add results to stats
stats_HCR <- res_df[res_df$scenario %in% 6515:6572, ]
stats_HCR <- cbind(stats_HCR, cor_stats)
### add short stock names
stats_HCR <- merge(stats_HCR, read.csv("input/names_short.csv", as.is = TRUE))

### plot contribution split by component for all stocks
df_plot <- gather(data = stats_HCR, key = "component", value = "value",
                  r, f, b)
ggplot(df_plot, aes(x = short, y = value, fill = component)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ fhist) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2)) +
  labs(y = "contribution of components to catch advice")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "correlations/catch_rule_contributions.png"),
       width = 20, height = 12, units = "cm", dpi = 300, type = "cairo-png")

### plot SSs
df_plot <- gather(data = stats_HCR, key = "component", value = "value",
                  SS_r, SS_b, SS_f)
ggplot(df_plot, aes(x = short, y = value, fill = component)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ fhist) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2)) +
  labs(y = "contribution of components to catch advice")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "correlations/catch_rule_contributions.png"),
       width = 20, height = 12, units = "cm", dpi = 300, type = "cairo-png")

### correlation between components and risk
ggplot(df_plot, aes(x = value, y = ssb_below_blim_total)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(fhist ~ component) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2)) +
  labs(x = "contribution of components to catch advice", y = "p(SSB<Blim)")
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                         "correlations/catch_rule_contributions_correlations.png"),
       width = 20, height = 12, units = "cm", dpi = 300, type = "cairo-png")

### check correlation
df_cors <- group_by(df_plot, fhist, component) %>% 
  summarise(r = cor.test(x = value, y = ssb_below_blim_total)$estimate,
            p = cor.test(x = value, y = ssb_below_blim_total)$p.value)
### plot again, with correlation test resuts
ggplot(df_plot, aes(x = value, y = ssb_below_blim_total)) +
  geom_point() +
  geom_text(data = df_cors, x = 0.45, y = 0.85,
            #aes(label = paste0("r=", signif(r, 2), "\np=", signif(p, 2)))) +
            # aes(label = paste0("italic(rho)==", signif(r, 2))), parse = TRUE) +
            aes(label = paste0("atop(italic(rho)==", signif(r, 2), 
                               ",italic(p)==", signif(p, 2),")")), parse = TRUE) + 
  facet_grid(fhist ~ component) +
  theme_bw() +
  # xlab(expression("p(SSB<B [lim] )"))
  labs(x = "contribution of components to catch advice",
       y = expression(italic(p)(SSB<B[lim])))
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/correlations/",
                         "catch_rule_contributions_correlations_stats.png"),
       width = 20, height = 12, units = "cm", dpi = 300, type = "cairo-png")

### ------------------------------------------------------------------------ ###
### 3.2.1 default - recreate full corrected time series ####
### plot stock dynamics and catch rule components
### ------------------------------------------------------------------------ ###

### functions for arithmetic operations with FLPar & FLQuant
`%FLQuant+FLPar%` <- function(e1, e2) {
  res <- e1
  for (yr in dimnames(e1)$year) res[, yr] <- res[, yr] + e2
  res
}
`%FLPar+FLQuant%` <- function(e1, e2) {
  res <- e2
  for (yr in dimnames(e2)$year) res[, yr] <- res[, yr] + e1
  res
}
`%//%` <- function(e1, e2) {
  res <- e1
  for (yr in dimnames(e1)$year) res[, yr] <- res[, yr] / e2
  res
}

### load stats & refpts
res_df <- readRDS("output/stats_scn_new.RDS") ### quants
#res_df_ <- res_df[res_df$scenario %in% 6515:6572, ]
### FLBRP reference points
refpts <- readRDS("input/refpts.rds")
### extract MSY length
refpts_OM <- readRDS("input/refpts_OM_det.rds")
Lmsy <- lapply(refpts_OM, function(x) median(c(x["LFeFmsy"])))

scenarios <- 1237:4484#c(6515:6572)#res_df_$scenario

res <- foreach(scenario = scenarios, .packages = "FLCore",
               .errorhandling = "pass") %dopar% {
  #browser()
  ### load data
  # stock <- readRDS(paste0("output/perfect_knowledge/combined/", scenario,
  #                         ".rds"))
  stock <- readRDS(paste0("/gpfs/afmcefas/simonf/output/combined/", scenario,
                          ".rds"))

  ### calculate component f for historical period
  tracking <- stock@tracking
  L_FM <- tracking["L_c"]
  ### calculate f
  L_FM <- (stock@lhpar["L_inf"] %FLPar+FLQuant% (2 * 1.5 * tracking["L_c"])) /
    (1 + 2 * 1.5)
  fac_f <- tracking["L_mean"] / L_FM
  ### insert into tracking
  tracking["HCR3.2.1f", ac(75:98)] <- fac_f[, ac(75:98)]
  
  ### find first F=5 for each iteration
  ### return numeric year or NA if F not maxed
  fmax <- apply(fbar(stock), 6, FUN = function(x){
    #browser()
    as.numeric(dimnames(x)$year[min(which(x >= 4.999999999))])
  })
  
  ### extract quantities to change
  stock.n <- stock.n(stock)
  catch.n <- catch.n(stock)
  harvest <- harvest(stock)
  ssb <- ssb(stock)
  rec <- rec(stock)
  fbar <- fbar(stock)
  catch <- catch(stock)
  #tracking <- stock@tracking
  
  ### assume collapse/extinction if SSB drops below 1 (0.1% of B0)
  for (i in which(is.finite(fmax))) {
    ### years to check
    years <- NA
    years <- seq(from = min(fmax[,,,,, i] + 1, 200),
                 to = dims(ssb(stock))$maxyear)
    ### find first year with SSB < 1
    min_year <- years[min(which(ssb(stock)[, ac(years),,,, i] < 1))]
    ### years to correct
    if (length(min_year) > 0) {
      if (!is.na(min_year)) {
        years_chg <- years[years >= min_year]
        
        ### change stock.n, catch.n and fishing mortality
        # stock.n[, ac(years_chg),,,, i] <- 0
        # catch.n[, ac(years_chg),,,, i] <- 0
        ### assume NA fishing mortality, as stock and catch both zero
        # harvest[, ac(years_chg),,,, i] <- NA
        ### tracking, set NAs after stock collapse
        tracking[, ac(years_chg),,,, i] <- NA
        ### derived quants
        ssb[, ac(years_chg),,,, i] <- 0
        catch[, ac(years_chg),,,, i] <- 0
        rec[, ac(years_chg),,,, i] <- 0
        fbar[, ac(years_chg),,,, i] <- 0 ### assume zero fishing after collapse
      }
    }
  }
  #browser()
  
  # ### insert "corrected" values
  # stock.n(stock) <- stock.n
  # catch.n(stock) <- catch.n
  # harvest(stock) <- harvest
  # 
  # ### update summary slots
  # stock(stock) <- computeStock(stock)
  # catch(stock) <- computeCatch(stock)
  
  ### extract/return usefull series
  res <- FLQuants(SSB = ssb, fbar = fbar, Catch = catch, Rec = rec, 
                  HCR3.2.1r = tracking["HCR3.2.1r"],
                  HCR3.2.1f = tracking["HCR3.2.1f"],
                  HCR3.2.1b = tracking["HCR3.2.1b"],
                  Lmean = tracking["L_mean"])
  res <- lapply(res, function(x) {
    dimnames(x)[1] <- "all"
    names(dimnames(x))[1] <- "quant"
    return(x)
  })
  
  ### save 
  #saveRDS(res, file = paste0("output/perfect_knowledge/combined/3.2.1_quants/",
  #                            scenario, ".rds"))
  saveRDS(res, file = paste0("/gpfs/afmcefas/simonf/output/combined/",
                             "3.2.1_quants/", scenario, ".rds"))

}
library(tidyr)
### plot stock dynamics and catch rule components
. <- foreach(scenario = scenarios, .packages = c("FLCore", "tidyr", "ggplot2"),
             .export = c("refpts", "Lmsy", "res_df")) %dopar% {
  ### load quants
  quants <- readRDS(paste0("output/perfect_knowledge/combined/3.2.1_quants/",
                           scenario, ".rds"))
  
  ### calculate percentiles
  quants_perc <- lapply(quants, function(x) {
    quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
  })
  
  ### coerce into data frame
  quants_df <- as.data.frame(quants_perc)
  
  ### format for plotting
  quants_df <- spread(data = quants_df, iter, data)
  
  ### reference lines
  stk_name <- ac(res_df$stock[res_df$scenario == scenario])
  df_refs <- data.frame(qname = c("HCR3.2.1r", "HCR3.2.1f", "HCR3.2.1b",
                                  "SSB", "fbar", "Catch", "Rec", "Lmean"),
                        value = c(rep(1, 3),
                                  refpts[[stk_name]]["msy", "ssb"],
                                  refpts[[stk_name]]["msy", "harvest"],
                                  refpts[[stk_name]]["msy", "yield"],
                                  refpts[[stk_name]]["msy", "rec"],
                                  Lmsy[[stk_name]]))
  
  p <- ggplot(data = quants_df, aes(x = year, y = `50%`)) +
    geom_vline(xintercept = 100, colour = "grey") +
    geom_hline(data = df_refs, aes(yintercept = value), lty = "dashed") +
    geom_ribbon(aes(ymin = `5%`, ymax = `95%`), alpha = 0.2) +
    geom_ribbon(aes(ymin = `25%`, ymax = `75%`), alpha = 0.3) +
    geom_line() +
    facet_wrap(~ qname, scales = "free_y") +
    theme_bw() +
    ylim(0, NA) +
    labs(y = "")
  
  ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/all/",
                           scenario, "_", stk_name, ".png"),
         width = 20, height = 12, units = "cm", dpi = 300, type = "cairo-png",
         plot = p)
  
  ### include individual iterations
  p2 <- p + geom_line(data = as.data.frame(quants), aes(y = data, group = iter),
                      alpha = 0.1, size = 0.2, colour = "red")
  ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/all/",
                           "iterations/", scenario, "_", stk_name, ".png"),
         width = 20, height = 12, units = "cm", dpi = 300, type = "cairo-png",
         plot = p2)
  
}


### plot 4plot with medians for all stocks
df_qts <- foreach(scenario = 6515:6572, .packages = c("FLCore"), 
             .export = c("res_df")) %dopar% {
  ### load quants
  qts <- readRDS(paste0("output/perfect_knowledge/combined/3.2.1_quants/",
                        scenario, ".rds"))[1:4]
  qts <- lapply(qts, iterMedians) ### median over all iterations
  ### coerce into data frame
  df_tmp <- as.data.frame(qts)
  ### add scenario definitions
  scn <- res_df[res_df$scenario == scenario, ]
  df_tmp <- cbind(df_tmp, fhist = scn$fhist, scenario = scenario,
                  stock = scn$stock)
  return(df_tmp)
}
df_qts <- do.call(rbind, df_qts)
### plot
ggplot(data = df_qts, aes(x = year, y = data, colour = stock)) +
  geom_line(size = 0.3) + geom_vline(xintercept = 100) + 
  theme_bw() +
  facet_grid(qname ~ fhist, scales = "free_y") +
  labs(y = "") +
  scale_colour_discrete(guide = FALSE)
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/correlations",
                         "/default_4plot.png"),
       width = 20, height = 12, units = "cm", dpi = 300, type = "cairo-png")

### plot SSB for all stocks
df_ssb <- foreach(scenario = 6515:6572, .packages = c("FLCore"), 
                  .export = c("res_df")) %dopar% {
  ### load quants
  qts <- readRDS(paste0("output/perfect_knowledge/combined/3.2.1_quants/",
                        scenario, ".rds"))[["SSB"]]
  qts <- iterMedians(qts) ### median over all iterations
  ### coerce into data frame
  df_tmp <- as.data.frame(qts)
  ### add scenario definitions
  scn <- res_df[res_df$scenario == scenario, ]
  df_tmp <- cbind(df_tmp, fhist = scn$fhist, scenario = scenario,
                  stock = scn$stock)
  return(df_tmp)
}
df_ssb <- do.call(rbind, df_ssb)
### retrieve BMSY
refpts_ssb <- lapply(refpts, function(x) {
  c(x["msy", "ssb"])
})
refpts_ssb <- as.data.frame(do.call(rbind, refpts_ssb))
refpts_ssb$stock <- row.names(refpts_ssb)
row.names(refpts_ssb) <- NULL
names(refpts_ssb)[1] <- "BMSY"
### plot
ggplot(data = df_ssb[df_ssb$fhist == "one-way", ], aes(x = year, y = data)) +
  geom_vline(xintercept = 100, colour = "grey") + geom_line(size = 0.3) +
  geom_hline(data = refpts_ssb, aes(yintercept = BMSY)) + 
  theme_bw() +
  facet_wrap("stock", scales = "free_y") +
  labs(y = "") +
  scale_colour_discrete(guide = FALSE)
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/correlations",
                         "/default_ssb_all_stocks.png"),
       width = 30, height = 18, units = "cm", dpi = 300, type = "cairo-png")

### plot SSB one-way in single plot
ggplot(data = df_ssb,
       aes(x = year, y = data, group = stock)) +
  geom_vline(xintercept = 100, colour = "grey") + geom_line(size = 0.3) + 
  theme_bw() + facet_wrap(~ fhist, nrow = 2) +
  labs(y = "Spawning Stock Biomass", x = "year") +
  scale_colour_discrete(guide = FALSE)
ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/correlations",
                         "/default_ssb_all_stocks_oneplot.png"),
       width = 10, height = 12, units = "cm", dpi = 300, type = "cairo-png")

# stk_tmp <- readRDS("output/perfect_knowledge/combined/6516.rds")
# quant_tmp <- readRDS("output/perfect_knowledge/combined/3.2.1_quants/6516.rds")
# plot(quant_tmp, iter = 1:5)



### ------------------------------------------------------------------------ ###
### 3.2.1 all stocks - constraints - 4plots ####
### ------------------------------------------------------------------------ ###
### load stats & refpts
res_df <- readRDS("output/stats_scn_new.RDS") ### quants
refpts <- readRDS("input/refpts.rds")
### extract MSY length
refpts_OM <- readRDS("input/refpts_OM_det.rds")
Lmsy <- lapply(refpts_OM, function(x) median(c(x["LFeFmsy"])))

### load scenario results
scenarios <- 1237:4484
df_tmp <- res_df[res_df$scenario %in% scenarios, ]

### loop through stocks
for (stk_pos in unique(df_plot$stk_pos)) {
  
  ### find required scenarios for current stocks
  scns_tmp <- df_tmp$scenario[df_tmp$stk_pos == stk_pos]
  ### load quants for all scenarios
  qts <- lapply(scns_tmp, function(scenario) {
    ### read
    qts_tmp <- readRDS(paste0("output/perfect_knowledge/combined/3.2.1_quants/",
                              scenario, ".rds"))
    ### extract percentiles
    qts_tmp <- lapply(qts_tmp, function(x) {
      quantile(x, c(0.05, 0.5, 0.95), na.rm = TRUE)
    })
    ### coerce into data frame
    cbind(as.data.frame(qts_tmp),
          scenario = scenario, 
          fhist = res_df$fhist[res_df$scenario == scenario],
          lower_limit = res_df$lower_constraint[res_df$scenario == scenario],
          upper_limit = res_df$upper_constraint[res_df$scenario == scenario]
          )
  })
  qts <- do.call(rbind, qts)
  qts$lower_limit <- as.factor(qts$lower_limit)
  qts$upper_limit <- as.factor(qts$upper_limit)
  qts$iter <- as.factor(qts$iter)
  
  ### reference points / lines
  ### reference lines
  stk_name <- ac(res_df$stock[res_df$scenario == scns_tmp[1]])
  df_refs <- data.frame(qname = c("HCR3.2.1r", "HCR3.2.1f", "HCR3.2.1b",
                                  "SSB", "fbar", "Catch", "Rec", "Lmean"),
                        value = c(rep(1, 3),
                                  refpts[[stk_name]]["msy", "ssb"],
                                  refpts[[stk_name]]["msy", "harvest"],
                                  refpts[[stk_name]]["msy", "yield"],
                                  refpts[[stk_name]]["msy", "rec"],
                                  Lmsy[[stk_name]]))
  
  ### plot with lower limit as basis
  . <- foreach(lower_limit = unique(qts$lower_limit), 
           .packages = "ggplot2") %dopar% {
    p <- ggplot(qts[qts$lower_limit == lower_limit, ],
           aes(x = year, y = data, colour = upper_limit, 
               linetype = iter, alpha = iter)) +
      geom_vline(xintercept = 100, colour = "grey") +
      geom_hline(data = df_refs, aes(yintercept = value), lty = "dashed") +
      geom_line() +
      facet_wrap(~ qname, scales = "free_y") +
      scale_linetype_manual("percentile",
                            values = c("dashed", "solid", "dashed")) +
      scale_alpha_manual("percentile", values = c(0.5, 1, 0.5)) +
      scale_colour_discrete("upper\nlimit") +
      theme_bw() +
      ylim(0, NA) +
      labs(x = "year", y = "", title = 
           paste0(df_tmp$fhist[df_tmp$stk_pos == stk_pos][1], ", ",
                 df_tmp$stock[df_tmp$stk_pos == stk_pos][1], ", ",
                 "lower limit = ", 
                 as.character(format(as.numeric(lower_limit)), nsmall = 2)))
    ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/constraints",
                             "/stock_plots/", 
                             df_tmp$fhist[df_tmp$stk_pos == stk_pos][1], "_",
                             df_tmp$stock[df_tmp$stk_pos == stk_pos][1], "_",
                             "lower", 
                             format(as.numeric(as.character(lower_limit)), 
                                    nsmall = 2), 
                             ".png"),
           width = 30, height = 18, units = "cm", dpi = 300, type = "cairo-png",
           plot = p)
  }
  ## plot with upper limit as basis
  . <- foreach(upper_limit = unique(qts$upper_limit),
               .packages = "ggplot2") %dopar% {
    p <- ggplot(qts[qts$upper_limit == upper_limit, ],
                aes(x = year, y = data, colour = lower_limit,
                    linetype = iter, alpha = iter)) +
      geom_vline(xintercept = 100, colour = "grey") +
      geom_hline(data = df_refs, aes(yintercept = value), lty = "dashed") +
      geom_line() +
      facet_wrap(~ qname, scales = "free_y") +
      scale_linetype_manual("percentile",
                            values = c("dashed", "solid", "dashed")) +
      scale_alpha_manual("percentile", values = c(0.5, 1, 0.5)) +
      scale_colour_discrete("lower\nlimit") +
      theme_bw() +
      ylim(0, NA) +
      labs(x = "year", y = "", title =
           paste0(df_tmp$fhist[df_tmp$stk_pos == stk_pos][1], ", ",
                  df_tmp$stock[df_tmp$stk_pos == stk_pos][1], ", ",
                  "upper limit = ",
                  format(as.numeric(as.character(upper_limit)), nsmall = 2)))
    ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/constraints",
                             "/stock_plots/",
                             df_tmp$fhist[df_tmp$stk_pos == stk_pos][1], "_",
                             df_tmp$stock[df_tmp$stk_pos == stk_pos][1], "_",
                             "upper",
                             format(as.numeric(as.character(upper_limit)),
                                    nsmall = 2),
                             ".png"),
           width = 30, height = 18, units = "cm", dpi = 300, type = "cairo-png",
           plot = p)
  }
}

### ------------------------------------------------------------------------ ###
### 3.2.1 all stocks - multiplier & exponent - 4plots ####
### ------------------------------------------------------------------------ ###

### load scenario results
scenarios <- 4485:6804
df_tmp <- res_df[res_df$scenario %in% scenarios, ]

### loop through stocks
for (stk_pos in unique(df_plot$stk_pos)) {
  
  ### find required scenarios for current stocks
  scns_tmp <- df_tmp$scenario[df_tmp$stk_pos == stk_pos]
  ### load quants for all scenarios
  qts <- lapply(scns_tmp, function(scenario) {
    ### read
    qts_tmp <- readRDS(paste0("output/perfect_knowledge/combined/3.2.1_quants/",
                              scenario, ".rds"))
    ### extract percentiles
    qts_tmp <- lapply(qts_tmp, function(x) {
      quantile(x, c(0.05, 0.5, 0.95), na.rm = TRUE)
    })
    ### coerce into data frame
    cbind(as.data.frame(qts_tmp),
          scenario = scenario, 
          fhist = res_df$fhist[res_df$scenario == scenario],
          HCRmult = res_df$HCRmult[res_df$scenario == scenario],
          b_z = res_df$b_z[res_df$scenario == scenario]
    )
  })
  qts <- do.call(rbind, qts)
  qts$HCRmult <- as.factor(qts$HCRmult)
  qts$b_z <- as.factor(qts$b_z)
  qts$iter <- as.factor(qts$iter)
  levels(qts$iter) <- c("95%", "50%", "5%")
  
  ### reference points / lines
  ### reference lines
  stk_name <- ac(res_df$stock[res_df$scenario == scns_tmp[1]])
  df_refs <- data.frame(qname = c("HCR3.2.1r", "HCR3.2.1f", "HCR3.2.1b",
                                  "SSB", "fbar", "Catch", "Rec", "Lmean"),
                        value = c(rep(1, 3),
                                  refpts[[stk_name]]["msy", "ssb"],
                                  refpts[[stk_name]]["msy", "harvest"],
                                  refpts[[stk_name]]["msy", "yield"],
                                  refpts[[stk_name]]["msy", "rec"],
                                  Lmsy[[stk_name]]))
  
  ### plot with HCRmult as basis
  . <- foreach(HCRmult = unique(qts$HCRmult), 
               .packages = "ggplot2") %dopar% {
  p <- ggplot(qts[qts$HCRmult == HCRmult, ],
             aes(x = year, y = data, colour = b_z, 
                 linetype = iter, alpha = iter)) +
   geom_vline(xintercept = 100, colour = "grey") +
   geom_hline(data = df_refs, aes(yintercept = value), lty = "dashed") +
   geom_line() +
   facet_wrap(~ qname, scales = "free_y") +
   scale_linetype_manual("percentile",
                         values = c("dashed", "solid", "dashed")) +
   scale_alpha_manual("percentile", values = c(0.5, 1, 0.5)) +
   scale_colour_discrete("b_z") +
   theme_bw() +
   ylim(0, NA) +
   labs(x = "year", y = "", title = 
          paste0(df_tmp$fhist[df_tmp$stk_pos == stk_pos][1], ", ",
                 df_tmp$stock[df_tmp$stk_pos == stk_pos][1], ", ",
                 "HCRmult = ", 
                 format(as.numeric(as.character(HCRmult)), nsmall = 2)))
  ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                          "multiplier_exponent/stock_plots/", 
                          df_tmp$fhist[df_tmp$stk_pos == stk_pos][1], "_",
                          df_tmp$stock[df_tmp$stk_pos == stk_pos][1], "_",
                          "HCRmult", 
                          format(as.numeric(as.character(HCRmult)), 
                                 nsmall = 2), 
                          ".png"),
        width = 30, height = 18, units = "cm", dpi = 300, type = "cairo-png",
        plot = p)
  }
  ## plot with b_z as basis
  . <- foreach(b_z = unique(qts$b_z), 
               .packages = "ggplot2") %dopar% {
    p <- ggplot(qts[qts$b_z == b_z, ],
               aes(x = year, y = data, colour = HCRmult, 
                   linetype = iter, alpha = iter)) +
     geom_vline(xintercept = 100, colour = "grey") +
     geom_hline(data = df_refs, aes(yintercept = value), lty = "dashed") +
     geom_line() +
     facet_wrap(~ qname, scales = "free_y") +
     scale_linetype_manual("percentile",
                           values = c("dashed", "solid", "dashed")) +
     scale_alpha_manual("percentile", values = c(0.5, 1, 0.5)) +
     scale_colour_discrete("HCRmult") +
     theme_bw() +
     ylim(0, NA) +
     labs(x = "year", y = "", title = 
            paste0(df_tmp$fhist[df_tmp$stk_pos == stk_pos][1], ", ",
                   df_tmp$stock[df_tmp$stk_pos == stk_pos][1], ", ",
                   "b_z = ", 
                   format(as.numeric(as.character(b_z)), nsmall = 2)))
    ggsave(filename = paste0("output/perfect_knowledge/plots/3.2.1/",
                            "multiplier_exponent/stock_plots/", 
                            df_tmp$fhist[df_tmp$stk_pos == stk_pos][1], "_",
                            df_tmp$stock[df_tmp$stk_pos == stk_pos][1], "_",
                            "b_z", 
                            format(as.numeric(as.character(b_z)), 
                                   nsmall = 2), 
                            ".png"),
          width = 30, height = 18, units = "cm", dpi = 300, type = "cairo-png",
          plot = p)
  }
}

### ------------------------------------------------------------------------ ###
### more stats ####
### yield vs yield achieved when fished at Fmsy
### inter-annual variation
### F/Fmsy, B/Bmsy
### ------------------------------------------------------------------------ ###
### reference points
refpts <- readRDS("input/refpts.rds")

### get sorted: find corresponding stk_pos and stock names
# stk_names <- names(readRDS("input/stock_list.rds"))
# stks_pos_names <- data.frame(stock = stk_names, stk_pos = 1:58)
# write.csv(x = stks_pos_names, file = "input/stock_names_pos.csv", 
#           row.names = FALSE)
stks_pos_names <- read.csv(file = "input/stock_names_pos.csv", as.is = TRUE)

### first fish stocks at Fmsy for entire simulation period & calculate catch
### go through all stocks
res <- foreach(stk_pos = stks_pos_names$stk_pos, 
               refpts_stk = refpts[stks_pos_names$stock],
               .packages = c("FLash"), .errorhandling = "pass") %dopar% {
                 # if (stk_pos > 2) stop()
  #browser()
  ### load OM
  load(paste0("input/stocks/observation_error/", stk_pos, ".RData"))
  
  ### fish at FMSY
  ctrl <- fwdControl(data.frame(year = 101:200, 
                               quantity = "f",
                               val = c(refpts_stk["msy", "harvest"])))
  stk_fwd <- fwd(stk, ctrl = ctrl, sr = sr.om, sr.residuals = sr.om.res,
                sr.residuals.mult = TRUE, maxF = 5)
  
  ### sum catch over projection period per iter
  catch <- apply(window(catch(stk_fwd), start = 101), 6, sum)
  ### sum up SSB
  ssb <- apply(window(ssb(stk_fwd), start = 101), 6, sum)
  
  return(list(catch = catch, ssb = ssb))

}
names(res) <- stks_pos_names$stock

### separate catch and ssb
res_catch <- lapply(res, "[[", "catch")
res_ssb <- lapply(res, "[[", "ssb")
### save them
saveRDS(object = res_catch, file = "input/all_stocks_catch_at_Fmsy.rds")
saveRDS(object = res_ssb, file = "input/all_stocks_ssb_at_Fmsy.rds")
res_catch <- readRDS("input/all_stocks_catch_at_Fmsy.rds")
res_ssb <- readRDS("input/all_stocks_ssb_at_Fmsy.rds")

### select scenarios to calculate yield
scns_stats <- 1237:6804

### load stats calculated earlier
res_df <- readRDS("output/stats_scn_new.RDS")

### life-history parameters & refpts
lhist <- readRDS("input/lhist_extended.rds")

res_df_tmp <- res_df[res_df$scenario %in% scns_stats, ]
#res_df_tmp <- res_df_tmp[1:100, ]

### loop through scenarios
res_tmp <- foreach(scenario = res_df_tmp$scenario, 
  pars = split(lhist[match(ac(res_df_tmp$stock), lhist$stock), ], 
                                1:nrow(res_df_tmp)),
  catch_msy = res_catch[res_df_tmp$stk_pos],
  ssb_msy = res_ssb[res_df_tmp$stk_pos],
                   .packages = "FLCore",
                   .errorhandling = "stop") %dopar% {
                     
  #browser()
  ### load quants
  # tmp <- readRDS(paste0("output/perfect_knowledge/combined/3.2.1_quants/",
  #                      scenario, ".rds"))
  tmp <- readRDS(paste0("/gpfs/afmcefas/simonf/output/combined/3.2.1_quants/",
                        scenario, ".rds"))
  
  ### load summary
  # stk <- readRDS(paste0("output/perfect_knowledge/combined/",
  #                      scenario, ".rds"))
  stk <- readRDS(paste0("/gpfs/afmcefas/simonf/output/combined/",
                        scenario, ".rds"))
  
  ### replace 0 catch in case SSB is 0, i.e. stock collapsed
  ### otherwise would imply stability
  catch <- tmp$Catch
  catch@.Data[which(tmp$SSB == 0)] <- NA
  
  variability <- iav(object = catch, period = 2, start = 99, 
                    summary_per_iter = mean, summary = NULL)
  
  ### F / Fmsy
  f_rel <- tmp$fbar / pars$fmsy
  ### B / Bmsy
  ssb_rel <- tmp$SSB / pars$bmsy
  
  ### average per iteration, median of that
  f_rel <- apply(f_rel, 6, mean, na.rm = TRUE)
    #apply(apply(f_rel, 6, mean, na.rm = TRUE), 1:2, median, na.rm = TRUE)
  
  ssb_rel <- apply(ssb_rel, 6, mean, na.rm = TRUE)
    #apply(apply(ssb_rel, 6, mean, na.rm = TRUE), 1:2, 
    #              median, na.rm = TRUE)
  
  ### catch & SSB relative to when fished at MSY
  ### load stock/quants from MSE simulation
  catch_MSE <- tmp$Catch
  catch_MSE <- apply(window(catch_MSE, start = 101), 6, sum)
  ssb_MSE <- tmp$SSB
  ssb_MSE <- apply(window(ssb_MSE, start = 101), 6, sum)
  
  ### calculate relative values
  catch_MSY_prop <- catch_MSE / catch_msy
  ssb_MSY_prop <- ssb_MSE / ssb_msy
  
  ### recalculate risks and keep iterations
  ssb <- window(tmp$SSB, start = 101)
  collapse <- apply(ssb, c(1, 6), function(x) {
    sum(x < 1)/100
  })
  
  return(list(iav = variability, f_rel = f_rel, ssb_rel = ssb_rel,
              catch_MSY_prop = catch_MSY_prop, ssb_MSY_prop = ssb_MSY_prop,
              collapse_risk = collapse))
  
}
names(res_tmp) <- res_df_tmp$scenario
### save results with iterations
saveRDS(object = res_tmp, file = "output/1237_6804_stats_iter.rds")
res_tmp <- readRDS(file = "output/1237_6804_stats_iter.rds")

### extract medians
res_tmp2 <- foreach(tmp = res_tmp, .packages = "FLCore",
                   .errorhandling = "pass") %dopar% {
  sapply(tmp, median, na.rm = TRUE)
}
res_tmp2 <- as.data.frame(do.call(rbind, res_tmp2))

### merge with stats table
res_tmp2$scenario <- res_df_tmp$scenario
res_df2 <- merge(res_df, res_tmp2, all = TRUE)
### save
saveRDS(object = res_df2, file = "output/stats_scn_new.RDS")
write.csv(x = res_df2, file = "output/stats_scn_new.csv", row.names = FALSE)



