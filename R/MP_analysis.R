### ------------------------------------------------------------------------ ###
### script for retrieving and analysing MSE results ####
### ------------------------------------------------------------------------ ###

OM_scn <- as.numeric(commandArgs(TRUE))
if (length(OM_scn) == 0) stop("OM_scn missing")
#OM_scn <- 29

library(mseDL)
library("FLife")
library("FLash")
library("dplyr")
library("tidyr")
library("ggplot2")
library("cowplot")

source("MP_functions.R")

### set up parallel computing environment
library(doParallel)
#cl <- makeCluster(parallel::detectCores() - 1)
cl <- makeCluster(24)
registerDoParallel(cl)

### generic stuff
refpts_all <- readRDS("input/refpts_paper.rds")
stk_pos <- read.csv("input/stock_names_pos.csv")
stocks <- read.csv("input/stock_list_full2.csv", stringsAsFactors = FALSE)
names_short <- read.csv("input/names_short.csv", stringsAsFactors = FALSE)
OM_scns <- read.csv("input/OM_scns.csv", stringsAsFactors = FALSE)

stock_list <- data.frame(stock = names_short$stock, 
                         paper = names_short$paper,
                         fhist = rep(c("one-way", "roller-coaster"), each = 29),
                         stk_pos = c(1:15, 31:44, 16:30, 45:58),
                         stringsAsFactors = FALSE)


### ------------------------------------------------------------------------ ###
### summarise results for specific scenario ####
### ------------------------------------------------------------------------ ###

#OM_scn <- 29

### file names
id <- OM_scns$id[OM_scn]
files <- paste0(stock_list$fhist, "/", stock_list$stock, ".rds")
names(files) <- stock_list$paper
dir.create(paste0("output/observation_error/", id, "/", "corrected/",
                  "one-way"), recursive = TRUE)
dir.create(paste0("output/observation_error/", id, "/", "corrected/",
                  "roller-coaster"), recursive = TRUE)
### load results
res_list <- lapply(files, function(x) 
  readRDS(paste0("output/observation_error/", id, "/", x))
)
### and reference points
refpts <- refpts_all[[id]]

### fish Fmsy
res_fmsy <- foreach(res_tmp = res_list,
                    name_i = stock_list$stock, 
                    .packages = c("FLCore", "FLash")) %dopar% {
  stk_tmp <- res_tmp@stock
  ### get reference points
  refpts_i <- refpts[[name_i]]
  ### fish at FMSY
  ctrl <- fwdControl(data.frame(year = 101:dims(stk_tmp)$maxyear, 
                                quantity = "f",
                                val = c(refpts_i["msy", "harvest"])))
  stk_fwd <- fwd(stk_tmp, ctrl = ctrl, sr = res_tmp@sr, 
                 sr.residuals = res_tmp@sr@residuals,
                 sr.residuals.mult = TRUE, maxF = 5)
  return(stk_fwd)
}
### save
# saveRDS(object = res_fmsy, file = paste0("output/observation_error/", id, 
#                                          "/stocks_at_fmsy.rds"))
### get catch and ssb over projection period
catch_MSY <- lapply(res_fmsy, catch)
saveRDS(object = catch_MSY, 
        file = paste0("output/observation_error/", id, 
                      "/stocks_at_fmsy_catch.rds"))
rm(res_fmsy)

### "correct" results
### i.e. make sure that once a stock collapsed, it does not recover afterwards
### and calculate SSB/F/catch relative to MSY
res_corrected <- foreach(name_i = stock_list$stock,
                         res_tmp = res_list,
                         fhist_i = stock_list$fhist,
                         catch_MSY_i = catch_MSY,
                         .packages = "FLCore") %dopar% {
  stk_tmp <- res_tmp@stock
  ### get reference points
  refpts_i <- refpts[[name_i]]
   
  ### find first F=5 for each iteration
  ### return numeric year or NA if F not maxed out
  fmax <- apply(fbar(stk_tmp), 6, FUN = function(x){
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
  catch_rel <- catch / refpts_i["msy", "yield"]
   
  ### save SSB and catch, both versions
  return_tmp <- list(ssb = SSB, catch = catch, fbar = fbar, 
                      fbar_rel = fbar_rel, ssb_rel = ssb_rel, 
                      catch_rel = catch_rel)
  saveRDS(object = return_tmp,
          file = paste0("output/observation_error/", id, "/corrected/",
                        fhist_i, "/", name_i, ".rds"))
  return(return_tmp)
}

### calculate Blim
### SSB where recruitment is impaired by 30%
BevHolt <- function(a, b, ssb) {
  return(a * ssb / (b + ssb))
}
### calculate Blim
ssbs = seq(0, 1000, 0.01)
res_blim <- lapply(res_list, function(x){
  rec <- BevHolt(a = c(x@sr@params["a", 1]),
                 b = c(x@sr@params["b", 1]),
                 ssb = ssbs)
  ssbs[tail(which(rec <= c(max(rec) * 0.7)), 1)]
})
names(res_blim) <- stock_list$paper
blim_old <- 162.79
saveRDS(res_blim, file = paste0("output/observation_error/", id, "/blim.rds"))

### calculate summary stats
stats <- foreach(seq_i = seq_along(stock_list$stock),
                 stk_id = stock_list$stk_pos,
                 name_i = stock_list$stock,
                 fhist_i = stock_list$fhist,
                 paper_i = stock_list$paper,
                 res_tmp = res_corrected,
                 blim_i = res_blim,
                 .packages = "FLCore") %dopar% {
  ### subset to simulation period
  quants <- window(FLQuants(res_tmp), start = 101)
  ### reference points
  refpts_i <- refpts[[name_i]]
  catch_msy_i <- catch_MSY[[seq_i]]
  ### list for storing results
  res <- list()
  ### number of iterations for current scenario
  n_iter <- dim(quants[[1]])[6]
  n_years <- dim(quants[[1]])[2]
  ### replace 0 catch in case SSB is 0, i.e. stock collapsed
  ### otherwise would imply stability
  catchNA <- res_tmp$catch
  catchNA@.Data[which(res_tmp$ssb == 0)] <- NA
  res$iav <- c(iav(object = window(catchNA, start = 99), period = 2, 
                    start = 99, summary_per_iter = mean, summary = median))
  res$iav_short <- c(iav(object = window(catchNA, start = 99, end = 120), 
                          period = 2, start = 99, 
                          summary_per_iter = mean, summary = median))
  res$iav_medium <- c(iav(object = window(catchNA, start = 99 + 20, end = 160), 
                          period = 2, start = 99 + 20, 
                          summary_per_iter = mean, summary = median))
   res$iav_long <- c(iav(object = window(catchNA, start = 99 + 60, end = 200), 
                        period = 2, start = 99 + 60, 
                        summary_per_iter = mean, summary = median))
  ### F / Fmsy
  res$f_rel <- median(apply(res_tmp$fbar_rel, 6, mean, na.rm = TRUE))
  res$f_rel_new <- median(apply(res_tmp$fbar_rel[, ac(101:200)], 6, mean, 
                                na.rm = TRUE))
  res$f_rel_short <- median(apply(res_tmp$fbar_rel[, ac(101:120)], 6, mean, 
                                   na.rm = TRUE))
  res$f_rel_medium <- median(apply(res_tmp$fbar_rel[, ac(121:160)], 6, mean, 
                                    na.rm = TRUE))
  res$f_rel_long <- median(apply(res_tmp$fbar_rel[, ac(161:200)], 6, mean, 
                                  na.rm = TRUE))
  ### B / Bmsy
  res$ssb_rel <- median(apply(res_tmp$ssb_rel, 6, mean, na.rm = TRUE))
  res$ssb_rel_new <- median(apply(res_tmp$ssb_rel[, ac(101:200)], 6, mean, 
                                   na.rm = TRUE))
  res$ssb_rel_short <- median(apply(res_tmp$ssb_rel[, ac(101:120)], 6, mean, 
                                    na.rm = TRUE))
  res$ssb_rel_medium <- median(apply(res_tmp$ssb_rel[, ac(121:160)], 6, mean, 
                                      na.rm = TRUE))
  res$ssb_rel_long <- median(apply(res_tmp$ssb_rel[, ac(161:200)], 6, mean, 
                                    na.rm = TRUE))
  ### catch & SSB relative to when fished at MSY
  ### load stock/quants from MSE simulation
  #catch_MSE <- apply(window(res_tmp$catch, start = 101), 6, sum)
  #ssb_MSE <- apply(window(res_tmp$ssb, start = 101), 6, sum)
  ### calculate relative values
  res$yield_rel_MSY_new <- median(apply(window(res_tmp$catch, start = 101, end = 200), 6, sum) / apply(window(catch_msy_i, start = 101, end = 200), 6, sum))
  res$yield_rel_MSY_short <- median(apply(window(res_tmp$catch, start = 101, end = 120), 6, sum) / apply(window(catch_msy_i, start = 101, end = 120), 6, sum))
  res$yield_rel_MSY_medium <- median(apply(window(res_tmp$catch, start = 121, end = 160), 6, sum) / apply(window(catch_msy_i, start = 121, end = 160), 6, sum))
  res$yield_rel_MSY_long <- median(apply(window(res_tmp$catch, start = 161, end = 200), 6, sum) / apply(window(catch_msy_i, start = 161, end = 200), 6, sum))
  
  ### collapse risk
  res$collapse_total <- sum(res_tmp$ssb < 1) / (n_iter*n_years)
  res$collapse_total_new <- mean(window(res_tmp$ssb, start = 101) < 1)
  res$collapse_total_short <- mean(window(res_tmp$ssb, start = 101, end = 120) < 1)
  res$collapse_total_medium <- mean(window(res_tmp$ssb, start = 121, end = 160) < 1)
  res$collapse_total_long <- mean(window(res_tmp$ssb, start = 161, end = 200) < 1)
  
  ### Blim risk
  res$ssb_below_blim <- sum(res_tmp$ssb < blim_i) / (n_iter*n_years)
  res$ssb_below_blim_new <- mean(window(res_tmp$ssb, start = 101) < blim_i)
  res$ssb_below_blim_short <- mean(window(res_tmp$ssb, start = 101, end = 120) < blim_i)
  res$ssb_below_blim_medium <- mean(window(res_tmp$ssb, start = 121, end = 160) < blim_i)
  res$ssb_below_blim_long <- mean(window(res_tmp$ssb, start = 161, end = 200) < blim_i)
   
  res$stk_pos <- stk_id
  res$stock <- name_i
  res$fhist <- fhist_i
  res$paper <- paper_i
  return(res)
}
stats <- do.call(rbind, lapply(stats, data.frame))
write.csv(stats, file = paste0("output/observation_error/", id, 
                               "/stats.csv"))

