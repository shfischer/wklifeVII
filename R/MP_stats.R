### ------------------------------------------------------------------------ ###
### script for retrieving and analysing MSE results ####
### ------------------------------------------------------------------------ ###

### arguments passed on to R
args <- commandArgs(TRUE)
print("arguments passed on to this script:")
print(args)

### evaluate arguments, if they are passed to R:
if (length(args) > 0) {
  
  ### extract arguments
  for (i in seq_along(args)) eval(parse(text = args[[i]]))
  if (!exists("id")) stop("id missing")
  if (!exists("tuning")) tuning <- FALSE
  if (!exists("correct")) correct <- FALSE
  if (!exists("calc_stats")) calc_stats <- TRUE
  
} else {
  
  stop("no argument passed to R")
  
}

#id <- 38

### load packages
library(mseDL)
library(FLife)
library(FLash)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

source("MP_functions.R")

### set up parallel computing environment
library(doParallel)
#cl <- makeCluster(parallel::detectCores() - 1)
cl <- makeCluster(24)
registerDoParallel(cl)

### generic stuff
refpts_all <- readRDS("input/refpts_paper.rds")
stocks <- read.csv("input/stock_list_full2.csv", stringsAsFactors = FALSE)
names_short <- read.csv("input/names_short.csv", stringsAsFactors = FALSE)
OM_scns <- read.csv("input/OM_scns.csv", stringsAsFactors = FALSE)
### subuset to current OM scenario definition
OM_scn <- OM_scns[OM_scns$idSEQ == id, ]

### observation error?
obs_error <- ifelse(OM_scn$obs_error, "observation_error",
                    "perfect_knowledge")

stock_list <- data.frame(stock = names_short$stock, 
                         paper = names_short$paper,
                         fhist = rep(c("one-way", "roller-coaster"), each = 29),
                         stk_pos = c(1:15, 31:44, 16:30, 45:58),
                         stringsAsFactors = FALSE)


### ------------------------------------------------------------------------ ###
### summarise results for specific scenario ####
### ------------------------------------------------------------------------ ###

### file names
path <- paste0("output/", obs_error, "/", OM_scn$id, "/")

### analyse default rule or tuning options?
if (isFALSE(tuning)) {
  lapply(paste0("output/", obs_error, "/", OM_scn$id, "/", "corrected/",
                  c("one-way", "roller-coaster")), dir.create, recursive = TRUE)
  tun_dirs <- ""
} else {
  tun_dirs <- list.dirs(path = paste0(path, "one-way"), full.names = FALSE)[-1]
  tun_opt <- strsplit(x = tun_dirs, split = "_")
  tun_opt <- lapply(tun_opt, function(x) {
    tmp <- strsplit(x, split = "-")
    names(tmp) <- sapply(tmp, "[[", 1)
    tmp <- sapply(tmp, "[[", 2)
    tmp <- sapply(tmp, as.numeric)
    return(tmp)
  })
  tun_opt <- do.call(bind_rows, tun_opt)
  tun_opt$dir <- tun_dirs
  
  tun_path <- apply(expand.grid(paste0("output/", obs_error, "/", OM_scn$id, 
                                       "/corrected"),
                           paste0(c("one-way", "roller-coaster")), tun_dirs), 1,
                    paste0, collapse = "/")
  lapply(tun_path, dir.create, recursive = TRUE)
}

### load reference points for OM
refpts <- refpts_all[[OM_scn$id]]


### go through scenarios
stats_scns <- foreach(i_scn = seq_along(tun_dirs)) %do% {
  
  ### ---------------------------------------------------------------------- ###
  ### calculate corrected quants ####
  ### ---------------------------------------------------------------------- ###
  
  ### file names with MSE results
  files <- paste0(stock_list$fhist, "/", tun_dirs[[i_scn]], "/", 
                  stock_list$stock, ".rds")
  names(files) <- stock_list$paper
  
  ### load results
  res_list <- lapply(files, function(x) 
    readRDS(paste0("output/", obs_error, "/", OM_scn$id, "/", x))
  )
  
  if (isTRUE(correct) & isFALSE(tuning)) {
      
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
            file = paste0("output/", obs_error, "/", OM_scn$id, 
                          "/stocks_at_fmsy_catch.rds"))
    rm(res_fmsy)
    
  } else {
    
    ### get catch at Fmsy
    catch_MSY <- readRDS(paste0("output/", obs_error, "/", OM_scn$id,
                                "/stocks_at_fmsy_catch.rds"))
    
  }
  
  ### "correct" results
  ### i.e. make sure that once a stock collapsed, it does not recover afterwards
  ### and calculate SSB/F/catch relative to MSY
  
  res_corrected <- foreach(name_i = stock_list$stock,
                           res_tmp = res_list,
                           fhist_i = stock_list$fhist,
                           catch_MSY_i = catch_MSY,
                           .packages = "FLCore") %dopar% {
    
    if (isTRUE(correct)) {
      
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
              file = paste0("output/", obs_error, "/", OM_scn$id, "/corrected/",
                            fhist_i, "/", tun_dirs[[i_scn]], "/", name_i, 
                            ".rds"))
      return(return_tmp)
      
    } else {
      
      return_tmp <- readRDS(paste0("output/", obs_error, "/", OM_scn$id, 
                                   "/corrected/", fhist_i, "/", 
                                   tun_dirs[[i_scn]], "/", name_i, ".rds"))
      return(return_tmp)
      
    }
    
  }
  
  if (isTRUE(correct) & isFALSE(tuning)) {
    
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
    saveRDS(res_blim, file = paste0("output/", obs_error, "/", OM_scn$id, 
                                    "/blim.rds"))
    
  }
  
  res_blim <- readRDS(paste0("output/", obs_error, "/", OM_scn$id, "/blim.rds"))
  
  ### ---------------------------------------------------------------------- ###
  ### calculate stats ####
  ### ---------------------------------------------------------------------- ###
  
  if (isTRUE(calc_stats)) {
  
    ### calculate summary stats
    stats <- foreach(seq_i = seq_along(stock_list$stock),
                     stk_id = stock_list$stk_pos,
                     name_i = stock_list$stock,
                     fhist_i = stock_list$fhist,
                     paper_i = stock_list$paper,
                     res_tmp = res_corrected,
                     blim_i = res_blim,
                     .packages = "FLCore") %dopar% {
      ### reference points
      catch_msy_i <- catch_MSY[[seq_i]]
      ### list for storing results
      res <- list()
      ### for ICV only: replace 0 catch in case SSB is 0, i.e. stock collapsed
      ### otherwise would imply stability
      catchNA <- res_tmp$catch
      catchNA@.Data[which(res_tmp$ssb == 0)] <- NA
      ### inter-annual catch variability
      res$iav_old <- c(iav(object = window(catchNA, start = 99), period = 2, 
                           from = 99, summary_per_iter = mean, 
                           summary_year = median))
      res$iav <- c(iav(object = catchNA, period = 2,
                       from = 99, summary_all = median))
      res$iav_short <- c(iav(object = catchNA, period = 2, from = 99, to = 120,
                              summary_all = median))
      res$iav_medium <- c(iav(object = catchNA, period = 2, from = 119, to = 160,
                              summary_all = median))
      res$iav_long <- c(iav(object = catchNA, period = 2, from = 159, to = 200,
                              summary_all = median))
      
      ### F / Fmsy
      res$f_rel_old <- median(apply(res_tmp$fbar_rel[, ac(101:200)], 6, mean, 
                                    na.rm = TRUE))
      res$f_rel <- median(res_tmp$fbar_rel[, ac(101:200)], na.rm = TRUE)
      res$f_rel_short <- median(res_tmp$fbar_rel[, ac(101:120)], na.rm = TRUE)
      res$f_rel_medium <- median(res_tmp$fbar_rel[, ac(121:160)], na.rm = TRUE)
      res$f_rel_long <- median(res_tmp$fbar_rel[, ac(161:200)], na.rm = TRUE)
  
      ### B / Bmsy
      res$ssb_rel_old <- median(apply(res_tmp$ssb_rel[, ac(101:200)], 6, mean, 
                                       na.rm = TRUE))
      res$ssb_rel <- median(res_tmp$ssb_rel[, ac(101:200)], na.rm = TRUE)
      res$ssb_rel_short <- median(res_tmp$ssb_rel[, ac(101:120)], na.rm = TRUE)
      res$ssb_rel_medium <- median(res_tmp$ssb_rel[, ac(121:160)], na.rm = TRUE)
      res$ssb_rel_long <- median(res_tmp$ssb_rel[, ac(161:200)], na.rm = TRUE)
      
      ### catch/MSY
      res$yield_rel_old <- median(apply(res_tmp$catch_rel[, ac(101:200)], 6, 
                                        mean, na.rm = TRUE))
        median(apply(window(res_tmp$catch, start = 101, 
                                               end = 200), 6, sum) / 
                                    apply(window(catch_msy_i, start = 101, 
                                                 end = 200), 6, sum))
      res$yield_rel <- median(res_tmp$catch_rel[, ac(101:200)], na.rm = TRUE)
      res$yield_rel_short <- median(res_tmp$catch_rel[, ac(101:120)], na.rm = TRUE)
      res$yield_rel_medium <- median(res_tmp$catch_rel[, ac(121:160)], 
                                     na.rm = TRUE)
      res$yield_rel_long <- median(res_tmp$catch_rel[, ac(161:200)], na.rm = TRUE)
      
      ### collapse risk
      res$risk_collapse <- mean(c(res_tmp$ssb[, ac(101:200)]) < 1, na.rm = TRUE)
      res$risk_collapse_short <- mean(c(res_tmp$ssb[, ac(101:120)]) < 1, 
                                      na.rm = TRUE)
      res$risk_collapse_medium <- mean(c(res_tmp$ssb[, ac(121:160)]) < 1,
                                       na.rm = TRUE)
      res$risk_collapse_long <- mean(c(res_tmp$ssb[, ac(161:200)]) < 1, 
                                     na.rm = TRUE)
      
      ### Blim risk
      res$risk_blim <- mean(c(res_tmp$ssb[, ac(101:200)]) < blim_i, na.rm = TRUE)
      res$risk_blim_short <- mean(c(res_tmp$ssb[, ac(101:120)]) < blim_i,
                                  na.rm = TRUE)
      res$risk_blim_medium <- mean(c(res_tmp$ssb[, ac(121:160)]) < blim_i,
                                   na.rm = TRUE)
      res$risk_blim_long <- mean(c(res_tmp$ssb[, ac(161:200)]) < blim_i,
                                 na.rm = TRUE)
       
      res$stk_pos <- stk_id
      res$stock <- name_i
      res$fhist <- fhist_i
      res$paper <- paper_i
      return(res)
    }
    stats <- do.call(rbind, lapply(stats, data.frame))
    if (isTRUE(tuning)) {
      stats <- bind_cols(stats, tun_opt[rep(i_scn, nrow(stats)), ])
    }
    
    write.csv(stats, file = paste0("output/", obs_error, "/", OM_scn$id, 
                                   "/stats", tun_dirs[[i_scn]], ".csv"))
    
  }
  
  return(stats)
  
}

### combine all tuning options
if (isTRUE(calc_stats) & isTRUE(tuning)) {
  
  stats_tuning <- do.call(bind_rows, stats_scns)
  write.csv(stats_tuning, file = paste0("output/", obs_error, "/", OM_scn$id, 
                                 "/stats_tuning.csv"))
  
}

