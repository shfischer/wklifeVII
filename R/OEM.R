### ------------------------------------------------------------------------ ###
### MSE for ICES WKLIFE VII 2017 ####
### create initial observations and add uncertainty
### ------------------------------------------------------------------------ ###

rm(list=ls())
### load packages
required_pckgs <- c("FLash", "FLAssess", "ggplotFL", "FLBRP", "data.table")
### save as object in order to avoid output to screen
. <- lapply(required_pckgs, function(x){
  suppressMessages(library(x, character.only = TRUE))
})

### load additional functions
source("functions/oFun.R")
source("functions/fFun.R")
set.seed(0)

### set up parallel computing environment
library(doParallel)
cl <- makeCluster(28)
registerDoParallel(cl)

### ------------------------------------------------------------------------ ###
### load stock objects
### ------------------------------------------------------------------------ ###
load("om_list.RData")
### add some details
OM_list <- lapply(seq_along(OM_list), function(x){
  c(OM_list[[x]], name = names(OM_list)[x], stk_pos = x)
})

### ------------------------------------------------------------------------ ###
### "loop" trough stocks
### ------------------------------------------------------------------------ ###

OM_list <- foreach(x = OM_list, .export = setdiff(ls(), "OM_list"),
                   .packages = required_pckgs) %dopar% {
  
  ### set seed 
  set.seed(1)
  
  ### ---------------------------------------------------------------------- ###
  ### create index, based on stock object
  ### ---------------------------------------------------------------------- ###
  
  ### get ages
  ages <- an(dimnames(stock.n(x$stk))[["age"]])
  
  ### define model for selectivity: logistic function
  ### inflection point of curve = 10% of max age
  q_model <- FLModelSim(model = ~max_q/(1+exp(-steepness*(age - age50))), 
                        params = FLPar(max_q = 1, steepness = 1, 
                                       age50 = max(ages)/10))
  ### model selectivity
  q_modeled <- predict(q_model, age = ages)
  
  ### create index template
  idx <- FLIndex(index = stock.n(x$stk))
  
  ### insert selectivity
  index.q(idx) <- c(q_modeled)
  
  ### add uncertainty
  ### no uncertainty for perfect knowledge scenarios
  ### log-normal noise, cv = 0.2
  #index.q(idx) <- index.q(idx) * rlnorm(n = length(index.q(idx)), sdlog = 0.2)
  
  ### calculate historical index values
  index(idx) <- index.q(idx) * stock.n(x$stk) * stock.wt(x$stk)
  
  # plot(FLQuants(idx = quantSums(index(idx)),
  #               tsb = quantSums(stock.n(x$stk)*stock.wt(x$stk)),
  #               ssb = ssb(x$stk)))
  
  ### save in observations object
  x$observations <- list(idx = FLIndices(idx = idx))
  
  ### ---------------------------------------------------------------------- ###
  ### estimate some required reference points
  ### ---------------------------------------------------------------------- ###
  
  ### I_lim = lowest observed index value
  I_lim <- apply(quantSums(index(x$observations[[1]][["idx"]])), 6, min, 
                 na.rm = TRUE)
  ### I_trigger not defined yet
  I_trigger <- NA
  
  ### add noise
  #I_lim <- I_lim * rlnorm(n = length(I_lim), sdlog = 0.05)
  #I_trigger <- I_trigger * rlnorm(n = length(I_trigger), sdlog = 0.05)
  
  ### add reference points to stk
  attr(x = x$stk, which = "refpts") <- FLPar(I_lim = I_lim,
                                           I_trigger = I_trigger)
  
  ### ---------------------------------------------------------------------- ###
  ### catch length frequencies
  ### ---------------------------------------------------------------------- ###
  
  ### create FLStockLen with catch length frequencies
  stk_len <- length_freq(x$stk, full_series = TRUE,
                         lhpar = attr(x$stk, "lhpar"),
                         len_noise_sd = 0.0, len_sd = 1,
                         len_sd_cut = 2)
  ### save catch.n as attribute in stk
  attr(x$stk, "catch_len") <- window(catch.n(stk_len), end = iy)
  rm(stk_len)
  
  ### try plotting
  # ggplot(as.data.frame(attr(x$stk, "catch_len")[, ac(75)]),
  #        aes(x = length, y = data, colour = as.factor(iter))) +
  #   geom_line(show.legend = FALSE) + geom_point(show.legend = FALSE) +
  #   facet_wrap(~ year)
  
  ### ---------------------------------------------------------------------- ###
  ### extract/rename required objects
  ### ---------------------------------------------------------------------- ###
  
  ### stock recruitment model
  names(x)[names(x) == "srbh"] <- "sr.om"
  ### recruitment residuals
  names(x)[names(x) == "srbh.res"] <- "sr.om.res"
  
  ### ---------------------------------------------------------------------- ###
  ### create reference points with perfect knowledge
  ### ---------------------------------------------------------------------- ###
  
  ### create stock without any uncertainy from brp
  stk_tmp <- as(x$brp, "FLStock")
  
  ### target FMSY for 100 years
  ctrl_tmp <- fwdControl(data.frame(year = 2:101, quantity = "f",
                      val = c(refpts(x$brp)["msy", "harvest"])))
  ### project
  stk_tmp <- fwd(stk_tmp, ctrl = ctrl_tmp, sr = x$sr.om)
  
  ### extract also C/I from here
  ### create perfect index value
  idx_tmp <- c(q_modeled) * (stock.wt(stk_tmp) * stock.n(stk_tmp))[, ac(100)]
  ### ratio catch/index
  I_F_proxy <- catch(stk_tmp)[, ac(100)] / quantSums(idx_tmp)
  
  ### create catch length frequencies, without uncertainty
  len_tmp <- length_freq(stk_tmp[, ac(100)], full_series = TRUE, 
                         lhpar = FLCore::iter(attr(x$stk, "lhpar"), 1),
                         len_noise_sd = 0.0, len_sd = 1,
                         len_sd_cut = 2)
  LFeFmsy <- calc_mean(catch.n(len_tmp), min = calc_Lc(catch.n(len_tmp)), 
                       include_min = FALSE)
  
  ### targt M
  ctrl_tmp <- fwdControl(data.frame(year = 2:101, quantity = "f",
                                    val = c(x$stk@lhpar["M", 1])))
  ### project
  stk_tmp <- fwd(stk_tmp, ctrl = ctrl_tmp, sr = x$sr.om)
  ### create catch length frequencies, without uncertainty, no spreading
  len_tmp <- length_freq(stk_tmp[, ac(100)], full_series = TRUE, 
                         lhpar = FLCore::iter(attr(x$stk, "lhpar"), 1),
                         len_noise_sd = 0.0, len_sd = 1,
                         len_sd_cut = 2)
  LFeM <- calc_mean(catch.n(len_tmp), min = calc_Lc(catch.n(len_tmp)), 
                       include_min = FALSE)
  
  ### save reference points in refpts attribute
  attr(x$stk, "refpts") <- rbind2(attr(x$stk, "refpts"), 
               FLPar(LFeFmsy = c(LFeFmsy), LFeM = c(LFeM),
                     F_MSY = refpts(x$brp)["msy", "harvest"],
                     F0.1 = refpts(x$brp)["f0.1", "harvest"],
                     I_F_proxy = I_F_proxy,
                     iter = dims(x$stk@lhpar)$iter))
  
  return(x)
  
}

### ---------------------------------------------------------------------- ###
### calculate F0.1 with YPR
### ---------------------------------------------------------------------- ###
OM_list <- foreach(x = OM_list, .packages = required_pckgs) %dopar% {
  
  ### needed for catch rule 3.2.1 from WKMSYCat34, option c for factor f
  
  ### run YPR for range of F's
  F_list <- seq(0, 3, 0.01)
  
  ### inverse von Bertalanffy growth function
  inv_vB <- function(L, L_inf, K, t0){
    res <- -log(1 - (L / L_inf)) / K + t0
    return(res)
  }
  
  ### length at first capture
  L_c <- yearMeans(calc_Lc(attr(x$stk, "catch_len")))
  ### convert into age  ### calculate t_c from L_c
  t_c <- inv_vB(L = L_c, 
                L_inf = attr(x$stk, "lhpar")["L_inf"], 
                K = attr(x$stk, "lhpar")["K"], t0 = attr(x$stk, "lhpar")["t0"])
  t_c <- round(t_c)  
  
  ### calculate F0.1 for 1 iteration
  F0.1 <- lapply((1:dims(x$stk)$iter)[1], function(y){
    
    ### conduct YPR for list of F's
    res_list <- lapply(F_list, function(f){
      
      YPR(a =	attr(x$stk, "lhpar")["a", y], b =	attr(x$stk, "lhpar")["b", y], 
          Linf =	attr(x$stk, "lhpar")["L_inf", y],
          K =	attr(x$stk, "lhpar")["K", y], t0 =	attr(x$stk, "lhpar")["t0", y], 
          tc =	c(FLCore::iter(t_c, y)),
          M =	attr(x$stk, "lhpar")["M", y], F =	f, 
          max_age = attr(x$stk, "lhpar")["max_age", y]) 
      
    })
    ### format
    res_list <- do.call(rbind, res_list)
    
    ### search for F0.1
    ### get slope (step size is equal between all elements)
    slope <- diff(res_list)
    
    ### find first F where slope is < 0.1 * initial slope
    res <- F_list[which(slope < 0.1 * slope[1])[1]]
    
    return(res)
    
  })
  F0.1 <- unlist(F0.1)
  
  ### save in refpts
  attr(x$stk, "refpts") <- rbind2(attr(x$stk, "refpts"),
                                  FLPar(F0.1YPR = F0.1, iter = dims(x$stk@lhpar)$iter))
  
  return(x)
  
}


### ------------------------------------------------------------------------ ###
### save
### ------------------------------------------------------------------------ ###

. <- foreach(x = OM_list) %dopar% {
  
  stk <- x$stk
  observations <- x$observations
  sr.om <- x$sr.om
  sr.om.res <- x$sr.om.res
  
  save(stk, observations, sr.om, sr.om.res, it, fy, y0, dy, iy, ny, nsqy, vy,
       file = paste0("input/stocks/perfect_knowledge/", x$stk_pos, ".RData"))
  
}

### ------------------------------------------------------------------------ ###
### observation error
### ------------------------------------------------------------------------ ###

OM_list <- foreach(x = OM_list, .packages = required_pckgs) %dopar% {
                     
  ### set seed 
  set.seed(1)
  
  ### ---------------------------------------------------------------------- ###
  ### add index uncertainty
  ### ---------------------------------------------------------------------- ###
  
  ### get index
  idx <- x$observations$idx$idx
  
  ### add uncertainty to catchability
  ### log-normal noise, cv = 0.2
  index.q(idx) <- index.q(idx) * rlnorm(n = length(index.q(idx)), sdlog = 0.2)
  
  ### update index values
  index(idx) <- index.q(idx) * stock.n(x$stk) * stock.wt(x$stk)
  
  ### save in observations object
  x$observations <- list(idx = FLIndices(idx = idx))
  
  ### ---------------------------------------------------------------------- ###
  ### update index reference points
  ### ---------------------------------------------------------------------- ###
  ### uncertainty already implemented in index values
  ### the new minimum values already include this uncertainty
  
  attr(x$stk, "refpts")["I_lim"] <- apply(quantSums(index(idx)), 6, min, 
                                          na.rm = TRUE)

  ### ---------------------------------------------------------------------- ###
  ### add uncertainty to life-history parameters and reference points
  ### ---------------------------------------------------------------------- ###
  
  ### reference points
  attr(x$stk, "refpts")["LFeFmsy"] <- attr(x$stk, "refpts")["LFeFmsy"] * 
    rlnorm(n = length(attr(x$stk, "refpts")["LFeFmsy"]), sdlog = 0.1)
  attr(x$stk, "refpts")["LFeM"] <- attr(x$stk, "refpts")["LFeM"] * 
    rlnorm(n = length(attr(x$stk, "refpts")["LFeM"]), sdlog = 0.1)
  attr(x$stk, "refpts")["F_MSY"] <- attr(x$stk, "refpts")["F_MSY"] * 
    rlnorm(n = length(attr(x$stk, "refpts")["F_MSY"]), sdlog = 0.1)
  attr(x$stk, "refpts")["F0.1"] <- attr(x$stk, "refpts")["F0.1"] * 
    rlnorm(n = length(attr(x$stk, "refpts")["F0.1"]), sdlog = 0.1)
  attr(x$stk, "refpts")["I_F_proxy"] <- attr(x$stk, "refpts")["I_F_proxy"] * 
    rlnorm(n = length(attr(x$stk, "refpts")["I_F_proxy"]), sdlog = 0.1)
  
  ### lhist
  attr(x$stk, "lhpar")["L_inf"] <- attr(x$stk, "lhpar")["L_inf"] * 
    rlnorm(n = length(attr(x$stk, "lhpar")["L_inf"]), sdlog = 0.1)
  attr(x$stk, "lhpar")["K"] <- attr(x$stk, "lhpar")["K"] * 
    rlnorm(n = length(attr(x$stk, "lhpar")["K"]), sdlog = 0.1)
  attr(x$stk, "lhpar")["t0"] <- attr(x$stk, "lhpar")["t0"] * 
    rlnorm(n = length(attr(x$stk, "lhpar")["t0"]), sdlog = 0.1)
  # attr(x$stk, "lhpar")["a"] <- attr(x$stk, "lhpar")["a"] * 
  #   rlnorm(n = length(attr(x$stk, "lhpar")["a"]), sdlog = 0.1)
  # attr(x$stk, "lhpar")["b"] <- attr(x$stk, "lhpar")["b"] * 
  #   rlnorm(n = length(attr(x$stk, "lhpar")["b"]), sdlog = 0.1)
  attr(x$stk, "lhpar")["M"] <- attr(x$stk, "lhpar")["M"] * 
    rlnorm(n = length(attr(x$stk, "lhpar")["M"]), sdlog = 0.1)
  
  ### ---------------------------------------------------------------------- ###
  ### catch length frequency uncertainty
  ### ---------------------------------------------------------------------- ###
  
  ### create FLStockLen with catch length frequencies
  stk_tmp <- x$stk
  attr(stk_tmp, "catch_len") <- NULL
  stk_len <- length_freq(stk_tmp, full_series = TRUE,
                         lhpar = attr(x$stk, "lhpar"),
                         len_noise_sd = 0.2, len_sd = 1,
                         len_sd_cut = 2)
  ### save catch.n as attribute in stk
  attr(x$stk, "catch_len") <- window(catch.n(stk_len), end = iy)
  
  return(x)
  
}

### ------------------------------------------------------------------------ ###
### save
### ------------------------------------------------------------------------ ###

. <- foreach(x = OM_list) %dopar% {
  
  stk <- x$stk
  observations <- x$observations
  sr.om <- x$sr.om
  sr.om.res <- x$sr.om.res
  
  save(stk, observations, sr.om, sr.om.res, it, fy, y0, dy, iy, ny, nsqy, vy,
       file = paste0("input/stocks/observation_error/", x$stk_pos, ".RData"))
  
}


stopCluster(cl)
