### ------------------------------------------------------------------------ ###
### set-up OMs for MSE simulation ####
### set additional parameters for simulation
### ------------------------------------------------------------------------ ###
### the R session needs to be restarted as FLash and FLasher environments
### interfere

### scenarios (OM sets) to run
scns <- 29

required_pckgs <- c("FLash", "FLAssess", "ggplotFL", "FLBRP", "data.table")
### save as object in order to avoid output to screen
. <- lapply(required_pckgs, function(x){
  suppressMessages(library(x, character.only = TRUE))
})
library(doParallel)

cl <- makeCluster(parallel::detectCores())
registerDoParallel(cl)

### load additional functions
source("MP_functions.R")

set.seed(0)

### ------------------------------------------------------------------------ ###
### load OMs
### ------------------------------------------------------------------------ ###

### OM scenarios (M, etc.)
OM_scns <- read.csv("input/OM_scns.csv", stringsAsFactors = FALSE)
OM_scns <- OM_scns[scns, ]
### stock list
stocks <- read.csv("input/stock_list_full2.csv", as.is = TRUE)

### create folders for storing OMs
for (id in OM_scns$id) {
  dir.create(path = paste0("input/stocks/perfect_knowledge/", id), 
             recursive = TRUE)
}

### ------------------------------------------------------------------------ ###
### specify dimensions of simulation ####
### ------------------------------------------------------------------------ ###

# it <- dims(OM_list[[1]]$stk)$iter # iterations
# fy <- dims(OM_list[[1]]$stk)$maxyear + 100 # final year
# y0 <- range(OM_list[[1]]$stk)[["minyear"]] # initial data year
# dy <- range(OM_list[[1]]$stk)[["maxyear"]] # final data year
# iy <- dims(OM_list[[1]]$stk)$maxyear # initial year of projection (also intermediate)
# ny <- fy - iy + 1 # number of years to project from intial year
nsqy <- 3 # number of years to compute status quo metrics
# vy <- ac(iy:fy) # vector of years to be projected


### ------------------------------------------------------------------------ ###
### "loop" through all stocks ####
### ------------------------------------------------------------------------ ###

OM_list <- foreach(OM_scn = split(OM_scns, 1:nrow(OM_scns)),
                   OM_scn_i = OM_scns$idSEQ,
                   OM_scn_id = OM_scns$id,
                   OM_iter = OM_scns$n_iter,
                   OM_n_years = OM_scns$n_years,
                   OM_I_trigger = OM_scns$I_trigger,
                   .errorhandling = "stop"
                   ) %:%
  foreach(x = rep(stocks$X, 2),
          lh_i = rep(split(stocks, 1:nrow(stocks)), 2),
          fhist_i = rep(c("one-way", "roller-coaster"), 
                        each = length(stocks$X)),
          #x = seq_along(OM_list[[OM_scn_i]]), 
          #OM_i = OM_list[[OM_scn_i]],
          .export = c("nsqy", "getCtrl"),
          .errorhandling = "stop", 
          .packages = c("FLash", "FLAssess", "ggplotFL", "FLBRP", 
                        "data.table")) %dopar% {
  #browser()
  ### load OM
  OM_i <- readRDS(paste0("input/OM1/", OM_scn$id, "/", fhist_i, "/",
                         lh_i$stock, ".rds"))
  set.seed(0)
  OM_iy <- dims(OM_i$stk)$maxyear
  ### simulation years to add
  OM_sim_yrs <- (OM_iy + OM_n_years) - OM_iy
  
  ### ---------------------------------------------------------------------- ###
  ### create residuals for SRR
  ### ---------------------------------------------------------------------- ###
  
  ### create residuals
  ### values from 2nd argument are used as exp(value), using 0 leads to values 
  ### distributed around 1
  #if (OM_scn_i %in% c(1, 17:19)) OM_rec_sd <- 0.3
  srbh.res <- FLife::rlnoise(n = OM_iter, 
    FLQuant(0, dimnames = list(year = ac(OM_iy:(OM_iy + OM_n_years)))), 
            sd = OM_scn$rec_sd_proj, b = OM_scn$rec_rho)
  OM_i$srbh.res <- srbh.res
  
  names(OM_i)[names(OM_i) == "sr"] <- "srbh"
  
  ### ---------------------------------------------------------------------- ###
  ### add life-history parameters to stock ####
  ### ---------------------------------------------------------------------- ###
  
  ### extract lhpar
  lhpar <- OM_i$lhpar
  ### add some more parameters
  ### natural mortality and max age
  add_pars <- FLPar(M = mean(m(OM_i$stk[, 1])),
                    max_age = range(OM_i$stk)[["max"]])
  ### add
  lhpar <- rbind2(lhpar, add_pars)
  ### change some names
  dimnames(lhpar)$params[dimnames(lhpar)$params %in% c("linf", "k")] <- 
   c("L_inf", "K")
  ### propagate
  lhpar <- propagate(lhpar, dims(OM_i$stk)$iter)
  ### save as attribute in stk
  attr(OM_i$stk, "lhpar") <- lhpar
  
  ### ---------------------------------------------------------------------- ###
  ### get reference points from FLBRP
  ### ---------------------------------------------------------------------- ###
  
  attr(OM_i$stk, "refpts") <- refpts(OM_i$brp)
  
  ### ---------------------------------------------------------------------- ###
  ### set up operating model: extend ####
  ### ---------------------------------------------------------------------- ###
  
  OM_i$stk <- stf(OM_i$stk, OM_sim_yrs, nsqy, nsqy)
  
  ### ---------------------------------------------------------------------- ###
  ### create initial observations ####
  ### ---------------------------------------------------------------------- ###
    
  ### set seed 
  set.seed(1)
  
  ### ---------------------------------------------------------------------- ###
  ### create index, based on stock object
  ### ---------------------------------------------------------------------- ###
  
  ### get ages
  ages <- an(dimnames(stock.n(OM_i$stk))[["age"]])
  
  ### define model for selectivity: logistic function
  ### inflection point of curve = 10% of max age
  q_model <- FLModelSim(model = ~max_q/(1+exp(-steepness*(age - age50))), 
                       params = FLPar(max_q = 1, steepness = 1, 
                                      age50 = max(ages)/10))
  ### model selectivity
  q_modeled <- predict(q_model, age = ages)
  
  ### create index template
  idx <- FLIndex(index = stock.n(OM_i$stk))
  
  ### insert selectivity
  index.q(idx) <- c(q_modeled)
  
  ### overwrite selectivity with maturity in case of perfect knowledge
  ### i.e. index is simply SSB
  if (isTRUE(OM_I_trigger == "B_trigger")) {
    index.q(idx) <- mat(OM_i$stk)
  }
  
  ### add uncertainty
  ### no uncertainty for perfect knowledge scenarios
  ### log-normal noise, cv = 0.2
  #index.q(idx) <- index.q(idx) * rlnorm(n = length(index.q(idx)), sdlog = 0.2)
  
  ### calculate historical index values

  ### if SSB index, reduce index by mortality
  if (isTRUE(OM_I_trigger == "B_trigger")) {
    index(idx) <- index.q(idx) * stock.n(OM_i$stk) * stock.wt(OM_i$stk) *
      exp(-(harvest(OM_i$stk) * harvest.spwn(OM_i$stk) +
                            m(OM_i$stk) * m.spwn(OM_i$stk)))
  } else {
  ### otherwise default calculation
    index(idx) <- index.q(idx) * stock.n(OM_i$stk) * stock.wt(OM_i$stk)
  }
  
  # plot(FLQuants(idx = quantSums(index(idx)),
  #               tsb = quantSums(stock.n(x$stk)*stock.wt(x$stk)),
  #               ssb = ssb(x$stk)))
  
  ### save in observations object
  OM_i$observations <- list(idx = FLIndices(idx = idx))
  
  ### ---------------------------------------------------------------------- ###
  ### estimate some required reference points
  ### ---------------------------------------------------------------------- ###
  
  ### I_lim = lowest observed index value
  I_lim <- apply(quantSums(index(OM_i$observations[[1]][["idx"]])), 6, min, 
                na.rm = TRUE)
  ### I_trigger not defined yet
  I_trigger <- NA
  ### define as 1/2Bmsy if requested
  if (isTRUE(OM_I_trigger == "B_trigger")) {
    I_trigger <- c(refpts(OM_i$brp)["msy", "ssb"])/2
  }
  
  ### add noise
  #I_lim <- I_lim * rlnorm(n = length(I_lim), sdlog = 0.05)
  #I_trigger <- I_trigger * rlnorm(n = length(I_trigger), sdlog = 0.05)
  
  ### add reference points to stk
  attr(x = OM_i$stk, which = "refpts") <- FLPar(I_lim = I_lim,
                                            I_trigger = I_trigger)
  
  ### ---------------------------------------------------------------------- ###
  ### catch length frequencies
  ### ---------------------------------------------------------------------- ###
  
  ### create FLStockLen with catch length frequencies
  stk_len <- length_freq(OM_i$stk, full_series = TRUE,
                        lhpar = attr(OM_i$stk, "lhpar"),
                        len_noise_sd = 0.0, len_sd = 1,
                        len_sd_cut = 2)
  ### save catch.n as attribute in stk
  attr(OM_i$stk, "catch_len") <- window(catch.n(stk_len), end = OM_iy)
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
  names(OM_i)[names(OM_i) == "srbh"] <- "sr.om"
  ### recruitment residuals
  names(OM_i)[names(OM_i) == "srbh.res"] <- "sr.om.res"
  
  ### ---------------------------------------------------------------------- ###
  ### create reference points with perfect knowledge
  ### ---------------------------------------------------------------------- ###
  
  ### create stock without any uncertainy from brp
  stk_tmp <- as(OM_i$brp, "FLStock")
  
  ### target FMSY for 100 years
  ctrl_tmp <- fwdControl(data.frame(year = 2:101, quantity = "f",
                                   val = c(refpts(OM_i$brp)["msy", "harvest"])))
  ### project
  stk_tmp <- fwd(stk_tmp, ctrl = ctrl_tmp, sr = OM_i$sr.om)
  
  ### extract also C/I from here
  ### create perfect index value
  idx_tmp <- c(q_modeled) * (stock.wt(stk_tmp) * stock.n(stk_tmp))[, ac(100)]
  ### ratio catch/index
  I_F_proxy <- catch(stk_tmp)[, ac(100)] / quantSums(idx_tmp)
  
  ### create catch length frequencies, without uncertainty
  len_tmp <- length_freq(stk_tmp[, ac(100)], full_series = TRUE, 
                        lhpar = FLCore::iter(attr(OM_i$stk, "lhpar"), 1),
                        len_noise_sd = 0.0, len_sd = 1,
                        len_sd_cut = 2)
  LFeFmsy <- calc_mean(catch.n(len_tmp), min = calc_Lc(catch.n(len_tmp)), 
                      include_min = FALSE)
  
  ### targt M
  ctrl_tmp <- fwdControl(data.frame(year = 2:101, quantity = "f",
                                   val = c(OM_i$stk@lhpar["M", 1])))
  ### project
  stk_tmp <- fwd(stk_tmp, ctrl = ctrl_tmp, sr = OM_i$sr.om)
  ### create catch length frequencies, without uncertainty, no spreading
  len_tmp <- length_freq(stk_tmp[, ac(100)], full_series = TRUE, 
                        lhpar = FLCore::iter(attr(OM_i$stk, "lhpar"), 1),
                        len_noise_sd = 0.0, len_sd = 1,
                        len_sd_cut = 2)
  LFeM <- calc_mean(catch.n(len_tmp), min = calc_Lc(catch.n(len_tmp)), 
                   include_min = FALSE)
  
  ### save reference points in refpts attribute
  attr(OM_i$stk, "refpts") <- rbind2(attr(OM_i$stk, "refpts"), 
  FLPar(LFeFmsy = c(LFeFmsy), LFeM = c(LFeM),
        F_MSY = refpts(OM_i$brp)["msy", "harvest"],
        F0.1 = refpts(OM_i$brp)["f0.1", "harvest"],
        I_F_proxy = I_F_proxy,
        iter = dims(OM_i$stk@lhpar)$iter))
  
  ### ---------------------------------------------------------------------- ###
  ### calculate F0.1 with YPR
  ### ---------------------------------------------------------------------- ###

  ### needed for catch rule 3.2.1 from WKMSYCat34, option c for factor f

  ### run YPR for range of F's
  F_list <- seq(0, 3, 0.01)

  ### inverse von Bertalanffy growth function
  inv_vB <- function(L, L_inf, K, t0){
    res <- -log(1 - (L / L_inf)) / K + t0
    return(res)
  }
   
  ### length at first capture
  L_c <- yearMeans(calc_Lc(attr(OM_i$stk, "catch_len")))
  ### convert into age  ### calculate t_c from L_c
  t_c <- inv_vB(L = L_c,
                L_inf = attr(OM_i$stk, "lhpar")["L_inf"],
                K = attr(OM_i$stk, "lhpar")["K"], 
                t0 = attr(OM_i$stk, "lhpar")["t0"])
  t_c <- round(t_c)
  
  ### calculate F0.1 for 1 iteration
  F0.1 <- lapply((1:dims(OM_i$stk)$iter)[1], function(y){

    ### conduct YPR for list of F's
    res_list <- lapply(F_list, function(f){

      YPR(a =	attr(OM_i$stk, "lhpar")["a", y], 
          b =	attr(OM_i$stk, "lhpar")["b", y],
          Linf =	attr(OM_i$stk, "lhpar")["L_inf", y],
          K =	attr(OM_i$stk, "lhpar")["K", y], 
          t0 =	attr(OM_i$stk, "lhpar")["t0", y],
          tc =	c(FLCore::iter(t_c, y)),
          M =	attr(OM_i$stk, "lhpar")["M", y], F =	f,
          max_age = attr(OM_i$stk, "lhpar")["max_age", y])

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
  attr(OM_i$stk, "refpts") <- rbind2(attr(OM_i$stk, "refpts"),
                                     FLPar(F0.1YPR = F0.1, 
                                           iter = dims(OM_i$stk@lhpar)$iter))

  ### ---------------------------------------------------------------------- ###
  ### stock numbers at begin of simulation ####
  ### ---------------------------------------------------------------------- ###
  ### needed if survey at beginning of advice year used
  
  ### target F value from last year
  ctrl <- getCtrl(values = fbar(OM_i$stk)[, ac(100)], quantity = "f", 
                  years = 100, it = dim(OM_i$stk)[6])
  
  ### forecast
  stk_fwd <- fwd(OM_i$stk, ctrl = ctrl, sr = OM_i$sr.om, 
                 sr.residuals = OM_i$sr.om.res,
                 sr.residuals.mult = TRUE, maxF = 5)[]
  
  ### insert stock numbers in year 101
  ### these numbers are the numbers at the begin of 101, i.e. based on year 100
  stock.n(OM_i$stk)[, ac(101)] <- stock.n(stk_fwd)[, ac(101)]
  
  
  ### ---------------------------------------------------------------------- ###
  ### save
  ### ---------------------------------------------------------------------- ###
    
  ### create local environment to change names
  local({
    
    stk <- OM_i$stk
    observations <- OM_i$observations
    sr.om <- OM_i$sr.om
    sr.om.res <- OM_i$sr.om.res
    it <- OM_iter
    fy <- dims(OM_i$stk)$maxyear
    y0 <- 75
    dy <- 100
    iy <- 100
    ny <- OM_n_years + 1
    nsqy <- nsqy
    vy <- ac(iy:fy)
    
    dir.create(paste0("input/OM2/perfect_knowledge/", OM_scn$id, "/", fhist_i, 
                      "/"), 
               recursive = TRUE)
    save(stk, observations, sr.om, sr.om.res, it, fy, y0, dy, iy, ny, nsqy, vy,
         file = paste0("input/OM2/perfect_knowledge/", OM_scn$id, "/", fhist_i, 
                       "/", lh_i$stock, ".RData"))
    saveRDS(list(stk = stk, observations = observations, sr.om = sr.om, 
                 sr.om.res = sr.om.res, it = it, fy = fy, y0 = y0, dy = dy, 
                 iy = iy, ny = ny, nsqy = nsqy, vy = vy),
            file = paste0("input/OM2/perfect_knowledge/", OM_scn$id, "/", 
                          fhist_i, "/", lh_i$stock, ".rds"))
  })
  
  #return(OM_i)
  
# }

### ------------------------------------------------------------------------ ###
### observation error ####
### ------------------------------------------------------------------------ ###

# OM_list <- foreach(OM_scn_i = OM_list[OM_scns$obs_error],
#                    OM_iter = OM_scns$n_iter[OM_scns$obs_error]) %:%
#   foreach(OM_i = OM_scn_i, stk_pos = seq_along(OM_scn_i),
#           .packages = required_pckgs, .export = c("nsqy")) %dopar% {
  
  ### set seed
  set.seed(1)
  
  ### ---------------------------------------------------------------------- ###
  ### add index uncertainty
  ### ---------------------------------------------------------------------- ###
  
  ### get index
  idx <- OM_i$observations$idx$idx
  
  ### add uncertainty to catchability
  ### log-normal noise, cv = 0.2
  q_error <- index(idx) %=% 1
  q_error[] <- rlnorm(n = length(index.q(idx)), sdlog = OM_scn$idxSD)
  
  ### add index uncertainty to each age independently (default)
  if (!isTRUE(OM_scn$idxSD_age == FALSE)) {
    
    ### insert age specific noise
    index.q(idx) <- index.q(idx) * q_error
  
  ### otherwise, add uncertainty to biomass at age
  } else {
    
    ### use only noise from first age and apply to all other ages
    q_error_bio <- q_error
    q_error_bio[] <- q_error[1, ]
    index.q(idx) <- index.q(idx) * q_error_bio
    
  }
  
  ### update index values
  index(idx) <- index.q(idx) * stock.n(OM_i$stk) * stock.wt(OM_i$stk)
  
  ### remove index value from projection period, if they exist
  index(idx)[, ac(101)] <- NA
  
  ### save in observations object
  OM_i$observations <- list(idx = FLIndices(idx = idx))
  
  ### ---------------------------------------------------------------------- ###
  ### update index reference points
  ### ---------------------------------------------------------------------- ###
  ### uncertainty already implemented in index values
  ### the new minimum values already include this uncertainty
  
  attr(OM_i$stk, "refpts")["I_lim"] <- apply(quantSums(index(idx)), 6, min,
                                             na.rm = TRUE)
  
  ### ---------------------------------------------------------------------- ###
  ### add uncertainty to life-history parameters and reference points
  ### ---------------------------------------------------------------------- ###
  
  ### reference points
  attr(OM_i$stk, "refpts")["LFeFmsy"] <- attr(OM_i$stk, "refpts")["LFeFmsy"] * 
    rlnorm(n = length(attr(OM_i$stk, "refpts")["LFeFmsy"]), sdlog = 0.1)
  attr(OM_i$stk, "refpts")["LFeM"] <- attr(OM_i$stk, "refpts")["LFeM"] * 
    rlnorm(n = length(attr(OM_i$stk, "refpts")["LFeM"]), sdlog = 0.1)
  attr(OM_i$stk, "refpts")["F_MSY"] <- attr(OM_i$stk, "refpts")["F_MSY"] * 
    rlnorm(n = length(attr(OM_i$stk, "refpts")["F_MSY"]), sdlog = 0.1)
  attr(OM_i$stk, "refpts")["F0.1"] <- attr(OM_i$stk, "refpts")["F0.1"] * 
    rlnorm(n = length(attr(OM_i$stk, "refpts")["F0.1"]), sdlog = 0.1)
  attr(OM_i$stk, "refpts")["I_F_proxy"] <- 
    attr(OM_i$stk, "refpts")["I_F_proxy"] * 
    rlnorm(n = length(attr(OM_i$stk, "refpts")["I_F_proxy"]), sdlog = 0.1)
  
  ### lhist
  attr(OM_i$stk, "lhpar")["L_inf"] <- attr(OM_i$stk, "lhpar")["L_inf"] * 
    rlnorm(n = length(attr(OM_i$stk, "lhpar")["L_inf"]), sdlog = 0.1)
  attr(OM_i$stk, "lhpar")["K"] <- attr(OM_i$stk, "lhpar")["K"] * 
    rlnorm(n = length(attr(OM_i$stk, "lhpar")["K"]), sdlog = 0.1)
  attr(OM_i$stk, "lhpar")["t0"] <- attr(OM_i$stk, "lhpar")["t0"] * 
    rlnorm(n = length(attr(OM_i$stk, "lhpar")["t0"]), sdlog = 0.1)
  # attr(OM_i$stk, "lhpar")["a"] <- attr(OM_i$stk, "lhpar")["a"] * 
  #   rlnorm(n = length(attr(OM_i$stk, "lhpar")["a"]), sdlog = 0.1)
  # attr(OM_i$stk, "lhpar")["b"] <- attr(OM_i$stk, "lhpar")["b"] * 
  #   rlnorm(n = length(attr(OM_i$stk, "lhpar")["b"]), sdlog = 0.1)
  attr(OM_i$stk, "lhpar")["M"] <- attr(OM_i$stk, "lhpar")["M"] * 
    rlnorm(n = length(attr(OM_i$stk, "lhpar")["M"]), sdlog = 0.1)
  
  ### ---------------------------------------------------------------------- ###
  ### catch length frequency uncertainty
  ### ---------------------------------------------------------------------- ###
  
  ### create FLStockLen with catch length frequencies
  stk_tmp <- OM_i$stk
  attr(stk_tmp, "catch_len") <- NULL
  stk_len <- length_freq(stk_tmp, full_series = TRUE,
                         lhpar = attr(OM_i$stk, "lhpar"),
                         len_noise_sd = OM_scn$lengthSD, len_sd = 1,
                         len_sd_cut = 2)
  ### save catch.n as attribute in stk
  attr(OM_i$stk, "catch_len") <- window(catch.n(stk_len), end = 100)

  ### --------------------------------------------------------------------   ###
  ### save
  ### ---------------------------------------------------------------------- ###
  
  ### create local environment to change object names
  local({
  
    stk <- OM_i$stk
    observations <- OM_i$observations
    sr.om <- OM_i$sr.om
    sr.om.res <- OM_i$sr.om.res
    it <- OM_iter
    fy <- dims(OM_i$stk)$maxyear
    y0 <- 75
    dy <- 100
    iy <- 100
    ny <- OM_n_years + 1
    nsqy <- nsqy
    vy <- ac(iy:fy)
    
    dir.create(paste0("input/OM2/observation_error/", OM_scn$id, "/", fhist_i, 
                      "/"), 
               recursive = TRUE)
    save(stk, observations, sr.om, sr.om.res, it, fy, y0, dy, iy, ny, nsqy, vy,
         file = paste0("input/OM2/observation_error/", OM_scn$id, "/", fhist_i, 
                       "/", lh_i$stock, ".RData"))
    saveRDS(list(stk = stk, observations = observations, sr.om = sr.om, 
                 sr.om.res = sr.om.res, it = it, fy = fy, y0 = y0, dy = dy, 
                 iy = iy, ny = ny, nsqy = nsqy, vy = vy),
            file = paste0("input/OM2/observation_error/", OM_scn$id, "/", 
                          fhist_i, "/", lh_i$stock, ".rds"))

  })
  
  return(NULL)
  
}


stopCluster(cl)

#quit(save = "no")
