### ------------------------------------------------------------------------ ###
### create biological stocks from life-history parameters ####
### as used for WKLIFE VII 2017
### ------------------------------------------------------------------------ ###
### author: Simon Fischer (Cefas), simon.fischer@cefas.co.uk
### FRAMEWORK based on the a4a standard MSE developed at JRC
### by Ernesto Jardim, Iago Mosqueira, Finlay Scott, et al.
### additional contributors:
### Karin Olsson
### ------------------------------------------------------------------------ ###
### created 08/2017
### last modifications:
### 2017-12 Simon Fischer
### ------------------------------------------------------------------------ ###


### load packages
library(FLife)
library(FLBRP)
library(FLasher)
library(foreach)
library(doParallel)

### set up cluster for parallel computing
cl <- makeCluster(28)
registerDoParallel(cl)

### load additional functions
### select R scripts from functions folder
load_files <- list.files("functions/")
load_files <- load_files[grepl(pattern = "*.R$", x = load_files)]
### source the scripts
invisible(lapply(paste0("functions/", load_files), source))

### ------------------------------------------------------------------------ ###
### load data ####
### ------------------------------------------------------------------------ ###

### extended list
stocks_lh <- read.csv("input/stock_list_full.csv")
### used at WKLIFE VII 2017:
#stocks_lh <- read.csv("input/wklife.csv")[, -1]
### subset to non-NA stocks
stocks_lh <- stocks_lh[!is.na(stocks_lh$a), ]
names(stocks_lh)[1] <- "name"

### use lmax as proxy for linf, if linf not provided but lmax exits
pos_lmax <- which(!is.na(stocks_lh$lmax) & is.na(stocks_lh$linf))
stocks_lh$linf[pos_lmax] <- stocks_lh$lmax[pos_lmax]

### set steepness to 0.75
stocks_lh$s <- 0.75

### change name for stocks that appear twice
stocks_lh$stock <- as.character(stocks_lh$stock)
if(nrow(stocks_lh) > 15){
  stocks_lh$stock[29] <- paste0(stocks_lh$stock[29], "_2")
}

### ------------------------------------------------------------------------ ###
### create OMs to find fbar range ####
### ------------------------------------------------------------------------ ###

# stks <- foreach(i = split(stocks_lh, 1:nrow(stocks_lh)),
#                 .errorhandling = "pass", 
#                 .packages = c("FLife", "FLasher", "FLBRP")) %dopar% {
#   
#   ### create brp
#   ### get lh params
#   lh_res <- c(dimnames(lhPar(FLPar(linf = 1)))$params, "l50")
#   lh_avail <- intersect(lh_res, names(i))
#   lh_pars <- i[, lh_avail]
#   lh_pars <- lh_pars[, !is.na(lh_pars)]
#   lh_pars <- as(lh_pars, "FLPar")
#   ### create missing pars
#   lh_pars <- lhPar(lh_pars)
#   
#   # Max age: age at l = 0.95 * linf
#   max_age <- ceiling(log(0.05)/(-c(lh_pars["k"]))+c(lh_pars["t0"]))
#   
#   ### create brp
#   brp <- lhEql(lh_pars, range = c(min = 1, max = max_age, 
#                                   minfbar = 1, maxfbar = max_age, 
#                                   plusgroup = max_age))
#   
#   ### coerce into FLStock
#   stk <- as(brp, "FLStock")
#   
#   ### create SRR
#   srr <- FLSR(params = params(brp), model = model(brp))
#   
#   ### harvest at Fmsy for 100 years
#   target <- fwdControl(data.frame(year = 3:101, 
#                                   value = c(refpts(brp)["msy", "harvest"]),
#                                   quant = "f"))
#   ### project
#   stk_fwd <- fwd(stk, control = target, sr = srr)
#   
#   return(stk_fwd)
# 
# }
# 
# #for(i in 1:length(stks)){
# (i <- i + 1) 
#   ggplot(data = as.data.frame((catch.n(stks[[i]]))[, ac(3:100)]),
#        aes(x = age, y = data, colour = as.factor(year))) +
#   geom_line() + geom_point()
#   ggplot(data = as.data.frame((catch.n(stks[[i]]) * catch.wt(stks[[i]]))[, ac(3:100)]),
#          aes(x = age, y = data, colour = as.factor(year))) +
#     geom_line() + geom_point()
# 
#   (catch <- (catch.n(stks[[i]]))[, ac(101)])
#   (catch_rel <- catch / sum(catch))
#   sum(catch_rel[4:20])
# #}
# 
# res <- lapply(1:29, function(x){
#   (catch <- (catch.n(stks[[x]]))[, ac(101)])
#   (catch_rel <- catch / sum(catch))
#   (catch2 <- (catch.n(stks[[x]]) * catch.wt(stks[[x]]))[, ac(101)])
#   (catch_rel2 <- catch2 / sum(catch2))
#   data.frame(n = sum(catch_rel[fbar_range[x,1]:fbar_range[x,2]]), 
#              b = sum(catch_rel2[fbar_range[x,1]:fbar_range[x,2]]))
# })
# res <- do.call(rbind, res)
# res

### ------------------------------------------------------------------------ ###
### create final brps ####
### ------------------------------------------------------------------------ ###

### set fbar range
stocks_lh$minfbar <- c(1,1,2,9,1,1,1,1,1,1,1,1,1,1,1,3,4,2,2,2,1,1,1,1,1,1,3,1,
                       4)[1:nrow(stocks_lh)]
stocks_lh$maxfbar <- c(3,6,15,23,5,3,4,5,2,5,4,5,9,8,4,13,11,11,15,13,3,3,5,5,
                       10,4,11,4,20)[1:nrow(stocks_lh)]

### create FLBRP objects from life-history data
brps <- foreach(i = split(stocks_lh, 1:nrow(stocks_lh)),
                .errorhandling = "pass", 
                .packages = c("FLife", "FLBRP")) %dopar% {
  
  ### create brp
  ### get lh params
  lh_res <- c(dimnames(lhPar(FLPar(linf = 1)))$params, "l50")
  lh_avail <- intersect(lh_res, names(i))
  lh_pars <- i[, lh_avail]
  lh_pars <- lh_pars[, !is.na(lh_pars)]
  lh_pars <- as(lh_pars, "FLPar")
  ### create missing pars
  lh_pars <- lhPar(lh_pars)
  
  # Max age: age at l = 0.95 * linf
  max_age <- ceiling(log(0.05)/(-c(lh_pars["k"]))+c(lh_pars["t0"]))
  
  ### create brp
  brp <- lhEql(lh_pars, range = c(min = 1, max = max_age, 
                                  minfbar = i$minfbar, maxfbar = i$maxfbar, 
                                  plusgroup = max_age))
  
  ### save life-history parameters in FLBRP
  attr(brp, "lhpar") <- lh_pars
  
  return(brp)

}
names(brps) <- stocks_lh$stock

### save brps
saveRDS(brps, file = "input/brps.rds")

### calculate M
stocks_lh$M <- unlist(lapply(brps, function(x){
  mean(m(x))
}))
### mature M (only mature proportion of stock considered)
stocks_lh$M_mat <- unlist(lapply(brps, function(x){
  weighted.mean(x = m(x), w = mat(x))
}))

### M/K
stocks_lh$MK <- with(stocks_lh, M_mat / k)

### get reference points
ref_pts <- lapply(brps, refpts)
names(ref_pts) <- stocks_lh$stock
saveRDS(ref_pts, "input/refpts.rds")

### save stock list
write.csv(file = "input/stock_list_full2.csv", x = stocks_lh)

### ------------------------------------------------------------------------ ###
### create FLStocks ####
### ------------------------------------------------------------------------ ###

### number of iterations
its <- 500

OMs <- foreach(i = brps, .errorhandling = "pass", 
                .packages = c("FLife", "FLasher", "FLBRP")) %dopar% {
  
  #### coerce FLBRP into FLStock
  ### keep only second year (first year with non zero catch)
  stk <- as(i, "FLStock")[, 2]
  ### name first year "1"
  stk <- qapply(stk, function(x) {
    dimnames(x)$year <- "1"; return(x)
  })
  
  ### extend object to year 100
  stk <- fwdWindow(stk, i, end = 100)
  
  ### propagate with requested number of iterations
  stk <- propagate(stk, its)
  
  ### create FLSR object
  stk_sr <- FLSR(params = params(i), model = model(i))
  
  ### create residuals for (historical) projection
  set.seed(0)
  residuals <- rlnoise(its, rec(stk) %=% 0, sd = 0.2, b = 0.3)
  
  ### project forward, targeting 0.5 F_MSY
  years_target <- ((dims(stk)$minyear + 1):75)
  stk <- fwd(stk, sr = stk_sr, 
             control = fwdControl(year = years_target, 
                                  value = refpts(i)['msy', 'harvest']*0.5, 
                                  quant = "f"),
             residuals = residuals[, ac(years_target)])
  
  ### project last 25 years with 2 scenarios: one-way & roller-coaster
  years_target <- 76:100
  ### one-way
  stk_one_way <- oneWayTrip(stk, sr = stk_sr, brp = i, years = years_target,
                            residuals = residuals[, ac(years_target)],
                            f0 = refpts(i)["msy", "harvest"]*0.5)
  
  ### roller-coaster
  stk_roller_coaster <- rollerCoaster(stk, sr = stk_sr, brp = i, 
                                      up = 0.2, down = 0.2, 
                                      years = years_target,
                                      residuals = residuals[, ac(years_target)],
                                      f0 = refpts(i)["msy", "harvest"]*0.5)
  
  ### extract life-history parameters
  lhpar <- attr(i, "lhpar")
  attr(i, "lhpar") <- NULL
  
  ### return list
  return(list(sr = stk_sr, brp = i, lhpar = lhpar, 
              stk_one_way = stk_one_way, 
              stk_roller_coaster = stk_roller_coaster))
  
}

### name the FLStock objects
OMs <- lapply(seq_along(OMs), function(x){
  name(OMs[[x]]$stk_roller_coaster) <- ac(stocks_lh$stock[x])
  name(OMs[[x]]$stk_one_way) <- ac(stocks_lh$stock[x])
  return(OMs[[x]])
})

### set names for list elements
names(OMs) <- stocks_lh$stock

### save list
saveRDS(OMs, "input/stock_list_full.rds")

### split into individual stocks
### and subset to last 26 years
OMs_one_way <- lapply(OMs, function(x){
  x$stk_roller_coaster <- NULL
  names(x)[length(names(x))] <- "stk"
  x$stk <- window(x$stk, start = 75)
  return(x)
})
OMs_roller_coaster <- lapply(OMs, function(x){
  x$stk_one_way <- NULL
  names(x)[length(names(x))] <- "stk"
  x$stk <- window(x$stk, start = 75)
  return(x)
})
### combine into single list
### sort, to get order used in WKLIFE VII
OMs <- c(OMs_one_way[1:15], OMs_roller_coaster[1:15],
         OMs_one_way[-c(1:15)], OMs_roller_coaster[-c(1:15)])

### save list
saveRDS(OMs, "input/stock_list.rds")
# 
# stk1 <- (function(){
#   load("input/stocks/perfect_knowledge/1.RData")
#   return(stk)
# })()

stopCluster(cl)

