### ------------------------------------------------------------------------ ###
### MSE for ICES WKLIFE VII 2017 ####
### load OMs and set additional parameters for simulation
### ------------------------------------------------------------------------ ###


rm(list=ls())
### load packages
required_pckgs <- c("FLash", "FLAssess", "ggplotFL", "FLBRP")
### save as object in order to avoid output to screen
. <- lapply(required_pckgs, function(x){
  suppressMessages(library(x, character.only = TRUE))
})
library(doParallel)

### load additional functions
source("functions/funs.R")

set.seed(0)

#==============================================================================
# Load the stock and files
#==============================================================================

### load list with OMs
OM_list <- readRDS("../../wklifeFL/exec/oms.rds")

### ------------------------------------------------------------------------ ###
### specify dimensions of simulation
### ------------------------------------------------------------------------ ###

it <- dims(OM_list[[1]][[1]])$iter # iterations
fy <- dims(OM_list[[1]][[1]])$maxyear + 100 # final year
y0 <- range(OM_list[[1]][[1]])[["minyear"]] # initial data year
dy <- range(OM_list[[1]][[1]])[["maxyear"]] # final data year
iy <- dims(OM_list[[1]][[1]])$maxyear # initial year of projection (also intermediate)
ny <- fy - iy + 1 # number of years to project from intial year
nsqy <- 3 # number of years to compute status quo metrics
vy <- ac(iy:fy) # vector of years to be projected


### ------------------------------------------------------------------------ ###
### "loop" through all stocks ####
### ------------------------------------------------------------------------ ###

names <- names(OM_list)
cl <- makeCluster(15)
registerDoParallel(cl)
OM_list <- foreach(x = seq_along(OM_list), .export = ls(),
                   .packages = required_pckgs) %dopar% {
  
  set.seed(0)
  ### ---------------------------------------------------------------------- ###
  ### create residuals for SRR
  ### ---------------------------------------------------------------------- ###
  
  ### create residuals
  ### values from 2nd argument are used as exp(value), using 0 leads to values 
  ### distributed around 1
  srbh.res <- FLife::rlnoise(n = it, FLQuant(0, dimnames = list(year = vy)), 
                             sd = 0.3, b = 0.2)
  OM_list[[x]]$srbh.res <- srbh.res
  
  names(OM_list[[x]])[names(OM_list[[x]]) == "sr"] <- "srbh"

  ### ---------------------------------------------------------------------- ###
  ### add life-history parameters to stock ####
  ### ---------------------------------------------------------------------- ###
  
  ### extract lhpar
  lhpar <- OM_list[[x]][["lhpar"]]
  ### add some more parameters
  ### natural mortality and max age
  add_pars <- FLPar(M = mean(m(OM_list[[x]][["stk"]][, 1])),
                    max_age = range(OM_list[[x]][["stk"]])[["max"]])
  ### add
  lhpar <- rbind2(lhpar, add_pars)
  ### change some names
  dimnames(lhpar)$params[dimnames(lhpar)$params %in% c("linf", "k")] <- 
    c("L_inf", "K")
  ### propagate
  lhpar <- propagate(lhpar, dims(OM_list[[1]][["stk"]])$iter)
  ### save as attribute in stk
  attr(OM_list[[x]][["stk"]], "lhpar") <- lhpar
  
  ### ---------------------------------------------------------------------- ###
  ### get reference points from FLBRP
  ### ---------------------------------------------------------------------- ###
  
  attr(OM_list[[x]]$stk, "refpts") <- refpts(OM_list[[x]]$brp)

  ### ---------------------------------------------------------------------- ###
  ### set up operating model: extend ####
  ### ---------------------------------------------------------------------- ###
  
  OM_list[[x]]$stk <- stf(OM_list[[x]]$stk, fy-dy, nsqy, nsqy)
  
  OM_list[[x]]

}
names(OM_list) <- names

### ------------------------------------------------------------------------ ###
### save
### ------------------------------------------------------------------------ ###
save(OM_list, it, fy, y0, dy, iy, ny, nsqy, vy, 
     file = "om_list.RData")

stopCluster(cl)
