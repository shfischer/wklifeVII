# SCRIPT TO CREATE OPERATING MODELS
# Karin Olsson
# 2017-08-30

# PACKAGES

library(FLife)
library(FLBRP)
library(FLasher)
library(foreach)
library(doParallel)
#library(FLa4a)

source("../R/oms.R")
source("../R/functions.R")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## LOADING WKLIFE DATA
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data(wklife)

# Change Nephrops data from (mature) female to mature male 
# (from Stock Annex for assesment for Norway Lobster in 9a))
wklife$sex <- as.character(wklife$sex)
wklife$sex[15] <- "M"; wklife$a[15] <- 0.00028; wklife$b[15] <- 3.229; 
wklife$linf[15] <- 70; wklife$l50[15] <- 28.4; wklife$k[15] <- 0.2

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## CREATE OPERATING MODELS
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

wkl <- rbind(cbind(wklife, scen="oneway", om=paste0(wklife$stock, "_oneway")),
             cbind(wklife, scen="rollercoaster", om=paste0(wklife$stock, "_rollercoaster")))
### force steepness of 0.75 for all stocks
wkl$s <- 0.75

# Values for minfbar, maxfbar estimated from inspecting plots of catch.n at stable F=0.5*Fmsy 
fbar_vals <- data.frame(stock=wklife$stock,
                        minfbar=c(1,1,2,9,1,1,1,1,1,1,1,1,1,1,1),
                        maxfbar=c(3,6,15,23,5,3,4,5,2,5,4,5,9,8,4))
fbar_vals <- rbind(fbar_vals, fbar_vals)

# Iterations
iters <- 500
  
registerDoParallel(28)

oms <- foreach (i=1:nrow(wkl)) %dopar% {
  
  set.seed(0)
  
  # Create the FLStock, FLBRP, FLSR
  om <- OM_stockmodel(data=wkl[i,], fbar_vals[i,], its=iters)
  
  # Create OM
  om$stk <- OM_scenario(stk=om$stk, sr=om$sr, brp=om$brp, resid.sd=0.2, resid.rho=0.3,
                      up=0.2, down=0.2, its=iters, scenario=wkl[i,"scen"])
  
  # Create plot of OM
  p <- plot(om$stk[, 75:100]) + ggtitle(sub("_", " ", wkl[i, "om"]))
  
  # ggsave(filename = paste0(wkl[i, "om"],".png"),
  #        plot = p, scale = 1, width = 29, height = 21, units = "cm", dpi = 300, 
  #        limitsize = TRUE)
  
  # Limit OM to last 25 years
  om$stk_full <- om$stk
  om$stk <- om$stk[, ac(75:100)]
  
  om

}

names(oms) <- wkl$om
saveRDS(oms, file="oms.rds", compress="xz")

#
# stks <- lapply(oms, function(x) x$stk[, 75:100])
# saveRDS(stks, file="stks.rds", compress="xz")
#

# plots <- lapply(oms, '[[', 'plot')
# saveRDS(plots, file="plots.rds", compress="xz")

