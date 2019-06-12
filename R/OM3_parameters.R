### ------------------------------------------------------------------------ ###
### derive more parameters for stocks ####
### ------------------------------------------------------------------------ ###

library(FLife)
library(plyr)
library(popbio)
library(dplyr)
library(tidyr)
library(ggplot2)
library(hclust)

### ------------------------------------------------------------------------ ###
### modify lhist as used in simulation ####
### ------------------------------------------------------------------------ ###
stocks_lh <- read.csv("input/stock_list_full.csv")
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
if (nrow(stocks_lh) > 15) {
  stocks_lh$stock[29] <- paste0(stocks_lh$stock[29], "_2")
}

### set default t0 to -0.1
stocks_lh$t0[is.na(stocks_lh$t0)] <- -0.1

### set fbar range
stocks_lh$minfbar <- c(1,1,2,9,1,1,1,1,1,1,1,1,1,1,1,3,4,2,2,2,1,1,1,1,1,1,3,1,
                       4)[1:nrow(stocks_lh)]
stocks_lh$maxfbar <- c(3,6,15,23,5,3,4,5,2,5,4,5,9,8,4,13,11,11,15,13,3,3,5,5,
                       10,4,11,4,20)[1:nrow(stocks_lh)]

### ------------------------------------------------------------------------ ###
### create BRPs and estimate parameters ####
### ------------------------------------------------------------------------ ###

stocks_lh_names <- split(stocks_lh, seq(nrow(stocks_lh)))
names(stocks_lh_names) <- stocks_lh$stock


growth <- vonB
res <- lapply(stocks_lh_names, function(object) {
  #browser()
  pars <- lhPar(object)

  ### max age
  max_age <- ceiling(log(0.05)/(-c(pars["k"])) + c(pars["t0"]))
    
  ### create BRP
  brp <- (lhEql(pars, range = c(min = 1, max = max_age, 
                               minfbar = pars["minfbar"],
                               maxfbar = pars["maxfbar"], 
                               plusgroup = max_age)))
  
  ### stock recruitment parameters
  srr_pars <- params(brp)
  dimnames(srr_pars)$params <- c("alpha", "beta")
  ### add to lhist
  pars <- rbind(pars, srr_pars)
  
  ### SPR0
  pars <- rbind(pars, FLPar(spr0 = c(spr0(brp))))
  
  ### MSY reference points: F, biomass and yield
  rfs <- FLPar(refpts(brp)["msy", c("harvest", "yield", "ssb"), drop = TRUE])
  dimnames(rfs)$params <- c("fmsy", "msy", "bmsy")
  ### add
  pars <- rbind(pars, rfs)
  
  ### Lopt (length at MSY)
  growth <- vonB
  lopt <- lopt(pars)
  pars <- rbind(pars,lopt)
  
  ### length at first capture
  ### calculation here does not make sense
  ### is length at a50 - ato95 and ato95 is always default=1
  # lc <- vonB(as.FLQuant(c(pars["a50"] - pars["ato95"]),
  #                       dimnames = list(iter = dimnames(pars)$iter)), pars)
  # lc <- FLPar(lc = (lc))
  # pars <- rbind(pars, lc)
  
  ### growth rate r
  r <- maply(seq(dims(brp)$iter), function(x) {
    log(lambda(leslie(iter(brp, x), 
                  fbar = c(refpts(brp)["crash", "harvest", x]))[drop = TRUE]))
  })
  r <- FLPar(r = array(c(r), c(1, length(c(r)))))
  pars <- rbind(pars, r)
  
  ### growth rate at MSY
  rc <- maply(seq(dims(brp)$iter), function(x) {
    log(lambda(leslie(iter(brp, x), 
                      fbar = c(refpts(brp)["msy", "harvest", x]))[drop = TRUE]))
  })
  rc <- FLPar(rc = array(c(rc), c(1, length(c(r)))))
  pars <- rbind(pars, rc)
  
  ### natural mortality, mature population
  m <- FLPar(m = weighted.mean(x = m(brp), w = mat(brp)))
  pars <- rbind(pars, m)
  
  ### M/k
  Mk <- FLPar(mk = pars["m"] / pars["k"])
  pars <- rbind(pars, Mk)
  
  ### Fmsy/m
  fmsym <- FLPar(fmsym = pars["fmsy"] / pars["m"])
  pars <- rbind(pars, fmsym)
  
  ### Bmsy/B0
  bmsyb0 <- FLPar(bmsyb0 = pars["bmsy"] / pars["v"])
  pars <- rbind(pars, bmsyb0)  
  
  return(pars)
  
})

res2 <- lapply(res, function(x) {
  data.frame(t(data.frame(x)))
})
res2 <- rbind.fill(res2)
res2$stock <- names(res)
rownames(res2) <- names(res)

### save
saveRDS(object = res2, file = "input/lhist_extended.rds")
write.csv(x = res2, file = "input/lhist_extended.csv", row.names = FALSE)


### ------------------------------------------------------------------------ ###
### correlations between parameters ####
### ------------------------------------------------------------------------ ###

### select parameters
res3 <- res2[, c("linf", "k", "t0", "a", "b", "a50", "s", "v", "l50", "minfbar", "maxfbar", "alpha", "beta", "spr0", "fmsy", "msy", "bmsy", "lopt", "r", "rc", "m", "mk", "fmsym", "bmsyb0", "lmax")]
### remove some (cause NAs)
res3 <- res3[, !names(res3) %in% c("s", "v", "l50", "lmax")]

### create correlation matrix
res_cor <- cor(res3)
### order with hclust
dd <- as.dist((1 - res_cor)/2)
hc <- hclust(dd)
res_cor <- res_cor[hc$order, hc$order]
### keep only upper triangle
res_cor[lower.tri(res_cor)] <- NA

### transform for nicer plotting
data_cor <- as.data.frame(res_cor)
data_cor$par <- rownames(data_cor)
data_cor <- data_cor %>% gather(key, value, -par)
data_cor$par <- factor(data_cor$par, levels = dimnames(res_cor)[[1]])
data_cor$key <- factor(data_cor$key, levels = dimnames(res_cor)[[1]])
data_cor <- data_cor[!is.na(data_cor$value), ]

data_cor %>% ggplot(aes(x = key, y = par, fill = value, 
                        label = round(value, 2))) +
  geom_tile(colour = "white") +
  geom_text(size = 2) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name = "Pearson\nCorrelation") +
  coord_fixed() +
  #theme_minimal()# + 
  theme_bw() +
  labs(x = "", y = "") +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(filename = "input/lhist_extended_correlations.png",
       width = 20, height = 15, units = "cm", dpi = 300, type = "cairo")
