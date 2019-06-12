### ------------------------------------------------------------------------ ###
### calculate stock numbers at beginning of first simulation year ####
### ------------------------------------------------------------------------ ###

required_pckgs <- c("FLash", "FLAssess")
### save as object in order to avoid output to screen
. <- lapply(required_pckgs, function(x){
  suppressMessages(library(x, character.only = TRUE))
})
library(doParallel)

cl <- makeCluster(4)
registerDoParallel(cl)

### load additional functions
### select R scripts from functions folder
load_files <- list.files("functions/")
load_files <- load_files[grepl(pattern = "*.R$", x = load_files)]
### source the scripts
invisible(lapply(paste0("functions/", load_files), source))

set.seed(0)



### ------------------------------------------------------------------------ ###
### find fishing pressure where Lmean = LFeM ####
### ------------------------------------------------------------------------ ###
library(data.table)

brps <- readRDS("input/brps.rds")

res <- foreach(brp = brps,
               .export = c("brps", "length_freq"),
               .packages = c("FLash", "FLAssess", "FLBRP", "foreach",
                             "data.table")) %dopar% {

  ### coerce BRP into FLStock
  stk <- as(brp, "FLStock")[, 2]
  ### name first year "1"
  stk <- qapply(stk, function(x) {
    dimnames(x)$year <- "1"; return(x)
  })
  
  ### extend to 100 years
  stk <- stf(stk, nyears = 99, wts.nyears = 1)
  
  ### recruitment
  stk_sr <- FLSR(params = params(brp), model = model(brp))
  
  ### set range of Fs to test
  fs <- seq(from = 0, to = c(refpts(brp)["crash", "harvest"]), by = 0.001)
  ### add Fmsy
  fs <- c(fs, c(refpts(brp)["msy", "harvest"]))
  
  ### loop through Fs
  stks_fwd <- foreach(Ftrgt = fs) %do% {
    ### ctrl object
    ctrl <- fwdControl(data.frame(year = 2:100, val = Ftrgt, quantity = "f"))
    ### project
    stk_fwd <- fwd(stk, ctrl = ctrl, sr = stk_sr)
    return(stk_fwd)
  }
  names(stks_fwd) <- fs
  
  ### get lhpar
  lhpar <- attr(brp, "lhpar")
  dimnames(lhpar)$params[1] <- "L_inf"
  
  ### calculate length frequencies of catch
  stks_lengths <- lapply(stks_fwd, function(x) {
    catch.n(length_freq(stk = x, lhpar = lhpar, full_series = TRUE))
  })

  ### calculate length at first capture
  L_c <- round(mean(sapply(lapply(stks_lengths, calc_Lc), mean, na.rm = TRUE),
                     na.rm = TRUE))

  ### calculate mean length above length of first capture
  L_mean <- lapply(stks_lengths, calc_mean, min = L_c)
  
  ### calculate LF=M reference length
  MK <- 1.5
  L_FM <- (lhpar["L_inf"] + 2 * MK * L_c) / (1 + 2 * MK)
  
  # blubb <- as.data.frame(FLQuants(L_mean))
  # blubb$qname <- as.numeric(as.character(blubb$qname))
  # ggplot(data = blubb, aes(x = year, y = data, colour = qname, group = qname)) +
  #   geom_line() +
  #   ylim(c(23.75, NA))
  ### lengths reasonable stable at end of 100 year period
  
  ### get last year 
  L_end <- sapply(L_mean, function(x) {
    c(x[, ac(100)])
  })
  
  
  # df <- data.frame(Lmean = L_end, fbar = fs)
  # ggplot(data = df[-1, ], aes(x = fbar, y = Lmean)) +
  #   geom_line() +
  #   geom_hline(yintercept = c(L_FM), colour = "red", show.legend = TRUE) +
  #   geom_vline(xintercept = c(refpts(brp)["msy", "harvest"]), col = "blue")

  ### find F where length is closest to LF=M
  #which.min(abs(L_end - c(L_FM)))
  
  #plot(stks_fwd[[which.min(abs(L_end - c(L_FM)))]])
  
  ### select stock
  stk_L <- stks_fwd[[which.min(abs(L_end - c(L_FM)))]]
  
  ### get ages
  ages <- an(dimnames(stock.n(stk_L))[["age"]])
  
  ### define model for selectivity: logistic function
  ### inflection point of curve = 10% of max age
  q_model <- FLModelSim(model = ~max_q/(1+exp(-steepness*(age - age50))), 
                        params = FLPar(max_q = 1, steepness = 1, 
                                       age50 = max(ages)/10))
  ### model selectivity
  q_modeled <- predict(q_model, age = ages)
  
  ### create perfect index value
  idx_tmp <- c(q_modeled) * tail(stock.wt(stk_L) * stock.n(stk_L))
  ### ratio catch/index
  FproxyMSY_MK <- c(tail(catch(stk_L)) / quantSums(idx_tmp))
  
  ### create index for stock fished at Fmsy
  idx_tmp_MSY <- c(q_modeled) * tail(stock.wt(stks_fwd[[length(stks_fwd)]]) * 
                                       stock.n(stks_fwd[[length(stks_fwd)]]))
  ### ratio
  FproxyMSY <- c(tail(catch(stks_fwd[[length(stks_fwd)]])) / quantSums(idx_tmp_MSY))
  
  # res <- nls(Lmean ~ exp(-fbar + b) + C,
  #            data = df[-1, ], 
  #            start = list(C = 0, b = 0))
  # res <- nls(Lmean ~ a^(-fbar+b)+c,
  #            data = df[-1, ], 
  #            start = list(a = 2, b = 0, c = 0), 
  #            lower = list(a = 0, b = -10, c = 0), 
  #            upper = list(a = 1000, b = 10, c = 1000),
  #            algorithm = "port")
  # 
  # plot((Lmean) ~ fbar, data = df[-1, ], cex = 0.5, pch = 1)
  # lines(y = predict(res), x = df$fbar[-1], col = "red")
  # 
  # 
  # res <- lm(log(Lmean) ~ fbar, data = df[-1, ])
  # lines(exp(predict(res)) ~ df$fbar[-1])
  # 
  # res <- loess(Lmean ~ fbar, data = df[-1, ])
  # plot((Lmean) ~ fbar, data = df[-1, ], cex = 0.5, pch = 1)
  # lines(y = predict(res, data.frame(fbar = seq(0, 0.358, 0.001))), 
  #        x = seq(0, 0.358, 0.001), col = "red")
  
  return(c("FproxyMSY_MK" = FproxyMSY_MK, "FproxyMSY" = FproxyMSY))
  
}
names(res) <- names(brps)

saveRDS(res, file = "input/brps_FproxyMSY.rds")
res <- readRDS("input/brps_FproxyMSY.rds")

### load stock positions
stock_names_pos <- read.csv("input/stock_names_pos.csv")

### go through OMs and add FproxyMSY values
om_files <- list.files("input/stocks/observation_error/", 
                       pattern = "[0-9]+\\.RData")

### create copies of original OMs
file.copy(from = paste0("input/stocks/observation_error/", om_files),
          to = gsub(x = paste0("input/stocks/observation_error/", om_files),
                    pattern = "\\.", replacement = "_old."))

. <- foreach(om_file = om_files,
             stk_pos = as.numeric(sapply(om_files, gsub, pattern = ".RData",
                                         replacement = "")),
             .packages = c("FLCore"), 
             .export = c("stock_names_pos", "res")) %do% {
  ### get stock name
  stk_name <- ac(stock_names_pos$stock[stock_names_pos$stk_pos == stk_pos])
  
  ### load OM
  load(paste0("input/stocks/observation_error/", om_file))
  
  ### get FproxyMSY
  pars_add <- propagate(FLPar(res[[stk_name]]), dims(stk)$iter)
  ### add and save
  attr(stk, "refpts") <- rbind2(attr(stk, "refpts"), pars_add)
  
  
  ### save OM again
  save(stk, observations, sr.om, sr.om.res, it, fy, y0, dy, iy, ny, nsqy, vy,
       file = paste0("input/stocks/observation_error/", om_file))
  
}

### check
# stk_old <- (function(){
#   load("input/stocks/observation_error/31_old.RData")
#   return(stk)
# })()
# stk_new <- (function(){
#   load("input/stocks/observation_error/31.RData")
#   return(stk)
# })()

### remove old OMs
file.remove(list.files("input/stocks/observation_error/", 
                       pattern = "[0-9]+\\_old.RData", full.names = TRUE))
