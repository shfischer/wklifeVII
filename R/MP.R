### ------------------------------------------------------------------------ ###
### run data-limited MSE with the FLR/mse package ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### get arguments passed on to R ####
### ------------------------------------------------------------------------ ###

### args
args <- commandArgs(TRUE)
print("arguments passed on to this script:")
print(args)

### evaluate arguments, if they are passed to R:
if (length(args) > 0) {
  
  ### extract arguments
  for (i in seq_along(args)) eval(parse(text = args[[i]]))
  if (!exists("stock")) stop("stock missing")
  
} else {
  
  stop("no argument passed to R")
  
}

### ------------------------------------------------------------------------ ###
### packages ####
### ------------------------------------------------------------------------ ###
required_pckgs <- c("FLash", "FLAssess", "ggplotFL",
                    "FLBRP", "data.table", "foreach", "mseDL", "doRNG")
### save as object in order to avoid output to screen
. <- lapply(required_pckgs, function(x){
  suppressMessages(library(x, character.only = TRUE))
})

### ------------------------------------------------------------------------ ###
### functions ####
### ------------------------------------------------------------------------ ###

source("MP_functions.R")

### -------------------------------------------------------------------- ###
### load scenario specifications and objects ####
### -------------------------------------------------------------------- ###

### load the objects for stk, recruitment, ...
OM_scns <- read.csv("input/OM_scns.csv", stringsAsFactors = FALSE)
stocks <- read.csv(file = "input/stock_list_full2.csv", as.is = TRUE)
stocks <- expand.grid(stock = stocks$stock, 
                      fhist = c("one-way", "roller-coaster"))

input <- readRDS(paste0("input/OM2/observation_error/",
                        OM_scns$id[id], "/", stocks$fhist[stock], "/",
                        stocks$stock[stock], ".rds"))

### ------------------------------------------------------------------------ ###
### prepare objects for MSE ####
### ------------------------------------------------------------------------ ###

input$sr.om <- propagate(input$sr.om, dims(input$stk)$iter, fill = TRUE)
input$sr.om@residuals <- input$sr.om.res
input$sr.om.res <- NULL

### operating model
om <- FLom(stock = input$stk, ### stock 
           sr = input$sr.om, ### stock recruitment and precompiled residuals
           ### j()
           fleetBehaviour = mseCtrl(method = ctrl_catch_workaround),
           projection = mseCtrl(method = fwd_attr))
tracking = c("Fperc", "convergence", "advice", 
             "Implementation",
             "IEM", "FleetDyn","OM.f", "OM.ssb", "OM.catch",
             "HCR3.2.1c", "HCR3.2.1r", "HCR3.2.1f",
             "HCR3.2.1b", "HCRmult", "L_c", "L_mean",
             "C_current", "I_current",
             "spict_f", "spict_b", "spict_fmsy", "spict_bmsy")
oem <- FLoem(method = obs_bio_len,
             observations = list(stk = input$stk, idx = input$observations$idx), 
             deviances = list(stk = FLQuants(catch_dev = FLQuant()), 
                              idx = FLQuants(idx_dev = FLQuant())),
             args = list(len_noise_sd = 0.2,
                         len_sd = 1,
                         len_sd_cut = 2))
ctrl.mp <- mpCtrl(list(
  ctrl.est = mseCtrl(method = wklife_3.2.1_f,
                     args = list(n_catch_yrs = 1, 
                                 option_f = "a",
                                 option_r = "a",
                                 option_b = "a",
                                 MK = 1.5,
                                 b_z = 1)),
  ctrl.phcr = mseCtrl(method = wklife_3.2.1_x,
                      args = list(n_catch_yrs = 1,
                                  multiplier = 1,
                                  interval = 2,
                                  start_year = 100)),
  ctrl.hcr = mseCtrl(method = wklife_3.2.1_h),
  ctrl.is = mseCtrl(method = wklife_3.2.1_k)
))

iem <- FLiem(method = noise.wrapper,
             args = list(fun = "rlnorm", mean = 0, sd = 0.1,
                         multiplicative = TRUE))

### genArgs
genArgs <- list(fy = 200, ### final simulation year
                y0 = 75, ### first data year
                iy = 100, ### first simulation (intermediate) year
                nsqy = 3, ### not used, but has to provided
                nblocks = 10, ### block for parallel processing
                seed = 1, ### random number seed before starting MSE
                seed_part = TRUE
)

# rm(list = obj_load)

### ------------------------------------------------------------------------ ###
### run MSE ####
### ------------------------------------------------------------------------ ###


### run in parallel
library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)
cl_length <- length(cl)

### load packages and functions into workers
. <- foreach(i = seq(cl_length)) %dopar% {
  #devtools::load_all("../mse/")
  . <- lapply(required_pckgs, function(x){
    suppressMessages(library(x, character.only = TRUE))
  })
  source("MP_functions.R")
}
### run MSE
res <- mpDL(om = om, oem = oem, iem = iem, ctrl.mp = ctrl.mp, genArgs = genArgs,
            scenario = "test", tracking = tracking, verbose = TRUE,
            cut_hist = FALSE)

### save output
path_out <- paste0("output/observation_error/", OM_scns$id[id], "/",
                   stocks$fhist[stock], "/")
file_out <- paste0(path_out, stocks$stock[stock], ".rds")
dir.create(path = path_out, recursive = TRUE)
saveRDS(object = res, file = file_out)

