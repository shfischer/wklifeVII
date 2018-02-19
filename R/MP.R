### ------------------------------------------------------------------------ ###
### MSE for ICES WKLIFE VII 2017
### testing catch rules from ICES WKMSYCat34 2017
### ------------------------------------------------------------------------ ###
### author: Simon Fischer (Cefas), simon.fischer@cefas.co.uk
### based on the a4a standard MSE developed at JRC
### by Ernesto Jardim, Iago Mosqueira, Finlay Scott, et al.
### ------------------------------------------------------------------------ ###
### created 08/2017
### last modifications:
### 2017-10 Simon Fischer
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### load arguments passed to R ####
### ------------------------------------------------------------------------ ###
### load arguments
args <- commandArgs(TRUE)
### extract them
for (i in 1:length(args)) eval(parse(text=args[[i]]))

### set default values
### number of cores, i.e. processes to spawn
if (!isTRUE(exists("n_cores"))){ 
  stop("n_cores need to be passed to R!")
} else {
  n_cores <- n_cores - 1 ### slaves, exluding master
}
### parallelization architecture
if (!isTRUE(exists("cluster_type"))) cluster_type <- 2
### split each scenario into n parts?
if (!isTRUE(exists("n_parts"))) n_parts <- 1
### scenarios to be simulated
if (!isTRUE(exists("scn_start")) | !isTRUE(exists("scn_end"))){
  scns <- TRUE
} else {
  scns <- scn_start:scn_end
}


### ------------------------------------------------------------------------ ###
### set up parallel computing environment ####
### ------------------------------------------------------------------------ ###
### needs to be load first, otherwise the MPI breaks down ...

### for local in-node/PC parallelization
if (cluster_type == 1) {
  
  library(doParallel)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  getDoParWorkers()
  getDoParName()
  
} else if (cluster_type == 2) {
  
  library(doMPI)
  cl <- startMPIcluster(count = n_cores) # one less than requested in the .job file
  registerDoMPI(cl)
  
}

### ------------------------------------------------------------------------ ###
### libraries ####
### ------------------------------------------------------------------------ ###
required_pckgs <- c("FLash", "FLAssess", "FLXSA", "ggplotFL",
                    "FLBRP", "data.table", "spict")
### save as object in order to avoid output to screen
. <- lapply(required_pckgs, function(x){
  suppressMessages(library(x, character.only = TRUE))
})

### ------------------------------------------------------------------------ ###
### functions ####
### ------------------------------------------------------------------------ ###

### select R scripts from functions folder
load_files <- list.files("functions/")
load_files <- load_files[grepl(pattern = "*.R$", x = load_files)]
### source the scripts
invisible(lapply(paste0("functions/", load_files), source))

### ------------------------------------------------------------------------ ###
### load scenario definitions ####
### ------------------------------------------------------------------------ ###
source("MP_scenarios.R")

### ------------------------------------------------------------------------ ###
### start simulations ####
### ------------------------------------------------------------------------ ###
### "loops" through scenarios
### nested foreach loop, splits scenarios into smaller junks, if requested
Sys.time()

### "loop" through scenarios and parts
res <- foreach(scn = seq_along(ctrl.mps)[scns], .packages = required_pckgs,
               .export = ls(), .errorhandling = "pass") %:%
         foreach(part = 1:n_parts, .errorhandling = "pass") %dopar% {
  ### load external functions
  invisible(lapply(paste0("functions/", load_files), source))
  
  ### set seed depending on part of scenario
  set.seed(part)
  
  ### -------------------------------------------------------------------- ###
  ### load scenario specifications and objects ####
  ### -------------------------------------------------------------------- ###
  ### load control object
  ctrl.mp <- ctrl.mps[[scn]]
  
  ### load the objects for stk, recruitment, ...
  load(paste0("input/stocks/", ctrl.mp$scn_desc$uncertainty, "/",
              ctrl.mp$ctrl.om$stk_pos, ".RData"))
  sr.om.res.mult <- ctrl.mp$ctrl.om$sr.res.mult
  
  ### ---------------------------------------------------------------------- ###
  ### split/subset objects, if required ####
  ### ---------------------------------------------------------------------- ###
  if (n_parts > 1) {
    
    ### check if splitting feasible
    if (it%%n_parts != 0) stop("'it' cannot be split into 'n_parts'!")
    ### get desired iterations
    it_part <- split(1:it, cut(1:it, n_parts))[[part]]
    it <- it/n_parts
    ### subset
    stk <- FLCore::iter(stk, it_part)
    if ("lhpar" %in% names(attributes(stk)))
      attr(stk, "lhpar") <- FLCore::iter(attr(stk, "lhpar"), it_part)
    if ("refpts" %in% names(attributes(stk)))
      attr(stk, "refpts") <- FLCore::iter(attr(stk, "refpts"), it_part)
    if ("catch_len" %in% names(attributes(stk)))
      attr(stk, "catch_len") <- FLCore::iter(attr(stk, "catch_len"), it_part)
    observations <- lapply(observations, function(x){
      lapply(x, function(y){
        FLCore::iter(y, it_part)
      })
    })
    #sr.om <- FLCore::iter(sr.om, it_part) ### already 1 iteration
    sr.om.res <- FLCore::iter(sr.om.res, it_part)
    
  }
  
  ### ---------------------------------------------------------------------- ###
  ### create object for tracking ####
  ### ---------------------------------------------------------------------- ###
  tracking <- FLQuant(NA,
    dimnames = list(metric = c("Fperc", "convergence", "advice", "Implementation",
                               "IEM", "FleetDyn","OM.f", "OM.ssb", "OM.catch",
                               "HCR3.2.1c", "HCR3.2.1r", "HCR3.2.1f",
                               "HCR3.2.1b", "HCRmult", "L_c", "L_mean",
                               "C_current", "I_current",
                               "spict_f", "spict_b", "spict_fmsy", "spict_bmsy"),
                    year = dimnames(catch(stk))$year,
                    iter = 1:it))
  ### start tracking
  tracking["Implementation", ac(iy)] <- catch(stk)[,ac(iy)]
  ### initialize advised catch
  tracking["advice", ac(98)] <- catch(stk)[,ac(99)]
  
  ### ---------------------------------------------------------------------- ###
  ### loop through simulation years ####
  ### ---------------------------------------------------------------------- ###
  for(i in vy[-length(vy)]){
    
    gc()
    ay <- an(i)
    cat(i, "> ")
    vy0 <- 1:(ay-y0) # data years (positions vector) - one less than current year
    sqy <- ac((ay-1):(ay-nsqy)) # years for status quo computations
    
    tracking["OM.f", ac(ay-1)] <- fbar(stk)[,ac(ay-1)]
    tracking["OM.ssb", ac(ay-1)] <- ssb(stk)[,ac(ay-1)]
    tracking["OM.catch", ac(ay-1)] <- catch(stk)[,ac(ay-1)]
    #tracking["advice", ac(ay-2)] <- catch(stk)[,ac(ay-1)]
    
    ### -------------------------------------------------------------------- ###
    ### OEM ####
    ### -------------------------------------------------------------------- ###
    
    ### -------------------------------------------------------------------- ###
    ### function o() - observations
    ctrl.oem <- ctrl.mp$ctrl.oem
    ctrl.oem$stk <- stk
    ctrl.oem$observations <- observations
    ctrl.oem$vy0 <- vy0
    ctrl.oem$ay <- ay
    ctrl.oem$tracking <- tracking
    if (!is.null(ctrl.mp$ctrl.oem)){
      o.out <- do.call("o", ctrl.oem)
    } else {
      ctrl.oem$method <- 'perfectInfo.wrapper'
      o.out <- do.call("o", ctrl.oem)
    }
    stk0 <- o.out$stk
    idx0 <- o.out$idx
    observations <- o.out$observations
    tracking <- o.out$tracking
    
    ### -------------------------------------------------------------------- ###
    ### MP ####
    ### -------------------------------------------------------------------- ###
    
    ### -------------------------------------------------------------------- ###
    ### Assessment/Estimator of stock statistics
    ### function f()
    if (!is.null(ctrl.mp$ctrl.f)){
      ctrl.f <- ctrl.mp$ctrl.f
      ctrl.f$stk <- stk0
      ctrl.f$idx <- idx0
      ctrl.f$tracking <- tracking
      out.assess <- do.call("f", ctrl.f)
      stk0 <- out.assess$stk
      tracking <- out.assess$tracking
    }
    tracking["Fperc",ac(ay)] <- fdy <- fbar(stk0)[,ac(ay-1)]
    
    
    ### -------------------------------------------------------------------- ###
    ### HCR parametrization
    ### function x()
    if (!is.null(ctrl.mp$ctrl.x)){
      ctrl.x <- ctrl.mp$ctrl.x
      ctrl.x$stk <- stk0
      ctrl.x$idx <- idx0
      ctrl.x$ay <- ay
      ctrl.x$iy <- iy
      ctrl.x$tracking <- tracking
      if(exists("hcrpars")) ctrl.x$hcrpars <- hcrpars
      out <- do.call("x", ctrl.x)
      hcrpars <- out$hcrpars
      tracking <- out$tracking
    }
    
    ### -------------------------------------------------------------------- ###
    ### HCR
    ### function h()
    if (!is.null(ctrl.mp$ctrl.h)){
      ctrl.h <- ctrl.mp$ctrl.h
      ctrl.h$stk <- stk0
      ctrl.h$ay <- ay
      ctrl.h$tracking <- tracking
      # If hcrpars was constructed by x() above, insert relevant parts into ctrl.h
      if(exists("hcrpars")){
        ctrl.h$hcrpars <- hcrpars
        hcrparnames <- names(hcrpars[names(hcrpars) %in% names(ctrl.h)])
        ctrl.h[hcrparnames] <- hcrpars[hcrparnames]
      }
      ctrl <- do.call("h", ctrl.h)
    } else {
      ctrl <- getCtrl(yearMeans(fbar(stk0)[,sqy]), "f", ay+1, it)
    }
    tracking["advice", ac(ay)] <- ctrl@trgtArray[ac(ay+1),"val",]
    
    ### -------------------------------------------------------------------- ###
    ### Management Implementation
    ### function k()
    if (!is.null(ctrl.mp$ctrl.k)){
      ctrl.k <- ctrl.mp$ctrl.k
      ctrl.k$ctrl <- ctrl
      ctrl.k$stk <- stk0
      ctrl.k$ay <- ay
      #ctrl.k$refCatch <- tracking["Implementation", ac(ay)]
      ctrl.k$tracking <- tracking
      out <- do.call("k", ctrl.k)
      ctrl <- out$ctrl
      tracking <- out$tracking
      tracking["Implementation", ac(ay)] <- ctrl@trgtArray[ac(ay+1),"val",]
    } else {
      tracking["Implementation", ac(ay)] <- tracking["advice", ac(ay+1)]
    }
    
    ### -------------------------------------------------------------------- ###
    ### Technical measures
    ### function w()
    if (!is.null(ctrl.mp$ctrl.w)){
      ctrl.w <- ctrl.mp$ctrl.w
      ctrl.w$stk <- stk0
      ctrl.w$sqy <- sqy
      ctrl.w$tracking <- tracking
      out <- do.call("w", ctrl.w)
      attr(ctrl, "snew") <- out$snew
      tracking <- out$tracking
    }
    
    ### -------------------------------------------------------------------- ###
    ### IEM ####
    ### -------------------------------------------------------------------- ###
    
    ### -------------------------------------------------------------------- ###
    ### implementation error
    ### function l()
    if (!is.null(ctrl.mp$ctrl.l)){
      ctrl.l <- ctrl.mp$ctrl.l
      ctrl.l$ctrl <- ctrl
      ctrl.l$tracking <- tracking
      out <- do.call("l", ctrl.l)
      ctrl <- out$ctrl
      tracking <- out$tracking
    }
    tracking["IEM",ac(ay)] <- ctrl@trgtArray[ac(ay+1), "val",]
    
    ### -------------------------------------------------------------------- ###
    ### OM
    ### -------------------------------------------------------------------- ###
    ### fleet dynamics/behaviour
    ### function j()
    
    if (!is.null(ctrl.mp$ctrl.j)){
      ctrl.j <- ctrl.mp$ctrl.j
      ctrl.j$ctrl <- ctrl
      ctrl.j$tracking <- tracking
      #if(exists("snew")) ctrl.l$snew <- snew
      out <- do.call("j", ctrl.j)
      ctrl <- out$ctrl
      tracking <- out$tracking
    }
    tracking["FleetDyn",ac(ay)] <- ctrl@trgtArray[ac(ay+1), "val",]
    
    ### -------------------------------------------------------------------- ###
    ### stock dynamics and OM projections
    ### function g()
    
    if(!is.null(attr(ctrl, "snew"))) harvest(stk)[,ac(ay+1)] <- attr(ctrl, "snew")
    ### project
    stk[] <- fwd(stk, ctrl = ctrl, sr = sr.om, sr.residuals = sr.om.res,
                 sr.residuals.mult = sr.om.res.mult, maxF = 5)[]
    
  } ### end of year loop
  cat("\n")
  
  ### save tracking in stk
  attr(stk, "tracking") <- tracking
  attr(stk, "ctrl.mp") <- ctrl.mp
  ### save stk to disk
  saveRDS(stk, file = paste0("/gpfs/afmcefas/simonf/output/", 
                             scn, "_", part, ".rds"))
  
} ### end of scenario loop

Sys.time()

### ------------------------------------------------------------------------ ###
### close down parallel workers
### ------------------------------------------------------------------------ ###

if (cluster_type == 1) {
  stopCluster(cl)
} else if (cluster_type == 2) {
  closeCluster(cl)
  mpi.quit()
}

### quit R, if not already closed by previous shutdown signals
q("yes")

