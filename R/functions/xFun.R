### ------------------------------------------------------------------------ ###
### x() HCR parametrization ####
### ------------------------------------------------------------------------ ###
#' prepare parameters for HCR
#' @param method name of the chosen HCR function
#' @param stk perceived stock
#' @param idx perceived index/indices
#' @param ay current (intermediate) year
#' @param tracking object for tracking
#' @param ... additional argument, passed on to function defined by method
#' @return list with two elements:
#'   \describe{
#'     \item{hcrpars}{parameters passed on to HCR}
#'     \item{tracking}{object for tracking, updated}
#'   }

x <- function(...) {
	args <- list(...)
	method <- args$method
	args$method <- NULL
	### Check inputs
	if (!is(args$stk,"FLS")) stop("stk argument must be an FLStock")
	### dispatch
	out <- do.call(method, args)
	### return
	out  
}


### ------------------------------------------------------------------------ ###
### WKLIFE VII
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### catch rule 3.2.1 from WKMSYCat34
### ------------------------------------------------------------------------ ###
### HCR 3.2.1 from WKMSYCat34 report
### so far handles only current catch and does not change anything else

wklife_3.2.1_x <- function(stk, tracking, ay, n_catch_yrs, interval = 1,
                           start_year, multiplier = NA, ...){
  
  ### last data year
  lst_yr <- range(stk)[["maxyear"]]
  
  hcrpars <- list()
  
  ### check if new advice requested
  if ((ay - start_year) %% interval == 0) {
    
    ### calculate current catch
    C_current <- yearMeans(catch(stk[, ac(seq(to = lst_yr, 
                                              length.out = n_catch_yrs))]))
    
    ### save in tracking object
    tracking["C_current", ac(lst_yr)] <- C_current
    
    ### add multiplier, if available
    tracking["HCRmult", ac(lst_yr)] <- multiplier
    
    ### extract required object from tracking for HCR
    hcrpars$factors <- tracking[c("HCR3.2.1r", "HCR3.2.1f", "HCR3.2.1b", 
                                  "C_current", "HCRmult"),
                                ac(lst_yr)]
    
  ### otherwise use values from last advice year
  } else {
    
    hcrpars$factors <- tracking[c("HCR3.2.1r", "HCR3.2.1f", "HCR3.2.1b", 
                                  "C_current", "HCRmult"),
                                ac(lst_yr - (interval - 1))]
    
  }
  
  ### return results
  return(list(tracking = tracking, hcrpars = hcrpars))
  
  
}

### ------------------------------------------------------------------------ ###
### catch rule 3.2.2 from WKMSYCat34
### ------------------------------------------------------------------------ ###
wklife_3.2.2_x <- function(stk, idx, tracking, ay, interval = 1, w = 1.4,
                           start_year, index_years = 1, multiplier = NA, 
                           Fproxy_type = "FproxyMSY",
                           ...) {

  ### last data year
  lst_yr <- range(stk)[["maxyear"]]
  
  hcrpars <- list()
  
  ### retrieve index reference values
  ### check if I_trigger available
  if (!all(is.na(attr(stk, "refpts")["I_trigger"]))) {
    I_trigger <- attr(stk, "refpts")["I_trigger"]
  } else {
    I_trigger <- attr(stk, "refpts")["I_lim"] * w
  }
  
  ### get index FproxyMSY value
  if (Fproxy_type == "FproxyMSY") {
    ### proxy when fished at MSY
    FproxyMSY <- attr(stk, "refpts")["FproxyMSY"]
  } else if (Fproxy_type == "FproxyMSY_MK") {
    ### proxy when fished so that Lmean is LF=M proxy
    FproxyMSY <- attr(stk, "refpts")["FproxyMSY_MK"]
  } else {
    stop("unknown Fproxy type")
  }
  
  ### check if new advice requested
  if ((ay - start_year) %% interval == 0) {
    
    ### get biomass index
    idx_tmp <- quantSums(index(idx[[1]]))
    ### subset to requested years
    idx_tmp <- idx_tmp[, tail(dimnames(idx_tmp)$year, index_years)]
    ### mean over years = current index value
    I_current <- quantMeans(idx_tmp)
    
    ### save in tracking object
    tracking["I_current", ac(lst_yr)] <- I_current
    
    ### add multiplier, if available
    if (!is.na(multiplier)) tracking["HCRmult", ac(lst_yr)] <- multiplier
    

    ### extract required object from tracking for HCR
    hcrpars$factors <- FLPar(I_current = I_current, HCR_mult = multiplier,
                             I_trigger = I_trigger, FproxyMSY = FproxyMSY)
    
    ### otherwise use values from last advice year
  } else {
    
    hcrpars$factors <- FLPar(I_current = tracking["I_current",
                                                  ac(lst_yr - (interval - 1))], 
                             HCR_mult = tracking["HCRmult",
                                                 ac(lst_yr - (interval - 1))],
                             I_trigger = I_trigger, 
                             FproxyMSY = FproxyMSY)
    
  }
  
  ### return results
  return(list(tracking = tracking, hcrpars = hcrpars))
  
}

### ------------------------------------------------------------------------ ###
### catch rule 3.1
### ------------------------------------------------------------------------ ###
### only extract values from tracking
wklife_3.1_x <- function(tracking, ay, ...) {
  
  ### extract advice
  hcrpars <- list()
  ### stored in year before (data year)
  hcrpars$factors <- tracking["advice", ac(ay - 1)]
  
  return(list(tracking = tracking, hcrpars = hcrpars))
  
}




#' Generate HCR parameters from the output of SPiCT assessment
#' 
#' The HCR parameters are generated from the SPiCT output which is attached to the perceieved stock as an attribute
#' The corresponding HCR is the standard ICES HCR.
#' For that the following parameters are generated:
#' fmin, ftrg, blim, bsafe, ssb_lag
spict_hcr_pars <- function(stk, tracking, ay, iy, ...){
  # hcrpars is a list
  niters <- length(stk@spict_fit)
  hcrpars <- list()
  hcrpars$fmin <- 0
  hcrpars$ftrg <-  unlist(lapply(stk@spict_fit, function(fit) get.par('logFmsy', fit, exp=TRUE)[,'est'])) # Fmsy
  hcrpars$blim <- unlist(lapply(stk@spict_fit, function(fit) min(get.par('logB', fit, exp=TRUE)[,'est']))) # lowest est B in timeseries
  hcrpars$bsafe <- unlist(lapply(stk@spict_fit, function(fit) get.par('logBmsy', fit, exp=TRUE)[,'est'])) # Bmsy
  hcrpars$ssb_lag <- 1
  return(list(hcrpars = hcrpars, tracking = tracking))
}



