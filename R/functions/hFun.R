# Harvest Control Rule function: h()
#' Evaluate the chosen HCR function
#'
#' Evaluate the chosen HCR function using the current stock perception and a control
#' @param method Name of the chosen HCR function.
#' @param stk The perceived stock.
#' @param ay The current year. The management control (e.g. F or effort) will be set in ay+1.
#' @param EFF Effort array (if effort management is being used).
#' @param EFF0 Tracking array.
#' @param control The control object for the chosen HCR function. A list of parameters.
h <- function(...){
  args <- list(...)
  method <- args$method
  args$method <- NULL
  # Check inputs
  if(!is(args$stk,"FLS")) stop("stk argument must be an FLStock")
  # dispatch
  ctrl <- do.call(method, args)
  # check outputs
  if(!is(ctrl, "fwdControl")) stop("The HCR must return and object of class fwdControl")	
  # return
  ctrl  
}

#' The typical HCR used by ICES
#'
#' The typical HCR used by ICES which sets a target F based on the SSB based on 4 parameters: blim, bsafe, fmin and ftrg.
#' F increases linearly between SSB = blim and SSB = bsafe, from F = fmin to F = ftrg.
#' If:
#' B < Blim, F = Fmin;
#' B > trigger, F = Fmsy;
#' B > Blim & B < trigger, F linear between Fbycatch and Fmsy;
#' F = ftrg is the maximum F, F = fmin is the minimum F.
#' F is set in year ay, based on SSB in year ay - ssb_lag.
#' The control argument is a list of parameters used by the HCR.
#' @param stk The perceived FLStock.
#' @param control A list with the elements fmin, ftrg, blim, bsafe and ssb_lag, all of which are numeric.
#' @param ay The year for which the target F is set, based on the SSB in year (ay - control$ssb_lag).
ices_hcr <- function(stk, fmin, ftrg, blim, bsafe, ssb_lag, ay){
  # rule
  ssb <- ssb(stk)[, ac(ay-ssb_lag)]
  fout <- FLQuant(fmin, dimnames=list(iter=dimnames(ssb)$iter))
  fout[ssb >= bsafe] <- ftrg
  inbetween <- (ssb < bsafe) & (ssb > blim)
  gradient <- (ftrg - fmin) / (bsafe - blim)
  fout[inbetween] <- (ssb[inbetween] - blim) * gradient + fmin
  
  # create control file
  ctrl <- getCtrl(c(fout), "f", ay+1, dim(fout)[6])
  # return
  ctrl
}

#' A fixed target f
#'
#' No matter what get F= Ftarget
#' The control argument is a list of parameters used by the HCR.
#' @param stk The perceived FLStock.
#' @param control A list with the element ftrg (numeric).
fxdFtrg_hcr <- function(stk, ftrg, ay){
	# rule 
	if(!is(ftrg, "FLQuant")) ftrg <- FLQuant(ftrg, dimnames=list(iter=dimnames(stk@catch)$iter))
	# create control file
	ctrl <- getCtrl(c(ftrg), "f", ay+1, dim(ftrg)[6])
	# return
	ctrl
}


movFtrg_hcr <- function(stk, hcrpars, ay){
	# rule 
	if(!is(hcrpars, "FLQuant")) hcrpars <- FLQuant(hcrpars, dimnames=list(iter=dimnames(stk@catch)$iter))
	# create control file
	ctrl <- getCtrl(c(hcrpars), "f", ay+1, dim(hcrpars)[6])
	# return
	ctrl
}


### ------------------------------------------------------------------------ ###
### catch rules for WKLIFE VII
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### catch rule 3.2.1 from WKMSYCat34
### ------------------------------------------------------------------------ ###

wklife_3.2.1_h <- function(hcrpars, ay, tracking, stk, 
                           upper_constraint = Inf, lower_constraint = 0, ...){

  ### parameters for HCR
  factors <- hcrpars$factors
  
  ### calculate catch advice by multiplying r, f, b
  adv <- apply(X = factors, MARGIN = 6, prod, na.rm = TRUE)
  
  ### apply TAC constraint, if requested
  if (!is.infinite(upper_constraint) | lower_constraint != 0){
    
    ### get last advice
    adv_last <- tracking["advice", ac(dims(stk)$maxyear - 1)]
    ### ratio of new advice/last advice
    adv_ratio <- adv/adv_last
    
    ### upper constraint
    if (!is.infinite(upper_constraint)){
      ### find positions
      pos_upper <- which(adv_ratio > upper_constraint)
      ### limit advice
      if (length(pos_upper) > 0){
        adv[,,,,, pos_upper] <- adv_last[,,,,, pos_upper] * upper_constraint
      }
    ### lower constraint
    }
    if (lower_constraint != 0){
      ### find positions
      pos_lower <- which(adv_ratio < lower_constraint)
      ### limit advice
      if (length(pos_lower) > 0){
        adv[,,,,, pos_lower] <- adv_last[,,,,, pos_lower] * lower_constraint
      }
    }
  }
  
  ### create fwdControl object
  ctrl <- getCtrl(values = c(adv), quantity = "catch", years = ay+1, 
                  it = dim(adv)[6])
  
  return(ctrl)
  
}
  
### ------------------------------------------------------------------------ ###
### catch rule 3.2.2 from WKMSYCat34
### ------------------------------------------------------------------------ ###

wklife_3.2.2_h <- function(hcrpars, ay, ...) {
  
  ### first part of advice: current index * proxy
  adv1 <- hcrpars$factors["I_current"] * hcrpars$factors["I_F_proxy"]
  
  ### ratio current index / trigger
  adv2 <- hcrpars$factors["I_current"] / hcrpars$factors["I_trigger"]
  ### capped at 1
  adv2 <- apply(adv2, 1:2, function(x){
    min(x, 1)
  })
  
  ### combine both parts
  adv <- adv1 * adv2
  
  ### multiply, if requested
  pos <- which(!is.na(hcrpars$factors["HCR_mult"]))
  if (length(pos) > 0) {
    adv[, pos] <- adv[, pos] * hcrpars$factors["HCR_mult", pos]
  }
  
  ### create fwdControl object
  ctrl <- getCtrl(values = c(adv), quantity = "catch", years = ay+1, 
                  it = dims(adv)$iter)
  
  return(ctrl)
  
}

### ------------------------------------------------------------------------ ###
### catch rule 3.1 from WKMSYCat34
### ------------------------------------------------------------------------ ###

wklife_3.1_h <- function(hcrpars, ay, ...){
  
  ### get advice
  adv <- hcrpars$factors
  
  ### create fwdControl object
  ctrl <- getCtrl(values = c(adv), quantity = "catch", years = ay+1, 
                  it = dim(adv)[6])
  
  return(ctrl)
  
}

