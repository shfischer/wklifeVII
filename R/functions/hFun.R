### ------------------------------------------------------------------------ ###
### h(): Harvest Control Rule ####
### ------------------------------------------------------------------------ ###
#' apply Harvest Control Rule (HCR)
#' @param method name of the chosen HCR function
#' @param stk perceived stock
#' @param ay current (intermediate) year
#' @param hcrpars parameters for HCR
#' @param tracking object for tracking
#' @param ... additional argument, passed on to function defined by method
#' @return fwdControl object used later in forecast
h <- function(...) {
  args <- list(...)
  method <- args$method
  args$method <- NULL
  ### Check inputs
  if (!is(args$stk,"FLS")) stop("stk argument must be an FLStock")
  ### dispatch
  ctrl <- do.call(method, args)
  # check outputs
  if (!is(ctrl, "fwdControl")) 
    stop("The HCR must return and object of class fwdControl")	
  ### return
  ctrl  
}

### ------------------------------------------------------------------------ ###
### catch rules for WKLIFE VII
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### catch rule 3.2.1 from WKMSYCat34
### ------------------------------------------------------------------------ ###

wklife_3.2.1_h <- function(hcrpars, ay, tracking, stk, 
                           upper_constraint = Inf, lower_constraint = 0, ...) {

  ### parameters for HCR
  factors <- hcrpars$factors
  
  ### calculate catch advice by multiplying r, f, b
  adv <- apply(X = factors, MARGIN = 6, prod, na.rm = TRUE)
  
  ### apply TAC constraint, if requested
  if (!is.infinite(upper_constraint) | lower_constraint != 0) {
    
    ### get last advice
    adv_last <- tracking["advice", ac(dims(stk)$maxyear - 1)]
    ### ratio of new advice/last advice
    adv_ratio <- adv/adv_last
    
    ### upper constraint
    if (!is.infinite(upper_constraint)) {
      ### find positions
      pos_upper <- which(adv_ratio > upper_constraint)
      ### limit advice
      if (length(pos_upper) > 0) {
        adv[,,,,, pos_upper] <- adv_last[,,,,, pos_upper] * upper_constraint
      }
    ### lower constraint
    }
    if (lower_constraint != 0) {
      ### find positions
      pos_lower <- which(adv_ratio < lower_constraint)
      ### limit advice
      if (length(pos_lower) > 0) {
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
  adv1 <- hcrpars$factors["I_current"] * hcrpars$factors["FproxyMSY"]
  
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
  ctrl <- getCtrl(values = c(adv), quantity = "catch", years = ay + 1, 
                  it = dims(adv)$iter)
  
  return(ctrl)
  
}

### ------------------------------------------------------------------------ ###
### catch rule 3.1 from WKMSYCat34
### ------------------------------------------------------------------------ ###

wklife_3.1_h <- function(hcrpars, ay, ...) {
  
  ### get advice
  adv <- hcrpars$factors
  
  ### create fwdControl object
  ctrl <- getCtrl(values = c(adv), quantity = "catch", years = ay + 1, 
                  it = dim(adv)[6])
  
  return(ctrl)
  
}

