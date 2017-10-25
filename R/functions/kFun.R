# Management implementation function: k()

#' Evaluate the chosen Management implementation function
#'
#' Evaluate the chosen HCR function using the current stock perception and a control.
#' For example, TAC control or effort control.
#' Returns a ctrl object for projecting the OM.
#' @param method Name of the chosen implementation function.
#' @param stk The perceived stock.
#' @param ay The current year. The management control (e.g. TAC or effort) will be set in ay+1.
#' @param EFF Effort array (only used if effort management is being used).
#' @param EFF0 Tracking array.
#' @param control The control object for the chosen management implementation function. A list of parameters.

k <- function(...){
  args <- list(...)
  method <- args$method
  args$method <- NULL
  # Check inputs
  if(!is(args$stk,"FLS")) stop("stk argument must be an FLStock")
  # dispatch
  out <- do.call(method, args)
  # check outputs
  if(!is(out$ctrl, "fwdControl")) stop("The HCR must return an object of class fwdControl")	
  # return
  out  
}

#' TAC implementation function
#'
#' Performs a short term forecast (STF) to hit the target F in year ay+1.
#' The resulting catch in year ay+1 is the TAC, i.e. the TAC that will result in Fbar = Ftarget.
#' The STF uses geometric mean recruitment. Fbar in the intermediate years (i.e. years between last data year and ay+1) are set as the mean of the last nsqy years.
#' The control argument is a list of parameters:
#' nsqy - number of years to average over to get Fbar for STF
#' delta_tac_min - constraint on the TAC
#' delta_tac_max - constraint on the TAC
#' @param stk The perceived FLStock.
#' @param imp_control A list with the elements: nsqy, delta_tac_min, delta_tac_max
#' @param ay The year for which the target F is set, based on the SSB in year (ay - control$ssb_lag).
#' @param EFF0 The tracking array
#' @param EFF Not used by this function but may be used by the other implementation functions

tac <- function(stk, ctrl, ay, nsqy=3, delta_tac_max=NA, delta_tac_min=NA, tracking){
  refCatch <- tracking["Implementation", ac(ay)]
  # Year range of perceived stock
  yrs <- as.numeric(dimnames(stock.n(stk))$year)
  last_data_yr <- yrs[length(yrs)]
  # Status quo years
  sqy <- (last_data_yr-nsqy+1):last_data_yr
  # Get the Fbar for the intermediate years
  fsq0 <- yearMeans(fbar(stk)[,ac(sqy)])
  # Number of intermediate years (between last data year and ay+1)
  ninter_yrs <- ay - last_data_yr
  # Control object for the STF
  ctrl <- getCtrl(c(rep(fsq0, times=ninter_yrs), ctrl@trgtArray[,"val",]), "f", (last_data_yr+1):(ay+1), dim(stock.n(stk))[6])
  # Number of projection years
  nproj_yrs <- (ay+1) - last_data_yr
  stkTmp <- stf(stk, nproj_yrs, wts.nyears=nsqy)
  # Set geomean sr relationship
  gmean_rec <- c(exp(yearMeans(log(rec(stk)[,ac(sqy)]))))
  # Project!
  stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=list(model="mean", params = FLPar(gmean_rec,iter=it)))
  # Get TAC for following year that results from hitting the F in STF
  TAC <- catch(stkTmp)[,ac(ay+1)]
  # Get TAC in previous year from EFF0
  #lastTAC <- (EFF0["Implementation",ac(ay)])
  # Apply limit to last TAC - doesn't matter if the deltas are NA
  #upper_limit <- lastTAC * delta_tac_max
  #lower_limit <- lastTAC * delta_tac_min
  upper_limit <- refCatch * delta_tac_max
  lower_limit <- refCatch * delta_tac_min
  TAC <- pmin(c(upper_limit), c(TAC), na.rm=TRUE)
  TAC <- pmax(c(lower_limit), c(TAC), na.rm=TRUE)
  # Update EFF0
  #EFF0["Implementation",ac(ay+1)] <- TAC
  # Make control for the OM projection
  #ctrl <- getCtrl(c(TAC), "catch", ay+1, it)
  #return(list(ctrl=ctrl, EFF0=EFF0))
  ctrl <- getCtrl(c(TAC), "catch", ay+1, it)
  list(ctrl = ctrl, tracking = tracking)
  
}

#' TAC implementation function using SPiCT output
#'
#' Uses SPiCT forecasting function to get F in ay+1.
#' Although TAC has already been set in year ay, we cannot forecast SPiCT with catch in ay and Fmsy in ay+1.
spict_tac <- function(stk, ctrl, ay, tracking){
  niters <- dim(catch(stk))[6]
  # Loop over iterations
  TAC <- rep(NA, niters)
  for (iter in 1:niters){
    # SPiCT projection sets F in managment years through a factor applied to F in final estimated year
    # Get F in last estimated year
    lastF <- get.par('logF', stk@spict_fit[[iter]], exp = TRUE)[ac(ay-1), "est"]
    # Calculate factor to hit target in management starts based on target F in control
    ffac <- ctrl@trgtArray[1,"val",iter] / lastF
    # Extract input from original fit
    proj_input <- stk@spict_fit[[iter]]$inp
    # Force managment to start in ay+1
    proj_input$manstart <- ay+1 # management starts
    proj <- prop.F(fac = ffac, inpin = proj_input, repin = stk@spict_fit[[iter]], dbg = 0)
    # Check F
    # get.par('logF', proj, exp = TRUE)[ac(2011:2015), "est"]
    # Get catch - should have set up timepredc in assessment function
    projC <- get.par('logCpred', proj, exp = TRUE)[, "est"]
    TAC[iter] <- projC[length(projC)]
  }
  ctrl <- getCtrl(c(TAC), "catch", ay+1, niters)
  return(list(ctrl = ctrl, tracking = tracking))
}

