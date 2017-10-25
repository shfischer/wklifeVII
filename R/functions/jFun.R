#--------------------------------------------------------------------
# fleet dynamics, mostly changes in f-at-age
j <- function(...){
	args <- list(...)
	method <- args$method
	args$method <- NULL
	# Check inputs
	if(!is(args$ctrl, "fwdControl")) stop("ctrl must be of class fwdControl")
	# dispatch
	out <- do.call(method, args)
	# check outputs
	if(!is(out$ctrl, "fwdControl")) stop("The HCR must return and object of class fwdControl")	
	# return
	out  
}

hyperstability.wrapper <- function(ctrl, beta=1, maxF=2, alpha=maxF^(1-beta), tracking){
	# Only operates on F targets - so nothing happens to TAC
	# This function creates a control file to be later used in the fwd()
	# function where two optional relations are established between
	# fishing effort and fishing mortality
	# Beta is in this MSE either 1 for a 1:1 linear relationship between
	# F and effort, if beta = 0.7, the relation is not linear and it can
	# mimick a hyperstability scenario.
	# alpha = maxF^(1-beta) # linear meets curve at maxF
	ctrl@trgtArray[ctrl@target[,"quantity"]=="f",,] <- alpha *
	ctrl@trgtArray[ctrl@target[,"quantity"]=="f",,]^beta
	list(ctrl=ctrl, tracking=tracking)
	
}

### function for avoiding the issue that the catch is higher than
### the targeted catch
ctrl_catch_workaround <- function(ctrl, tracking, ...) {
  
  ### duplicate target
  ctrl@target <- rbind(ctrl@target, ctrl@target, ctrl@target)
  ### replace catch in second row with landings
  ctrl@target$quantity[1] <- "landings"
  ctrl@target$quantity[3] <- "catch"
  
  ### extract target values
  val_temp <- ctrl@trgtArray[, "val", ]
  
  ### extend trgtArray
  ### extract dim and dimnames
  dim_temp <- dim(ctrl@trgtArray)
  dimnames_temp <- dimnames(ctrl@trgtArray)
  ### duplicate years
  dim_temp[["year"]] <- dim_temp[["year"]] * 3
  dimnames_temp$year <- rep(dimnames_temp$year, 3)
  
  ### create new empty array
  trgtArray <- array(data = NA, dim = dim_temp, dimnames = dimnames_temp)
  
  ### fill with values
  ### first as target
  trgtArray[1, "val", ] <- val_temp
  ### then again, but as max
  trgtArray[2, "max", ] <- val_temp
  ### min F
  trgtArray[3, "max", ] <- val_temp
  
  ### insert into ctrl object
  ctrl@trgtArray <- trgtArray
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}  
