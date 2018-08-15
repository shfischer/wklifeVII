### ------------------------------------------------------------------------ ###
### w(): Technical measures ####
### ------------------------------------------------------------------------ ###
#' implement technical measures
#' @param method name of the chosen HCR function
#' @param stk perceived stock
#' @param tracking object for tracking
#' @param ... additional argument, passed on to function defined by method
#' @return ...


w <- function(...){
	args <- list(...)
	method <- args$method
	args$method <- NULL
	# Check inputs
	if (!is(args$stk,"FLS")) stop("stk argument must be an FLStock")
	# dispatch
	out <- do.call(method, args)
	# check outputs
	if (!is(out$snew, "FLQuant")) stop("The technical measures must return and object of class FLQuant")	
	# return
	out  
}

### not used for data limited MSE

