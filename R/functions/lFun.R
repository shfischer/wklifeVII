### ------------------------------------------------------------------------ ###
### l(): Implementation Error ####
### ------------------------------------------------------------------------ ###
#' add noise to targets
#' @param method name of the chosen HCR function
#' @param ctrl fwdControl object with targets
#' @param tracking object for tracking
#' @param ... additional argument, passed on to function defined by method
#' @return list with two elements:
#'   \describe{
#'     \item{ctrl}{fwdControl object with modified targets}
#'     \item{tracking}{object for tracking, updated}
#'   }

l <- function(...) {
	args <- list(...)
	method <- args$method
	args$method <- NULL
	### Check inputs
	if (!is(args$ctrl, "fwdControl")) stop("ctrl must be of class fwdControl")
	### dispatch
	out <- do.call(method, args)
	# check outputs
	if (!is(out$ctrl, "fwdControl")) 
	  stop("The HCR must return and object of class fwdControl")	
	### return
	out  
}

### add random noise to targets
noise.wrapper <- function(ctrl, fun = "rlnorm", mean = 0, sd = 0.1, 
                          multiplicative = TRUE, tracking) {
  
  ### define arguments for random number generation
	iem <- list(mean = mean, sd = sd, n = sum(!is.na(ctrl@trgtArray)))
	### random numbers created for all targets, including min/max/val
	### if NA, they stay NA
	### multiply targets or add them
	if (multiplicative) {
		ctrl@trgtArray <- do.call(fun, iem)*ctrl@trgtArray
	} else {
		ctrl@trgtArray <- do.call(fun, iem) + ctrl@trgtArray
	}
	### return control object and tracking
	list(ctrl = ctrl, tracking = tracking)
	
}

