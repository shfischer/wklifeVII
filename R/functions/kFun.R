### ------------------------------------------------------------------------ ###
### k(): Management implementation ####
### ------------------------------------------------------------------------ ###
#' implement management from HCR
#' @param method name of the chosen HCR function
#' @param stk perceived stock
#' @param ctrl fwdControl object with targets
#' @param ay current (intermediate) year
#' @param tracking object for tracking
#' @param ... additional argument, passed on to function defined by method
#' @return list with two elements:
#'   \describe{
#'     \item{ctrl}{fwdControl object used later in forecast}
#'     \item{tracking}{object for tracking, updated}
#'   }

k <- function(...){
  args <- list(...)
  method <- args$method
  args$method <- NULL
  # Check inputs
  if (!is(args$stk,"FLS")) stop("stk argument must be an FLStock")
  # dispatch
  out <- do.call(method, args)
  # check outputs
  if (!is(out$ctrl, "fwdControl")) 
    stop("The HCR must return an object of class fwdControl")	
  # return
  out  
}

### not used for data limited MSE

