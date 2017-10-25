# functions.R - DESC
# functions.R

# Copyright 2003-2012 FLR Team. Distributed under the GPL 2 or later
# Maintainer: Iago Mosqueira, JRC
# $Id: $

# oneWayTrip {{{
oneWayTrip <- function(stk, sr, brp, fmax=refpts(brp)['crash', 'harvest'] * 0.80,
	years=seq(dims(stk)$minyear + 1, dims(stk)$maxyear),
  residuals=FLQuant(1, dimnames=dimnames(rec(stk))),f0=NULL) {

	# limits
  if(is.null(f0)){ 
    f0 <- c(fbar(stk)[,1])
    } else{
      f0=as.vector(f0)
      }

	fmax <- c(fmax)
	rate <- exp((log(fmax) - log(f0)) / (length(years)))

	# linear trend: f <- rate ^ (seq(0, length(years))) * f0
	fs <- (matrix(rate, nrow=length(years) + 1, ncol=length(rate)) ^
  	matrix(seq(0, length(years)), nrow=length(years) + 1, ncol=length(rate)) *
	  matrix(f0, nrow=length(years) + 1, ncol=length(rate)))[-1,]

	# fwdControl
  ftar <- FLQuant(c(fs), dimnames=list(year=years, iter=seq(length(rate))), quant='age')
  
	
  # fwd
	res <- fwd(stk, control=as(FLQuants(f=ftar), "fwdControl"), sr=sr, residuals=residuals)
	
	return(res)
} # }}}

# rollerCoaster {{{
rollerCoaster <- function(stk, sr, brp, fmax=refpts(brp)['crash', 'harvest']*0.75,
	fmsy=refpts(brp)['msy', 'harvest'], up=0.06, down=0.04,
	years=seq(dims(stk)$minyear + 1, dims(stk)$maxyear),
	residuals=FLQuant(1, dimnames=dimnames(rec(stk))),f0=NULL) {

	# F
	f <- rep(NA, length(years))

	# limits
	if(is.null(f0)){
	  f0 <- c(fbar(stk)[,1])[1]
	  } else{
	  f0=as.vector(f0)
	}
	fmax <- c(fmax)

	# linear trend
	rateup <- log(fmax/f0) / log(1 + up)
	fup <- f0 * ((1 + up) ^ (0:ceiling(rateup)))
	lfup <- length(fup)
	f[1:lfup] <- fup

	# at the top
	f[lfup:(lfup+5)] <- fup[lfup]

	# coming down!
	ratedo <- c(log(fmsy/f[length(fup)+5]) / log(1 + down))
	lfdo <- length(f) - (lfup +6) + 1
	fdo <- f[lfup + 5] * ((1 + down) ^ seq(0, ceiling(ratedo), length=lfdo))
	f[(lfup+6):length(f)] <- fdo[1:lfdo]
	
	# fwdControl
	fctl <- fwdControl(data.frame(year=years, quant='f', value=f))

	# fwd
	stk <- fwd(stk, control=fctl, sr=sr, residuals=residuals)
	
	return(stk)
} # }}}
