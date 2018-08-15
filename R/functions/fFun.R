### ------------------------------------------------------------------------ ###
### f() Assessment/Estimator of stock statistics ####
### ------------------------------------------------------------------------ ###
#' generate observations: observed index/indices and stock
#' @param method name of the chosen HCR function
#' @param stk observed stock
#' @param idx observed index/indices
#' @param tracking object for tracking
#' @param ... additional argument, passed on to function defined by method
#' @return list with two elements:
#'   \describe{
#'     \item{stk}{observed stock, after assessment/estimation}
#'     \item{tracking}{object for tracking, updated}
#'   }

f <- function(...) {
	args <- list(...)
	method <- args$method
	args$method <- NULL
	### checks 
	if (!is(args$stk, "FLS")) stop("stk must be of class FLStock")
	if (!is(args$idx, "FLIndices")) stop("idx must be of class FLIndices")
	### dispatch
	out <- do.call(method, args)
	if (!is(out$stk, "FLS")) stop("stk must be of class FLStock")
  ### return
	out
}


### ------------------------------------------------------------------------ ###
### WKLIFE VII
### ------------------------------------------------------------------------ ###
### calculate length at first capture and mean length
### and factors from WKMSYCat34 report for catch rule 3.2.1

wklife_3.2.1_f <- function(stk, tracking, n_catch_yrs, 
                           option_f = 0, option_r = 0, option_b = 0, ...) {
  
  ### search for empty years and only compute values if no values available
  ### find year position without data
  pos <- unlist(lapply(dimnames(tracking)$year, function(x){
    all(is.na(c(tracking["L_c", x])))
  }))
  ### corresponding years
  pos_years <- dimnames(tracking)$year[pos]
  ### years coinciding with range in stk
  yrs_req <- intersect(pos_years, dims(stk)$minyear:dims(stk)$maxyear)
  
  ### calculate length at first capture
  L_c <- calc_Lc(object = catch.n(stk)[, ac(yrs_req)])
  
  ### calculate mean length above length of first capture
  L_mean <- calc_mean(object = catch.n(stk)[, ac(yrs_req)], 
                      min = L_c[, ac(yrs_req)])
  
  ### save in tracking object
  tracking["L_c", ac(yrs_req)] <- L_c
  tracking["L_mean", ac(yrs_req)] <- L_mean
  
  ### continue with calculations, if advice required for this year
  
  ### calculate current mean length
  L_current <- yearMeans(L_mean[, tail(dimnames(L_mean)$year, n_catch_yrs)])
  
  ### extract entire time series of mean lengths, required for some options
  L_mean <- tracking["L_mean", dimnames(catch(stk))$year]
  ### compute factor f, depending on requested choice
  ### arguments required
  # args <- list(lhpar = attr(stk, "lhpar"),
  #              L_c = L_c, L_mean = L_mean,
  #              ...)
  
  args <- list(lhpar = attr(stk, "lhpar"), refpts = attr(stk, "refpts"),
               L_c = L_c[, dimnames(L_current)$year], L_current = L_current,
               L_mean = L_mean, n_catch_yrs = n_catch_yrs,
               option_f = option_f, option_r = option_r, option_b = option_b,
               ...)
  
  ### calculate factor f
  ### proxy for ratio FMSY/current exploitation
  factor_f <- do.call(paste0("wklife_3.2.1_f_", args$option_f), args)
  ### avoid negative values
  if (is(factor_f, "FLQuant")) factor_f <- apply(factor_f, 1:6, function(x) {
    max(x, 0)
  })
  
  ### save in tracking
  tracking["HCR3.2.1f", ac(tail(yrs_req, 1))] <- factor_f
  
  ### calculate factor r
  ### trend in stock biomass
  factor_r <- do.call(paste0("wklife_3.2.1_r_", args$option_r), args)
  
  ### save in tracking
  tracking["HCR3.2.1r", ac(tail(yrs_req, 1))] <- factor_r
  
  ### add tracking to args
  args$tracking <- tracking
  
  ### calculate factor b
  ### ~ proxy for ratio (current stock size)/MSYBtrigger
  factor_b <- do.call(paste0("wklife_3.2.1_b_", args$option_b), args)
  
  ### save in tracking
  tracking["HCR3.2.1b", ac(tail(yrs_req, 1))] <- factor_b
  
  ### return results
  return(list(stk = stk, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### options for r
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### r, option a
### ------------------------------------------------------------------------ ###
### 2 over 3 rule
### i.e. mean of recent 2 years divided by mean of 3 preceding years
wklife_3.2.1_r_a <- function(idx, n_1 = 2, n_2 = 3, ...) {
  
  ### extract index (use only first element of FLIndices)
  idx_temp <- quantSums(index(idx[[1]]))
  
  ### get years for calculation of ratio
  yrs_a <- tail(dimnames(idx_temp)$year, n_1)
  yrs_b <- tail(setdiff(dimnames(idx_temp)$year, yrs_a), n_2)
  
  ### calculate ratio of means
  factor_r <- yearMeans(idx_temp[, yrs_a]) / yearMeans(idx_temp[, yrs_b])
  
  return(factor_r)
  
}
### ------------------------------------------------------------------------ ###
### r, option b
### ------------------------------------------------------------------------ ###
### based on  r = exp(w * slope of lm of log stock size index)
wklife_3.2.1_r_b <- function(idx, r_w = 1, ...) {
  
  ### extract index (use only first element of FLIndices)
  idx_temp <- quantSums(index(idx[[1]]))
  
  ### get last 5 years
  yrs_used <- tail(dimnames(idx_temp)$year, 5)
  ### subset index to this period
  idx_temp <- idx_temp[, yrs_used]
  
  ### apply linear model (regression to all iterations)
  slopes <- apply(X = idx_temp, MARGIN = c(1, 6), FUN = function(x){
    ### extract only slope
    lm(log(x) ~ an(dimnames(x)$year))$coefficients[[2]]
  })
  
  ### calculate r
  factor_r <- exp(r_w * slopes)
  
  return(factor_r)
  
}

### ------------------------------------------------------------------------ ###
### r, option 0
### ------------------------------------------------------------------------ ###
wklife_3.2.1_r_0 <- function(...) {
  return(NA)
}

### ------------------------------------------------------------------------ ###
### options for f
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### f, option a
### ------------------------------------------------------------------------ ###
wklife_3.2.1_f_a <- function(L_c, L_current, lhpar, refpts, 
                             perfect_knowledge = FALSE, MK = NULL, ...) {
  ### retrieve or calculate LF=M
  if (isTRUE(perfect_knowledge)) {
    L_FM <- refpts["LFeFmsy"]
  } else {
    ### use shortcut for M/K, if requested
    if (!is.null(MK)) {
      L_FM <- (lhpar["L_inf"] + 2 * MK * L_c) / 
        (1 + 2 * MK)
    } else {
      L_FM <- (lhpar["L_inf"] + 2 * lhpar["M"]/lhpar["K"] * L_c) / 
        (1 + 2 * lhpar["M"]/lhpar["K"])
    }
  }
  ### calculate f
  return(L_current / L_FM)
}
### ------------------------------------------------------------------------ ###
### f, option b
### ------------------------------------------------------------------------ ###
wklife_3.2.1_f_b <- function(L_c, L_current, lhpar, refpts, 
                             perfect_knowledge = FALSE, ...) {
  ### calculate current Z
  Z_current <- lhpar["K"] * (lhpar["L_inf"] - L_current) / (L_current - L_c)
  ### calculate f
  if (isTRUE(perfect_knowledge)) {
    return(refpts["F_MSY"] / (Z_current - lhpar["M"]))
  } else {
    return(lhpar["M"] / (Z_current - lhpar["M"]))
  }
}
### ------------------------------------------------------------------------ ###
### f, option c
### ------------------------------------------------------------------------ ###
### f = F0.1 / (Z_current - M) where Z_current is calculated from Gedamke
### and Hoenig (2006) and F0.1 from length-based YPR

wklife_3.2.1_f_c <- function(lhpar, L_mean, L_c, n_catch_yrs, 
                             perfect_knowledge = FALSE, refpts, ...) {
  
  ### estimate Z from Gedamke and Hoenig (2006) method
  
  ### start value for optimization
  Z_start <- 0.5
  
  ### calculate/optimize to get Z
  # Z <- lapply(seq(to = dims(L_mean)[["maxyear"]], length.out = n_catch_yrs), 
  #             function(year_z){
  
  ### get years to be used
  yrs_used <- tail(dimnames(L_mean)$year, n_catch_yrs)
    
  #dbg_tmp <- NULL
  Z_iter <- lapply(dimnames(L_mean)$iter, function(x) {
    #dbg_tmp <<- x
    ### optimize
    res <- optim(par = Z_start, fn = GH, method = "BFGS", 
                 control = list(maxit = 1e+6, abstol = 1e-7), hessian = FALSE,
                 L_current_i = L_mean[, yrs_used,,,, x], 
                 L_c_i = L_c[, yrs_used,,,, x], 
                 L_inf_i = lhpar["L_inf", x], K_i = lhpar["K", x])
    
    ### alternatively, faster & only positive, but cannot handle NAs/NaN/Inf
    # res <- optim(par = Z_start, fn = GH3, method = "L-BFGS-B",
    #              control = list(maxit = 1e+6, pgtol = 1e-7), hessian = FALSE,
    #              L_current_i = window(L_mean[,,,,, x], end = last_year),
    #              L_c_i = L_c[, ac(last_year),,,, x],
    #              L_inf_i = L_inf[, x], K_i = K[, x],
    #              lower = 0)
    
    ### return parameter
    #return(res$par)
    return(data.frame(year = tail(yrs_used, 1), iter = x, data = res$par))
    
  })
  
  ### convert into FLQuant
  #Z_iter <- FLQuant(Z_iter, dimnames = list(year = year_z, iter = iters))
  #Z_iter <- do.call(rbind, Z_iter)
    
  # })
  
  Z <- do.call(rbind, Z_iter)
  Z$iter <- as.numeric(as.character(Z$iter))
  Z <- as(Z, "FLQuant")
  
  ### mean over years
  Z <- yearMeans(Z)
  #dimnames(Z)$year <- last_year
  

  ### calculate f
  if (isTRUE(perfect_knowledge)) {
    f <- c(refpts["F0.1"] / (Z - lhpar["M"]))
  } else {
    f <- c(refpts["F0.1YPR"] / (Z - lhpar["M"]))
  }
  
  ### coerce into FLQuant
  f <- FLQuant(f, dimnames = list(year = tail(dimnames(L_mean)$year, 1),
                                  iter = dimnames(L_mean)$iter))
  
  return(f)
  
}
### ------------------------------------------------------------------------ ###
### f, option 0
### ------------------------------------------------------------------------ ###
### return 1
wklife_3.2.1_f_0 <- function(...) {
  return(NA)
}

### ------------------------------------------------------------------------ ###
### options for b
### ------------------------------------------------------------------------ ###
### ~ proxy for ratio current stock size/MSYBtrigger
### ------------------------------------------------------------------------ ###
### b: option a
### ------------------------------------------------------------------------ ###
### b = min(1, I_current / I_trigger)
### with I_trigger = w * I_lim where w >= 1, default = 1.4
### use lowest observed value as default for I_lim
wklife_3.2.1_b_a <- function(idx, b_w = 1.4, refpts = NULL, tracking, 
                             b_z = NULL, ...) {
  
  ### extract index (use only first element of FLIndices)
  idx_temp <- quantSums(index(idx[[1]]))
  
  ### calculate B_trigger, if not provided
  ### use min observed index value * w
  if (is.null(refpts)) {
    I_trigger <- apply(X = idx_temp, MARGIN = 6, FUN = min, na.rm = TRUE) * b_w
  } else if (all(is.na(refpts["I_trigger"]))) {
    I_trigger <- refpts["I_lim"] * b_w
  } else {
    I_trigger <- refpts["I_trigger"]
  }
  
  ### calculate ratio I_current / I_lim
  I_ratio <- idx_temp[, tail(dimnames(idx_temp)$year, 1)] / I_trigger
  
  ### calculate b
  factor_b <- apply(X = I_ratio, MARGIN = 6, FUN = function(x) {
    min(1, x)
  })
  
  ### modify b, if requested
  if (!is.null(b_z)) {
    factor_b <- factor_b^b_z
  }
  
  return(factor_b)
  
}
### ------------------------------------------------------------------------ ###
### b: option b
### ------------------------------------------------------------------------ ###
### ~ PA buffer
wklife_3.2.1_b_b <- function(buffer_size = 20, buffer_interval = 4, tracking,
                             idx, ...) {
  
  ### get last data year
  lst_yr <- range(idx)[["maxyear"]]
  
  ### convert buffer_size into factor b
  b_size <- 1 - buffer_size/100
  
  ### find years/position of last application and apply buffer, if required
  factor_b <- apply(X = window(tracking["HCR3.2.1b"], end = lst_yr), 
                    MARGIN = 6, FUN = function(x) {
    
    ### find non-NAs
    nas <- which(!is.na(x))
    
    ### if only NAs, return TRUE
    if (length(nas) == 0) {
      return(b_size)
    ### otherwise find interval since last application
    } else {
      ### return buffer, if time since last application >= requested interval
      if ((length(x) - max(nas)) >= buffer_interval) {
        return(b_size)
      ### otherwise return NA
      } else {
        return(NA)
      }
    }
      
    
  })
  
  return(factor_b)
  
}
### ------------------------------------------------------------------------ ###
### b: option 0
### ------------------------------------------------------------------------ ###
### return NA
wklife_3.2.1_b_0 <- function(...) {
  return(NA)
}
### ------------------------------------------------------------------------ ###
### calculate length at first capture
### ------------------------------------------------------------------------ ###
### retrieve first length where abundance/catch is above 50% of max 
calc_Lc <- function(object ### FLQuant object with lengths
){
  
  ### apply over entire FLQuant object
  L_c <- apply(object, MARGIN = c(2, 6), function(x){
    
    ### return NA if only NAs in data
    if (all(is.na(x))) {
      return(NA)
    }
    
    ### length with max catch
    N_max <- max(x, na.rm = TRUE)
    ### length positions with catch > N_max/2
    N_max_2 <- which(x >= N_max/2)
    ### min position in there
    N_max_2_min <- head(N_max_2, 1)
    ### get corresponding length
    L_c <- an(dimnames(object)[[1]][N_max_2_min])
    ### return value
    return(L_c)
    
  })
  
  return(L_c)
  
}

### ------------------------------------------------------------------------ ###
### calculate mean length/age
### ------------------------------------------------------------------------ ###
### mean length/age, weighted by abundance at length/age

calc_mean <- function(object, ### FLQuant with lengths/ages
                      min = NULL, ### minimum class to be used
                      include_min = FALSE ### include min or use only above
){
  
  ### if min is numeric, coerce into FLQuant
  if (is.numeric(min) & !is(min, "FLQuant")) {
    
    min <- FLQuant(min, dimnames = list(year = dimnames(object)$year, 
                                        iter = dimnames(object)$iter))
    
  }
  
  ### remove lengths below min
  if (!is.null(min)) {
    
    ### get available length classes
    lclass <- dimnames(object)$length
    
    ### go through object and min
    ### loop through years
    for (year_temp in dimnames(object)$year) {
      
      ### go to next year, if no min value available
      if (all(is.na(min[, year_temp]))) {
        next()
      }
      
      ### loop through iters
      for (iter_temp in dimnames(object)$iter) {
        
        ### get min value
        min_temp <- c(min[, year_temp,,,, iter_temp])
        
        ### go to next iteration, if no min value available
        if (is.na(min_temp)) {
          next()
        }
        
        ### replace numbers with NAs
        ### keep all lengths above min
        if (!isTRUE(include_min)) {
          object[lclass[an(lclass) <= min_temp], year_temp,,,, iter_temp] <- NA
        } else {
          ### keep lengths above and including min
          object[lclass[an(lclass) < min_temp], year_temp,,,, iter_temp] <- NA
        }
        
      }
      
    }
    
  }
  
  ### calculate mean length/age over all years and iters
  apply(X = object, MARGIN = c(2, 6), FUN = function(x) {
    
    ### calculate
    res <- weighted.mean(x = an(dimnames(x)$length), 
                         w = ifelse(is.na(x), 0, x), na.rm = TRUE)
    
    ### check if result obtained
    ### if all cathc at all lengths = 0, return 0 as mean length
    if (is.nan(res)) {
      if (all(ifelse(is.na(x), 0, x) == 0)) {
        res[] <- 0
      }
    }
    
    return(res)
    
  })
  
  
}

### ------------------------------------------------------------------------ ###
### simple YPR model
### ------------------------------------------------------------------------ ###

### definition of the YPR function
YPR <- function(a, b, ### age length conversion
                tm = 0, ### age at first/50% maturity
                Linf, K, t0, ### von Bertalanffy growth parameters
                tc, ### age at first capture
                F, ### fishing mortality
                M, ### natural mortality
                max_age, ### max age
                n_start = 1, ### number of recruits at age 0
                report = FALSE ### report YPR value or full table
){
  
  ### template
  YPR_df <- data.frame(age = 0:max_age,
                       length = NA,
                       weight = NA,
                       maturity = NA,
                       selectivity = NA,
                       survivors = NA,
                       ssb = NA,
                       yield_no = NA,
                       yield_bio = NA)
  
  ### mean length
  ### calculate with von Bertalanffy growth function
  YPR_df$length <- Linf * (1 - exp(-K * (YPR_df$age - t0)))
  
  ### mean weight
  ### calculate with W = a*L^b
  YPR_df$weight <- a * YPR_df$length ^ b
  
  ### maturity
  ### knife edged maturity
  YPR_df$maturity <- ifelse(tm > YPR_df$age, 0, 1)
  
  ### knife edge selectivity
  YPR_df$selectivity <- ifelse(tc > YPR_df$age, 0, 1)
  
  ### calculate surivors at age
  ### for ages < tc, only natural mortality occurs
  survivors <- n_start * exp(-M * YPR_df$age)
  ### fishing starts at age tc
  survivors <- survivors[YPR_df$age == (tc)] * 
    exp(-(M + YPR_df$selectivity * F) * (YPR_df$age - tc))
  ### save
  YPR_df$survivors <- survivors
  
  ### spawning biomass
  YPR_df$ssb <- YPR_df$maturity * YPR_df$weight * YPR_df$survivors
  
  ### yield in numbers
  YPR_df$yield_no <- (F * YPR_df$selectivity / (F * YPR_df$selectivity + M)) * 
    YPR_df$survivors * (1 - exp(-YPR_df$selectivity * F - M))
  
  ### yield biomass
  YPR_df$yield_bio <- YPR_df$yield_no * YPR_df$weight
  
  ### calculate YPR
  YPR <- sum(YPR_df$yield_bio)
  
  ### return result
  if (!isTRUE(report)) {
    return(YPR)
  } else {
    return(YPR_df)
  }
  
}

### ------------------------------------------------------------------------ ###
### Gedamke and Hoenig model
### ------------------------------------------------------------------------ ###
### without different mortality periods
### code adapted from ICES methods: https://github.com/ices-tools-dev/ICES_MSY

GH <- function(Z, L_current_i, L_c_i, L_inf_i, K_i){
  
  ### predict L with given parameters
  L <- L_inf_i * (1 - ((1 - exp(-Z)) / (1 - exp(-Z - K_i))) * (1 - L_c_i/L_inf_i))
  
  ### remove NAs
  L_current_i <- L_current_i[!is.na(L_current_i)]
  
  ### difference to mean length
  L_diff <- L_current_i - c(L)
  
  ### sigma
  sigma <- sqrt(sum(1 * L_diff^2) / length(L_current_i))
  
  ### log-likelihood
  LL <- -length(L_current_i) * log(sigma) - sum(1 * L_diff^2) / 
    (2 * sigma^2)
  
  
  ### mean squared difference
  #diff <- sqrt(sum(L_diff^2) / length(L_diff))
  
  ### LL ???
  #LL <- -length(L_diff) * log(diff) - sum(diff^2) / (2 * diff^2)
  
  return(-LL)
  
}

### ------------------------------------------------------------------------ ###
### SPiCT assessment and forecast
### ------------------------------------------------------------------------ ###
#stk_bckp <- stk
wklife_3.1_f <- function(stk, tracking, idx, interval, start_year, ...){

  ### available years
  data_years <- dims(stk)$minyear:dims(stk)$maxyear
  
  ### check if assessment required
  if (((tail(data_years, 1) + 1) - start_year) %% interval == 0) {
    
    ### create biomass index
    idx_biomass <- quantSums(index(idx$idx))
    
    ### -------------------------------------------------------------------- ###
    ### SPiCT in loop ####
    ### -------------------------------------------------------------------- ###
    res <- lapply(1:dim(stk)[6], function(iter_i) {
      
      
      ### create input object for spict
      input <- list(obsC  = c(FLCore::iter(catch(stk), iter_i)),
                    timeC = data_years,
                    obsI  = c(FLCore::iter(idx_biomass, iter_i)),
                    timeI = data_years,
                    dteuler = 1/16,
                    manstart = tail(data_years, 1) + 2, ### start of management
                    timepredc = tail(data_years, 1) + 2,
                    timepredi = tail(data_years, 1) + 2,
                    dtpredc = 1,
                    ffac = 1,
              ### turn of covariance reporting, saves memory and computing time!
                    getReportCovariance = FALSE
      )
      
      ### ------------------------------------------------------------------ ###
      ### fit model
      ### ------------------------------------------------------------------ ###
      
      ### fit model
      ### printout to screen is captured
      ### warnings and error message are ignored
      fitted <- NULL ### default state of fitted if spict does not work
      . <- capture.output(try(fitted <- fit.spict(input)))
      
      ### define default return if model failed
      return_vals <- c(converged = 1, b = NA, f = NA, bmsy = NA, fmsy = NA,
                       advice = NA)
      
      ### stop if model failed entirely
      if (is.null(fitted)) return(return_vals)
      
      ### -------------------------------------------------------------- ###
      ### continue: perform SPiCT forecast
      ### -------------------------------------------------------------- ###
      
      ### select forecast year pointers
      maninds <- which(fitted$inp$time >= fitted$inp$manstart)
      
      ### create new input object
      input_temp <- fitted$inp[c("dteuler", "timeC", "obsC", "timeI", "obsI", 
                                 "timepredc", "timepredi", "manstart",
                                 "getReportCovariance")]
      
      ### get Fmsy
      Fmsy <- get.par('logFmsy', fitted, exp = TRUE)[, "est"]
      ### get last F value
      Flast <- get.par('logF', fitted, exp = TRUE)[(fitted$inp$indpred[1]-1), "est"]
      ### get last estimated biomass
      Blast <- get.par('logB', fitted, exp = TRUE)[(fitted$inp$indpred[1]-1), "est"]
      ### get Bmsy
      Bmsy <- get.par('logBmsy', fitted, exp = TRUE)[, "est"]
      ### MSYBtrigger = 1/2 Bmsy
      MSYBtrigger <- Bmsy/2
      
      ### calculate target F
      ### Fmsy as factor of current F, 
      ### corrected by B/MSYBtrigger (max 1)
      f_target <- (Fmsy / Flast) * min(1, Blast/MSYBtrigger)
      
      ### stop forecast if target does not exist
      ### convergence value 2: missing target
      if (isTRUE(is.na(f_target)) | is.null(f_target) | length(f_target) == 0) {
        return_vals[[1]] <- 2
        return(return_vals)
      }
      
      ### forecast
      fitted_fwd <- NULL
      try(fitted_fwd <- (prop.F(fac = f_target, inpin = input_temp, repin = fitted, 
                       maninds = maninds, dbg = 0)))
      
      ### stop if forecast failed
      ### convergence value 3: forecast failed
      if (is.null(fitted_fwd)) {
        return_vals[[1]] <- 3
        return(return_vals)
      }
      
      ### ------------------------------------------------------------------ ###
      ### return results
      ### ------------------------------------------------------------------ ### 
      
      return(c(converged = fitted$opt$convergence,
               b = sumspict.states(fitted)[1, "estimate"],
               f = sumspict.states(fitted)[2, "estimate"],
               bmsy = sumspict.srefpoints(fitted)[1, "estimate"],
               fmsy = sumspict.srefpoints(fitted)[2, "estimate"],
               advice = sumspict.predictions(fitted_fwd)[5, "prediction"] ))
      
    })
    
    
    ### -------------------------------------------------------------------- ###
    ### extract results from model fitting as vectors and save in summary ####
    ### -------------------------------------------------------------------- ###
    
    for (iter_i in 1:dim(stk)[6]) {
      tracking[c("convergence", "spict_b", "spict_f", "spict_bmsy", "spict_fmsy",
                 "advice"), ac(tail(data_years, 1)),,,, iter_i] <- res[[iter_i]]
    }
    
    return(list(stk = stk, tracking = tracking))
    
  ### otherwise copy advice from last year
  } else {
    
    tracking["advice", ac(tail(data_years, 1))] <- c(tracking["advice", ac(tail(data_years, 1) - 1)])
    
    return(list(stk = stk, tracking = tracking))
    
  }

}


################################################################################
### Finlay's approach

spict.wrapper <- function(stk, idx, ...){
  args <- list(...)
  idx_data_years <- range(idx[[1]])["minyear"]:range(idx[[1]])["maxyear"]
  catch_data_years <- range(stk)["minyear"]:range(stk)["maxyear"]
  # index = n * q
  # biomass index = sum(index * wt)
  biomass_index <- quantSums(index(idx[[1]]))
  if(!(is.null(args$dteuler))){
    dteuler <- args$dteuler
  }
  else {
    dteuler <- 1/16
  }
  niters <- dim(catch(stk))[6]
  spict_fit <- vector("list", length = niters)
  # Loop over iterations
  for (iter in 1:niters){
    # Organise data for SPiCT
    spict_input <- list(obsC = c(iter(catch(stk),iter)),
                        timeC = catch_data_years,
                        obsI = c(iter(biomass_index,iter)),
                        timeI = idx_data_years,
                        dteuler = dteuler,
                        timepredc = max(catch_data_years)+2,
                        dtpredc = 1,
                        timepredi =  max(catch_data_years)+2,
                        manstart = max(catch_data_years)+2,
                        ffac <- 1,
                        getReportCovariance = FALSE
    )
    # What do we do if the fit does not work
    
    ### fit model
    ### printout to screen is captured
    ### warnings and error message are ignored
    fitted <- NULL ### default state of fitted if spict does not work
    . <- capture.output(try(fitted <- fit.spict(spict_input)))
    #try(fitted <- fit.spict(spict_input))
    
    # If fit is OK store predicted catch, estimated biomass and the entire fit
    if(!inherits(fitted, "try-error")) {
      FLCore::iter(stk@catch,iter) <- FLQuant(get.par('logCpred', fitted, exp = TRUE)[1:(length(catch_data_years)-1),'est'], dimnames=list(year=catch_data_years))
      FLCore::iter(stk@stock,iter) <- FLQuant(get.par('logB', fitted, exp = TRUE)[ac(catch_data_years),'est'], dimnames=list(year=catch_data_years))
      spict_fit[[iter]] <- fitted
    }
    tracking["convergence",ac(range(stk)["maxyear"]+1),1,1,1,iter] <- fitted$opt$convergence
  }
  attr(stk, "spict_fit") <- spict_fit
  return(list(stk=stk, tracking=tracking))
}




# sep.wrapper <- function(stk, idx, control){
#   if(!"fmod" %in% names(control)){
#     print("Using default fmod sca settings")
#     fmod <- getFmod(stk, dfm=c(2/3, 1/2), model='separable')
#   } else {
#     fmod <- control$fmod
#   }
#   if(!"qmod" %in% names(control)){
#     print("Using default qmod sca settings")
#     qmod <- getQmod(idx)
#   } else {
#     qmod <- control$qmod
#   }
#   fit <- sca(stk, idx, fmodel=fmod, qmodel=qmod)
#   ## convergence diagnostic (quick and dirty)
#   maxgrad <- fit@fitSumm["maxgrad",]
#   stk <- stk + fit
#   return(list(stk = stk, converge = maxgrad))
# }
# 
# ma.wrapper <- function(stk, idx, xsa.control=list(), sca.control=list()){
#   
#   # XSA
#   if(is.null(xsa.control)){
#     print("Using default xsa control settings")
#     xsa.control  <- FLXSA.control()
#   }
#   # Fit XSA
#   fit <- FLXSA(stk, idx, xsa.control)
#   # convergence diagnostic (quick and dirty)
#   maxit <- fit@control@maxit
#   # Update stk0
#   stk.xsa <- stk+fit
# 
#   # SCA
#   fit <- sca(stk, idx)
#   stk.sca1 <- stk+fit  
# 
#   # SCA
#   if(!"fmod" %in% names(sca.control)){
#     print("Using default fmod sca settings")
#     fmod <- getFmod(stk, dfm=c(2/3, 1/2))
#   } else {
#     fmod <- control$fmod
#   }
#   if(!"qmod" %in% names(sca.control)){
#     print("Using default qmod sca settings")
#     qmod <- getQmod(idx)
#   } else {
#     qmod <- control$qmod
#   }
#   fit <- sca(stk, idx, fmodel=fmod, qmodel=qmod)
#   stk.sca <- stk+fit  
# 
#   # SCA separable
#   if(!"fmod" %in% names(sca.control)){
#     print("Using default fmod sca settings")
#     fmod <- getFmod(stk, dfm=c(2/3, 1/2), model='separable')
#   } else {
#     fmod <- control$fmod
#   }
#   if(!"qmod" %in% names(sca.control)){
#     print("Using default qmod sca settings")
#     qmod <- getQmod(idx)
#   } else {
#     qmod <- control$qmod
#   }
#   fit <- sca(stk, idx, fmodel=fmod, qmodel=qmod)
#   stk.sep <- stk+fit  
#   
#   catch.n(stk) <- (catch.n(stk.sca1)+catch.n(stk.sca)+catch.n(stk.xsa)+catch.n(stk.sep))/4
#   stock.n(stk) <- (stock.n(stk.sca1)+stock.n(stk.sca)+stock.n(stk.xsa)+stock.n(stk.sep))/4
#   harvest(stk) <- (harvest(stk.sca1)+harvest(stk.sca)+harvest(stk.xsa)+harvest(stk.sep))/4
#   list(stk=stk, convergence=NA)
# }
# 
# 
# ivw.wrapper <- function(stk, idx, xsa.control=list(), sca.control=list()){
#   
# 
#   # SCA
#   if(!"fmod" %in% names(sca.control)){
#     fmod <- getFmod(stk, dfm=c(2/3, 1/2))
#   } else {
#     fmod <- control$fmod
#   }
#   if(!"qmod" %in% names(sca.control)){
#     qmod <- getQmod(idx)
#   } else {
#     qmod <- control$qmod
#   }
#   fit <- sca(stk, idx, fmodel=fmod, qmodel=qmod)
#   stk.sca <- stk+fit  
# 	v1 <- sum((residuals(fit, stk, idx)$catch.n)^2)
# 
#   # SCA separable
#   if(!"fmod" %in% names(sca.control)){
#     print("Using default fmod sca settings")
#     fmod <- getFmod(stk, dfm=c(2/3, 1/2), model='separable')
#   } else {
#     fmod <- control$fmod
#   }
#   if(!"qmod" %in% names(sca.control)){
#     print("Using default qmod sca settings")
#     qmod <- getQmod(idx)
#   } else {
#     qmod <- control$qmod
#   }
#   fit <- sca(stk, idx, fmodel=fmod, qmodel=qmod)
#   stk.sep <- stk+fit  
# 	v2 <- sum((residuals(fit, stk, idx)$catch.n)^2)
#   
#   w1 <- (1/v1)/(1/v1+1/v2)
#   
#   catch.n(stk) <- w1*catch.n(stk.sca)+(1-w1)*catch.n(stk.sep)
#   stock.n(stk) <- w1*stock.n(stk.sca)+(1-w1)*stock.n(stk.sep)
#   harvest(stk) <- w1*harvest(stk.sca)+(1-w1)*harvest(stk.sep)
#   list(stk=stk, convergence=NA)
# }
# 
# 
# eqw.wrapper <- function(stk, idx, xsa.control=list(), sca.control=list()){
# 
#   # SCA
#   if(!"fmod" %in% names(sca.control)){
#     fmod <- getFmod(stk, dfm=c(2/3, 1/2))
#   } else {
#     fmod <- control$fmod
#   }
#   if(!"qmod" %in% names(sca.control)){
#     qmod <- getQmod(idx)
#   } else {
#     qmod <- control$qmod
#   }
#   fit <- sca(stk, idx, fmodel=fmod, qmodel=qmod)
#   stk.sca <- stk+fit  
# 
#   # SCA separable
#   if(!"fmod" %in% names(sca.control)){
#     print("Using default fmod sca settings")
#     fmod <- getFmod(stk, dfm=c(2/3, 1/2), model='separable')
#   } else {
#     fmod <- control$fmod
#   }
#   if(!"qmod" %in% names(sca.control)){
#     print("Using default qmod sca settings")
#     qmod <- getQmod(idx)
#   } else {
#     qmod <- control$qmod
#   }
#   fit <- sca(stk, idx, fmodel=fmod, qmodel=qmod)
#   stk.sep <- stk+fit  
# 
#   w1 <- 0.5
#   
#   catch.n(stk) <- w1*catch.n(stk.sca)+(1-w1)*catch.n(stk.sep)
#   stock.n(stk) <- w1*stock.n(stk.sca)+(1-w1)*stock.n(stk.sep)
#   harvest(stk) <- w1*harvest(stk.sca)+(1-w1)*harvest(stk.sep)
#   list(stk=stk, convergence=NA)
# }
