### ------------------------------------------------------------------------ ###
### wrapper for creating biomass index and length frequencies
### ------------------------------------------------------------------------ ###
obs_bio_len <- function(...){
  #browser()
  args <- list(...)
  args$method <- NULL
  ### create biomass index
  res <- do.call(idx_bio, args)
  
  ### use only data years for stock data
  args$stk <- res$stk
  
  ### add life-history parameters
  args$lhpar <- attr(args$stk, "lhpar")
  
  ### create catch lengths
  res$stk <- do.call(length_freq, args)
  ### extract history
  # catch_len_hist <- attr(res$stk, "catch_len")
  # ### insert new values
  # catch_len_hist[, dimnames(catch_len)$year] <- catch_len
  # ### cut of at last data year
  # catch_len_hist <- window(catch_len_hist,
  #                          end = an(last(dimnames(catch_len)$year)))
  # ### save
  # attr(res$stk, "catch_len") <- catch_len_hist
  
  return(res)
  
}

### ------------------------------------------------------------------------ ###
### create biomass at age index
### ------------------------------------------------------------------------ ###

idx_bio <- function(stk, observations, genArgs, tracking,
                    lst_catch = -1, ### last catch/idx year, relative to 
                    lst_idx = -1,   ### intermediate year (ay)
                    ssb_idx = NULL, ### use SSB as index
                    ...){
  ay <- genArgs$ay
  ### stock: subset to requested years
  ### default: lst_catch = -1, i.e. data up to the year before ay
  stk0 <- window(stk, end = ay + lst_catch)
  ### assume perfect knowledge of stock characteristics (catch etc.)
  
  ### next line commented, interferes dramatically with catch numbers at age
  ### used in creation of length frequencies
  #catch.n(stk0) <- (catch.n(stk0) + 1) # avoid zeros
  
  ### get indices
  idx <- observations$idx
  
  ### get years for which index calculation is required
  ### always calculate index in year ay
  ### add more years, if requested
  yrs_new <- ay
  if (lst_idx > 0) yrs_new <- seq(from = ay, to = ay + lst_idx)
  
  ### add/calculate index of current year (ay)
  for (idx_count in 1:length(idx)) {
    index(idx[[idx_count]])[, ac(yrs_new)] <- stock.n(stk)[, ac(yrs_new)] *
      index.q(idx[[idx_count]])[, ac(yrs_new)]*stock.wt(stk)[, ac(yrs_new)]
  }
  
  ### if SSB index, reduce index by mortality
  if (isTRUE(ssb_idx)) {
    for (idx_count in 1:length(idx)) {
      index(idx[[idx_count]])[, ac(yrs_new)] <- 
        index(idx[[idx_count]])[, ac(yrs_new)] * 
        exp(-(harvest(stk)[, ac(yrs_new)] * harvest.spwn(stk)[, ac(yrs_new)] +
                m(stk)[, ac(yrs_new)] * m.spwn(stk)[, ac(yrs_new)]))
    }
  }
  
  ### subset all indices to requested years -> observations
  idx0 <- lapply(idx, window, end = ay + lst_idx)
  
  ### return observations
  list(stk = stk0, idx = idx0, observations = list(idx = idx), 
       tracking = tracking)
}


### ------------------------------------------------------------------------ ###
### length frequencies
### ------------------------------------------------------------------------ ###
length_freq <- function(stk, genArgs = list(),
                        lst_catch = -1,
                        lengths = NULL, ### FLQuant with lengths
                        len_src = "catch", ### slot for numbers at age
                        full_series = FALSE, ### years to be used
                        lhpar = FLPar(a = 0.00001, b = 3,     ### length-weight
                                      L_inf = 100), ### relationship
                        len_sd = 1, ### standard deviation for spreading lngth
                        len_sd_cut = 2, ### cut of length spread
                        len_dist = "normal", ### normal or lognormal
                        len_noise_sd = 0, ### observation error
                        ...) {
  
  ay <- genArgs$ay
  ### available years
  yrs_avail <- range(stk)[["minyear"]]:range(stk)[["maxyear"]]
  
  ### use all years, if not otherwise specified
  if (isTRUE(full_series)) {
    yrs_new <- yrs_avail
  } else {
    ### otherwise use last data year only
    yrs_new <- ay + lst_catch 
  }
  
  ### extract length-weight parameters
  if ("lhpar" %in% names(attributes(stk))) lhpar <- attr(stk, "lhpar")
  a <- lhpar["a"]
  b <- lhpar["b"]
  L_inf <- round(mean(lhpar["L_inf"]))
  ### replicate for correct multiplication
  a <- rep(c(a), each = dims(stk)$age * length(yrs_new))
  b <- rep(c(b), each = dims(stk)$age * length(yrs_new))
  
  ### create FLStockLen template
  ### with all years from stk, but only 1 iteration to save memory,
  ### gets extended later, where needed
  stk_len <- FLStockLen(FLQuant(NA, dimnames = list(length = 1:L_inf, 
                                                    year = yrs_avail)))
  ### set harvest unit
  units(stk_len)$harvest <- "f"
  ### insert catch
  catch_temp <- catch(stk)
  names(dimnames(catch_temp))[1] <- "length"
  catch(stk_len) <- catch_temp
  
  ### extend iter dimension for catch lengths
  catch.n(stk_len) <- propagate(catch.n(stk_len), dims(stk)$iter)
  
  ### now create the length frequencies
  
  ### extract numbers
  numbers <- as.data.frame(get(paste0(len_src, ".n"))(stk)[, ac(yrs_new)])
  names(numbers)[names(numbers) == "data"] <- "numbers"
  
  ### calculate length from age with a & b
  weights <- get(paste0(len_src, ".wt"))(stk)[,ac(yrs_new)]
  lengths <- as.data.frame((weights / c(a))^(1/c(b)))
  names(lengths)[names(lengths) == "data"] <- "length"
  
  ### merge
  res <- merge(numbers, lengths)
  
  ### get available lengths
  lengths_avail <- unique(lengths$length)
  ### remove NAs
  lengths_avail <- lengths_avail[!is.na(lengths_avail)]
  
  ### create list with spread lengths
  spread_lst <- lapply(lengths_avail, function(x){
    
    ### get new spread lengths
    new_lengths <- seq(from = round(x - len_sd_cut*len_sd),
                       to = round(x + len_sd_cut*len_sd))
    ### keep only positive lengths and < max length
    new_lengths <- new_lengths[new_lengths > 0 &
                                 new_lengths <= dims(stk_len)$max]
    
    ### if only lengths above max length, use max length as plusgroup
    ### occurs only if stock is at or close to virgin biomass
    if (length(new_lengths) == 0) {
      new_lengths <- dims(stk_len)$max
    }
    
    ### calculate probabilities at length
    ### normal distribution
    if (len_dist == "normal") {
      probs <- dnorm(x = new_lengths, mean = x, sd = len_sd)
    } else if (len_dist == "lognormal") {
      probs <- dlnorm(x = new_lengths, meanlog = log(x), sd = log(len_sd))
    } else if (len_dist == "uniform") {
      probs <- dunif(x = new_lengths, min = min(new_lengths),
                     max = max(new_lengths))
    } else stop("unknown distribution")
    
    ### re-normalize
    probs <- probs / sum(probs)
    
    ### workaround to avoid zero probabilites for all lengths
    ### can happen if sd=0 or very low and rounded length is outside
    ### distribution
    if (sum(probs, na.rm = TRUE) == 0) {
      probs <- 1 / length(probs)
    }
    
    ### return spread lengths
    cbind(length = x, new_lengths, probs)
    
  })
  ### combine
  #browser()
  spread_lst <- do.call(rbind, spread_lst)
  
  ### merge with data
  res <- merge(x = res, y = spread_lst, by = "length", all = TRUE)
  
  ### calculate new numbers at length by spreading them
  res$numbers <- res$numbers * res$probs
  
  ### overwrite old lengths with new ones
  res$length <- res$new_lengths
  
  ### remove columns not required anymore
  res <- res[, !colnames(res) %in% c("new_lengths", "probs")]
  
  ### rename numbers
  names(res)[names(res) == "numbers"] <- "data"
  
  ### aggregate per length 
  ### use data.table for efficiency
  res <- data.table(res)
  res <- res[, list(data = sum(data)), by = c("year","iter","length")]
  
  ### add noise
  res$data <- res$data * rlnorm(n = length(res$data), sdlog = len_noise_sd)
  
  ### add missing length classes
  res <- merge(x = data.frame(length = 1:L_inf),
               y = res, all = TRUE)
  
  ### coerce iter into numeric
  res$iter <- as.numeric(as.character(res$iter))
  
  ### convert into FLQuant
  res_FLQuant <- as(res, "FLQuant")
  
  ### insert into stk_len
  catch.n(stk_len)[, ac(yrs_new)] <- res_FLQuant
  
  ### insert historical length frequencies, if they exist
  if ("catch_len" %in% names(attributes(stk))) {
    ### find years in both objects
    yrs_hist <- intersect(dimnames(attr(stk, "catch_len"))$year, 
                          dimnames(catch.n(stk_len))$year)
    lengths_insert <- intersect(dimnames(attr(stk, "catch_len"))$length, 
                                dimnames(catch.n(stk_len))$length)
    catch.n(stk_len)[lengths_insert, yrs_hist] <- attr(stk, "catch_len")[lengths_insert, yrs_hist]
  }
  
  # ### cut off length frequencies at last data year
  # stk_len@catch.n <- window(stk_len@catch.n, end = ay + lst_lngth)
  
  ### save refpts and lhpar, if they exist
  if ("refpts" %in% names(attributes(stk))) {
    attr(stk_len, "refpts") <- attr(stk, "refpts")
  }
  if ("lhpar" %in% names(attributes(stk))) {
    attr(stk_len, "lhpar") <- attr(stk, "lhpar")
  }
  
  return(window(stk_len, end = tail(yrs_new, 1)))
  
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
    if (b_z != 1) {
      factor_b <- factor_b^b_z
    }
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


### ------------------------------------------------------------------------ ###
### x()
### ------------------------------------------------------------------------ ###
### HCR 3.2.1 from WKMSYCat34 report
### so far handles only current catch and does not change anything else

wklife_3.2.1_x <- function(stk, tracking, genArgs, n_catch_yrs, interval = 1,
                           start_year, multiplier = NA, ...){
  
  ay <- genArgs$ay
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
### h()
### ------------------------------------------------------------------------ ###

wklife_3.2.1_h <- function(hcrpars, genArgs, tracking, stk, 
                           upper_constraint = Inf, lower_constraint = 0, ...) {
  
  ay <- genArgs$ay
  
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
  ctrl <- getCtrl(values = c(adv), quantity = "catch", years = ay + 1, 
                  it = dim(adv)[6])
  
  tracking["advice", ac(ay)] <- c(adv)
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### k() ####
### ------------------------------------------------------------------------ ###
wklife_3.2.1_k <- function(ctrl, tracking, genArgs, ...) {
  
  ay <- genArgs$ay
  
  tracking["Implementation", ac(ay)] <- ctrl@trgtArray[ac(ay + 1),"val", ]
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}




### ------------------------------------------------------------------------ ###
### IEM ####
### ------------------------------------------------------------------------ ###

### add random noise to targets
noise.wrapper <- function(ctrl, genArgs, fun = "rlnorm", mean = 0, sd = 0.1, 
                          multiplicative = TRUE, tracking, ...) {
  
  ay <- genArgs$ay
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
  
  tracking["IEM", ac(ay)] <- ctrl@trgtArray[ac(ay + 1), "val", ]
  
  ### return control object and tracking
  list(ctrl = ctrl, tracking = tracking)
  
}


### ------------------------------------------------------------------------ ###
### j()
### ------------------------------------------------------------------------ ###

### function for avoiding the issue that the catch is higher than
### the targeted catch
### can happen due to bug in FLash if >1 iteration provided
### sometimes, FLash struggles to get estimates and then uses F estimate from
### previous iteration
### workaround: target same value several times and force FLash to try again
ctrl_catch_workaround <- function(ctrl, tracking, genArgs, ...) {
  
  ay <- genArgs$ay
  tracking["FleetDyn", ac(ay)] <- ctrl@trgtArray[ac(ay + 1), "val", ]
  
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

### ------------------------------------------------------------------------ ###
### g() projection
### ------------------------------------------------------------------------ ###
fwd_attr <- function(stk, ctrl,
                      sr, ### stock recruitment model
                      sr.residuals, ### recruitment residuals
                      sr.residuals.mult = TRUE, ### are res multiplicative?
                      maxF = 5, ### maximum allowed Fbar
                      ...) {
  
  ### project forward with FLash::fwd
  stk[] <- fwd(object = stk, control = ctrl, sr = sr, 
               sr.residuals = sr.residuals, 
               sr.residuals.mult = sr.residuals.mult,
               maxF = maxF)
  
  ### return stock
  return(list(object = stk))
  
}




### ------------------------------------------------------------------------ ###
### adapt accessor functions for FLStockLen ####
### ------------------------------------------------------------------------ ###
setMethod("fbar", signature(object = "FLStockLen"), 
  function(object, ...){
    selectMethod("fbar", signature = "FLStock")(object)
})
setMethod("ssb", signature(object = "FLStockLen"), 
  function(object, ...){
    selectMethod("ssb", signature = "FLStock")(object)
})

### ------------------------------------------------------------------------ ###
### iter subset  ####
### ------------------------------------------------------------------------ ###

iter_attr <- function(object, iters, subset_attributes = TRUE) {
  
  ### subset object to iter
  res <- FLCore::iter(object, iters)

  if (isTRUE(subset_attributes)) {

    ### get default attributes of object class
    attr_def <- names(attributes(new(Class = class(object))))
    
    ### get additional attributes
    attr_new <- setdiff(names(attributes(object)), attr_def)
    
    ### subset attributes
    for (attr_i in attr_new) {
      attr(res, attr_i) <- FLCore::iter(attr(res, attr_i), iters)
    }
    
  }
  
  return(res)

}


### default attributes of 
# attr_def <- c("range", "catch", "catch.n", "catch.wt", "discards", "discards.n",
#               "discards.wt", "landings", "landings.n", "landings.wt", "stock",
#               "stock.n", "stock.wt", "m", "mat", "harvest", "harvest.spwn", 
#               "m.spwn", "name", "desc", "class")
# 
# is(stk.om)


### ------------------------------------------------------------------------ ###
### combine FLComp, including FLQuant attributes ####
### ------------------------------------------------------------------------ ###
### modified version from 
### https://github.com/flr/FLCore/blob/15a0240e8a689c06c3d18b0166598533d0f16e1a/R/FLComp.R
# setMethod('combine', signature(x='FLComp', y='FLComp'),
#           function(x, y, ..., check=FALSE) {}
  
# combine_attr <- function(x, y, ..., check = FALSE) {
#   
#   args <- c(list(x, y), list(...))
#   
#   # CHECK input classes match exactly
#   if (length(unique(lapply(args, is))) > 1)
#     stop("combine can only operate on objects of identical class")
#   
#   ds <- lapply(args, dims)
#   
#   # CHECK dimnames but iter
#   if (check) {
#     
#     idi <- names(ds[[1]]) !="iter"
#     
#     # COMPARE dims(x)[-'iter')]
#     if (length(unique(lapply(ds, "[", idi))) > 1)
#       stop("Object dimensions but iter must match")
#   }
#   
#   # CALCULATE iters
#   its <- sum(unlist(lapply(ds, "[[", "iter")))
#   
#   # PROPAGATE object
#   res <- propagate(x[,,,,, 1], its)
#   
#   # KEEP iter dimnames if unique
#   itns <- unlist(lapply(args, function(x) dimnames(x)$iter))
#   
#   # CHECK iter dimanmes are unique
#   if (length(itns) > length(unique(itns)))
#     itns <- ac(seq(1, its))
#   
#   # GET iter limits
#   itE <- cumsum(unlist(lapply(ds, "[", "iter")))
#   itS <- itE - unlist(lapply(ds, "[", "iter")) + 1
#   
#   for (i in seq(length(itS)))
#     res[,,,,, seq(itS[i], itE[i])] <- args[[i]]
#   
#   dimnames(res) <- list(iter = itns)
#   
#   ### combine additional FLQuant attributes
#   ### find additional attributes
#   attr_add <- setdiff(names(attributes(res)), 
#                       c(slotNames(res), "class"))
#   if (length(attr_add) > 0) {
#     
#     ### FLQuants
#     attr_FLQuant <- attr_add[sapply(lapply(attr_add, function(x) attr(res, x)), 
#                                     is, "FLQuant")]
#     if (length(attr_FLQuant) > 0) {
#       ### go through FLQuant attributes
#       for (attr_i in attr_FLQuant) {
#         ### expand dimensions
#         attr(res, attr_i) <- propagate(attr(res, attr_i)[,,,,, 1], its)
#         ### fill
#         for (i in seq(length(itS))) {
#           attr(res, attr_i)[,,,,, seq(itS[i], itE[i])] <- attr(args[[i]], attr_i)
#         }
#         ### dimnames iter
#         dimnames(attr(res, attr_i))$iter <- itns
#       }
#     }
#     ### FLPars
#     attr_FLPar <- attr_add[sapply(lapply(attr_add, function(x) attr(res, x)), 
#                                     is, "FLPar")]
#     if (length(attr_FLPar) > 0) {
#       ### go through FLQuant attributes
#       for (attr_i in attr_FLPar) {
#         ### expand dimensions
#         attr(res, attr_i) <- propagate(attr(res, attr_i)[, 1], its)
#         ### fill
#         for (i in seq(length(itS))) {
#           attr(res, attr_i)[, seq(itS[i], itE[i])] <- attr(args[[i]], attr_i)
#         }
#         ### dimnames iter
#         dimnames(attr(res, attr_i))$iter <- itns
#       }
#     }
#     
#   }
#   
#   return(res)
# }


# )

### ------------------------------------------------------------------------ ###
### inter-annual variability ####
### ------------------------------------------------------------------------ ###
#' calculate inter-annual variability of FLQuant
#'
#' This function calculates survey indices from the numbers at age of an 
#' FLStock object
#'
#' @param object Object of class \linkS4class{FLQuant} with values.
#' @param period Select every n-th year, e.g. biennial (optional).
#' @param from,to Optional year range for analysis.
#' @param summary_per_iter Function for summarising per iter. Defaults to mean.
#' @param summary Function for summarising over iter. Defaults to mean.
#' @return An object of class \code{FLQuant} with inter-annual variability.
#'
#' @export
#' 

setGeneric("iav", function(object, period, from, to, summary_per_iter, 
                           summary_year, summary_all) {
  standardGeneric("iav")
})

### object = FLQuant
#' @rdname iav
setMethod(f = "iav",
  signature = signature(object = "FLQuant"),
  definition = function(object, 
                        period, ### periodicity, e.g. use every 2nd value 
                        from, to,### year range
                        summary_per_iter, ### summarise values per iteration
                        summary_year,
                        summary_all) {
  
  ### subset years
  if (!missing(from)) object <- window(object, start = from)
  if (!missing(to)) object <- window(object, end = to)
  
  ### get years in object
  yrs <- dimnames(object)$year
  
  ### select every n-th value, if requested
  if (!missing(period)) {
    yrs <- yrs[seq(from = 1, to = length(yrs), by = period)]
  }
  
  ### reference years
  yrs_ref <- yrs[-length(yrs)]
  ### years to compare
  yrs_comp <- yrs[-1]
  
  ### calculate variation (absolute values, ignore pos/neg)
  res <- abs(1 - object[, yrs_comp] / object[, yrs_ref])
  
  ### replace Inf with NA (compared to 0 catch)
  res <- ifelse(is.finite(res), res, NA)
  
  ### summarise per iteration
  if (!missing(summary_per_iter)) {
    res <- apply(res, 6, summary_per_iter, na.rm = TRUE)
  }
  
  ### summarise per year
  if (!missing(summary_year)) {
    res <- apply(res, 1:5, summary_year, na.rm = TRUE)
  }
  
  ### summarise over everything
  if (!missing(summary_all)) {
    
    res <- summary_all(res, na.rm = TRUE)
    
  }
  
  return(res)

})




### ------------------------------------------------------------------------ ###
### more functions from funs.R ####
### ------------------------------------------------------------------------ ###

#--------------------------------------------------------------------
setGeneric("rz", function(object, ...) standardGeneric("rz"))
setMethod("rz", "FLIndex", function(object){
  flq <- index(object)
  flq[flq==0] <- min(flq[flq!=0])
  index(object) <- flq
  object
})

setMethod("rz", "FLIndices", function(object){
  lapply(object, rz)
})

#--------------------------------------------------------------------
yearCumsum <- function(x, ...){
  x[] <- apply(x, c(1,3:6), cumsum)
  return(x)
}

#--------------------------------------------------------------------
yearDiffPerc <- function(x, ...){
  #x[,-1] <- x[,-1]/x[,-ncol(x)]-1
  x[,-1] <- x[,-1]/c(x[,1])-1
  x[,1] <- 0
  return(x)
}

#--------------------------------------------------------------------
as.table.FLQuants <- function(x){
  x <- mcf(x)
  df0 <- as.data.frame(do.call("cbind", lapply(x, c)))
  row.names(df0) <- dimnames(x)[[2]]
  df0
}

#--------------------------------------------------------------------
iterQuantiles <- function(x, prob=0.5, ...){
  return(apply(x, c(1:5), quantile, prob=prob, na.rm = FALSE))
}

#--------------------------------------------------------------------
an <- function(x, ...) as.numeric(x, ...)

cbind.FLQuant <- function(x, y){
  lst <- mcf(list(x, y))
  res <- lst[[1]]
  res[,dimnames(y)[[2]]] <- y
  res
}

#--------------------------------------------------------------------
getFmod <- function(stk, model="interaction", dfm=c(1/2, 1/2)){
  dis <- dims(stk)
  KY=floor(dfm[1] * dis$year)
  KA=ceiling(dfm[2] *dis$age)
  if(model=="interaction"){
    frm <- substitute(~ te(age, year, k=c(KA, KY)) + s(age, k=KA) + s(year, k=KY), list(KA=KA, KY=KY))
  } else if (model=="separable"){
    frm <- substitute(~ s(age, k = KA) + s(year, k = KY), list(KA=KA, KY=KY))
  }
  as.formula(frm)
}

#--------------------------------------------------------------------
getQmod <- function(idx, dfm=3/4){
  lds <- lapply(idx, dims)
  lds <- lapply(lds, function(x){
    if(x$age<=3){
      frm <- ~factor(age)
    } else {
      frm <- substitute(~s(age, k=KA), list(KA=ceiling(dfm * x$age)))
    }
    as.formula(frm)
  })
  lds
}

#--------------------------------------------------------------------
getCtrl <- function(values, quantity, years, it){
  dnms <- list(iter=1:it, year=years, c("min", "val", "max"))
  arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
  arr0[,,"val"] <- unlist(values)
  arr0 <- aperm(arr0, c(2,3,1))
  ctrl <- fwdControl(data.frame(year=years, quantity=quantity, val=NA))
  ctrl@trgtArray <- arr0
  ctrl
}

#--------------------------------------------------------------------

#iem <- function(iem, ctrl){
#	iem0 <- iem
#  # strip out fun and multiplicative from list of parameters,
#  # to use do.call() later
#	iem0$fun <- iem0$multiplicative <- NULL
#  # eh? number of non NA values in target. But if there are NAs then we need to
#  # know their position for the *ctrl@trgtArray later
#	iem0$n <- sum(!is.na(ctrl@trgtArray))
#	if(iem$multiplicative){
#		ctrl@trgtArray <- do.call(iem$fun, iem0)*ctrl@trgtArray
#	} else {
#		ctrl@trgtArray <- do.call(iem$fun, iem0) + ctrl@trgtArray
#	}
#	ctrl
#}


#--------------------------------------------------------------------
getF <- function(trgt, fc, mxy, ay, object, multiplicative=TRUE, correction=TRUE){
  # F trajectory to Ftrg
  # note f is estimated for ay-1
  # ay is assessment year, also the intermediate year. For this year the
  # decisions were takenon ay-1 but we don't have observations of its
  # implementation
  # ry is the decrease needed in ay
  # ry1 is the decrease needed in ay+1
  # object has to be a FLQuant
  
  # work in medians and add variability in the end
  object0 <- object
  object <- iterMedians(object)
  fc <- iterMedians(fc)
  
  if(multiplicative){
    ry1 <- (trgt/object[,ac(ay)])^(1/(mxy-ay))
    object[,ac(ay+1)] <- ifelse((ay+1) <= mxy, object[,ac(ay)]*ry1, trgt)
    if(correction){
      cny <- mxy-ay+1
      ry <- ifelse((ay+1) < mxy, (trgt/fc)^(1/(cny)), 1)
      object[,ac(ay+1)] <- object[,ac(ay+1)]*(object[,ac(ay+1)]/(ry*fc))
    }
  } else {
    stop("Only additive implemented\n")
    #  ry <- (trgt-fc)/(mxy-ay+1)
    #  ry1 <- 2*ry
    #  object[,ac(ay+1)] <- fc+ry1
    #  if(correction){
    #      object[,ac(ay+1)] <- object[,ac(ay+1)]-(object[,ac(ay)]-(fc+ry))
    #  }
  }
  object0[,ac(ay+1)] <- object[,ac(ay+1)]/object[,ac(ay)]*object0[,ac(ay)]
  
  object0
}


#--------------------------------------------------------------------
# ar1rlnorm {{{
ar1rlnorm <- function(rho, years, iters=1, margSD=0.6) {
  n <- length(years)
  rhosq <- rho ^ 2
  
  res <- matrix(rnorm(n*iters, mean=0, sd=margSD), nrow=n, ncol=iters)
  res <- apply(res, 2, function(x) {
    for(i in 2:n)
      x[i] <- sqrt(rhosq) * x[i-1] + sqrt(1-rhosq) * x[i]
    return(exp(x))
  }
  )
  return(FLQuant(array(res, dim=c(1,n,1,1,1,iters)),
                 dimnames=list(year=years, iter=seq(1, iters))))
} # }}}


# add uncertainty using RW or correlation matrix
setGeneric("genFLQuant", function(object, cv, method, niter) standardGeneric("genFLQuant"))
setMethod("genFLQuant", c("FLQuant"),
          function(object, cv = 0.2, method = "rw", niter = 250) {
            # use log transform, to be expanded on later versions
            mu <- log(object)
            if(method == "ac") {
              Rho <- cor(t(mu[drop = TRUE]))
              flq <- mvrnorm(niter * dim(mu)[2], rep(0, nrow(Rho)), cv^2 * Rho)
              mu <- propagate(mu, niter)
              flq <- FLQuant(c(t(flq)), dimnames = dimnames(mu))
              flq <- exp(mu + flq)
            }
            if(method == "rw") {
              n.ages  <- dim(mu)[1]
              n.years <- dim(mu)[2]
              n <- n.ages * n.years
              # set up lcs to extract posterior means
              B = diag(n)
              B[B==0] <- NA
              lcs = inla.make.lincombs(Predictor = B)
              # treat mu as a GMRF model -
              # an independant random walk for each age (same variance accross ages)
              form <- x ~ f(years, model = 'rw1', replicate = ages)
              data <- list(x = c(mu), years = rep(1:n.years, each = n.ages), ages  = rep(1:n.ages, n.years))
              result <- inla(form, data = data, control.predictor = list(compute = TRUE),
                             lincomb = lcs, control.inla = list(lincomb.derived.correlation.matrix = TRUE))
              # the covariance of the fitted RW
              RW.cov <- result $ misc $ lincomb.derived.correlation.matrix
              # two options for the mean:
              #  1) use the mean estimate of RW process from INLA
              #     - this is potentially very smooth and lacking in strucure
              #	mu.hat <- result $ summary.linear.predictor $ mean
              #	flq <- mvrnorm(niter, mu.hat, cv^2 * RW.cov)
              #  2) use the original data and add the noise to that
              #  2 is more consistent with ac method and always maintains data structure
              flq <- exp(mvrnorm(niter, c(mu), cv^2 * RW.cov))
              flq <- FLQuant(c(t(flq)), dimnames = dimnames(propagate(mu, niter)))
            }
            units(flq) <- units(object)
            return(flq)
          })


#--------------------------------------------------------------------
# add uncertainty to FLStock using quants
setGeneric("au", function(object, R, C, F, ...) standardGeneric("au"))

setMethod("au", c("FLStock", "FLQuant", "missing", "FLQuant"),
          function(object, R, C, F, ...){
            # requires checking dimensions
            if(!identical(dim(catch.n(object))[-c(1,6)], dim(R)[-c(1,6)]))
              stop("Recruitment vector must have consistent dimensions with the stock object")
            if(!identical(dim(catch.n(object))[-6]    , dim(F)[-6]))
              stop("Harvest matrix must have consistent dimensions with the stock object")
            if(dim(R)[6]!=dim(R)[6]) stop("R and F must have the same number of iterations")
            # get dims and set flq
            dms <- dims(object)
            nages <- dms$age
            nyrs <- dms$year
            niters <- dim(R)[6]
            flq <- FLQuant(dimnames=dimnames(F))
            
            # compute cumulative Z
            Z <- FLCohort(F + m(object))
            Z[] <- apply(Z, c(2:6), cumsum)
            
            # expand variability into [N] by R*[F]
            #Ns <- FLCohort(R[rep(1,nages)])
            #Ns <- Ns*exp(-Z)
            Ns <- R[rep(1,nages)]
            # stupid hack to get the right dimensions
            Z[,dimnames(Ns)[[2]]] <- Ns*exp(-Z[,dimnames(Ns)[[2]]])
            Ns <- as(Z, "FLQuant")
            
            # Update object
            stock.n(object) <- flq
            # [R]
            stock.n(object)[1] <- R
            # [N]
            stock.n(object)[-1,-1] <- Ns[-nages,-nyrs]
            # plus group
            stock.n(object)[nages,-1] <- Ns[nages,-nyrs] + stock.n(object)[nages,-1]
            stock(object) <- computeStock(object)
            # [F]
            harvest(object) <- F
            # [C]
            Z <- harvest(object) + m(object)
            Cs <- harvest(object)/Z*(1-exp(-Z))*stock.n(object)
            catch.n(object) <- Cs
            catch(object) <- computeCatch(object)
            # [L] & [D] rebuilt from C
            # Ds=D/(D+L)*Cs where Cs is the simulated catch
            D <- discards.n(object)
            L <- landings.n(object)
            discards.n(object) <- D/(D+L)*Cs
            discards(object) <- computeDiscards(object)
            landings.n(object) <- Cs - discards.n(object)
            landings(object) <- computeLandings(object)
            # out
            object
          })

#--------------------------------------------------------------------
# fake fbar for signature FLStockLen
setMethod("fbar", signature(object = "FLStockLen"), 
          function(object, ...){
            selectMethod("fbar", signature = "FLStock")(object)
          })

#--------------------------------------------------------------------
# functions for stock creation

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

### ------------------------------------------------------------------------ ###
### combine severals parts of FLStock in one object
### ------------------------------------------------------------------------ ###
### function that combines list of object over iter dimension
ibind <- function(object) {
  
  if (!is.list(object)) stop("object has to be a list")
  
  ### get iterations per list element
  iters_list <- unlist(lapply(object, function(x){dims(x)$iter}))
  
  ### create list with iterations for each part
  iter_end <- cumsum(iters_list)
  iter_start <- iter_end - iters_list + 1
  iters <- lapply(seq_along(iter_start), function(x) {
    seq(iter_start[x], iter_end[x])
  })
  
  ### create template with correct iter dimension
  template <- propagate(FLCore::iter(object[[1]], 1), sum(iters_list))
  
  ### fill with values from list
  for (part in seq_along(iters)) {
    
    FLCore::iter(template, iters[[part]]) <- object[[part]]
    
  }
  
  return(template)
  
}

stock_ibind <- function(...) {
  
  stk_list <- list(...)
  
  if (!all(sapply(stk_list, is, "FLS"))) stop("all elements have to be FLS")
  
  ### combine stocks
  stk_tmp <- ibind(stk_list)
  
  ### find additional attributes
  attr_add <- setdiff(names(attributes(stk_tmp)), 
                      c(slotNames(stk_tmp), "class", "ctrl.mp"))
  
  ### loop through them
  for (attr_i in attr_add) {
    attr(stk_tmp, attr_i) <- ibind(lapply(stk_list, attr, attr_i))
  }
  
  return(stk_tmp)
}
