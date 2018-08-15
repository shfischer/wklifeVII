### ------------------------------------------------------------------------ ###
### o() generate observations ####
### ------------------------------------------------------------------------ ###
#' generate observations: observed index/indices and stock
#' @param method name of the chosen HCR function
#' @param stk operating model (stock)
#' @param observations list with index/indices
#' @param ay current (intermediate) year
#' @param tracking object for tracking
#' @param ... additional argument, passed on to function defined by method
#' @return list with four elements:
#'   \describe{
#'     \item{stk}{observed stock}
#'     \item{idx}{observed index/indices}
#'     \item{observations}{full observations, including full index/indices}
#'     \item{tracking}{object for tracking, updated}
#'   }

o <- function(...) {
	args <- list(...)
	method <- args$method
	args$method <- NULL
	### checks 
	if (!is(args$stk, "FLS")) stop("stk must be of class FLStock")
	### dispatch
	out <- do.call(method, args)
	if (!is(out$stk, "FLS")) stop("stk must be of class FLStock")
	if (!is(out$idx, "FLIndices")) stop("idx must be of class FLIndices")
	out
}


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
  args$lhpar <- attr(stk, "lhpar")
  
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

idx_bio <- function(stk, observations, ay, tracking,
                    lst_catch = -1, ### last catch/idx year, relative to 
                    lst_idx = -1,   ### intermediate year (ay)
                    ...){
  
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
  
  ### subset all indices to requested years -> observations
  idx0 <- lapply(idx, window, end = ay + lst_idx)
  
  ### return observations
  list(stk = stk0, idx = idx0, observations = list(idx = idx), 
       tracking = tracking)
}


### ------------------------------------------------------------------------ ###
### length frequencies
### ------------------------------------------------------------------------ ###
length_freq <- function(stk, ay = NULL,
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
