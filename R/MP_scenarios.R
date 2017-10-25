### ------------------------------------------------------------------------ ###
### definition of scenarios ####
### ------------------------------------------------------------------------ ###
### the output is a list ctrl.mps
### each list element is one scenario


### ------------------------------------------------------------------------ ###
### catch rule 3.2.1: perfect knowledge scenarios
### ------------------------------------------------------------------------ ###
### no uncertainty except for recruitment uncertainty
### no noise for observation, lhist, refpts, implementation...
### factor f:
### option a) deterministic LFeM
### option b) M from OM in numerator
### option c) F0.1 from YPR from OM

### template with default options
ctrl.mp_3.2.1 <- list(scns = data.frame(id = "wklife_test",
                                  stringsAsFactors = FALSE),
                stock_name = NA,
                ctrl.om = list(stk_pos = NA, stk = "stk", sr = "srbh", sr.res = "srbh.res",
                               observations = "observations", sr.res.mult = TRUE),
                ctrl.oem = list(method = "obs_bio_len", len_noise_sd = 0.0, 
                                len_sd = 1, len_sd_cut = 2),
                ctrl.f = list(method = "wklife_3.2.1_f", n_catch_yrs = 1),
                ctrl.x = list(method = "wklife_3.2.1_x", n_catch_yrs = 1, multiplier = NA,
                              interval = 2, start_year = 100),
                ctrl.h = list(method = "wklife_3.2.1_h"),
                ctrl.k = NULL,
                ctrl.w = NULL,
                ctrl.j = list(method = "ctrl_catch_workaround"),
                ctrl.l = NULL,
                scn_desc = list(catch_rule = "3.2.1", 
                                uncertainty = "perfect_knowledge")
)

### create list with control objects
### list with additional options
### test all factor individually = 6 scenarios
add_f_options <- list(list(option_r = "a"), list(option_r = "b"),
                      list(option_f = "a"), list(option_f = "b"),
                      list(option_f = "c"),
                      list(option_b = "a"))

### create list with these 6 scenarios
ctrl.mps <- lapply(seq_along(add_f_options), function(x){
  ### create ctrl.mp object, based on template
  temp <- ctrl.mp_3.2.1
  ### add ctrl.f option
  temp$ctrl.f <- c(temp$ctrl.f, unlist(add_f_options[x]))
  ### create name
  name_temp <- NULL
  for(i in 1:length(unlist(add_f_options[x]))){
    name_temp <- append(name_temp, 
                        c(names(unlist(add_f_options[x]))[i], 
                          unlist(add_f_options[x])[i]))
  }
  temp$ID <- paste(temp$ID, paste(name_temp, collapse = "_"), sep = "_", collapse = "_")
  return(temp)
})

### apply these 6 scenarios for all 30 stocks
ctrl.mps <- lapply(1:30, function(x){
  ### template
  temp <- ctrl.mps
  ### insert stock position
  temp <- lapply(temp, function(y){
    y$ctrl.om$stk_pos <- x
    y$ID <- paste(y$ID, "stock", x, sep = "_", collapse = "_")
    y
  })
  temp
})
### create flat list structure
ctrl.mps <- unlist(ctrl.mps, recursive = FALSE)

### ------------------------------------------------------------------------ ###
### catch rule 3.2.2
### ------------------------------------------------------------------------ ###

### template
ctrl.mp_3.2.2 <- list(scns = data.frame(id = "wklife_test",
                                  stringsAsFactors = FALSE),
                stock_name = NA,
                ctrl.om = list(stk_pos = NA, stk = "stk", sr = "srbh", sr.res = "srbh.res",
                               observations = "observations", sr.res.mult = TRUE),
                ctrl.oem = list(method = "idx_bio"),
                ctrl.f = NULL,
                ctrl.x = list(method = "wklife_3.2.2_x", multiplier = NA,
                              interval = 2, start_year = 100, index_years = 1),
                ctrl.h = list(method = "wklife_3.2.2_h"),
                ctrl.k = NULL,
                ctrl.w = NULL,
                ctrl.j = list(method="ctrl_catch_workaround"),
                ctrl.l = NULL,
                ID = "perfect_knowledge",
                scn_desc = list(catch_rule = "3.2.2", 
                                uncertainty = "perfect_knowledge")
)
### add stocks 1-30
ctrl.mp_3.2.2_list <- lapply(1:30, function(x){
  ### load template
  res <- ctrl.mp_3.2.2
  ### add stock position
  res$ctrl.om$stk_pos <- x
  return(res)
})

### add to ctrl list
ctrl.mps <- c(ctrl.mps, ctrl.mp_3.2.2_list)

### ------------------------------------------------------------------------ ###
### catch rule 3.2.1
### ------------------------------------------------------------------------ ###
### perfect knowledge
### factor f:
### option a) replace LFeM with length when fished at FMSY
### option b) replace M in numerator with FMSY
### option c) replace F0.1 from YPR with F0.1 from OM

### find scenario definitions
f_a <- seq(3, 180, by = 6)
f_b <- seq(4, 180, by = 6)
f_c <- seq(5, 180, by = 6)
fs <- sort(c(f_a, f_b, f_c))
rm(f_a, f_b, f_c)
### copy them
ctrl.mp_3.2.1_perfect_knowledge <- ctrl.mps[fs]
### add the perfect knowledge flag
ctrl.mp_3.2.1_perfect_knowledge <- lapply(ctrl.mp_3.2.1_perfect_knowledge,
                                          function(x){
  x$ctrl.f$perfect_knowledge <- TRUE
  return(x)
})
### add to list
ctrl.mps <- c(ctrl.mps, ctrl.mp_3.2.1_perfect_knowledge)

### ------------------------------------------------------------------------ ###
### catch rule 3.2.1
### ------------------------------------------------------------------------ ###
### perfect knowledge
### but MK = 1.5 shortcut

f_a <- seq(3, 180, by = 6)
ctrl.mp_3.2.1_MK <- ctrl.mps[f_a]
### add the perfect knowledge flag
ctrl.mp_3.2.1_MK <- lapply(ctrl.mp_3.2.1_MK,
                          function(x){
                            x$ctrl.f$MK <- 1.5
                            return(x)
                          })
### add to list
ctrl.mps <- c(ctrl.mps, ctrl.mp_3.2.1_MK)

### ------------------------------------------------------------------------ ###
### catch rule 3.2.2
### with annual TAC
### ------------------------------------------------------------------------ ###

### load template
ctrl.mp_3.2.2_list_1 <- ctrl.mp_3.2.2_list

### change to annual TAC
ctrl.mp_3.2.2_list_1 <- lapply(ctrl.mp_3.2.2_list_1, function(x){
  x$ctrl.x$interval <- 1
  return(x)
})

### add to list
ctrl.mps <- c(ctrl.mps, ctrl.mp_3.2.2_list_1)

### ------------------------------------------------------------------------ ###
### catch rule 3.2.2
### with observation error
### ------------------------------------------------------------------------ ###

### load template
ctrl.mp_3.2.2_list_2 <- ctrl.mp_3.2.2_list

### add error
ctrl.mp_3.2.2_list_2 <- lapply(ctrl.mp_3.2.2_list_2, function(x){
  x$scn_desc$uncertainty <- "observation_error"
  x$ctrl.l <- list(method = "noise.wrapper", fun = "rlnorm", mean = 0, sd = 0.1,
                   multiplicative=TRUE)
  return(x)
})

### add to list
ctrl.mps <- c(ctrl.mps, ctrl.mp_3.2.2_list_2)

### ------------------------------------------------------------------------ ###
### catch rule 3.2.1
### ------------------------------------------------------------------------ ###
### perfect knowledge & combinations
### factor f:
### option a) replace LFeM with length when fished at FMSY
### option b) replace M in numerator with FMSY
### combinations f:a/b with r:a/b and b:a

### find template scenarios
# subset(scn_df, catch_rule == "3.2.1" & perfect_knowledge == TRUE &
#        options == "option_f:a perfect_knowledge:TRUE")

### find scenario definitions
pos_f_a <- seq(211, 300, by = 3)
pos_f_b <- seq(212, 300, by = 3)

### copy them
ctrl.mp_3.2.1_combs_fa <- ctrl.mps[pos_f_a]
ctrl.mp_3.2.1_combs_fb <- ctrl.mps[pos_f_b]

### add r:a to f:a
ctrl.mp_3.2.1_combs_fa_ra <- lapply(ctrl.mp_3.2.1_combs_fa, function(x){
  x$ctrl.f$option_r <- "a"
  return(x)
})
### add r:b to f:a
ctrl.mp_3.2.1_combs_fa_rb <- lapply(ctrl.mp_3.2.1_combs_fa, function(x){
  x$ctrl.f$option_r <- "b"
  return(x)
})
### add r:a to f:a
ctrl.mp_3.2.1_combs_fb_ra <- lapply(ctrl.mp_3.2.1_combs_fb, function(x){
  x$ctrl.f$option_r <- "a"
  return(x)
})
### add r:b to f:a
ctrl.mp_3.2.1_combs_fb_rb <- lapply(ctrl.mp_3.2.1_combs_fb, function(x){
  x$ctrl.f$option_r <- "b"
  return(x)
})

### combine them
ctrl.mp_3.2.1_combs <- c(ctrl.mp_3.2.1_combs_fa_ra, ctrl.mp_3.2.1_combs_fa_rb,
                         ctrl.mp_3.2.1_combs_fb_ra, ctrl.mp_3.2.1_combs_fb_rb)

### add b:a to all combinations
ctrl.mp_3.2.1_combs <- lapply(ctrl.mp_3.2.1_combs, function(x){
  x$ctrl.f$option_b <- "a"
  return(x)
})

### add to list
ctrl.mps <- c(ctrl.mps, ctrl.mp_3.2.1_combs)

### ------------------------------------------------------------------------ ###
### 3.2.2 with lower FproxyMSY
### ------------------------------------------------------------------------ ###
### HCR multiplier (multiply adviced catch from HCR)
### 0.5 for all stocks and fishing histories
### 1 - 0.5 in 0.05 steps for herring and sandeel


### load template (with observation error)
ctrl.mp_3.2.2_list_3 <- ctrl.mp_3.2.2_list_2

### add multiplier for all stocks
ctrl.mp_3.2.2_list_3 <- lapply(ctrl.mp_3.2.2_list_3, function(x){
  x$ctrl.x$multiplier <- 0.5
  return(x)
})

### extract pelagic stocks
ctrl.mp_3.2.2_list_4 <- ctrl.mp_3.2.2_list_3[c(1, 9, 16, 24)]

### add multipliers
### (0.5 and 1 already defined in earlier scenarios)
mults <- seq(0.55, 0.95, by = 0.05)
ctrl.mp_3.2.2_list_4 <- lapply(mults, function(x){
  ### loop through pelagic stocks
  res <- lapply(ctrl.mp_3.2.2_list_4, function(y){
    y$ctrl.x$multiplier <- x
    return(y)
  })
})
### return as plain list
ctrl.mp_3.2.2_list_4 <- unlist(ctrl.mp_3.2.2_list_4, recursive = FALSE)

### add to list
ctrl.mps <- c(ctrl.mps, ctrl.mp_3.2.2_list_3, ctrl.mp_3.2.2_list_4)

### ------------------------------------------------------------------------ ###
### 3.2.2 with w != 1.4
### ------------------------------------------------------------------------ ###
### for herring and sandeel

### load template (with observation error)
ctrl.mp_3.2.2_list_5 <- ctrl.mps[c(361, 369, 376, 384)]

# scn_df[scn_df$catch_rule == "3.2.2" & scn_df$stock %in% c("her-nis", "san-ns4") & is.na(scn_df$HCRmult) & scn_df$TAC == 2 & scn_df$uncertainty == "observation_error",]

### add various w values
ws <- c(2, 3, 4, 5, 10)
ctrl.mp_3.2.2_list_5 <- lapply(ws, function(x){
  ### loop through pelagic stocks
  res <- lapply(ctrl.mp_3.2.2_list_5, function(y){
    y$ctrl.x$w <- x
    return(y)
  })
})
### return as plain list
ctrl.mp_3.2.2_list_5 <- unlist(ctrl.mp_3.2.2_list_5, recursive = FALSE)

### add to list
ctrl.mps <- c(ctrl.mps, ctrl.mp_3.2.2_list_5)

### ------------------------------------------------------------------------ ###
### 3.1: SPiCT
### ------------------------------------------------------------------------ ###

### template
ctrl.mp_spict <- list(scns = data.frame(id = "wklife_test",
                                  stringsAsFactors = FALSE),
                stock_name = NA,
                ctrl.om = list(stk_pos = NA, stk = "stk", sr = "srbh", sr.res = "srbh.res",
                               observations = "observations", sr.res.mult = TRUE),
                ctrl.oem = list(method = "idx_bio"),
                ctrl.f = list(method = "wklife_3.1_f", interval = 2,
                              start_year = 100),
                ctrl.x = list(method = "wklife_3.1_x", interval = 2),
                ctrl.h = list(method = "wklife_3.1_h"),
                ctrl.l = list(method = "noise.wrapper", fun = "rlnorm", 
                              mean = 0, sd = 0.1, multiplicative=TRUE),
                ctrl.j = list(method="ctrl_catch_workaround"),
                scn_desc = list(catch_rule = "3.1", 
                                uncertainty = "observation_error")
)
### add all 30 stocks
ctrl.mp_spict_list <- lapply(1:30, function(x){
  ### template
  tmp <- ctrl.mp_spict
  ### add stock
  tmp$ctrl.om$stk_pos <- x
  return(tmp)
  
})

### add to list
ctrl.mps <- c(ctrl.mps, ctrl.mp_spict_list)

### ------------------------------------------------------------------------ ###
###  3.2.1 combinations with error
### ------------------------------------------------------------------------ ###
# scn_df[scn_df$catch_rule == "3.2.1" &
#          scn_df$options %in% c("option_f:a perfect_knowledge:TRUE option_r:a option_b:a",
#                                "option_f:a perfect_knowledge:TRUE option_r:b option_b:a"), ]

ctrl.mps_3.2.1_combs_error <- ctrl.mps[391:450]

ctrl.mps_3.2.1_combs_error <- lapply(ctrl.mps_3.2.1_combs_error, function(x){
  
  ### length frequency error
  x$ctrl.oem$len_noise_sd = 0.2
  ### select OMs with implemented noise (parameters and index)
  x$scn_desc$uncertainty <- "observation_error"
  x$ctrl.f$perfect_knowledge <- NULL
  x$ctrl.f$MK <- 1.5
  ### implementation error
  x$ctrl.l <- list(method = "noise.wrapper", fun = "rlnorm", mean = 0, sd = 0.1,
                   multiplicative=TRUE)
  x$ID <- "3.2.1_combination_noise"

  return(x)

})

### add to list
ctrl.mps <- c(ctrl.mps, ctrl.mps_3.2.1_combs_error)

### ------------------------------------------------------------------------ ###
###  3.2.1 combinations with error & multiplier for pol-nsea
### ------------------------------------------------------------------------ ###
# scn_df[scn_df$catch_rule == "3.2.1" &
#        scn_df$options %in% c("option_f:a option_r:a option_b:a MK:1.5",
#                              "option_f:a option_r:b option_b:a MK:1.5") &
#        scn_df$stock == "pol-nsea", ]

### extract templates for pol-nsea
ctrl.mps_3.2.1_combs_mult <- ctrl.mps[c(628, 643, 658, 673)]

### add multipliers
mults <- c(seq(from = 0.5, to = 0.95, by = 0.05))
ctrl.mps_3.2.1_combs_mult <- lapply(ctrl.mps_3.2.1_combs_mult, function(x){
  ### "loop" through multipliers
  res_x <- lapply(mults, function(y){
    ### load template
    res_y <- x
    ### add multiplier
    res_y$ctrl.x$multiplier <- y
    return(res_y)
  })
  return(res_x)
})
ctrl.mps_3.2.1_combs_mult <- unlist(ctrl.mps_3.2.1_combs_mult, recursive = FALSE)

### add to list
ctrl.mps <- c(ctrl.mps, ctrl.mps_3.2.1_combs_mult)

### ------------------------------------------------------------------------ ###
###  3.2.1 combinations with error & b_w for pol-nsea
### ------------------------------------------------------------------------ ###
### extract templates for pol-nsea
ctrl.mps_3.2.1_combs_w <- ctrl.mps[c(628, 643, 658, 673)]

### add w
mults <- c(1.4, 1.6, 1.8, 2, 3, 4, 5)
ctrl.mps_3.2.1_combs_w <- lapply(ctrl.mps_3.2.1_combs_w, function(x){
  ### "loop" through multipliers
  res_x <- lapply(mults, function(y){
    ### load template
    res_y <- x
    ### add multiplier
    res_y$ctrl.f$b_w <- y
    return(res_y)
  })
  return(res_x)
})
ctrl.mps_3.2.1_combs_w <- unlist(ctrl.mps_3.2.1_combs_w, recursive = FALSE)

### add to list
ctrl.mps <- c(ctrl.mps, ctrl.mps_3.2.1_combs_w)

### ------------------------------------------------------------------------ ###
###  3.2.1 combs & noise & meeting changes
### ------------------------------------------------------------------------ ###
### z
# scn_df[scn_df$catch_rule == "3.2.1" & scn_df$uncertainty == "observation_error" &
#          is.na(scn_df$HCRmult) & is.na(scn_df$w) & is.na(scn_df$b_w) &
#          scn_df$options == "option_f:a option_r:a option_b:a MK:1.5" &
#          scn_df$stk_pos2 %in% c(2,6,11,15), ]$scenario
### template
ctrl.mps_red <- ctrl.mps[c(628, 632, 637, 641, 643, 647, 652, 656)]
### z values
ctrl.mps_z <- lapply(ctrl.mps_red, function(ctrl){
  ### go through z list
  res <- lapply(c(1, 1.1, 1.2, 1.5, 2, 3, 5), function(z){
    ctrl_tmp <- ctrl
    ctrl_tmp$ctrl.f$b_z <- z
    ctrl_tmp
  })
})
ctrl.mps_z <- unlist(ctrl.mps_z, recursive = FALSE)

### multiplier
# scn_df[scn_df$catch_rule == "3.2.1" & scn_df$uncertainty == "observation_error" &
#          is.na(scn_df$HCRmult) & is.na(scn_df$w) & is.na(scn_df$b_w) &
#          scn_df$options == "option_f:a option_r:a option_b:a MK:1.5" &
#          scn_df$stk_pos2 %in% c(6,11,15), ]$scenario
### template
ctrl.mps_red2 <- ctrl.mps[c(632, 637, 641, 647, 652, 656)]
### z values
ctrl.mps_mult <- lapply(ctrl.mps_red2, function(ctrl){
  ### go through z list
  res <- lapply(seq(0.5, 0.9, 0.1), function(mult){
    ctrl_tmp <- ctrl
    ctrl_tmp$ctrl.x$multiplier <- mult
    ctrl_tmp
  })
})
ctrl.mps_mult <- unlist(ctrl.mps_mult, recursive = FALSE)

### upper constraints
### template
ctrl.mps_red
ctrl.mps_upper <- lapply(ctrl.mps_red, function(ctrl){
  ### go through z list
  res <- lapply(c(1.1, 1.15, 1.2, 1.25, 1.3, 1.4, 1.5, 2, 3), function(upper){
    ctrl_tmp <- ctrl
    ctrl_tmp$ctrl.h$upper_constraint <- upper
    ctrl_tmp
  })
})
ctrl.mps_upper <- unlist(ctrl.mps_upper, recursive = FALSE)

### add
ctrl.mps <- c(ctrl.mps, ctrl.mps_z, ctrl.mps_mult, ctrl.mps_upper)

### ------------------------------------------------------------------------ ###
###  3.2.1 combs & noise & meeting changes: upper/lower constraints combinations
### ------------------------------------------------------------------------ ###

ctrl.mps_upper_lower <- lapply(ctrl.mps_red, function(ctrl){
  ### go through upper constraints
  res <- lapply(c(1.1, 1.2, 1.3), function(upper){
    ### go through lower constraints
    res_tmp <- lapply(seq(0.5, 0.9, 0.1), function(lower){
      ctrl_tmp <- ctrl
      ctrl_tmp$ctrl.h$upper_constraint <- upper
      ctrl_tmp$ctrl.h$lower_constraint <- lower
      return(ctrl_tmp)
    })
  })
  return(unlist(res, recursive = FALSE))
})
ctrl.mps_upper_lower <- unlist(ctrl.mps_upper_lower, recursive = FALSE)
### add to list
ctrl.mps <- c(ctrl.mps, ctrl.mps_upper_lower)
### 913-1032

### ------------------------------------------------------------------------ ###
###  3.2.1 combs & noise & meeting changes: mult & b exponent
### ------------------------------------------------------------------------ ###

ctrl.mps_mult_z <- lapply(ctrl.mps_red, function(ctrl){
  ### go through upper constraints
  res <- lapply(seq(0.8, 1, 0.05), function(mult){
    ### go through lower constraints
    res_tmp <- lapply(seq(1, 3, 1), function(z){
      ctrl_tmp <- ctrl
      ctrl_tmp$ctrl.x$multiplier <- mult
      ctrl_tmp$ctrl.f$b_z <- z
      return(ctrl_tmp)
    })
  })
  return(unlist(res, recursive = FALSE))
})
ctrl.mps_mult_z <- unlist(ctrl.mps_mult_z, recursive = FALSE)

### add to list
ctrl.mps <- c(ctrl.mps, ctrl.mps_mult_z)
### 1033-1152



### ------------------------------------------------------------------------ ###
### create table with specifications
### ------------------------------------------------------------------------ ###

### load wklife stock data
wklife <- read.csv("input/wklife.csv")

scn_df <- data.frame(scenario = seq_along(ctrl.mps))

### stock position
scn_df$stk_pos <- unlist(lapply(ctrl.mps, function(x){
  x$ctrl.om$stk_pos
}))
### position in list of 15 stocks
scn_df$stk_pos2 <- ifelse(scn_df$stk_pos <= 15, scn_df$stk_pos, scn_df$stk_pos - 15)

### set stock name
scn_df <- merge(x = scn_df, y = wklife[, c("X", "stock")], 
                by.x = "stk_pos2", by.y = "X", all = TRUE)
### sort
scn_df <- scn_df[order(scn_df$scenario), ]
scn_df <- scn_df[, c("scenario", "stk_pos", "stk_pos2", "stock")]

### fishing history
scn_df$fhist <- ifelse(scn_df$stk_pos <= 15, "one-way", "roller-coaster")

### catch rule
scn_df$catch_rule <- unlist(lapply(ctrl.mps, function(x){
  x$scn_desc$catch_rule
}))

### TAC periodicity
scn_df$TAC <- unlist(lapply(ctrl.mps, function(x){
  ifelse(is.null(x$ctrl.x$interval), NA, x$ctrl.x$interval)
}))

### catch rule options
scn_df$options <- unlist(lapply(ctrl.mps, function(x){
  opts <- unlist(lapply(names(x$ctrl.f[3:length(x$ctrl.f)]), function(y){
    y
    #res_temp <- strsplit(y, split = "_")[[1]]
    #ifelse(length(res_temp) == 1, res_temp[1], res_temp[2])
  }))
  pars <- unlist(x$ctrl.f[3:length(x$ctrl.f)])
  paste(opts, pars, sep = ":", collapse = " ")
}))

### uncertainty
scn_df$uncertainty <- unlist(lapply(ctrl.mps, function(x){
  x$scn_desc$uncertainty
}))

### perfect knowledge
scn_df$perfect_knowledge <- unlist(lapply(ctrl.mps, function(x){
  res <- x$ctrl.f$perfect_knowledge
  if(is.null(res)) return(FALSE)
  return(res)
}))

### HCR multiplier
scn_df$HCRmult <- unlist(lapply(ctrl.mps, function(x){
  ifelse(is.null(x$ctrl.x$multiplier), NA, x$ctrl.x$multiplier)
}))

### w
scn_df$w <- unlist(lapply(ctrl.mps, function(x){
  ifelse(is.null(x$ctrl.x$w), NA, x$ctrl.x$w)
}))

### w
scn_df$b_w <- unlist(lapply(ctrl.mps, function(x){
  ifelse(is.null(x$ctrl.f$b_w), NA, x$ctrl.f$b_w)
}))

### z (b exponent)
scn_df$b_z <- unlist(lapply(ctrl.mps, function(x){
  ifelse(is.null(x$ctrl.f$b_z), NA, x$ctrl.f$b_z)
}))

### constraints
scn_df$upper_constraint <- unlist(lapply(ctrl.mps, function(x){
  ifelse(is.null(x$ctrl.h$upper_constraint), NA, x$ctrl.h$upper_constraint)
}))
scn_df$lower_constraint <- unlist(lapply(ctrl.mps, function(x){
  ifelse(is.null(x$ctrl.h$lower_constraint), NA, x$ctrl.h$lower_constraint)
}))

### get rid of row names
row.names(scn_df) <- NULL

rm(ctrl.mp_3.2.1, ctrl.mp_3.2.1_combs_fa, 
   ctrl.mp_3.2.1_combs_fa_ra, ctrl.mp_3.2.1_combs_fa_rb, ctrl.mp_3.2.1_combs_fb,
   ctrl.mp_3.2.1_combs_fb_ra, ctrl.mp_3.2.1_combs_fb_rb, ctrl.mp_3.2.1_MK,
   ctrl.mp_3.2.1_perfect_knowledge, ctrl.mp_3.2.2, ctrl.mp_3.2.2_list,
   ctrl.mp_3.2.1_combs,
   ctrl.mps_upper_lower, ctrl.mps_mult_z, ctrl.mps_red)

