####
## Determine which of so many alternative models to fit
####

if (!fit_only_mehg) {
if (sing_pop) {
 which_model_fit <- read.csv("which_model_fit.csv", header = T)
 which_model_fit %<>% left_join(
   ., data.frame(
      pop_spec = unique(data.all$pop_spec) %>% as.character()
    , pop_spec_num = seq(n_distinct(data.all$pop_spec))
   )
 )
} else {
  
if (!multi_spec) {
  if (fit_ind_mehg) {
    if (red_p_model) {
      which_stan_file <- "CMR_multiple_populations_mehg_ssp_alt_p_len"
    } else {
      which_stan_file <- "CMR_multiple_populations_mehg_ssp"
    }
  } else {
    if (red_p_model) {
      which_stan_file <- "CMR_multiple_populations_ssp_alt_p_len"
    } else {
      which_stan_file <- "CMR_multiple_populations_ssp"
    }
  }

} else {
  
  if (fit_ind_mehg) {
   which_stan_file      <- "CMR_multiple_populations_mehg"
  } else {
    if (multi_spec_red) {
      which_stan_file   <- "CMR_multiple_populations_red"
    } else {
      if (red_p_model) {
        which_stan_file <- "CMR_multiple_populations_alt_p_len"
      } else {
       which_stan_file  <- "CMR_multiple_populations"
      }
      
    }
  }
  
}
  
}
} else {
  which_stan_file <- "CMR_mehg"
}
