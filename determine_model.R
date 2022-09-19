####
## Determine which of so many alternative models to fit
####

# unique(data.all$pop_spec)

## Viable pops to include in multi-pop MeHg individual-level model
 ## Total: 4, 5, 6, 8, 9, 15, 16, 17, 18, 19, 21
 ## ANBO: Blackrock H, Jones Pond, Sonoma Mt
 ## RANA: Dilman Meadows, Fox Creek, Jones Pond, Lost Horse, San Fransisquito, Three Creeks
 ## BCF: Lily Pond, Matthews Pond

if (!fit_only_mehg) {
if (sing_pop) {
 which_model_fit <- read.csv("which_model_fit.csv")
 which_stan_file <- which_model_fit[which(which_model_fit[, 1] == unique(data.all$pop_spec) %>% as.character()), ][2]
} else {
  
if (!multi_spec) {
  if (fit_ind_mehg) {
   which_stan_file <- "CMR_multiple_populations_mehg_ssp"
  } else {
   which_stan_file <- "CMR_multiple_populations_ssp"
  }

} else {
  
  if (fit_ind_mehg) {
   which_stan_file     <- "CMR_multiple_populations_mehg"
  } else {
    if (multi_spec_red) {
      which_stan_file  <- "CMR_multiple_populations_red"
    } else {
      if (red_p_model) {
      # which_stan_file <- "CMR_multiple_populations_alt_p"
        which_stan_file <- "CMR_multiple_populations_alt_p_len"
      } else {
       which_stan_file <- "CMR_multiple_populations"
      }
      
    }
  }
  
}
  
}
} else {
  which_stan_file <- "CMR_mehg"
}

model_name      <- paste("fits/", paste(which_stan_file, "_", sep = ""), Sys.Date(), ".Rds", sep = "")
which_stan_file <- paste("stan_current/", which_stan_file, ".stan", sep = "")
