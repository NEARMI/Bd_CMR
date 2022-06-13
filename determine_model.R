####
## Determine which of so many alternative models to fit
####

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
   which_stan_file <- "CMR_multiple_populations_mehg"
  } else {
    if (multi_spec_red) {
      which_stan_file <- "CMR_multiple_populations_red"
    } else {
      which_stan_file <- "CMR_multiple_populations"
    }
  }
  
}
  
}

model_name      <- paste("fits/", paste(which_stan_file, "_", sep = ""), Sys.Date(), ".Rds", sep = "")
which_stan_file <- paste("stan_current/", which_stan_file, ".stan", sep = "")
