#####################################
## Fit CMR model to amphibian data ##
#####################################

source("packages_functions.R")

stan.iter     <- 1500
stan.burn     <- 500
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin

red_ind    <- FALSE
single_pop <- TRUE

for (this_pop_now in 1:22) {

source("data_load.R")
  
if (single_pop) {
 which.dataset <- unique(data.all$pop_spec)[this_pop_now]
 data.all      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
 sampling      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
}

source("data_manip.R")
source("data_stan.R")
source("data_covariates.R")
source("stan_fit_single.R")

rm(stan.fit)
rm(stan_data)
rm(capt_history)
rm(capt_history.p)
rm(capt_history.phi)
gc()

}

