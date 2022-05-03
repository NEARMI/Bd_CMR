########################################################################
## Fit only a single population at a time, utilizing the command line ##
## to pick which population (for submitting multiple jobs at once)    ##
########################################################################

source("packages_functions.R")

stan.iter     <- 3500
stan.burn     <- 500
stan.thin     <- 3
stan.length   <- (stan.iter - stan.burn) / stan.thin
red_ind       <- FALSE
some_pops     <- TRUE

source("data_load.R")

if (some_pops) {
which.dataset <- unique(data.all$pop_spec)[-c(10:14)] %>% droplevels()
data.all      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
sampling      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
}
  
source("data_manip.R")
source("data_stan.R")
source("data_covariates.R")
source("stan_indices.R")
source("establishing_mm.R")
source("stan_fit_mm.R") 
source("plotting_multipop.R")
