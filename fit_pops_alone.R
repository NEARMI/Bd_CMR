#####################################
## Fit CMR model to amphibian data ##
#####################################

source("packages_functions.R")
source("../ggplot_theme.R")

stan.iter     <- 1500
stan.burn     <- 500
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin

red_ind    <- FALSE
single_pop <- TRUE

which_model_fit <- read.csv("which_model_fit.csv", header = FALSE)

for (this_pop_now in 1:nrow(which_model_fit)) {

source("data_load.R")

if (single_pop) {
 which.dataset <- which_model_fit[this_pop_now, 1]
 data.all      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
 sampling      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
 this_model_fit <- which_model_fit[this_pop_now, 2] %>%
   paste("stan_current/", ., sep = "") %>% paste(., ".stan", sep = "")
}
  
source("data_manip.R")
source("data_stan.R")
source("data_covariates.R")
source("stan_indices.R")
source("dataset_notes.R")
source("stan_fit_single.R")
source("plotting.R")

rm(stan.fit)
rm(stan_data)
rm(capt_history)
rm(capt_history.p)
rm(capt_history.phi)
gc()

}

