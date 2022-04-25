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

which_model_fit <- read.csv("which_model_fit.csv")

# for (this_pop_now in c(1:4, 7:14, 16:19, 21, 6, 5, 20, 15)) {
for (this_pop_now in 1:nrow(which_model_fit)) {

source("data_load.R")
  
if (single_pop) {
 which.dataset <- unique(data.all$pop_spec)[this_pop_now]
 data.all      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
 sampling      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
 this_model_fit <- which_model_fit[which(which_model_fit$pop_spec == as.character(which.dataset)), ]$model %>%
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

