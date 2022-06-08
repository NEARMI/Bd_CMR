#####################################
## Fit CMR model to amphibian data ##
#####################################

#### Notes June 8, 2022 ---- 

## 0.1) Couldn't get anything to work to fix intercept uncertainty apart from just shifting entirely to fixed effects.
 ## Will still want to use _some_ random effects (for individuals etc.) but it does seems like a general issue that because there is so
  ## much variability among the populations random effects may not be so viable at a population-level

## Fixed model likely will become the default one (i.e., other models moved to /dev, but leaving all
 ## stan models in the main folder at least for now)

## 1) To do:
 ## A) [x] Update code to accommodate the single species multi pop fixed effects model
 ## B) [x] Update folder of stan models, update repo
   ##   -- Done, but still need to compare model fits to older model fits to check everything is as it should be
 ## C) [ ] Run new single species jobs with fixed effects for intercepts and slopes on Yeti
 ## D) [ ] Make similar updates to the multi-species model
 ## E) [ ] Run the multi-species model
 ## F) [ ] Update all of the helper scripts such as the plotting scripts
 ## G) [ ] Figure out what to do with the Newt MA and PA models
 ## H) [ ] Run the whole model

#### Code ----

#### NOTE: In this file and all other files search *** for current choices that could potentially change

## Packages and Functions
source("packages_functions.R")
source("../ggplot_theme.R")

## Read in data
source("data_load.R")

## For dev and debug purposes pick a subset of locations
some_pops  <- TRUE

# data.all %>% group_by(Year, pop_spec, Mark) %>% filter(BdSample == "Y") %>% summarize(nswab = n()) %>% arrange(desc(nswab)) %>% as.data.frame()       

if (some_pops) {
#which.dataset <- unique(data.all$pop_spec)[c(1:9, 11, 13, 15:21)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(15:21)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(3)] %>% droplevels()
which.dataset <- unique(data.all$pop_spec)[c(3:7)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 18)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(5, 6, 15, 16, 17, 18, 21)] %>% droplevels()
#which.dataset  <- unique(data.all$pop_spec)[c(15:21)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[-c(10:14)] %>% droplevels()
#which.dataset  <- unique(data.all$pop_spec)[12] %>% droplevels() 
data.all      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
sampling      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
}

## For dev and debug purposes also can subset total number of individuals 
 ## (done randomly though a seed is set in packages_functions.R)
red_ind    <- FALSE
if (red_ind) {
num_ind    <- 200
}

## Create the capture history scaffold from the raw data
source("data_manip.R")

## Create the indexing vectors and capture history structure needed for the stan model
source("data_stan.R")

## Deal with all of the individual level and population level covariates 
source("data_covariates.R")

## Processing of indices for the stan model to reduce looping for increasing computational speeds
source("stan_indices.R")

## And finally, created all of the necessary model matrices for the various linear predictors inside the model
if (length(which.dataset) != 1) {
source("establishing_mm.R")
}

## Quick look at a given population
#source("capt_plot.R")
#source("capt_plot_multi.R")

## And finally run the stan model
stan.iter     <- 800
stan.burn     <- 400
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin
stan.chains   <- 1
stan.cores    <- 1
stan.refresh  <- 10
if (length(which.dataset) == 1) {
source("stan_fit_single.R")
} else {
source("stan_fit_mm.R") 
}

## And some diagnostics and such
if (length(which.dataset) == 1) {
source("plotting.R")
} else {
source("multipop_plotting.R")
}
