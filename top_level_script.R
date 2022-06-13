#####################################
## Fit CMR model to amphibian data ##
#####################################

#### Notes June 13, 2022 ---- 

## 1) At least the reduced multi-spec fit seems to be fine (reduced meaning no species-specific fixed effects).
## 2) Unclear about model with species-specific fixed effects, model running now

## 3) So the ToDo list has changed again
 ## A) Check model with species-specific fixed effects
 ## B) Adjust other remaining model to whichever seems better (fixed effects or not)
 ## C) Figure out what to do with Newt populations
 ## D) Run all final models
 ## E) Update publication style overleaf

#### Code ----

#### NOTE: In this file and all other files search *** for current choices that could potentially change

## Packages and Functions
source("packages_functions.R")
source("../ggplot_theme.R")

## Read in data
source("data_load.R")

## Some choices to determine what model will be fit. Can't be perfectly dynamic because populations are manually chosen
 ## and some choices won't work with certain models, so will need to double check. Mismatches will lead to errors

 ## 0) Single population?
sing_pop       <- FALSE
 ## 1) Multiple species?
multi_spec     <- TRUE
 ## 1.2) If multiple species, fit a reduced model with no species-specific fixed effects?
multi_spec_red <- TRUE
 ## 2) Not all populations?
some_pops      <- TRUE
 ## 3) Fit individual-level MeHg?
fit_ind_mehg   <- FALSE

## From these choices find the model to fit
source("determine_model.R")
print(paste("Model to fit:  ", which_stan_file, sep = ""))

## If a subset of populations, pick which ones
if (some_pops) {
which.dataset <- unique(data.all$pop_spec)[c(1:9, 11, 13, 15:21)] %>% droplevels()
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
if (multi_spec) {
source("establishing_mm.R")
}

## Print some details about the dataset being fit, if its a single population
if (sing_pop) {
source("dataset_notes.R")
}

## Quick look at a given population
#source("capt_plot.R")
#source("capt_plot_multi.R")

## And finally run the stan model
stan.iter     <- 1000
stan.burn     <- 500
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
