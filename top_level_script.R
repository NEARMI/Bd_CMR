#####################################
## Fit CMR model to amphibian data ##
#####################################

#### Notes May 4 (continuation of To Do list from April 28) ---- 

## 0.1) [Over Time] Fit the next trial models and compare to the standard population
  ##      -- [ ] CMR_multiple_populations_mehg_ssp
  ##      -- [ ] CMR_multiple_populations_mehg
  ##      -- [ ] CMR_multiple_populations_ssp

## 0.2) Priorities prior to Point Pelee
 ## -- [Thurs] Send all non-continuous fits individually
 ## -- [Thurs] Send two Newt populations with continuous fits
 ## -- [Thurs] Send a bigger multi-pop model
 ## -- [Fri]   Clean up code and repo and make plan

## 0.3) Priorities for the following week (Point Pelee)
 ## -- Debug fits
 ## -- Upload new figures to overleaf
 ## -- Update overleaf writing
 ## -- Finish next models (E below)
 ## -- Start sketching more complicated disease/CMR model

## 1) Moving forward that next critical steps are to:
 ## -- A) [ ] Run all of the individual populations alone on Yeti
  ##           -- [ ] update individual_model_runs
  ##           -- [ ] upload new plot fits to overleaf
 ## -- B) [ ] Run a multi pop fit
 ## -- C) [ ] Create a number of subsidiary models from the master models for other fits
  ##           -- [x] Create a single_population model that has Bd load in detection to see if it can solve some of the crazy patterns in RANA
    ##              -- Fits fine, doesn't help at all
  ##           -- [-in progress-] A multi-pop model that has only one species
    ##              -- [ ] Compiles, but going to be slow, need to run on other computer for debugging
  ##           -- [-in progress-] Multi-pop MeHg model that can support all populations in which MeHg can be estimated
    ##              -- [ ] Compiles, but going to be slow, need to run on other computer for debugging 
  ##           -- [-in progress-] Multi-pop single species MeHg model (combination of the above 2)
    ##              -- [ ] Compiles, but going to be slow, need to run on other computer for debugging
  ##           -- [-done for now-] Refine the ``Simple'' continuous time Bd model
    ##              -- could conceivably modify to use some form of temp, but probably will save that for the model below
  ##           -- [ ] {Probably an entirely different paper} Write an actual disease model, which is informed by CMR style data
    ##              -- 

#### Code ----

#### NOTE: In this file and all other files search *** for current choices that could potentially change

## Packages and Functions
source("packages_functions.R")
source("../ggplot_theme.R")

## Read in data
source("data_load.R")

## For dev and debug purposes pick a subset of locations
some_pops  <- TRUE

if (some_pops) {
#which.dataset <- unique(data.all$pop_spec)[c(3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 18)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(5, 6, 15, 16, 17, 18, 21)] %>% droplevels()
which.dataset  <- unique(data.all$pop_spec)[c(15:21)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[-c(10:14)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[4] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(1, 2, 13)] %>% droplevels()
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
stan.iter     <- 1000
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
