#####################################
## Fit CMR model to amphibian data ##
#####################################

#### Notes May 5 (continuation of To Do list from April 28) ---- 

## New models seem good, unclear of their utility but definitely function
## Code and repo cleaned up 

## 0.1) Priorities prior to Point Pelee
 ## -- [Thurs] [ ] Create a Number spreadsheet of all desired fits and what jobs are needed
 ## -- [Thurs] [ ] Clean up code and repo 
 ## -- [Fri]   [ ] Send all non-continuous fits individually
 ## -- [Fri]   [ ] Send two Newt populations with continuous fits
 ## -- [Fri]   [ ] Send a series of multi-pop models

## 0.2) Priorities for the following week (Point Pelee)
 ## -- Debug fits
 ## -- Upload new figures to overleaf
 ## -- Update overleaf writing
 ## -- Start sketching more complicated disease/CMR model

## 1) Moving forward, the next overall critical steps are to:
 ## -- A) [ ] Run all of the desired runs on Yeti
 ## -- B) [ ] Make sense of the output and choose and upload some new plots to overleaf
 ## -- C) [ ] Form a plan for the direction of this paper
 ## -- D) [ ] Start on a more expansive disease CMR model for the newt data

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
which.dataset <- unique(data.all$pop_spec)[c(3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 18)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(5, 6, 15, 16, 17, 18, 21)] %>% droplevels()
#which.dataset  <- unique(data.all$pop_spec)[c(15:21)] %>% droplevels()
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
