#####################################
## Fit CMR model to amphibian data ##
#####################################

## ** Ended the day working on (C) figuring out model matrix strategy for the phi and p linear predictors
## ** Next steps:
 ## A) Get this model matrix setup working and check output
 ## B) Get the new version of the length imputation integrated into all models that do length imputation
 ## C) Code and repo cleaning for updated models (currently lots needed)

#### Notes April 26 ---- 

## 1) Another day of working on speeding up code. Mutli-population model is still really slow. 
 ##   Aim today is to figure out what the culprits are and fix them
 ## A) [x] MeHg part of the model is fine and runs very quickly
 ## B) [x] Length part of the model was not ok but is fine now after fixing model matrix issues. Best to just use normal model
 ##    [ ] Model updated (see len_trial_normal2.stan) but these update still need to be applied to all of the actual models
 ## C) [ ] Lots of problems remain with the specification of the phi and p linear predictors. 
 ##      -- Both need to be converted to model matrix form, which is difficult for the mix of categorical * continuous covaraites that have
 ##         both fixed and random effects. I think it will be possible to not use a _full_ (and thus completely opaque) random effect specification
 ##         but instead break this up by int, and each covaraite. Currently in progress

## 2) **Older note because I am not sure what is going to happen once the model matrix form is written out much better** But leaving here for now
 ##   Moving forward to fill all models: I think the most viable strategy is going to be to fit a few different models that accommodate different populations.
 ## A) The first of these models would utilize most of the multipop structure, but drop the
  ##    `by species' fixed effects, and just use random effects for variation among the populations, so
  ##    in essence all that is really needed is to collapse the multipop model a bit.
  ##    The goal then would be to still fit all of the populations together, but simply with less stuff to fit to speed things up
  ##    [ ] remove all `by species' fixed effects
  ##    [ ] check
 ## B) The second option would be to use the above model but fit different species by themselves (at least RANA and ANBO)
  ##    [ ] try
 ## C) The second actually different model would be an adjusted Bd-MeHg interaction model that would incorporate
  ##    all populations in which there are enough MeHg measured individuals to impute
  ##    [ ] Expand the single population model with a single fixed effect and multiple random effects
 ## D) The third model is a simple regression-form Bd over time for NOVI
  ##    [ ] Nearly finished but needs a bit of debugging
 ## E) The fourth, which could be a completely separate paper would be a better disease model
  ##    [ ] Much work to be done here


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
which.dataset  <- unique(data.all$pop_spec)[c(1, 2, 4, 5, 8, 9, 15, 16)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[-c(10:14)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[17] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(1, 2, 13)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[17] %>% droplevels()
data.all      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
sampling      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
}

## For dev and debug purposes also can subset total number of individuals 
 ## (done randomly though a seed is set in packages_functions.R)
red_ind    <- TRUE
# red_ind_PA_debug <- FALSE ## Temp debug switch, will integrate if model works
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
#source("establishing_mm.R")

## Quick look at a given population
#source("capt_plot.R")
#source("capt_plot_multi.R")

## And finally run the stan model
stan.iter     <- 700
stan.burn     <- 300
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin
if (length(which.dataset) == 1) {
source("stan_fit_single.R")
} else {
source("stan_fit.R") 
}

## And some diagnostics and such
#source("diagnostics.R")
