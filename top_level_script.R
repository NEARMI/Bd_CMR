#####################################
## Fit CMR model to amphibian data ##
#####################################

#### Notes April 25 ---- 

## 1) Day of working on speeding up code.
 ## -- Very successful in getting the model without lengths sped up with vectors. 
 ## -- Single population model done, but many remain unchecked
 ## -- Multiple population model compiling, but not yet run, hard to debug because of how long it takes. Will attempt overnight


## 2) Moving forward to fill all models: I think the most viable strategy is going to be to fit a few different models that accommodate different populations.
 ## A) The first of these models would utilize most of the multipop structure, but drop the
  ##    `by species' fixed effects, and just use random effects for variation among the populations, so
  ##    in essence all that is really needed is to collapse the multipop model a bit.
  ##    ANBO, RANA, [and some NOVI?] can be fit in this way
 ## B) The second model would be an adjusted Bd-MeHg interaction model that would incorporate
  ##    all populations in which there are enough MeHg measured individuals to impute
 ## C) The third model is a simple regression-form Bd over time for NOVI
 ## D) The fourth, which could be a completely separate paper would be a better disease model


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
