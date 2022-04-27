#####################################
## Fit CMR model to amphibian data ##
#####################################

### To do rest of this week and next

## April 28 (Thursday)
 ## 1) D, E, F below

## April 29 (Friday)
 ## 1) G, H below
 ## 2) Read through reviews

## Next week:
 ## 1) Update overleaf with some new model thoughts
 ## 2) By end of week get overleaf to a point for sharing
 ## 3) Write clear steps for individual tasks for week in Point Pelee

#### Notes April 27 ---- 

## 1) Model matrix and no-loop version of multi-pop model seems to be about 8x faster than the old version
 ##    -- comparing the predictions from the individual model for RANA.SummitMeadow and the estimates from the 
 ##    -- multi-population are extremely similar, suggesting that the multi-pop model is built correctly

## 2) Also made some decent progress on code and repo cleaning, but still have some to do

## 3) Moving forward that next critical steps are to:

 ## -- D) [ ] clean up all of the fit_pops_XXXX.R scripts 
 ## -- E) [ ] add back in a daily average detection probability from which to calculate population size
 ## -- F) [ ] update all of the plotting scripts for the new model structure. Backup the old plotting scripts?

 ## -- G) [ ] run the individual populations alone and update individual_model_runs and upload new scripts to overleaf
 ## -- H) [ ] start making the rest of the models (see 3 below)



## 3) **Older note because I am not sure what is going to happen once the model matrix form is written out much better** But leaving here for now
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
some_pops  <- FALSE

if (some_pops) {
#which.dataset  <- unique(data.all$pop_spec)[c(4, 5, 8, 9, 11, 15, 16, 20)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[-c(10:14)] %>% droplevels()
which.dataset <- unique(data.all$pop_spec)[20] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(1, 2, 13)] %>% droplevels()
data.all      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
sampling      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
}

## For dev and debug purposes also can subset total number of individuals 
 ## (done randomly though a seed is set in packages_functions.R)
red_ind    <- FALSE
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
if (length(which.dataset) == 1) {
source("establishing_mm.R")
}

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
source("stan_fit_mm.R") 
}

## And some diagnostics and such
#source("diagnostics.R")
