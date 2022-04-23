#####################################
## Fit CMR model to amphibian data ##
#####################################

#### Notes April 22 ---- 

## Monday morning continue with the speedups. See all of the untracked stan models with XXXX_speedup.stan and also
  ## continue with **** -->> listed below. First step will be to fully remove the segment() from the 1 ~ bernoulli statements
    ## by creating index vectors in R so that there is no need for any looping in the stan model

## ^^ The best way to continue with this will be to create a new R script in the pipeline after data_stan.R that 
 ## creates the indices 

## 1) Running a few different detection models shows that sex in detection is pretty important to be able to ground 
 ##   background survival estimates
  ##   BUT STILL NEED TO  [Comment at Lunchtime Friday April 22]:
   ## -- [ ] adjust all of the CMR_single_populations models appropriately for sex in p
   ## -- [ ] update stan_readme to correctly describe the current files

## 2) In the midst of working to speed up the model, starting with the simplest model first without any imputation or anything.
 ## Lots of extra models and code currently being generated that will need to be cleaned up.
  ## Keeping notes on my laptop desktop for now.
  ## Waiting to push models until I figure out what is going on.
   ##  -- At the very least, it seems like doing some indexing in R to reduce the loops over phi, p, and bd will save maybe 30%
   ##  -- Semi-vectorizing the likelihood section (1 ~ bernoulli for phi and p), while still looping over n_ind and using the segment()
   ##     function really seems to slow things down.
   ##  **** -->> I think instead just coming up with all of the phi, p, and chi indices in R that go into each should, I think 

 ## ^^^ NEED TO TAKE REALLY NICE NOTES AND MAKE A README FOR ALL ATTEMPTED SPEEDUPS because this will be valuable information
  ## for working on the much larger models. Also need to make sure to name all models sensibly

## 3) I think the most viable strategy is going to be to fit a few different models that accommodate different populations.
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
#which.dataset <- unique(data.all$pop_spec)[-c(10:14)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[17] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(1, 2, 13)] %>% droplevels()
which.dataset <- unique(data.all$pop_spec)[17] %>% droplevels()
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
