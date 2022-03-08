#####################################
## Fit CMR model to amphibian data ##
#####################################

####
## Notes as of March 8:
####

## Step 1 for tomorrow is to double check the few missing swab values in the adams data set
 ## and to figure out what to do with the dpulicate swabs

## Code seems to be working, simple model seems to be working pretty well for individual populations (at least RALU so far)
 ## Next critical piece is to update the script the fits each and every population individually and then run it.
  ## -- The steps to get this running will be: 
   ## 1) Fill out a lookup table for what covaraties can be fit for different populations
   ## 2) Dynamically build the stan model for each population based on this lookup table

## Summary of notes from the past few days:
 ## DATA:
  ## -- Many unresolved questions still with Brian T
  ## -- Jill working on Mercury samples for 2019 and 2020
  ## -- Injury covariate and drop dead individuals

 ## MODEL:
  ## -- Most importantly it is still a pretty big question of what to do with sub-populations
   ##   (opens a can of worms about meta-population, super-population, movement between subsites etc.)
   ##   The major problems are about bias in detection because of assuming individuals are potentially found when they really are not
   ##     (because they are in a different subpopulation)
   ##   SEE OLDER NOTES FROM A PREVIOUS COMMIT FOR A LONGER DISCUSSION
  ## -- And a slightly less major question about how to deal with primary periods and secondary periods vs continuous times between captures
  ## -- Still many open questions about what covariates to use and how to merge discrete covariates
  ## -- Potential mercury latent model
  ## -- Multiple imputation for unobserved covariates


## Packages and Functions
source("packages_functions.R")

## Read in data
source("data_load.R")

## Construct modified data frame of recapture histories for each individual in each population
 ## For single species debug purposes pick a single data set
single_pop <- FALSE

if (single_pop) {
which.dataset <- unique(data.all$pop_spec)[12]
data.all      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
sampling      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
}

red_ind    <- FALSE
if (red_ind) {
num_ind    <- 150
}

source("data_manip.R")

## Create the indexing vectors and capture history structure needed for the stan model
source("data_stan.R")

## Deal with all of the individual level and population level covariates 
source("data_covariates.R")

## Quick look at a given population
source("capt_plot.R")
#source("capt_plot_multi.R")

## And finally run the stan model
stan.iter     <- 1500
stan.burn     <- 500
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin
if (single_pop) {
source("stan_fit_single.R")
} else {
source("stan_fit.R") 
}

## And some diagnostics and such
source("diagnostics.R")
