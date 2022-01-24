#####################################
## Fit CMR model to amphibian data ##
#####################################

####
## Notes as of Jan 21:
####

## With the current model population variation in between season survival as a function of Bd load is
 ## not being recovered. It may be sensible to try and fit a fixed effect of species and a random
  ## effect of location as a joint random effect of species:location isn't cutting it

## 


## Weirdly, when fitting individual populations one at a time the same species in different places sometimes
 ## produces pretty different responses (e.g., ANBO sometimes is estimated to have increasing "survival" with higher bd
  ## load and sometimes decreasing "survival")
 ## -- The solution to this is likely to be to play around with either a random effect or fixed effect for species and
  ## a random effect for location (though it may be hard to pull these effects apart)


### --- Some holdover notes from before --- ###

## Some model cleanup:
 ## some minor detection model intercept issues with covaraite for year and month present as well

## Next steps overall
 ## 1) Continue poking around the detection model and survival within and between season with the frogs and the newts. 

## Some next steps for this script (not in chronological order):
 ## 1) With newts basically not working when dropping temporal component and frogs having no chance of fitting
  ## a temporal component the first logical next step is to try and fit a combined model
 ## 2) Need to try multi-pop model with more pops that don't have a temporal component
  ## -- basically not possible until I see much more of the data however
 ## 3) Possibly try conditional infection given that some individuals have 0 load when measured

## Packages and Functions
source("packages_functions.R")

## Read in data
source("data_load.R")

## Construct modified data frame of recapture histories for each individual in each population
 ## For single species debug purposes pick a single data set
single_pop <- TRUE

if (single_pop) {
which.dataset <- unique(data.all$pop_spec)[4]
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

