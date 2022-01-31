#####################################
## Fit CMR model to amphibian data ##
#####################################

####
## Notes as of Jan 31:
####

## Some definite craziness in the multi-species model:
 ## 1) Am not going to be able to have a species fixed effect and a location random effect -- too correlated given that most species
  ## only are measured in one or two locations. --- At least the first attempt gave really wide CI for each species and population CI
   ## were more variable -- thus it seems that so far (without more location specific covariates) the population and not the species
    ## soaks up a lot of the variation. 
     ## --- Possibly with more location-specific covaraites this will turn out better, but will have to check back in later

## HOWEVER: looking at the fits of the populations individually, there is definitely some weird stuff happening with detection
 ## 1) Possibly putting too much emphasis on bd and detection, maybe need some other covaraites

## ------ There are definitely plenty of open questions about model complexity ------
 ## 1) In theory a fixed effect could be fit for each species for each coefficient (forms of detection, forms of survival), but
   ## its unclear if the model will actually fit
 ## 2) Also pretty uncertain about when to use species vs pop_spec vs location...

## Probably just need to set this down until I get the rest of the covariate data



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
