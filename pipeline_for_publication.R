#################################################################
## Fit CMR model to amphibian data                             ##
## Specific pipeline for just the two fits for the publication ##
#################################################################

## Packages and Functions
source("packages_functions.R")
source("ggplot_theme.R")

## Read in data
source("complete_data.R") 

## Choices for what model is fit
sing_pop         <- FALSE ## Single population?
multi_spec       <- TRUE  ## Multiple species?
multi_spec_red   <- FALSE ## If multiple species, fit a reduced model with no species-specific fixed effects?
some_pops        <- TRUE  ## Not all populations?

## ONLY THING THAT CHANGES BETWEEN THE TWO PUB FITS
fit_ind_mehg     <- FALSE ## Fit individual-level MeHg? Be careful what populations to choose

fit_only_mehg    <- FALSE ## Fit a model only predicting MeHg? (no survival) Setting as TRUE invalidates many other options pertaining to survival
red_p_model      <- TRUE  ## Reduced detection model? (if FALSE fits a random effect level for every day in every population)
red_ind          <- FALSE

if (fit_ind_mehg) {
  which_stan_file <- "CMR_multiple_populations_mehg"
  pops_to_fit     <- c(4, 5, 12, 13, 14, 15, 16, 17, 18, 19)
} else {
  which_stan_file <- "CMR_multiple_populations_alt_p_len"
  pops_to_fit     <- -c(4, 5, 12, 13, 14, 15, 16, 17, 18, 19)
}

which.dataset <- unique(data.all$pop_spec)[pops_to_fit] %>% droplevels()
data.all      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
sampling      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()


## Create the capture history scaffold from the raw data
source("data_manip.R")

## Create the indexing vectors and capture history structure needed for the stan model
source("data_stan.R")

## Deal with all of the individual level and population level covariates 
source("data_covariates.R")

## Processing of indices for the stan model to reduce looping for increasing computational speeds
source("stan_indices.R")

## And finally, created all of the necessary model matrices for the various linear predictors inside the model
source("establishing_mm.R")

## And finally run the stan model
stan.iter     <- 3000
stan.burn     <- 1000
stan.thin     <- 2
stan.length   <- (stan.iter - stan.burn) / stan.thin
stan.chains   <- 1
stan.cores    <- 1
stan.refresh  <- 10
source("stan_fit_mm.R")

## Summary of model fits, figures, etc.
## NOTES!!!
## 1) Probably wont work top to bottom with source, so open and run what you need/want
## 2) Runs from a saved model output (loads it in as a .Rds -- should be reasonably self explanatory)
## "figures_functions.R"

