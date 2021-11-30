#####################################
## Fit CMR model to amphibian data ##
#####################################

## Packages and Functions
source("packages_functions.R")

## Read in data
source("data_load.R")

## Construct modified data frame of recapture histories for each individual in each population
source("data_manip.R")

## Create the indexing vectors and capture history structure needed for the stan model
source("data_stan.R")

## Deal with all of the individual level and population level covariates 
source("data_covariates.R")

## And finally run the stan model
stan.iter     <- 2000
stan.burn     <- 500
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin
source("stan_fit.R")

## And some diagnostics and such
source("diagnostics.R")




