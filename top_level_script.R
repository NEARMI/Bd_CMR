#####################################
## Fit CMR model to amphibian data ##
#####################################

## Packages and Functions
source("packages_functions.R")

## Read in data for Jones Pond RALU CORT Analysis
source("JP_data.R")

## Further data cleaning and organization for stan model
source("data_manip.R")

## Convert this master data frame created in data_manip.R into the various
 ## data objects needed to fit the stan model
source("data_stan.R")

## More covariate prep
source("data_covariates_JP_cort.R")

## Finalizing indices for long-form Stan model
source("stan_indices.R")

## And finally run the stan model
source("stan_fit_single_JP_cort.R")

## Model fit exploration
 ## Figures and Tables for manuscript 
  ## NOTE: You will be better served opening this script and run code line-by-line
   ## instead of sourcing
source("explore_model_fit.R")
