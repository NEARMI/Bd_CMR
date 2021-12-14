#####################################
## Fit CMR model to amphibian data ##
#####################################

####
## Notes as of Dec 10:
####

## Received some new data Dec 8. Adjusting code a bit and trying fits with each data set to explore the variation
 ## in parameter values across the data sets

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
which.dataset <- 1
data.all %<>% filter(dataset == which.dataset)
source("data_manip.R")

## Create the indexing vectors and capture history structure needed for the stan model
source("data_stan.R")

## Deal with all of the individual level and population level covariates 
source("data_covariates.R")

## Quick look at a given population
source("capt_plot.R")

## And finally run the stan model
stan.iter     <- 2000
stan.burn     <- 500
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin
source("stan_fit.R")

## And some diagnostics and such
source("diagnostics.R")




