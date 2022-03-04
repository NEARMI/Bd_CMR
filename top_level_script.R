#####################################
## Fit CMR model to amphibian data ##
#####################################

####
## Notes as of March 4:
####

## 1) Now have all of the data. Still have a number of questions about it, but it seems close to clean
 ## A) Update and strip the site level data 
 ## B) Need to send to Brian T for feedback and QA/QC on the whole dataset

## 2) The very first step with the new data is to update the code for the new [final?] structure to load
 ## and join all of the data files (AS WELL AS the site covaraite data).
  ## A) In short, make sure I can make it all the way to running a model (for now not worrying to much about
   ## how specifically sensible that model is)

## A few new things:
 ## -- dead
 ## -- injured

  ## B)
   ## -- After this need to extract the flagged rows and send those to Brian T
   ## -- Note: there is a list of questions for Brian T and Evan/Dave in Mac Notes

## 3) For trying the new data, the first step is to try and fit a few individual population models: 
 ## -- Best to use RALU which is well behaved and then go from there. A sensible population to start with will be
  ## LostHorse (LH_RALU)
 ## For a mid-sized population try the three strategies of fitting within vs between season survival
  ## 1. Continuous time survival
  ## 2. Many Primary Periods in a year
  ## 3. Primary Period just means year

####
## Comments from the previous few commits:
####

## Received data from Brian on site-level covariates and MeHg. Need to link these and bring them into the model

## 1) Absolutely no chance of fitting a pop * year random effect for survival _apart_
 ## from the effect of bd. Going to need to simplify the model quite a bit

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
