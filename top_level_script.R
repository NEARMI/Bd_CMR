#####################################
## Fit CMR model to amphibian data ##
#####################################

#### Notes June 6, 2022 ---- 

## Returning to the project after a two week break to work on some other things. 

## 0.1) Today played with a simple simulation to find that the pattern that the intercept is estimated with
  ## increasing error with higher variance in the random effects holds with glmer. Not quite sure what to do about this
## 0.2) Also reformulated the model to use Cholesky. Fitting now. Seems marginally faster but need to compare coefficient estimates
## 0.3) Did some research on multiple random effects with different grouping variables, but it seems that my current strategy is fine

## 1) Main issue now is to try and figure out a way to specify the joint model that avoids the uncertain 
 ## intercept problems I am facing.
  ## -- Some thoughts in model_transformation.txt 

## 2) Other issues:
 ## -- Figure out what to do with Scotia Barrens

#### Code ----

#### NOTE: In this file and all other files search *** for current choices that could potentially change

## Packages and Functions
source("packages_functions.R")
source("../ggplot_theme.R")

## Read in data
source("data_load.R")

## For dev and debug purposes pick a subset of locations
some_pops  <- TRUE

# data.all %>% group_by(Year, pop_spec, Mark) %>% filter(BdSample == "Y") %>% summarize(nswab = n()) %>% arrange(desc(nswab)) %>% as.data.frame()       

if (some_pops) {
#which.dataset <- unique(data.all$pop_spec)[c(1:9, 11, 13, 15:21)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(15:21)] %>% droplevels()
which.dataset <- unique(data.all$pop_spec)[c(3:7)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 18)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(5, 6, 15, 16, 17, 18, 21)] %>% droplevels()
#which.dataset  <- unique(data.all$pop_spec)[c(15:21)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[-c(10:14)] %>% droplevels()
#which.dataset  <- unique(data.all$pop_spec)[12] %>% droplevels() 
data.all      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
sampling      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
}

## For dev and debug purposes also can subset total number of individuals 
 ## (done randomly though a seed is set in packages_functions.R)
red_ind    <- FALSE
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
if (length(which.dataset) != 1) {
source("establishing_mm.R")
}

## Quick look at a given population
#source("capt_plot.R")
#source("capt_plot_multi.R")

## And finally run the stan model
stan.iter     <- 800
stan.burn     <- 400
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin
stan.chains   <- 1
stan.cores    <- 1
stan.refresh  <- 10
if (length(which.dataset) == 1) {
source("stan_fit_single.R")
} else {
source("stan_fit_mm.R") 
}

## And some diagnostics and such
if (length(which.dataset) == 1) {
source("plotting.R")
} else {
source("multipop_plotting.R")
}
