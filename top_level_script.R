#####################################
## Fit CMR model to amphibian data ##
#####################################

#### Notes June 9, 2022 ---- 

## 1) Ahhh... Backtracking a bit and now have some messy code. It really is quite unsatisfying to be using fixed effects and possibly just 
 ## wrong (maybe cheating in a way)
  ## A) So now have two scripts for establishing_mm and stan_fit_mm.
   ##    --- these will eventually need to be collapsed with some if statements

## 2) Found an error in a model matrix setup that caused much data to be unused in the fitting of the multi-species model
 ## -- While I don't expect miracles (I think fixed effects at the species level are still going to be a problem), it does feel worth it
  ##   to rerun this model corrected before moving forward to far
   ## -- This is esspecially true because the model I have been banging my head aganist with the random vs fixed effects specification is not
   ##    a 'final' model anway (as it includes only a single species)

## 3) So the ToDo list has changed a bit from yesterday
 ## A) Get the multi-population model fit and diagnose issues
 ## B) From there update the Overleaf and decide what the next steps should be
 ## C) Still do need to figure out what to do with the continuous Newt populations though...

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
which.dataset <- unique(data.all$pop_spec)[c(1:9, 11, 13, 15:21)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(15:21)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(3)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(3:7)] %>% droplevels()
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
# source("establishing_mm.R")
# source("establishing_mm_rand.R")
source("establishing_mm_rand_red.R")
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
