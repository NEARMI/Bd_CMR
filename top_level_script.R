#####################################
## Fit CMR model to amphibian data ##
#####################################

#### Notes April 20 ---- 

## 1) Model compiling for every population on its own. Next step is to clean up the Yeti batch script and send jobs
 ## A) A few covaraites possibly subject to change, inclduing:
  ##    -- Sex in detection
  ##    -- Bd indetection
  ##    -- Injury
 ## B) Potential desire to increase the complexity of the continuous time Bd model
  ##    -- Some form of temperature instead of date
  ##    -- Sex in detection
  ##    -- Bd in detection
  ##    -- Hurdle model?

## 2) Still need to make a few final decisions about the joint population model in order to get it to run
 ## A) biggest issue is what to do with Newts and non-newts
 ## B) covariates to simplify model
 ## C) speed increases?
 ## D) non-centered random effect issue with multiple "intercepts"

## 3) Run models, debug, and share Overleaf

#### Code ----

#### NOTE: In this file and all other files search *** for current choices that could potentially change

## Packages and Functions
source("packages_functions.R")
source("../ggplot_theme.R")

## Read in data
source("data_load.R")

## For dev and debug purposes pick a subset of locations
some_pops  <- FALSE

if (some_pops) {
#which.dataset <- unique(data.all$pop_spec)[-c(10:14)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[17] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(1, 2, 13)] %>% droplevels()
which.dataset <- unique(data.all$pop_spec)[4] %>% droplevels()
data.all      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
sampling      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
}

## For dev and debug purposes also can subset total number of individuals 
 ## (done randomly though a seed is set in packages_functions.R)
red_ind          <- FALSE
# red_ind_PA_debug <- FALSE ## Temp debug switch, will integrate if model works
if (red_ind) {
num_ind    <- 200
}

## Create the capture history scaffold from the raw data
source("data_manip.R")

## Create the indexing vectors and capture history structure needed for the stan model
source("data_stan.R")

## Deal with all of the individual level and population level covariates 
source("data_covariates.R")

## Quick look at a given population
#source("capt_plot.R")
#source("capt_plot_multi.R")

## And finally run the stan model
stan.iter     <- 1000
stan.burn     <- 400
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin
if (length(which.dataset) == 1) {
source("stan_fit_single.R")
} else {
source("stan_fit.R") 
}

## And some diagnostics and such
#source("diagnostics.R")
