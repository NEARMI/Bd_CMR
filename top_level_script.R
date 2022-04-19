#####################################
## Fit CMR model to amphibian data ##
#####################################

#### Notes April 18 ---- 

## 1) After substantial exploration of Bd across Julian dates and temperatures, it seems pretty clear that any control
 ## for temporal variation in sampling to model Bd over time isn't going to happen. This is ok for the main analysis generally,
  ## but does raise a few potential issues:
 ## A) If Scotia Barrens (and to a lesser extent SPR) gets collapsed (time), there is a chance for bias in individual predictions.
  ##   not exactly sure what to do about this...
   ##   -- The model is already too complicated so wanting to drop covariates and populations anyway...

## 1) Unfortunately fits timed out after three days... With only 200 individuals per population and no Eastern Newts...
 ## Looks like with 0.92 adapt delta and 13 treedepth will take about 3.5 days to fit.
  ## Which makes a full model close to impossible with current computational resources. Will need to make some sacrifices
   ## in terms of model complications (probably in terms of random effects)

## 2) Next step (Tuesday April 19 will be to insure that the model with 200 individuals "worked") and then figure out
 ## how to simplify the model to make it tractable

## 3) While that model is running need to progress with the individual variation over time model (PA) by adjusting the indices
 ## to model continuous Bd and summarize that Bd into max for estimates of between-season survival.
  ## -- Initial model finished and compiles. Need to run and debug

#### Next stuff to do on this project ----

# x = partial; X = done for now

## 0) [x] [Tuesday] code, model, and repo cleaning
  ##        -- got a reasonable start, some more will be needed when debugging continues
## 1) [x] [Tuesday] Figure out what to do with FL and update the spreadsheet
  ##        -- some solid progress here but need to check in with Evan and Dave to make sure
## 1) [ ] [Later] Work through debugging of the various model attempts
  ##        -- would have liked to on Wednesday, but the models are still running, and it looks like they may time out
## 2) [ ] [After debugging (1 above)] Make some final decisions for model runs
## 3) [ ] [After final decisions on model structure (2 above)] Send individual (with .sh file) and multi-pop jobs
## 4) [ ] [Eventually] Debug the long fits
## 5) [ ] [Eventually] Start expanding the current non-MeHg individual model for Bd trajectories over time for PA newts
## 6) [ ] [Eventually] Start cleaning Overleaf and creating some fitting figures in prep for sending to PIs


#### Various other ToDos for later ----

## 1) Model modifications / adjustments
 ## A) [ ] Possible effect of "injury" on survival not yet added
 ## B) [ ] Possible use of sex in detection
 ## C) [ ] Finer grained temperature data??

## 2) Potential larger modifications
 ## A) [x] Florida sampling scheme being fundamentally different to the other populations
   ## -- Think this gets fixed with defining a reasonable gap between samples
 ## B) [ ] What individuals were actually marked for NOVI for the Wisconsin populations
   ## -- Seems like all? But still need to email PI


#### Code ----

#### NOTE: In this file and all other files search *** for current choices that could potentially change

## Packages and Functions
source("packages_functions.R")
source("../ggplot_theme.R")

## Read in data
source("data_load.R")

## For dev and debug purposes pick a subset of locations
some_pops  <- TRUE

if (some_pops) {
#which.dataset <- unique(data.all$pop_spec)[-c(10:14)] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[17] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(1, 2, 13)] %>% droplevels()
which.dataset <- unique(data.all$pop_spec)[12] %>% droplevels()
data.all      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
sampling      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
}

## For dev and debug purposes also can subset total number of individuals 
 ## (done randomly though a seed is set in packages_functions.R)
red_ind          <- TRUE
red_ind_PA_debug <- TRUE ## Temp debug switch, will integrate if model works
if (red_ind) {
num_ind    <- 250
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
stan.iter     <- 550
stan.burn     <- 250
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin
if (length(which.dataset) == 1) {
source("stan_fit_single.R")
} else {
source("stan_fit.R") 
}

## And some diagnostics and such
#source("diagnostics.R")
