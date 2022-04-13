#####################################
## Fit CMR model to amphibian data ##
#####################################

#### Notes April 13 ---- 
# (note: moving series of older notes potentially interesting for Methods or Results sections to note_dump.txt)

## 1) Few coding and organization updates:
 ## A) Plotting script updated for various models that will be fit in individual populations
 ## B) Spreadsheet updated to list the individual models to fit by individual population
 ## C) Script for running all populations individually on Yeti complete, but waiting to test because of low
  ##   priority running the other jobs
 ## D) In the end, the model-matrix and index-version for the MeHg-Bd interaction model took an identical amount of time to run
  ##   both fit identically and all indexing checked out

## 2) Yeti runs with model-matrix and index-version moving glacially slow. May not even finish in the 
 ##   requested time of three days. Will debug when complete, but may need to resend

## 3) Given the glacial pace of the multi-population fits on Yeti, at least can check the indexing of both model versions running with
 ##   basically no individuals in each population


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
which.dataset <- unique(data.all$pop_spec)[17] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(1, 2, 13)] %>% droplevels()
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
