#####################################
## Fit CMR model to amphibian data ##
#####################################

#### Notes April 11 ---- 
# (note: moving series of older notes potentially interesting for Methods or Results sections to note_dump.txt)

## 1) Length imputation as a gamma regression with sex indeed working better than the alternative -- adopting for all models that
 ## can support it (not counting populations with no lengths for one sex)

## 2) Sex seems to help reasonably well for survival.
 ## -- I do wonder though if it will be sensible to add sex to detection (seems bias possible with only in survival -- big effects
  ## at present)

## 3) Testing with Blackrock complex, there seems to be essentially no difference in using an index-based or model-matrix based 
 ## method for categorical variables. 
  ## -- However, it seems that the correlation is lower between the levels, so probably will be helpful when there are many
  ## -- Yet, it seems that the model-matrix version is substantially slower. Jury is still out on why.

## 4) Running 200 individuals for each non-newt population on Yeti to make sure all of the indexing is still aligned with the
 ## model matrix version 
 ## ** [ ] Priority number one for tomorrow (April 12)

## 5) Running Jones Pond RALU with the model-matrix style MeHg-Bd interaction model 
 

#### Next stuff to do on this project ----

## 0) [ ] [Tuesday] code, model, and repo cleaning
## 1) [ ] [Wednesday Hopefully] Work through debugging of the various model attempts
## 2) [ ] [Wednesday Hopefully] Make some final decisions for model runs
## 4) [ ] [By Wednesday Night] Send 5 individual model and two big model (sans newts and with newts) fits to Yeti
## 5) [ ] [Friday] Debug some of the individual fits
## 6) [ ] [Friday] Build the script to run all populations individually
## 7) [ ] [Monday] Send the individual model script .sh to Yeti
## 8) [ ] [Monday] Debug the long fits
## 9) [ ] [Tuesday +] Start expanding the current non-MeHg individual model for Bd trajectories over time for PA newts
## 10) [ ] [Tuesday +] Start cleaning Overleaf and creating some fitting figures in prep for sending to PIs


#### Various other ToDos for later ----

## 1) Model modifications / adjustments
 ## A) [ ] Possible effect of "injury" on survival not yet added
 ## B) [ ] Possible use of sex in detection

## 2) Potential larger modifications
 ## A) [ ] Florida sampling scheme being fundamentally different to the other populations
   ## -- Current plan to maybe fit a completely different model
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
