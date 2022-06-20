#####################################
## Fit CMR model to amphibian data ##
#####################################

#### Notes June 20, 2022 ---- 

## 1) Job below failed because of ram issues. But ran about a full day faster than expected.
 ##   -- Submitted two new jobs with a subset of parameters that shouldn't have a ram issue

## 2) Spent most of the day on the Introduction.
 ## -- To do list for Overleaf as follows:
  ## A) Adjust detection model description to match the new random effects structure
  ## B) Drop the description of the continuous time model, will save that for a different paper
  ## C) Expand out the Results with some of the placeholder figures, can still write about these before the full runs are done
  ## D) expand the list of caveats in the Discussion

## -- Sent a job only making the small change of the reduced form P model. Could be enough to get the model to run,
 ## but will have to see. 
  ## -- An initial slow fit shows that 
   ##  -- 18 pops without PA and MA:
     ##     -- 1600 Seconds per chain:
      ##       -- 31k, 14k --> 78, 28 seconds per sample -- Ratio of 0.04875, 0.0175 relative to seconds estimate
      ##       -- 33k, 23k --> 83, 46 seconds per sample -- Ratio of 0.051875, 0.02875 relative to seconds estimate
      ##       -- 38k, 26k --> 95, 52 seconds per sample -- Ratio of 0.059375, 0.0325 relative to seconds estimate
   ##  -- 20 pops with these:
     ##     -- 10000 Seconds per chain:
      ##       -- Time estimate of about [at worst?]:
       ##          -- 0.059375 * 10000 per, so with 400 warmup -- 66h, with desired 500   -- 82h
       ##          -- 0.0325   * 10000 per, so with 500 samples -- 45h, with desired 3000 -- 270h
      ##    In total expect 4.625 days, so finishing ~ Monday at Midnight?
      ##      -- If on pace, expect warmup to be done ~ Sunday at about 10pm. Can check this
      ##    For the final complete model ~ 15 days, 11 days with 2000 samples and 2 thin
  ## -- Actually ended up fitting in 78h, but failed to allocate enough memory (23.5gb)

## 2) So the full ToDo list is now basically:
 ## A) Figure out what to do with Newt populations
 ## B) Run all final models
 ## C) Update publication style overleaf

#### Code ----

#### NOTE: In this file and all other files search *** for current choices that could potentially change

## Packages and Functions
source("packages_functions.R")
source("../ggplot_theme.R")

## Read in data
source("data_load.R")

## Some choices to determine what model will be fit. Can't be perfectly dynamic because populations are manually chosen
 ## and some choices won't work with certain models, so will need to double check. Mismatches will lead to errors
  ## List of populations that can be fit with individual-level mercury listed in "determine_model.R"

 ## 0) Single population?
sing_pop       <- FALSE
 ## 1) Multiple species?
multi_spec     <- TRUE
 ## 1.2) If multiple species, fit a reduced model with no species-specific fixed effects?
multi_spec_red <- FALSE
 ## 2) Not all populations?
some_pops      <- TRUE
 ## 3) Fit individual-level MeHg?
fit_ind_mehg   <- FALSE
 ## 4) Reduced detection model? (if FALSE fits a random effect level for every day in every population)
red_p_model    <- TRUE

## From these choices find the model to fit
source("determine_model.R")
print(paste("Model to fit:  ", which_stan_file, sep = ""))

## If a subset of populations, pick which ones
if (some_pops) {
which.dataset <- unique(data.all$pop_spec)[-10] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(1:9, 11, 13, 15:21)] %>% droplevels()
# which.dataset <- unique(data.all$pop_spec)[c(4, 5, 6, 15, 16, 17, 18, 19, 21)] %>% droplevels()
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
if (multi_spec) {
source("establishing_mm.R")
}

## Print some details about the dataset being fit, if its a single population
if (sing_pop) {
source("dataset_notes.R")
}

## Quick look at a given population
#source("capt_plot.R")
#source("capt_plot_multi.R")

## And finally run the stan model
stan.iter     <- 900
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
