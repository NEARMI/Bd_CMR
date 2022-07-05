#####################################
## Fit CMR model to amphibian data ##
#####################################

#### Notes July 5, 2022 ---- 

## 1) Worked today on plotting of MeHg fit and writing Results on Overleaf

## Next To Do:

## 1) Clean up Overleaf Methods section
 ##   -- [ ] Finish off section about population sizes. Make the extra supplemental file
 ##   -- [ ] Go over all of the supplemental figures and their captions
## 3) Overleaf Discussion section (mostly caveats and potential issues)
## 4) Some code and repo cleaning (mostly stan folders and stan files)

## Some potentially remaining issues to be cleaned up:

## 1) Maybe want to pursue a better Bd load level (e.g. a zero inflated regression for example)
  ## --> Mainly because scaled Bd makes it really weird to estimate survival when an individual is uninfected, because
   ##    with the current strategy "uninfected" doesn't really have a defined meaning as Bd is simply a continuous state
   ##    and anything greater than 2 sd from the mean in this continuous load isn't very sensible.
   ##   --> Doing a zero-inflated model would lead to more sensible estimates of what "uninfected" means
## 2) May want to put length into detection because of the funny length effect for Rana survival
## 3) Need to double check how I am calculating population sizes given the issue with 0 captures
 ## -- > May want to estimate less frequently than every day, maybe population size using average captures at the level of the random effect
   ##    (intersection of primary period and population)?
## 4) May want to try and fit the MeHg model with the other populations as well. Just because a low proportion of individuals get measured doesn't 
 ##      immediately mean it isn't enough to estimate the effect

########
#### Code ----
########

#### NOTE: In this file and all other files search *** for current choices that could potentially change

####
## Some choices of how this script will be used. These are parameters that must be adjusted when this script is run
####

## Can be false if just plotting of output is desired, which still requires the data cleaning, or plotting can also be false
 ## if saving a fit is the only desire
fit_model  <- FALSE
plot_model <- TRUE
## Flag to determine how plotting will proceed (after fitting, or from a saved model, or from extracted chains)
if (plot_model) {
plot_from <- {
  if (fit_model) {
    "fit"
  } else {
    "saved_model"
   # "saved_samples"
  }
}
}

## Print statements throughout to track progress if running from the command line or on an external server
if (fit_model) {
  print("A model will be fit")
} 
if (plot_model) {
  print("Plotting of model fit will occur")
}

## Name of the saved samples or saved model
if (plot_model) {
if (plot_from == "saved_model" | plot_from == "saved_samples") {
  if (plot_from == "saved_model") {
      saved_model <- "fits/stan_fit_multipop_all_full_2022-06-23.Rds" 
    # saved_model <- "fits/stan_fit_multipop_mehg_2022-06-30.Rds"
   if (file.exists(saved_model)) {
     print(paste("Plotting will occur using the saved model:", saved_model, sep = " "))
   } else {
     print("Wrongly named or missing saved fitted model")
     break
   }
  } else if (plot_from == "saved_samples") {
      saved_samples <- "samples/stan_multipop_samples_all.Rds" 
    # saved_samples <- "samples/stan_multipop_mehg_cleaned.Rds"
   print(paste("Plotting will occur using the saved chains:", saved_samples, sep = " "))
  } else {
   print("Plotting will not occur because of an unknown command or a missing file")
   break
  }
} else if (plot_from == "fit") {
  print("Plotting will be conducted from the model fit in this session")
}
}

## Packages and Functions
source("packages_functions.R")
source("../ggplot_theme.R")

## Read in data
source("data_load.R")

## Some choices to determine what model will be fit. 
 ## Can't be perfectly dynamic because populations are manually chosen
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
 ## 3) Fit individual-level MeHg? Be careful what populations to choose
fit_ind_mehg   <- FALSE
 ## 4) Reduced detection model? (if FALSE fits a random effect level for every day in every population)
red_p_model    <- TRUE

print("Fitting choices are:")
print(paste("sing_pop =", sing_pop, sep = " "))
print(paste("multi_spec =", sing_pop, sep = " "))
print(paste("multi_spec_red =", sing_pop, sep = " "))
print(paste("some_pops =", sing_pop, sep = " "))
print(paste("fit_ind_mehg =", sing_pop, sep = " "))
print(paste("red_p_model =", sing_pop, sep = " "))

## From these choices find the model to fit
source("determine_model.R")

## Some population choices:
 ## -10                                -- Full fit
 ## c(4, 5, 6, 15, 16, 17, 18, 19, 21) -- MeHg fit
 ## c(1:9, 11, 13, 15:21)              -- Full fit without Springfield and Scotia Barrens

## If a subset of populations, pick which ones
if (some_pops) {
which.dataset <- unique(data.all$pop_spec)[-10] %>% droplevels()
# which.dataset <- unique(data.all$pop_spec)[c(1:9, 11, 13, 15:21)] %>% droplevels()
# which.dataset <- unique(data.all$pop_spec)[c(4, 5, 6, 8, 9, 15, 16, 17, 18, 19, 21)] %>% droplevels()
data.all      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
sampling      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
}

## For dev and debug purposes also can subset total number of individuals 
 ## (done randomly though a seed is set in packages_functions.R)
red_ind    <- FALSE
if (red_ind) {
num_ind    <- 200
}

print(paste("From these choices, the following model will be [has been] fit:  ", which_stan_file, sep = ""))
print(paste("Using the following populations:"
  , as.character(which.dataset) %>% paste(collapse = " -- ")
  , sep = " "))
print(paste("With red_ind =", red_ind, sep = " "))

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

## Print more in-depth details about the dataset being fit, if its a single population, enough
 ## already printed to follow what is being fit if multi-pop
if (sing_pop) {
source("dataset_notes.R")
}

## Quick look at a given population
#source("capt_plot.R")
#source("capt_plot_multi.R")

## And finally run the stan model
stan.iter     <- 3500
stan.burn     <- 500
stan.thin     <- 3
stan.length   <- (stan.iter - stan.burn) / stan.thin
stan.chains   <- 1
stan.cores    <- 1
stan.refresh  <- 10
if (fit_model) {
if (length(which.dataset) == 1) {
source("stan_fit_single.R")
} else {
source("stan_fit_mm.R") 
}
}

## And some diagnostics and such
if (plot_model) {
if (length(which.dataset) == 1) {
source("plotting.R")
} else {
source("multipop_plotting.R")
}
}

