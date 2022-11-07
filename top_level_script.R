#####################################
## Fit CMR model to amphibian data ##
#####################################

####
## Some Project Notes
####

## 1) In this file and all other files search *** for current model assumptions
## 2) For some notes on the current model see "Project_Notes.R"

####
## Some choices of how this script will be used. These are parameters that must be adjusted when this script is run
####

fit_model  <- TRUE   ## Fit a model?           (can summarize and plot previous model fit)
plot_model <- TRUE   ## Plot the model output? (summarizing can take a bit, can just fit and store fit)

## Flag to determine how plotting will proceed (after fitting, or from a saved model, or from extracted chains)
if (plot_model) {
plot_from <- {
  if (fit_model) {
   "fit"             ## plot from a model fit in this session
  } else {
   "saved_model"     ## plot from the raw fitted model
 # "saved_samples"   ## plot from summarized model output
  }
 }
}

## Print statements here and throughout to track progress if running from the command line or on an external server
if (fit_model) {
  print("A model will be fit")
} 
if (plot_model) {
  print("Plotting of model fit will occur")
}

## Name of the saved samples or saved model
 ## Adjust these if plotting from a model that was previously run
if (plot_model) {
if (plot_from == "saved_model" | plot_from == "saved_samples") {
  if (plot_from == "saved_model") {
      saved_model <- "fits/stan_fit_multipop_all_full_2022-06-23.Rds" 
   if (file.exists(saved_model)) {
     print(paste("Plotting will occur using the saved model:", saved_model, sep = " "))
   } else {
     print("Wrongly named or missing saved fitted model")
     break
   }
  } else if (plot_from == "saved_samples") {
      saved_samples <- "samples/stan_multipop_mehg_cleaned.Rds"
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
source("complete_data.R") 

## Some choices to determine what model will be fit. 
 ## Can't be perfectly dynamic because populations are manually chosen
  ## and some choices won't work with certain models, so will need to double check. Mismatches will lead to errors
   ## List of populations that can be fit with individual-level mercury listed in "determine_model.R"

 ## 0) Single population?
sing_pop       <- FALSE
 ## 1) Multiple species?
multi_spec     <- FALSE
 ## 1.2) If multiple species, fit a reduced model with no species-specific fixed effects?
multi_spec_red <- FALSE
 ## 2) Not all populations?
some_pops      <- TRUE
 ## 3.1) Fit individual-level MeHg? Be careful what populations to choose
fit_ind_mehg   <- FALSE
 ## 3.2) Fit a model only predicting MeHg? (no survival) Setting as TRUE invalidates many other options pertaining to survival
fit_only_mehg  <- FALSE
 ## 4) Reduced detection model? (if FALSE fits a random effect level for every day in every population)
red_p_model    <- TRUE

## More printing to better track specific fit when using slurm
print("Fitting choices are:")
print(paste("sing_pop =", sing_pop, sep = " "))
print(paste("multi_spec =", multi_spec, sep = " "))
print(paste("multi_spec_red =", multi_spec_red, sep = " "))
print(paste("some_pops =", some_pops, sep = " "))
print(paste("fit_ind_mehg =", fit_ind_mehg, sep = " "))
print(paste("fit_only_mehg = ", fit_only_mehg, sep = " "))
print(paste("red_p_model =", red_p_model, sep = " "))

## From these choices find the model to fit
source("determine_model.R")

## Pop*spec
 # 1  - AMCI.SMNWR_E
 # 2  - AMCI.SMNWR_W
 # 3  - ANBO.Blackrock  
 # 4  - ANBO.JonesPond
 # 5  - ANBO.SonomaMountain
 # 6  - ANBO.TwoMedicine
 # 7  - NOVI.KettleMoraine
 # 8  - NOVI.MudLake
 # 9  - NOVI.ScotiaBarrens
 # 10 - NOVI.SMNWR_W
 # 11 - NOVI.Springfield
 # 12 - PSMA.LilyPond
 # 13 - PSMA.MatthewsPond
 # 14 - RANA.DilmanMeadows
 # 15 - RANA.FoxCreek
 # 16 - RANA.JonesPond
 # 17 - RANA.LittleThreeCreeks
 # 18 - RANA.LostHorse
 # 19 - RANA.SanFrancisquito
 # 20 - RANA.SummitMeadow

## Some population choices:
 ## Full fit without the biggest populations (KettleMoraine, Springfield, and Scotia Barrens)
   ## -c(7, 9, 11) 
 ## MeHg fit (with survival)
   ## c(4, 5, 6, 14:19)
 ## MeHg --> Bd (no survival)
   ## Very best data
     ## c(4, 5, 14, 16, 17, 18) 
   ## All but the very worst populations and the very large populations
     ## c(3, 4, 5, 6, 14:19)

## If a subset of populations, pick which ones
if (some_pops) {
# which.dataset <- unique(data.all$pop_spec)[10] %>% droplevels() 
# which.dataset <- unique(data.all$pop_spec)[2] %>% droplevels() 
# which.dataset <- unique(data.all$pop_spec)[11] %>% droplevels()
 which.dataset <- unique(data.all$pop_spec)[grep("ANBO", unique(data.all$pop_spec))] %>% droplevels()
# which.dataset <- unique(data.all$pop_spec)[grep("RANA", unique(data.all$pop_spec))] %>% droplevels()
# which.dataset <- unique(data.all$pop_spec)[c(3, 13, 16)] %>% droplevels() 
# which.dataset <- unique(data.all$pop_spec)[c(4, 5, 12, 13, 14, 15, 16, 17, 18, 19)] %>% droplevels()
  data.all      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
  sampling      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
  
if (sing_pop) {
  which_stan_file <- (which_model_fit %>% filter(pop_spec == which.dataset %>% as.character()))$model
  this_model_fit  <- paste("stan_current/", which_stan_file, ".stan", sep = "")
}
  
} else {
  which.dataset <- unique(data.all$pop_spec) %>% droplevels()
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
if (!sing_pop) {
 source("establishing_mm.R")
}

## Print more in-depth details about the dataset being fit, if its a single population, enough
 ## already printed to follow what is being fit if multi-pop
if (sing_pop) {
 source("dataset_notes.R")
}

## Quick look at a given population
# source("capt_plot.R")
# source("capt_plot_multi.R")

## And finally run the stan model
stan.iter     <- 1500
stan.burn     <- 500
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin
stan.chains   <- 1
stan.cores    <- 1
stan.refresh  <- 10
if (fit_model) {
 if (!fit_only_mehg) {
  if (length(which.dataset) == 1) {
   source("stan_fit_single.R")
  } else {
   source("stan_fit_mm.R") 
  }
 } else {
  source("stan_fit_mehg_only.R")
 }
}

## And some diagnostics and such
if (plot_model) {
 if (length(which.dataset) == 1) {
  source("plotting.R")
 } else {
  source("plotting_multipop.R")
 }
}

