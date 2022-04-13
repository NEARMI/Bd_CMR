########################################################################
## Fit only a single population at a time, utilizing the command line ##
## to pick which population (for submitting multiple jobs at once)    ##
########################################################################

source("packages_functions.R")

stan.iter     <- 3500
stan.burn     <- 500
stan.thin     <- 2
stan.length   <- (stan.iter - stan.burn) / stan.thin
red_ind       <- FALSE
single_pop    <- TRUE

which_model_fit <- read.csv("stan_current/which_model_fit.csv")

## Use parameter values from the command line
args           <- commandArgs(TRUE) 
which.dataset  <- args[1]

this_model_fit <- args[2]
this_model_fit %<>% paste("stan_current/", ., sep = "") %>% paste(., ".stan", sep = "")

source("data_load.R")
  
data.all      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
sampling      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()

source("data_manip.R")
source("data_stan.R")
source("data_covariates.R")
source("dataset_notes.R")
source("stan_fit_single.R")
source("plotting.R")
