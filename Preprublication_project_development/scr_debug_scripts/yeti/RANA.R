#####################################
## Fit CMR model to amphibian data ##
#####################################

## All RANA populations

fit_model      <- TRUE
plot_model     <- FALSE
source("packages_functions.R")
source("complete_data.R") 
sing_pop       <- FALSE
multi_spec     <- FALSE
multi_spec_red <- FALSE
some_pops      <- TRUE
fit_ind_mehg   <- TRUE
fit_only_mehg  <- FALSE
red_p_model    <- TRUE
red_ind        <- FALSE

print("Fitting choices are:")
print(paste("sing_pop =", sing_pop, sep = " "))
print(paste("multi_spec =", sing_pop, sep = " "))
print(paste("multi_spec_red =", sing_pop, sep = " "))
print(paste("some_pops =", sing_pop, sep = " "))
print(paste("fit_ind_mehg =", sing_pop, sep = " "))
print(paste("fit_only_mehg = ", fit_only_mehg, sep = " "))
print(paste("red_p_model =", sing_pop, sep = " "))

source("determine_model.R")

model_name      <- paste("fits/", paste(
  paste(which_stan_file, "RANA", sep = "_")
  , "_", sep = ""), Sys.Date(), ".Rds", sep = "")
which_stan_file <- paste("stan_current/", which_stan_file, ".stan", sep = "")

which.dataset <- unique(data.all$pop_spec)[grep("RANA", unique(data.all$pop_spec))] %>% droplevels()
data.all      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
sampling      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()

print(paste("From these choices, the following model will be [has been] fit:  ", which_stan_file, sep = ""))
print(paste("Using the following populations:"
  , as.character(which.dataset) %>% paste(collapse = " -- ")
  , sep = " "))
print(paste("With red_ind =", red_ind, sep = " "))

source("data_manip.R")
source("data_stan.R")
source("data_covariates.R")
source("stan_indices.R")
source("establishing_mm.R")

stan.iter     <- 1500
stan.burn     <- 500
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin
stan.chains   <- 3
stan.cores    <- 3
stan.refresh  <- 10

source("stan_fit_mm.R") 

