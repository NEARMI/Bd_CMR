###############################################################
## Simulate data for and fit an individual-CMR model in Stan ##
###############################################################

## Script whose predecessor is individual_CMR_expanding.R -- converting that script into functions for easier simulation of multiple
 ## populations for the more complicated CMR model

####
## Notes as of Dec 14:
####

## For extensive notes on the simulation model itself see the file without "_yeti" attached

## This script loops over parameter values and stores output for a power analysis of sorts

####
## Packages and misc
####
needed_packages <- c("magrittr", "dplyr", "tidyr", "lme4", "ggplot2", "rstan")
lapply(needed_packages, require, character.only = TRUE)
source("../../ggplot_theme.R")
set.seed(10006)

yeti_run <- TRUE  ## loop over many parameter combinations and fit models on cluster

yeti_vals <- expand.grid(
  pop_size = 250
, times    = c(2, 5, 8)
, bd_eff   = c(-0.05, -0.15, -0.25, -0.35)
, ind_var  = c(0, 1, 2)
)

for (yeti_set in 1:nrow(yeti_vals)) {

yeti_val <- yeti_vals[yeti_set, ]
  
####
## Parameters 
####
source("CMR_parameters_yeti.R")  

## alternative model parameterization that better matches to classic CMR models
use_prim_sec    <- TRUE

####
## Functions for simulation
####
source("CMR_functions.R")

####
## Run the sim to create the data
####
source("CMR_datasim.R")

####
## Clean up simulated data and build structure for stan model
####

## Simplified model that only fits 
collapse.mod <- TRUE 

if (use_prim_sec) {
  source("CMR_dataclean_collapsed.R")
} else {
  source("CMR_dataclean.R")
}

if (length(unique(expdat.all$ind)) < 300) {
expdat.all %>% {
  ggplot(., aes(times, ind, fill = as.factor(detected))) + 
    geom_tile(aes(alpha = sampling_days)) +
    geom_point(data = 
        expdat.all %>% 
        filter(bd_swabbed == 1)
      , aes(x = times, y = ind, z = NULL), lwd = 0.7) +
    geom_line(data = 
        expdat.all %>% 
        filter(dead == 1)
      , aes(
         group = ind
        , x = times, y = ind, z = NULL), lwd = 0.5, alpha = 0.5) +
    facet_grid(pop~periods) +
    xlab("Week of the year") +
    ylab("Individual") +
    scale_fill_manual(
        values = c("dodgerblue4", "firebrick4")
      , name   = "Detected?"
      , labels = c("No", "Yes")) +
    scale_x_continuous(breaks = seq(1, max(expdat.all$times), by = 2)) +
    guides(alpha = FALSE) +
    theme(
      axis.text.y = element_text(size = 6)
    , legend.text = element_text(size = 12)
    , legend.key.size = unit(.55, "cm")
    ) 
  }
}

####
## Run the stan model
####

if (use_prim_sec) {
  source("CMR_fit_collapsed.R")
} else {
  source("CMR_fit.R")
}

####
## Model diagnostics (best to open the script and run line by line because of all of the plots)
####
source("CMR_diagnostics_yeti.R")

}
