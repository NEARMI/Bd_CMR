###############################################################
## Simulate data for and fit an individual-CMR model in Stan ##
###############################################################

## Script whose predecessor is individual_CMR_expanding.R -- converting that script into functions for easier simulation of multiple
 ## populations for the more complicated CMR model

####
## Notes as of Dec 2:
####

## -- General progress update -- ##

 ## 1) Single population simulation working reasonably well for both full and collapsed simulation.
  ## A) As expected collapsed model doesn't capture within season process well at all given that the whole
   ## temporal dynamics of bd are ignored
  ## B) Streamlined quite a bit of code so do need to check a few population sizes etc. to make sure expected
   ## patterns are returned (i.e., better estimates with large populations and more sampling)

 ## 2) Full list of next steps below, but the first important thing to do is to get the multi-pop model working
  ## which has two steps really: 1) checking to make sure the code as is runs and gives sensible results; 2) adjusting
   ## the random effects to the between season process to make the collapsed model more useful

## -- Next to do -- ##

 ## 1) There is now some redundancy in .stan models that have _sim_ attached and those that don't 
  ## (specifically CMR_full_pr.stan and CMR_full_sim_pr.stan). The idea is that these will begin to diverge
   ## once more empirical data becomes available (like the non pr models has started to do)

 ## 2) The current formulation of random effects is rather non-sensible for capturing the between population
  ## processes that we really care about for the empirical data. Thus, need to adjust the random effects and 
   ## simulation so that between season survival (and more of a minor point--detection) are the processes
    ## that vary among populations (and not within season survival as it is now)

 ## 3) There is still a minor issue with survival from the last time time point through the
  ## offseason because survival from the last point to the end of the season is ignored
   ## (that is, for offseason == 1 timegaps isn't factored in to also scale survival -- which coould
    ## lead to a minorly biased estimate of between season survival (higher mortality than reality))

####
## Simulation Notes as of Dec 2, 2021
####

### -- Parameters and general simulation info -- ###

## 1) Temperature is simulated as being quadratic over time, then bd is simulated as being a linear
  ## function of temperature. The model estimates bd as a quadratic function over time
  ## with a linear term for temperature, which is a tad different than the simulation but nearly equivalent
## 2) Survival is estimated in two processes. Survival is first estimated over single time steps solely as
  ## a function of simulated bd. Survival isn't explicitly modeled as a function of time, it is just emergent
  ## that survival will differ depending on the time between samples. Survival between seasons is a separate 
  ## process where survival is estimated between seasons with a baseline and as a function of max bd that year
   ## Thus, things like primary and secondary periods (over which the population is assumed to be closed) can
   ## be assessed for things like bias (as it isn't realistic at all that the population truly is closed --
   ## though it is important to note the time frame between reality (e.g. secondary periods separated by days
   ## and primary periods separated by months isn't perfectly represented in this simulation model. Need
   ## some updates to the simulation to deal with this variable time frame))
## 3) Detection is directly estimated on each sampling occasion as a function of some baseline and bd load.
  ## This is something to consider expanding in complexity in the near future

### -- Models that can be fit -- ###

## Currently there are two supported models to be fit from this simulation:
 ## 1) A "full" model that estimates the quadratic bd load and week to week survival as a function of bd.
  ## This model requires relatively extensive sampling in order to fit
 ## 2) A "collapsed" model that does not estimate the time dependence of bd, but instead just estimates the max
  ## that each individual reaches (i.e., as if bd is an individual-specific trait). 
  ## In this model secondary periods are defined with a user-defined space between sampling events that are
  ## considered to be close enough to assume a closed population (potentially something interesting to explore)

####
## Packages and misc
####
needed_packages <- c("magrittr", "dplyr", "tidyr", "lme4", "ggplot2", "rstan")
lapply(needed_packages, require, character.only = TRUE)
source("../../ggplot_theme.R")
set.seed(10002)

####
## Parameters 
####
source("CMR_parameters.R")
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
collapse.mod <- FALSE

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
## NOTE: Diagnostics are not yet set up to debug the multi-population simulation well
 ## (or the collapsed model as well either, but that can be the next step as there is less to change)
source("CMR_diagnostics.R")


