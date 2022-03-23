#####################################
## Fit CMR model to amphibian data ##
#####################################

####
## Notes as of March 24:
####

## Ended yesterday in the midst of working on breaking up individual length imputation by species

## TO DO today:

 ## 1) Debug the ANBO 5 pop fit from yesterday
 ## 2) Check why length imputation isn't working
 ## 3) Covaraite TODOs (see book) for detection and MeHg

## TO DO later:
 
 ## 1) Reach out to various PIs about fits
 ## 2) Reach out to Brian T for what he is doing with MeHg and survival

####
## Notes as of March 23:
####

## 1) Returned to a summary of all populations today to try and figure out what the full model will look like. A few things stand out that I will
 ## need to deal with:
  ## --- A) Huge overlap between species and populations, and many species are rarely sampled, making species as a fixed effect probably a no-go[?]
    ##        ^^ Evan's suggestion is to collapse all of the Rana species, as this will allow for a fixed effect of species and a random effect of location 
  ## --- B) One species doesn't have any MeHg and the location isn't shared. Will have to try and use a fixed effect + pop_spec random effect?
    ##        ^^ The hope here would be to use an intercept and a species and random location so at least the intercept can inform

## Some older notes that are still relevant

## 1) Multi-pop ANBO is a reasonable success. Lots of debugging show that the estimates match the single population
 ##   fits very well. Some things I learned:
  ## --- A) The small populations can teach us basically nothing. Because of this using single intercepts for "species" is likely
  ##        going to be pretty dangerous. The strategy is likely going to have to be to have species as a random effect -- but then we
  ##        are back to the original problem of how to specify these random effects (as location and species are super correlated)
  ## --- B) ^^ Continuing this thought, having something like one effect of "size" gets washed out across populations when so many populations
  ##        cant help resolve this relationship. Will want to have population-unique deviates for basically all covariates

## 2) Next steps are to:
  ## --- A) [X DONE, debug fit in progress] Add random effect of size impact on survival
  ## --- B) Get some other site-level covariates into the model. But where in the model and which covariates?
  ## --- C) [X DONE, debug fit in progress] Figure out what to do with MeHg
  ## --- D) Try and fit with two more populations of RALU and see if a fixed effect of species shows any hope of working, or what
  ##        combination of fixed and random effects will work, which will first require:
       ## -- i) Length imputation broken up by species

####
## Older notes about what I learned about the single population model that I am keeping around for now:
####

 ## A) Estimating the effect of Bd on survival is hard... Most estimates overlap 0

 ## B) Collapsing subsites to site and using a day-level detection random effect and periods of population closure
  ##   seems to be the most viable strategy moving forward (See C and D below)

 ## C) It seems like estimating survival in very narrow time windows just isn't feasible
  ##   ^^ Which means using periods of open vs closed populations is likely the way forward
  ##   ^^ But what then is the threshold with so many different sampling schemes?

 ## D) Adding individual and day random effects for detection costs little
  ##    (though individual random effect for within-season survival is unidentifiable)
  ##   Mixing still fast, estimates very similar to them not included except for detection which improves

 ## E) With three processes (between season survival, within season survival, and detection) it is still an open
  ##   question of what covariates to include in which processes

 ## F) I have indexing columns that modification of the model structure is now easy. There is probably less need
  ##   to have so many .stan files in stan_current

##################

## Packages and Functions
source("packages_functions.R")

## Read in data
source("data_load.R")

## Construct modified data frame of recapture histories for each individual in each population
 ## For single species debug purposes pick a single data set
single_pop <- FALSE
if (!single_pop) {
some_pops  <- TRUE
} else {
some_pops  <- FALSE 
}

if (some_pops) {
which_spec    <- c("ANBO", "RALU")
data.all      %<>% filter(Species %in% which_spec) %>% droplevels()
sampling      %<>% filter(Species %in% which_spec) %>% droplevels()
}

if (single_pop) {
which.dataset <- unique(data.all$pop_spec)[12]
data.all      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
sampling      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
}

red_ind    <- TRUE
if (red_ind) {
num_ind    <- 100
}

source("data_manip.R")

## Create the indexing vectors and capture history structure needed for the stan model
source("data_stan.R")

## Deal with all of the individual level and population level covariates 
source("data_covariates.R")

## Quick look at a given population
source("capt_plot.R")
#source("capt_plot_multi.R")

## And finally run the stan model
stan.iter     <- 800
stan.burn     <- 300
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin
if (single_pop) {
source("stan_fit_single.R")
} else {
source("stan_fit.R") 
}

## And some diagnostics and such
source("diagnostics.R")
