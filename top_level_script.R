#####################################
## Fit CMR model to amphibian data ##
#####################################

####
## Notes as of March 17:
####

########## ---- 0) Some nice progress 

## A) Fitting and plotting scripts automated and working [?]


########## ---- 1) Things I have learned about the single population model

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

########## ---- 2) Thoughts on how to deal with a multi-population model

 ## A) My -hope- will be to have species as a fixed effect and population as a random effect
  ##    ^^ Though this may be impossible because of the overlap of these two factors

 ## B) Will need lots of nested random effects (different variance by population -- for Bd, p, and phi)

 ## C) My aim will be to figure out some way to impute MeHg, but this could get confusing if it is correlated with
  ##   both Bd and length. 
   ##   ^^ It seems _feasible_ to write down a joint latent / multiple imputation problem, but it is unclear if this will just become soup

 ## D) It is still an open question of what covaraites to use and in which processes.
  ##   As well as how to average categorical covariate 
   ##   ^^ But maybe this averaging isn't necessary if the individual capture on each day is traceable to the subpopulation and the covaraties there

########## ---- 3) The next steps

## 1) Make a final decision for a reasonable structure for the single population model
## 2) Fit to all individual populations
## 3) Make necessary adjustments to multi-population model, STARTING with just ANBO
## 4) Proceed from there

##################

## Packages and Functions
source("packages_functions.R")

## Read in data
source("data_load.R")

## Construct modified data frame of recapture histories for each individual in each population
 ## For single species debug purposes pick a single data set
single_pop <- TRUE

if (single_pop) {
which.dataset <- unique(data.all$pop_spec)[4]
data.all      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
sampling      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
}

red_ind    <- FALSE
if (red_ind) {
num_ind    <- 150
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
