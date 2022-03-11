#####################################
## Fit CMR model to amphibian data ##
#####################################

####
## Notes as of March 11:
####

## --- There are still a few issues with the single population model to be resolved ---

 ## 1) phi_ones designating a closed population still in the model with length of time between time gaps -- confounded
        ## ------ [ ]
        ## There are some populations were the sub-sites are sampled at dates pretty far from one another
        ## And occasionally one subsite will be sampled multiple times in advance of a different subsite
        ## This makes it hard to think how to deal with date, closed vs open, and effort
        ## Blackrock-C is a good example of this, where it is pretty hard to designate Primary and Secondary periods
        ##   because it isn't clear at all if these should get designated at a subsite level or a site level...
        ## ** Just a note that the current SecNumConsec are pretty wack for this location

 ## 2) gamma still scaling p but in a funny way if we are assuming an open population. Need to do something
  ##   more sensible with gamma to try and scale p (as the individual could have been present before they were sampled
  ##   which could impact estimates of bd-detection and bd-surival. 
  ##   BUT: if Bd being removed from the detection model, maybe we don't need gamma at all -- this will simplify the model
        ## ------ [ ]

 ## 3) "effort" is behaving a bit funny. Need to seriously consider adjusting dates a day forward or backward so each 
  ## "date" is ONE sampling event of each sub-site
        ## ------ [ ]

 ## 4) model for MeHg -- potentially estimating site-level mean [and sd?] from the individual distribution
        ## ------ [ ]

## --- Other list of unresolved issues ---

 ## DATA:
  ## -- Some confusing notes in the FL data set remain 
   ##     [Will email soon]
  ## -- DAYMET data for 2021 to be available soon
   ##     [Brian T working on this]

 ## MODEL:
  ## -- Most importantly it is still a pretty big question of what to do with sub-populations
   ##   (opens a can of worms about meta-population, super-population, movement between subsites etc.)
   ##   The major problems are about bias in detection because of assuming individuals are potentially found when they really are not
   ##     (because they are in a different subpopulation)
   ##   SEE OLDER NOTES FROM A PREVIOUS COMMIT FOR A LONGER DISCUSSION
    ##    [Really not sure what to do about this yet....]
    ##    [Continue for now with collapsing and having an "effort" covariate of the number of subsites sampled per day]
     ##    [OR -- may want to collapse consecutive days into the same "period" if different subsites were sampled each day]
      ##    [This is annoying because it may be site-by-site dependent and hard to do dynamically, but may have to go this route]
  
  ## -- And a slightly less major question about how to deal with primary periods and secondary periods vs continuous times between captures
   ##     [Try a covariate for N days between sampling events]
   ##     [But may want to revert to figuring out how to write out primary and secondary periods for all populations]
  ## -- Still many open questions about what covariates to use and how to merge discrete covariates
    ##    [ ]
  ## -- Potential mercury latent model
   ##     [Use distribution of mercury for a site-level covariate]
    ##    [But first look better at the mercury data to make sure there is a reasonable amount of among-population variation]
  ## -- Multiple imputation for unobserved covariates
   ##     [Do multiple imputation for most species] 
   ##     [For the really unmeasured species use a strong prior based on that species]
  ## -- What to do with multiple swabs
   ##     [Collapse as I already have done]
  ## -- How to deal with length measured sometimes and size measured sometimes
   ##     [Just use length]

## Packages and Functions
source("packages_functions.R")

## Read in data
source("data_load.R")

## Construct modified data frame of recapture histories for each individual in each population
 ## For single species debug purposes pick a single data set
single_pop <- TRUE

if (single_pop) {
which.dataset <- unique(data.all$pop_spec)[12]
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
stan.iter     <- 600
stan.burn     <- 200
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin
if (single_pop) {
source("stan_fit_single.R")
} else {
source("stan_fit.R") 
}

## And some diagnostics and such
source("diagnostics.R")
