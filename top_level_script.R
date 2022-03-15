#####################################
## Fit CMR model to amphibian data ##
#####################################

####
## Notes as of March 15:
####

## -- A note about model structure -- ##

## There are enough indexing columns created at this point that modifications to the largest model
 ## are no longer difficult to make. For this reason I am moving most of the stan models to the
  ## dev folder from here forward

## -- Still to be cleaned up -- ## 

## 0) Looking back at the model I am now worried a bit about a few things:
 ## A) Animals that were only ever captured on the last day can't contribute to the likelihood. 
  ##    ^^ What is currently going on with these animals in my model???

## -- Some progress -- ## 

## Made some nice progress on things by adding day-level and individual-level random effects
 ## for detection and survival 

## Stepping back to periods where the population is assumed to be closed and other periods where the population
 ## is assumed to be open is providing some reasonable fits and seems like a way forward
  ## ^^ It seems like estimating survival in very narrow time windows just isn't feasible

## Saved a nice plot showing little effect of different random effects at the day level and individual
 ## level on survival estimates, but which seem to help clean up detection
  ## ^^ Still need to try this on Blackrock

## -- Still open questions -- ##

## Honestly not too much remaining for individual populations before trying to fit every population
 ## individually, BUT, still need to figure out:
   ## 1) a threshold for all populations to assume open
   ## 2) what covariates to include for within vs between season survival
   ## 3) if to include any covariates in survival

## -- The next steps -- ##

## 1) Make a final decision for a reasonable structure
## 2) Fit to all individual populations
## 3) Make necessary adjustments to multi-population model
  ##   STARTING with just ANBO
## 4) Proceed from there


####
## Notes as of March 11:
####

## --- There are still a few issues with the single population model to be resolved ---

 ## 1) What to do with Primary and Secondary periods / allowing for open populations when each population is sampled so differently?

  ## CONJECTURE: This problem mostly arises because of the need to collapse subsites, which is needed because:
   ## -- some populations don't have subsites while some have many
   ## -- some are very poorly measured so estimating subsite level deviates will be hard
   ## --  most importantly in many cases animals move between subsites, so to not collapse subsites would likely require extra model layers of connectivity and
   ##     the like, which is too much for these data (specifically the rarity of recaptures and bd swabbing)
     ## -- There are some populations were some sub-sites are sampled at dates pretty far from one another but other sub-sites are sampled more commonly
     ## -- Some populations are measured for many months, others only for a single month, which makes defining "closed" confusing [e.g., by year seems odd]
   ## Blackrock-C is a good example of these issues ^^

  ## All taken together means that using a "closed" population over some primary period in which the population is sampled in secondary periods difficult to do well
  
  ## CONJECTURE: The alternative is to just allow for the population to be open (and collapse to site)
    ## -- Which requires some choice for how to control for the gaps between sampling events
     ## [a simple covariate controlling survival by time gap has produced counter-intuitive results]
    ## -- Which also requires a choice for how to scale detection to deal with variable sampling effort by day
     ## [a simple covaraite controlling detection by "effort" defined by the number of subsites sampled has produced counter-intuitive results]
    
     ## !!! a second option is to completely rework Secondary Periods into blocks of time explicitly with Site in mind doing as best as possible to
      ## group sampling of multiple SubSites on consecutive days into defined time blocks. 
       ## -- While this is sensible it becomes difficult for a lot of the sampling, e.g., Blackrock-C where individual sites were sampled much more
        ## frequently than others 
         
    ## CONJECTURE: This style of sampling is a real problem for this model (when there are sites with pretty different forms of sampling) for a few reasons:
     ## A) This will bias detection on certain animals that reside more often in that subsite vs other subsites
     ## B) Which could affect the survival estimates because detection and survival are already hard to pull apart
      ## ** A wish would be to have individual-level detection parameters, but that isn't feasible because so many individuals are only captured once,
       ## which will make it really hard to separate that fine-grained level of detection from covaraite effects on survival

 ## 2) What to do with detection in periods before the animal was caught?
   ## A) The simple answer is to just not inform detection in these periods.
   ## B) However, if the animal really was in the population (and we think Bd effects detection for example), estimates of Bd survival could get biased down
    ## by an overly high detection probability 
     ## -- I don't think we really can estimate if we think an individual was present before we captured it though without a extra and pretty intricate 
      ## layres of latent processes, so my inclination is to just avoid trying to deal with this

 ## 3) What to do with MeHg?
   ## A) Try to impute an individual's latent MeHg and use that to predict survival
   ## B) Estimate site-level mean [and sd?] from the individual distribution and use MeHg as a site-level covariate

## --- Other unresolved issues ---

 ## DATA:
  ## i) DATA: Some confusing notes in the FL data set remain 
   ##     [Will email soon]
  ## ii) DATA: DAYMET data for 2021 to be available soon
   ##     [Brian T working on this]
  ## iii) MODEL: Still many open questions about what covariates to use and how to merge discrete covariates
   ##     TBD

## Packages and Functions
source("packages_functions.R")

## Read in data
source("data_load.R")

## Construct modified data frame of recapture histories for each individual in each population
 ## For single species debug purposes pick a single data set
single_pop <- TRUE

if (single_pop) {
which.dataset <- unique(data.all$pop_spec)[11]
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
