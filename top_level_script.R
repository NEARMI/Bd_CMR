#####################################
## Fit CMR model to amphibian data ##
#####################################

#### Notes as of April 7 ----

## 1) Finished writing methods today. It actually was pretty hard to write the text given all of the difficulties with
 ## indices and such (especially what entries for each individual inform which processes etc and the fact that phi is one
  ## shorter than p even though both have their limitations (p cant be informed on the first occasion and phi on the last)
 ## -- **After doing some double checking, *I THINK* that all of my writing matches my indices and that all of my indices are still
  ## correct, though it may be worth a triple check.
 ## -- Also, may need to write a little bit about the length of phi being T - 1? Or maybe that can be explained away in the model and
  ## doesn't have to enter the main paper...

## 2) The next modeling steps include:
 ##  [Multi] [Single no MeHg] [Single MeHg]
 ## B) [x]     [X]                [X] Length to predicted Bd at the individual level
   ## ---- [ ] Most likely worth adding. Compiles, but need to test still
 ## C) [-]     [-]                [X] Interaction term for Bd and MeHg
   ## ---- [ ] Most likely only worth adding for the few populations where individual-MeHg values can be imputed. Compiles, but need to test still
 ## D) [X]     [X]                [X] Sex in survival
   ## ---- [ ] Probably worth adding, but unclear how much it gains us. Compiles but need to test still
 ## E) [ ]     [ ]                [ ] Effect of "injury" on survival
   ## ---- [ ] Still to do. Few actual individuals fall into this category so probably a small effect...
 ## F) [ ]     [X]                [ ] Sex in length imputation
   ## ---- [ ] Looks needed given U for sex often is a stand-in for juveniles which are smaller
    ##         Will work basically in all cases apart from a case where all of one sex wasn't measured for length, but I think this shouldn't exist...
 ## G) [ ]     [ ]                [ ] If sticking with this fixed effects specification, convert from using pop_spec for the random effect to site
   ## ---- [ ] Still to do, just code, no model changes I don't think

## 3) Data checks and potential larger modifications
 ## A) Florida sampling scheme being fundamentally different to the other populations
 ## B) What individuals were actually marked for NOVI for the Wisconsin populations

## 4) Longer term strategy:
 ## A) Big model
 ## B) Bd-MeHg interaction model for the best sampled populations
 ## C) Seasonal variation model for the newts

##### Older notes from the fits run the week of March 31 ----

## 1) The model run with all of the individuals from all populations apart from NOVI ran with no divergent transitions.
 ## A) though 400 warmup and 600 samples will be far too little for the full model. Probably something like 500 warmup and 2500
  ##   samples could be enough
 ## B) Even with all of the individuals, the fixed effect coefficient estimates are still pretty garbage
    ##   ^^ Most overlapping zero with wide CI
    ##   ^^ RANA having a super strong length effect ** Need to add sex?
 ## C) The model has changed a lot so need to work back through all of the priors to make sure they are still sensible 

##### Some older notes that are still relevant ----

#### ---- Multi-pop model
 ## A) Because of difficulty with the overlap of species and location I am collapsing
  ##   all of the Rana species, as this will allow for a fixed effect of species and a random effect of location 
 ## B) One species doesn't have any MeHg and the location isn't shared. Will have to try and use a fixed effect + pop_spec random effect?
  ##   The hope here would be to use an intercept and a species and random location so at least the intercept can inform

#### ---- Some things I learned from playing with the multi-pop model:
  ## --- A) The small populations can teach us basically nothing. Because of this using single intercepts for "species" is likely
  ##        going to be pretty dangerous. The strategy is likely going to have to be to have species as a random effect -- but then we
  ##        are back to the original problem of how to specify these random effects (as location and species are super correlated)
  ## --- B) ^^ Continuing this thought, having something like one effect of "size" gets washed out across populations when so many populations
  ##        cant help resolve this relationship. Will want to have population-unique deviates for basically all covariates

#### ---- Single-pop model

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

##### Code -----

#### NOTE: In this file and all other files search *** for current choices that could potentially change

## Packages and Functions
source("packages_functions.R")
source("../ggplot_theme.R")

## Read in data
source("data_load.R")

## For dev and debug purposes pick a subset of locations
some_pops  <- TRUE

if (some_pops) {
# which.dataset <- unique(data.all$pop_spec)[c(3:7, 17, 18)] %>% droplevels()
which.dataset <- unique(data.all$pop_spec)[-c(10:14)] %>% droplevels()
# which.dataset <- unique(data.all$pop_spec)[20] %>% droplevels()
data.all      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
sampling      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
}

## For dev and debug purposes also can subset total number of individuals 
 ## (done randomly though a seed is set in packages_functions.R)
red_ind    <- FALSE
if (red_ind) {
num_ind    <- 200
}

## Create the capture history scaffold from the raw data
source("data_manip.R")

## Create the indexing vectors and capture history structure needed for the stan model
source("data_stan.R")

## Deal with all of the individual level and population level covariates 
source("data_covariates.R")

## Quick look at a given population
#source("capt_plot.R")
#source("capt_plot_multi.R")

## And finally run the stan model
stan.iter     <- 1000
stan.burn     <- 400
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin
if (length(which.dataset) == 1) {
source("stan_fit_single.R")
} else {
source("stan_fit.R") 
}

## And some diagnostics and such
#source("diagnostics.R")
