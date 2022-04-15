#####################################
## Fit CMR model to amphibian data ##
#####################################

#### Notes April 15 ---- 

## IMPORTANT: Before continuing too far with modeling, need to:
 ## A) make a graph of different ways of looking at days, cumulative temperatures, and degree days
  ## in thermal optima against Bd. 
   ## -- Upload these to Overleaf to show the huge amount of noise present across space, time, and
    ## among individuals that seems unexplained by covariates.
 ## B) Clean up the Overleaf text to refer to these figures and essentially say this noise and sparse
  ## sampling over days makes it difficult to do better than average Bd in a year...
 ## C) YET, we really do need to control for when a swab was taken to be able to estimate an individuals
  ## load experienced in that year, so we have to do SOMETHING to control for time in the season

## 1) Unfortunately fits timed out after three days... With only 200 individuals per population.
 ## This is not a good sign. Resent with lower adapt delta and lower tree depth. Hopefully will fit this time
  ## Very likely I will need to reduce the size of the model, potential candidates include:
   ## A) some of the categorical variables in between season survival

## 2) With finer-grained temperature data it is potentially conceivable to be able to "control"
 ## for when each individual was swabbed and trapped to get a better estimate of individual
  ## bd loads -- this could _maybe_ lead to continuous variation in Bd loads over time
## A) There are a lot of ifs here, but it could possibly work out.
## B) Today I pulled together the temperature data and put together the desired Bd model. 
 ## [ ] The first step is to run this on the PA population and go from there -> but the real
  ## strength is to use the shared temperature data across populations to help inform bd (i.e, if there
   ## were far more favorable days in PA than WY, this scaled temperature covariate could help to inform
    ## the Bd levels seen in WY whenever that population was sampled)
 ## C) And finally, this finer-grained temp could maybe allow for a temp*bd interaction in survival
  ## as was found in past research, but I am still somewhat unclear on that

## ^^ Continuing the above thought, an individual's Bd on Day X in Pop Y now becomes a scaled
 ## value of degree days until Day X, as well as species and population deviates. 
  ## ^^ From there, survival becomes either a function of cumulative load throughout the season
   ## or just max Bd in the season (wouldn't have to scale if max, would have to scale if cumulative)



#### Next stuff to do on this project ----

# x = partial; X = done for now

## 0) [x] [Tuesday] code, model, and repo cleaning
  ##        -- got a reasonable start, some more will be needed when debugging continues
## 1) [x] [Tuesday] Figure out what to do with FL and update the spreadsheet
  ##        -- some solid progress here but need to check in with Evan and Dave to make sure
## 1) [ ] [Later] Work through debugging of the various model attempts
  ##        -- would have liked to on Wednesday, but the models are still running, and it looks like they may time out
## 2) [ ] [After debugging (1 above)] Make some final decisions for model runs
## 3) [ ] [After final decisions on model structure (2 above)] Send individual (with .sh file) and multi-pop jobs
## 4) [ ] [Eventually] Debug the long fits
## 5) [ ] [Eventually] Start expanding the current non-MeHg individual model for Bd trajectories over time for PA newts
## 6) [ ] [Eventually] Start cleaning Overleaf and creating some fitting figures in prep for sending to PIs


#### Various other ToDos for later ----

## 1) Model modifications / adjustments
 ## A) [ ] Possible effect of "injury" on survival not yet added
 ## B) [ ] Possible use of sex in detection
 ## C) [ ] Finer grained temperature data

## 2) Potential larger modifications
 ## A) [x] Florida sampling scheme being fundamentally different to the other populations
   ## -- Think this gets fixed with defining a reasonable gap between samples
 ## B) [ ] What individuals were actually marked for NOVI for the Wisconsin populations
   ## -- Seems like all? But still need to email PI


#### Code ----

#### NOTE: In this file and all other files search *** for current choices that could potentially change

## Packages and Functions
source("packages_functions.R")
source("../ggplot_theme.R")

## Read in data
source("data_load.R")

## For dev and debug purposes pick a subset of locations
some_pops  <- FALSE

if (some_pops) {
#which.dataset <- unique(data.all$pop_spec)[-c(10:14)] %>% droplevels()
which.dataset <- unique(data.all$pop_spec)[17] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(1, 2, 13)] %>% droplevels()
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
