#####################################
## Fit CMR model to amphibian data ##
#####################################

#### Notes April 19 ---- 

## 1) A little work with the PA population sugests[?] that it may be sensible to try and pursue a hurdle model
 ##     i.e., is an animal infected, if so, what is its bd load?
  
 ## A) There is a higher chance that this is needed when modeling Bd over time given that a Normal distribution will split the middle of the 0s and non-0s
  ##   and have an inflated error variance

 ## B) I am a little less clear if this will be necessary for the model that estimates only a mean per year, as the individual deviate
  ## seems to be doing a pretty good job, but it probably would be better...

 ## C) This does raise a larger point though: that collapsing the Newt populations to fit in with the other populations and
  ## estimating only a mean will _likely_ lead to some problems / potential bias given the large time span over which individuals were sampled
   ## And trying to estimate this time period by population is not possible, and using the same time period for all populations is also
    ## not a good idea. Makes it hard to think how to include all populations at once...

 ## D) Finally, with already a very slow running model, the hope is to simplify the model anyway
  ##   -- given that the model is taking like 4 days with only 200 individuals per population...

## 2) It is possible that I have some random effect specification problems in a few places:
 ##    -- Definitely need an intercept for the non-centered parameterization
 ##    -- [Simulations show that no intercept definitely behaves the worst]
  ##     -- Which potentially raises some confusion when multiple random intercepts are used
   ##      -- Which means I may be specifying the random effects a bit wrong in the presence of multiple categorical
   ##         covariates, which _may_ explain why the multi pop model is so slow at present
    ##          -- It is possible that I may need to rewrite the random effect specifications to multivariate normal?

## 3) A few needed updates to the multi pop model
 ## A) Can we do better than yearly average temp? How about average temperature from the start of the year through
  ##   to when sampling started?
 ## B) I still want to get Bd and maybe sex into detection

## 4) Next step (today hopefully? will be to insure that the model with 200 individuals "worked") and then figure out
 ## how to simplify the model to make it tractable

## 5) Initial PA model looking reasonable, but unclear until the full model is fit, which will take a very long time with
 ## the current structure. Need to come back to this after deciding on some methods to speed up the model generally

#### Next stuff to do on this project ----

# x = partial; X = done for now

## 2) [ ]  Make some final decisions for model runs
## 3) [ ]  Send individual (with .sh file) and multi-pop jobs
## 4) [ ]  Debug the long fits
## 5) [ ]  Start expanding the current non-MeHg individual model for Bd trajectories over time for PA newts
## 6) [ ]  Start cleaning Overleaf and creating some fitting figures in prep for sending to PIs


#### Various other ToDos for later ----

## 1) Model modifications / adjustments
 ## A) [ ] Possible effect of "injury" on survival not yet added
 ## B) [ ] Possible use of sex in detection
 ## C) [ ] Finer grained temperature data??

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
#which.dataset <- unique(data.all$pop_spec)[17] %>% droplevels()
#which.dataset <- unique(data.all$pop_spec)[c(1, 2, 13)] %>% droplevels()
which.dataset <- unique(data.all$pop_spec)[4] %>% droplevels()
data.all      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
sampling      %<>% filter(pop_spec %in% which.dataset) %>% droplevels()
}

## For dev and debug purposes also can subset total number of individuals 
 ## (done randomly though a seed is set in packages_functions.R)
red_ind          <- FALSE
# red_ind_PA_debug <- FALSE ## Temp debug switch, will integrate if model works
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
