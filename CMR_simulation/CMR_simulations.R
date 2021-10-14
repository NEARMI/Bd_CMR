###############################################################
## Simulate data for and fit an individual-CMR model in Stan ##
###############################################################

## Script whose predecessor is individual_CMR_expanding.R -- converting that script into functions for easier simulation of multiple
 ## populations for the more complicated CMR model

####
## Packages and functions
####
needed_packages <- c("magrittr", "dplyr", "tidyr", "lme4", "ggplot2", "rstan")
lapply(needed_packages, require, character.only = TRUE)
source("../../ggplot_theme.R")
set.seed(10002)
'%notin%' <- Negate('%in%')

####
## "Design" parameters
####
nsim      <- 1                    ## number of simulations (1 to check model, could be > 1 for some sort of power analysis or something)
ind       <- 60                   ## number of individuals in the population being modeled
periods   <- 3                    ## number of primary periods (years in most cases)
new_ind   <- rep(10, periods - 1) ## individuals added in each new period
inbetween <- seq(1.5, periods, by = 1)
all_ind   <- ind + sum(new_ind)   ## number of individuals ever to exist in the population
times     <- 20                   ## number of time periods (in the real data probably weeks; e.g., May-Sep or so)
samp      <- 10                   ## number of sampling events occurring over 'times' (e.g., subset of 'times' weeks when sampling occurred)
if (periods > 1) {
samp <- rep(samp, periods)      ## for now assume same number of periods per year, but this model allows variable sampling dates by season
between_season_duration <- 10   ## number of time periods that elapse between the on-season
}
when_samp <- "random"           ## random = sampling occurs on a random subset of possible days

####
## bd_parameters
####
bd_beta <- c(
    1             ## Intercept
  , 0.1           ## Time effect
  , 0.3           ## Linear effect of temp on bd
)
bd_sigma  <- 2  ## observation noise
bd_theta  <- 5  ## random effect variance covariance
bd_mort   <- c(decay = -0.2, offset = 6)         ## logistic response coefficients for mortality across log(bd_load)
bd_detect <- c(decay = 0.1, offset = -0.5)       ## logistic response coefficients for detection across log(bd_load)

## bd sampling sampling scheme
bd_swabs  <- "PAT"  ## ALL = assume all captured individuals have their Bd swabbed on every capture
                    ## IND = assume specific captured individuals always have their Bd swabbed while others are never swabbed
                    ## PAT = assume patchy bd swabbing among all captured individuals 

       if (bd_swabs == "IND") {
bd_drop <- 20       ## number of individuals we will assume didn't have their Bd measured
} else if (bd_swabs == "PAT") {
bd_perc <- .50      ## proportion of all captures with bd swabs taken
bd_drop <- ind * samp * (1 - bd_perc)
}

####
## other parameters
####

## mortality probability in-between periods
p_mort_type <- "max" ## con = single value; max = based on max bd; cum = based on cumulative bd
if (p_mort_type == "con") {
p_mort    <- c(0.30) 
} else if (p_mort_type == "max") {
p_mort    <- 0.02
} else if (p_mort_type == "cum") {
p_mort    <- 0.0001  
}

## observation noise in bd (** SEPT 21: yes, on the log scale so this is weird to have the same var across log bd load, will clean this up later)
obs_noise <- 0.25   

## Stan model parameters
stan.iter     <- 1500
stan.burn     <- 500
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin


## Simulate bd for a population
bd.simulate <- function (
  periods, times, all_ind
, bd_beta, bd_sigma, bd_theta
) {

## Simulate data using lme4 mixed model structure
expdat <- expand.grid(
  periods = seq(periods)
, times   = seq(times)     
, ind     = factor(seq(all_ind))
  )

## ** Make dynamic later
expdat %<>% mutate(
  temp = rlnorm(n()
  ,   (scale(times, center = T, scale = F)[, 1] * .05) - 
    .05 * scale(times, center = T, scale = F)[, 1]^2 + 3.5
  , .2)
)

expdat %<>% mutate(
  bd_load   = simulate(~times + temp + (1 | ind)
  , nsim    = nsim
  , family  = Gamma(link = "log")
  , newdata = expdat
  , newparams = list(
      beta  = bd_beta
    , sigma = bd_sigma
    , theta = bd_theta
    )
  )$sim_1
) %>% mutate(
  bd_load     = round(bd_load, digits = 0)
, log_bd_load = round(log(bd_load), digits = 0)
) %>% mutate(
  log_bd_load = ifelse(is.infinite(log_bd_load), 0, log_bd_load)
)

return(expdat)
  
}



## List of parameters

