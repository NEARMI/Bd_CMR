####
## Parameters for the CMR simulation
####

### "Design" parameters, one entry for each population. 

## Written in a way for ease of entry. A bit lengthy looking but should be easier for many pops

n_pop      <- 3

## Note: all in matrix form based on n_pop for easy population of parameter lists

## number of individuals in the population being modeled
ind       <- matrix(data = rep(30, n_pop)
  , ncol = 1, nrow = n_pop)        

## number of primary periods (years in most cases). AS OF OCT 14 must be the same for all populations. To be updated later
periods   <- matrix(data = rep(3, n_pop)
  , ncol = 1, nrow = n_pop)

## individuals added in each new period.
new_ind   <- matrix(data = rep(rep(5, n_pop), periods - 1)
  , ncol = periods - 1, nrow = n_pop, byrow = T)

## number of individuals ever to exist in each population
all_ind   <- matrix(data = ind + rowSums(new_ind)
  , ncol = 1, nrow = n_pop)

## number of time periods (in the real data probably will use weeks; e.g., May-Sep or so)
 ## For now (and probably forever[?] given the structure of the bd sumbodel[?] equal by population)
times     <- matrix(data = rep(20, n_pop)
  , ncol = 1, nrow = n_pop)

## number of sampling events occurring over 'times'
 ## for now assume same number of periods per year, but this model allows variable sampling dates by season
samp      <- matrix(data = rep(10, n_pop)
  , ncol = 1, nrow = n_pop)
samp      <- mapply(rep, samp, periods) %>% t()

## number of time periods that elapse between the on-season
between_season_duration <- matrix(data = rep(20, n_pop)
  , ncol = 1, nrow = n_pop)   

## random = sampling occurs on a random subset of possible days
when_samp <- matrix(data = rep("random", n_pop)
  , ncol = 1, nrow = n_pop)    

## vector to designate the offseasons
inbetween <- apply(periods, 1, FUN = function(x) seq(1.5, x, by = 1)) %>% t()

####
## bd_parameters
####

## intercept and slope coefficients
bd_beta <- matrix(
  data = c(
    sample(seq(-3, 3, length = n_pop), n_pop)       ## Intercept
  , rep(-0.2, n_pop) ## Time effect
  , rep(0.3, n_pop)  ## Linear effect of temp on bd
), nrow = n_pop, ncol = 3, byrow = F)

## error
bd_sigma  <- matrix(data = rep(20, n_pop)
  , ncol = 1, nrow = n_pop)  

## random effect variance covariance
bd_theta  <- matrix(
  data = c(
    rep(1, n_pop)    ## Intercept
  , rep(0.3, n_pop) ## Time effect
  , rep(0.2, n_pop)  ## Linear effect of temp on bd
), nrow = n_pop, ncol = 3, byrow = F)

## logistic response coefficients for mortality across log(bd_load)
bd_mort <- matrix(
  data = c(
    sample(seq(-0.5, -0.05, length = n_pop), n_pop)      ## logistic slope
  , rep(6, n_pop)        ## intercept
), nrow = n_pop, ncol = 2, byrow = F)

## logistic response coefficients for detection across log(bd_load)
bd_detect <- matrix(
  data = c(
    rep(0.1, n_pop)   ## logistic slope
  , rep(-0.5, n_pop)  ## intercept
), nrow = n_pop, ncol = 2, byrow = F)

bd_noinf <- matrix(
  data = c(
    rep(0.1, n_pop)
  ), nrow = n_pop, ncol = 1, byrow = F)

## OCT 14 NOTE: bd sampling sampling scheme (see Individual_CMR_expanding.R for various sampling schemes,
 ## for now just using PAT here (assume patchy bd swabbing among all captured individuals))

## proportion of all captures with bd swabs taken
bd_perc <- matrix(data = rep(.50, n_pop)
  , nrow = n_pop, ncol = 1)

bd_drop <- sweep(samp, 1, ind, "*") %>% sweep(., 1, 1 - bd_perc, "*")

## Observation noise in bd
obs_noise <- matrix(data = rep(3, n_pop)
  , nrow = n_pop, ncol = 1)  

####
## other parameters needed for the simulation function
####

## mortality probability in-between periods
 ## OCT 14 NOTE: for other options see Individual_CMR_expanding.R, for now just proceeding with "max"
  ## will add back other options later
p_mort    <- matrix(data = rep(.02, n_pop)
  , nrow = n_pop, ncol = 1)

####
## other parameters for outside the function
####

## Stan model parameters
stan.iter     <- 1500
stan.burn     <- 500
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin
