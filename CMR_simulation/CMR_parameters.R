####
## Parameters for the CMR simulation
####

### "Design" parameters, one entry for each population. 

## Written in a way for ease of entry. A bit lengthy looking but should be easier for many pops

n_pop      <- 1

## Note: all in matrix form based on n_pop for easy population of parameter lists

## number of individuals in the population being modeled
ind       <- matrix(data = rpois(n_pop, 50) # rep(50, n_pop)
  , ncol = 1, nrow = n_pop)        

## number of primary periods (years in most cases).
periods   <- matrix(data = rep(3, n_pop)
  , ncol = 1, nrow = n_pop)
## Play with different periods by year
if (n_pop > 1) {
  periods[n_pop] <- 2
}

## individuals added in each new period. A little weird to set it up this way, can
 ## maybe try a different method later
new_ind         <- vector("list", n_pop)
new_ind_per_pop <- rnbinom(n_pop, 10, .6)
new_ind_per_pop <- ifelse(new_ind_per_pop == 0, 1, new_ind_per_pop)
for (i in 1:n_pop) {
 new_ind[[i]] <- rep(new_ind_per_pop[i], periods[i, 1] - 1)
}

## number of individuals ever to exist in each population
all_ind   <- matrix(data = ind + lapply(new_ind, sum) %>% unlist()
  , ncol = 1, nrow = n_pop)

## number of time periods (in the real data probably will use weeks; e.g., May-Sep or so)
times     <- matrix(data = rpois(n_pop, 15) #rep(20, n_pop)
  , ncol = 1, nrow = n_pop)
## Play with different times per period
if (n_pop > 1) {
  times[n_pop] <- 15
}

## number of sampling events occurring over 'times'
samp         <- vector("list", n_pop)
samp_prop    <- 0.5 ## on average, what proportion of times are sampling events

for (i in 1:n_pop) {
 samp[[i]] <- rbinom(periods[i, 1], times[i, 1], samp_prop)
}

## number of time periods that elapse between the on-season
between_season_duration <- matrix(data = rpois(n_pop, 30) #rep(40, n_pop)
  , ncol = 1, nrow = n_pop)   

## random = sampling occurs on a random subset of possible days
when_samp <- matrix(data = rep("random", n_pop)
  , ncol = 1, nrow = n_pop)    

## vector to designate the offseasons
inbetween    <- vector("list", n_pop)
for (i in 1:n_pop) {
 inbetween[[i]] <- seq(1.5, periods[i, 1], by = 1)
}

####
## bd_parameters
####

## intercept and slope coefficients
bd_beta <- matrix(
  data = c(
    sample(seq(-3, 3, length = n_pop), n_pop)       ## Intercept
  , rep(-0.1, n_pop) ## Time effect
  , rep(0.3, n_pop)  ## Linear effect of temp on bd
), nrow = n_pop, ncol = 3, byrow = F)

## For debugging of a single population
if (n_pop == 1) {
 bd_beta[1, 1] <- 0
}

## error
bd_sigma  <- matrix(data = rep(20, n_pop)
  , ncol = 1, nrow = n_pop)  

## random effect variance covariance
bd_theta  <- matrix(
  data = c(
    rep(0.5, n_pop)  ## Intercept
  , rep(0.3, n_pop)  ## Time effect
  , rep(0.2, n_pop)  ## Linear effect of temp on bd
), nrow = n_pop, ncol = 3, byrow = F)

## logistic response coefficients for mortality across log(bd_load)
bd_mort <- matrix(
  data = c(
    sample(seq(-0.05, -0.15, length = n_pop), n_pop)      ## logistic slope
  , rep(4, n_pop)        ## intercept
), nrow = n_pop, ncol = 2, byrow = F)

## For debugging of a single population
if (n_pop == 1) {
 bd_mort[1, 1] <- -0.1
}

## logistic response coefficients for detection across log(bd_load)
bd_detect <- matrix(
  data = c(
    rep(0.2, n_pop)   ## logistic slope
  , sample(seq(0, -1, length = n_pop), n_pop) ## intercept 
), nrow = n_pop, ncol = 2, byrow = F)

bd_noinf <- matrix(
  data = c(
    rep(0.1, n_pop)
  ), nrow = n_pop, ncol = 1, byrow = F)

## proportion of all captures with bd swabs taken
bd_perc <- vector("list", n_pop)

for (k in 1:n_pop) {
 bd_perc[[k]] <- rep(.50, length(samp[[k]]))
}

## Observation noise in bd
obs_noise <- matrix(data = rep(1, n_pop)
  , nrow = n_pop, ncol = 1)  

####
## other parameters needed for the simulation function
####

## mortality probability in-between periods
p_mort          <- matrix(data = rlnorm(n_pop, -2.5, 0.5) * -1
  , nrow = n_pop, ncol = 1)

background_mort <- matrix(data = rep(1, n_pop)
  , nrow = n_pop, ncol = 1)

## Values that are used if the sim is set up to try and better approximate primary and secondary periods
use_prim_sec    <- TRUE

p_mort_within   <- matrix(data = plogis(rnorm(n_pop, 3, 0.5))
  , nrow = n_pop, ncol = 1)

## time gaps over which we assume the population is closed
sec_per_time    <- matrix(data = 1, nrow = n_pop, ncol = 1)

####
## other parameters for outside the function
####

## Stan model parameters
stan.iter     <- 1500
stan.burn     <- 500
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin
