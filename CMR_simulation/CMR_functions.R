#################################
## Functions for CMR simulator ##
#################################

'%notin%' <- Negate('%in%')

## Simulate bd for a population
bd.simulate <- function (
  periods, times, all_ind
, bd_beta, bd_sigma, bd_theta
, bd_noinf, obs_noise
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
  bd_load   = simulate(~times + temp + (1 + times | ind)
  , nsim    = 1
  , family  = Gamma(link = "log")
  , newdata = expdat
  , newparams = list(
      beta  = bd_beta
    , sigma = bd_sigma
    , theta = bd_theta
    )
  )$sim_1
) %>% mutate(
  bd_load     = rlnorm(n(), log(bd_load), obs_noise)
, bd_load     = round(bd_load, digits = 0)
)

## pick a random subset of individuals to never have gotten infected 
nev_inf <- sample(seq(all_ind), round(bd_noinf * all_ind))

if (length(nev_inf) > 0) {
 expdat[expdat$ind %in% nev_inf, ]$bd_load <- 0
}

expdat %<>% mutate(
  log_bd_load = log(bd_load)
, log_bd_load = round(log_bd_load, digits = 1)
) %>% mutate(
  log_bd_load = ifelse(is.infinite(log_bd_load), 0, log_bd_load)
)

return(expdat)
  
}

## Simulate the whole population CMR with bd sampling and all of the extra pieces needed to run a model
 ## for a single population
bd.sampling <- function (
  expdat
, all_ind, new_ind, times, periods, when_samp, samp, bd_perc
, inbetween, between_season_duration
, bd_mort, bd_detect, p_mort, background_mort, pop_ind
) {

## drop new_ind random individuals from each prior period to simulate individuals that migrated into the population
 ## after the first period
which_new_ind <- vector("list", length(new_ind))
temp_ind      <- seq(all_ind)
for (i in seq_along(new_ind)) {
  temp_samp          <- sort(sample(temp_ind, new_ind[i]))
  which_new_ind[[i]] <- temp_samp
  temp_ind           <- temp_ind[-which(temp_ind %in% temp_samp)]
}

which_new_ind            <- reshape2::melt(which_new_ind)
names(which_new_ind)     <- c("ind", "ind_gained")
which_new_ind$ind_gained <- which_new_ind$ind_gained + 1
which_new_ind            <- rbind(
  data.frame(
   ind        = temp_ind
 , ind_gained = 1
)
, which_new_ind
)

## Not dynamic yet, need to update for periods > 2
expdat %<>% 
  mutate(ind = as.numeric(ind)) %>% 
  left_join(., which_new_ind) 

## On which days is sampling occurring?
 ## NOTE: Oct 15: No other sampling scheme is currently supported
if (when_samp == "random") {
  
 ## Allow different sampling days by period
 sampling_days <- expand.grid(
  times   = seq(times)
, periods = seq(periods)
   ) %>% group_by(periods) %>% 
   mutate(
     sampling_days = ifelse(times %in% sample(times, samp[periods[1]]), 1, 0)
       ) %>% 
   mutate(all_times = interaction(times, periods)) %>% 
   mutate(all_times = as.numeric(all_times))
 
 ## Determine the number of time periods that elapse between back to back samples
 time_gaps   <- sampling_days %>% filter(sampling_days == 1) %>%
   ungroup() %>% mutate(time_gaps = (all_times - lag(all_times, 1))) %>%
   group_by(periods) %>%  
   mutate(last_sample = min(times)) %>%
   filter(!is.na(time_gaps))
 
 time_gaps[(time_gaps$times == time_gaps$last_sample) & (time_gaps$periods != min(time_gaps$periods)), ]$time_gaps <- 
   time_gaps[(time_gaps$times == time_gaps$last_sample) & (time_gaps$periods != min(time_gaps$periods)), ]$time_gaps + between_season_duration
   
  time_gaps <- time_gaps$time_gaps
  
  if (length(time_gaps) != (sum(samp) - 1)) {
    print("time_gaps is an incorrect length, check parameters"); break
  }
 
 ## Also need a vector of when sampling occurs
 sampling_times_all <- (sampling_days %>% filter(sampling_days == 1))$all_times
 
}

## grab the log of the range of Bd load
bd_range <- round(with(expdat, seq(min(log_bd_load), max(log_bd_load), by = .1)), 1)

## establish "true" relationship between load and mortality probability and detection probability
 ## Note: here "mort" is actually survival probability from time t to t+1
bd_probs <- data.frame(
  log_bd_load  = bd_range
, mort         = exp(bd_mort[2] + bd_mort[1] * bd_range) /
    (1 + exp(bd_mort[2] + bd_mort[1] * bd_range))
) %>% left_join(.
  , data.frame(
  log_bd_load  = bd_range
, detect       = exp(bd_detect[2] + bd_detect[1] * bd_range) /
    (1 + exp(bd_detect[2] + bd_detect[1] * bd_range))
)
  )

## Add the offseason survival probability to the simulated data for the individuals
 ## (All we really need here is the mort column, but need the other ones to add it to the
  ## other data frame to simulate cumulative survival)
expdat %<>% 
  left_join(., bd_probs) %>%
  group_by(ind) %>%
  arrange(periods, ind, times)

## calculate a summary of the bd load within a period for mortality between periods
expdat.bd_sum <- expdat %>% 
  group_by(periods, ind) %>%
  summarize(
    max_bd = max(log_bd_load)
  , cum_bd = sum(log_bd_load)
  )

 off_season <- expand.grid(
     periods     = inbetween
   , times       = 1
   , ind         = seq(all_ind)
   , temp        = 0
   , bd_load     = 0
   , log_bd_load = 0
   , mort        = 1 - p_mort
   , detect      = 0
    )
 
off_season %<>% left_join(., which_new_ind)

for (i in 1:nrow(off_season)) {
  temp_per <- off_season[i, ]$periods
  temp_ind <- off_season[i, ]$ind
  
  temp_row <- expdat.bd_sum[
    expdat.bd_sum$periods == temp_per - 0.5 &
    expdat.bd_sum$ind     == temp_ind
  , 
  ]
  
## if (p_mort_type == "max")
 ## NOTE: OCT 15: only max supported for now 
  off_season[i, ]$log_bd_load <- temp_row$max_bd
  off_season[i, ]$mort        <- plogis(background_mort + p_mort * off_season[i, ]$log_bd_load)

}

expdat %<>% rbind(off_season, .)

## For all individuals in periods prior to when they were gained set their detection to 0
expdat[(expdat$ind_gained > expdat$periods), ]$detect <- 0
expdat[(expdat$ind_gained > expdat$periods), ]$mort   <- 1

expdat %<>% arrange(ind, periods, times)
expdat %<>% ungroup() %>% group_by(ind)
expdat %<>% mutate(cum_surv = cumprod(mort)) 

####
## create the true state of the population from the simulated bd values
####

expdat %<>%
  mutate(dead = rbinom(n(), 1, 1 - mort)) %>%
  group_by(ind) %>% 
  mutate(
    dead = cumsum(dead)
  , dead = ifelse(dead > 1, 1, dead)) 

expdat %<>% filter(periods %notin% inbetween)

## On sampling days check for detection
expdat %<>% left_join(., sampling_days) %>%
  mutate(detected     = ifelse(dead == 0 & sampling_days == 1, rbinom(n(), 1, detect), 0))

####
## Create the sampling data from the simulated population 
####

## Drop all individuals that were never caught and update parameters
never_detected <- expdat %>% 
  group_by(ind) %>% 
  summarize(total_detection = sum(detected)) %>% 
  filter(total_detection == 0)

expdat %<>% dplyr::filter(ind %notin% never_detected$ind) %>% droplevels() %>%
  mutate(ind = as.factor(ind)) %>% mutate(ind = as.numeric(ind))

## total number of individuals ever captured
all_ind     <- length(unique(expdat$ind))

expdat.capt <- expdat %>% filter(sampling_days == 1)

## Also check which individuals were seen in each sampling period
ind_seen <- expdat.capt %>% 
  group_by(periods, ind) %>% 
  summarize(seen = sum(detected)) %>%
  mutate(seen = ifelse(seen > 0, 1, seen)) %>% 
  filter(seen == 1) %>%
  group_by(ind) %>%
  slice(1)

present <- matrix(
  data = 0
, ncol = periods
, nrow = all_ind)

recruited_inds <- (expdat %>% 
  group_by(ind) %>% 
  summarize(new_ind = mean(ind_gained)))$new_ind

for (i in seq_along(recruited_inds)) {
 present[i, recruited_inds[i]:periods] <- 1
}

## Create the capture array in the correct structure for stan.
capture_matrix <- matrix(
  data = (expdat.capt %>% arrange(ind, all_times))$detected
, nrow = all_ind
, ncol = sum(samp)
, byrow = T
)

## becuase of ind as a factor things could get messed up. Check to make sure order was retained
if (sum(capture_matrix[2, ] == 
    (expdat.capt %>% arrange(ind, all_times) %>% 
        filter(ind == unique(expdat.capt$ind)[2]))$detected) != sum(samp)) {
  print("capture_matrix filled in incorrectly"); break
}

## First create an index vector for the offseason and then drop the offseason from expdat
 ## Note: offseason_vec is associated with the phi values, so offseason is the last entry of each on-season
offseason_vec   <- rep(0, sum(samp))[-sum(samp)]
which_offseason <- cumsum(samp)
which_offseason <- which_offseason[-length(which_offseason)]
offseason_vec[which_offseason] <- 1

## Across all sampling_days across all seasons determine when each individual was first and last captured
capture_range <- expdat %>% 
  group_by(ind) %>%
  filter(sampling_days == 1) %>%  
  summarize(
    first = min(which(detected == 1))
  , final = max(which(detected == 1))
    ) %>%
  mutate(
    first = ifelse(is.infinite(first), 0, first)
  , final = ifelse(is.infinite(final), 0, final)
    )

## Will be used eventually for generated quantities
capture_total <- expdat %>% 
  group_by(all_times) %>% 
  filter(sampling_days == 1) %>%
  summarize(total_capt = sum(detected)) %>%
  arrange(periods
    , all_times
    )

####
## Set up the bd sampling data
####

measured_bd <- matrix(
  data = (expdat.capt %>% arrange(ind, all_times))$log_bd_load
, nrow = all_ind
, ncol = sum(samp)
, byrow = T
)

## Establish which sampling days included bd swabs
expdat %<>% left_join(.
  , {
    expdat %>% 
      filter(sampling_days == 1, detected == 1) %>% 
      ungroup() %>% 
      mutate(bd_swabbed = rbinom(n(), 1, bd_perc[periods]))
  }) 
  
expdat %<>% mutate(
  bd_swabbed = ifelse(is.na(bd_swabbed), 0, bd_swabbed)
)

expdat %<>% mutate(pop = pop_ind)
  
expdat.swabbed <- expdat %>% filter(sampling_days == 1)

## And finally, the individuals on the sampling days that were swabbed for bd
bd_measured <- matrix(
  data = (expdat.swabbed %>% arrange(ind, all_times))$bd_swabbed
, nrow = all_ind
, ncol = sum(samp)
, byrow = T
)

## Another check of another matrix
if (sum(bd_measured[2, ] == 
    (expdat.swabbed %>% arrange(ind, all_times) %>% 
        filter(ind == unique(expdat.swabbed$ind)[2]))$bd_swabbed) != sum(samp)) {
  print("measured_bd filled in incorrectly"); break
}

## And finally, finally, vectors to designate periods for times and sampling occasions
periods_time <- (expdat %>% filter(ind == unique(expdat$ind)[1]))$periods
periods_occ  <- (expdat %>% filter(ind == unique(expdat$ind)[1]
  , sampling_days == 1))$periods
 
return(list(
  expdat             = expdat
, all_ind            = all_ind
, sampling_times_all = sampling_times_all
, offseason_vec      = offseason_vec
, periods_time       = periods_time
, periods_occ        = periods_occ
, measured_bd        = measured_bd
, bd_measured        = bd_measured
, time_gaps          = time_gaps
, capture_matrix     = capture_matrix
, capture_range      = capture_range
, present            = present
, bd_probs           = bd_probs
, off_season         = off_season
))

}

## Organize the simulated population into long form for the stan model
bd.stan_org <- function (
  one_pop
, pop_ind
, times
, periods
, samp
) {
  
## Indices to determine which of the latent bd values can be informed by measured bd
X_bd.m <- reshape2::melt(one_pop$measured_bd)
names(X_bd.m) <- c("ind", "times", "bd")
X_bd.m %<>% mutate(times = plyr::mapvalues(times, from = seq(sum(samp)), to = one_pop$sampling_times_all))
bd_m.m <- reshape2::melt(one_pop$bd_measured)
names(bd_m.m) <- c("ind", "times", "meas")
bd_m.m %<>% mutate(times = plyr::mapvalues(times, from = seq(sum(samp)), to = one_pop$sampling_times_all))
all_bd <- expand.grid(
  ind   = seq(one_pop$all_ind)
, times = seq(times * periods)
)
all_bd %<>% left_join(., bd_m.m)
samp_indices <- which(all_bd$meas == 1)

X_bd.m %<>% left_join(., bd_m.m)
X_bd.m %<>% filter(meas == 1)
X_bd.m %<>% arrange(ind)

## Indices to determine which entries of phi need to be set to 0
y.m    <- reshape2::melt(one_pop$capture_matrix) %>% arrange(Var1)
temp   <- (one_pop$expdat %>% group_by(all_times) %>% summarize(mean_temp = mean(temp)))$mean_temp

## Indices for which entries of phi must be 0
phi_zeros <- matrix(data = 0, nrow = one_pop$all_ind, ncol = sum(samp) - 1)
for (i in 1:one_pop$all_ind) {
  phi_zeros[i, ] <- c(
    rep(1, one_pop$capture_range$first[i] - 1)
  , rep(0, ncol(phi_zeros) - (one_pop$capture_range$first[i] - 1))
    )
}
phi_zeros   <- (phi_zeros %>% reshape2::melt() %>% arrange(Var1))$value

## Index of which entries of phi are each individual
indiv_index <- rep(sum(samp), each = one_pop$all_ind)

## Index for whether an individual was known to be present in the population at each time
p_zeros <- matrix(data = 0, nrow = one_pop$all_ind, ncol = sum(samp))
for (i in 1:one_pop$all_ind) {
  p_zeros[i, ] <- rep(one_pop$present[i, ], samp)
}
p_zeros   <- (p_zeros %>% reshape2::melt() %>% arrange(Var1))$value

ind_time <-  data.frame(
   ind_time_rep  = rep(seq(one_pop$all_ind), each = times * periods)
 , time_rep      = rep(seq(times), one_pop$all_ind * periods)
)

## housing the necessary pieces for the long-form stan model for convenience in exporting
ind_occ_phi <- data.frame(
   ind_occ_min1_rep    = rep(seq(one_pop$all_ind), each = (sum(samp) - 1))
 , sampling_events_phi = rep(one_pop$sampling_times_all[-sum(samp)], one_pop$all_ind)
 , offseason           = rep(one_pop$offseason_vec, one_pop$all_ind)
 , time_gaps           = rep(one_pop$time_gaps, one_pop$all_ind)
 , phi_zeros           = phi_zeros
 , pop                 = rep(pop_ind, length(phi_zeros))
) %>% mutate(pop = pop_ind, ind = interaction(ind_occ_min1_rep, pop))

ind_occ_p <- data.frame(
   ind_occ_rep       = rep(seq(one_pop$all_ind), each = (sum(samp)))
 , sampling_events_p = rep(one_pop$sampling_times_all, one_pop$all_ind)
 , periods_occ       = rep(one_pop$periods_occ, one_pop$all_ind)
 , captures          = y.m$value
 , p_zeros           = p_zeros
 , pop               = rep(pop_ind, length(p_zeros))
) %>% mutate(pop = pop_ind, ind = interaction(ind_occ_rep, pop))

one_pop$capture_range %<>% mutate(pop = 1, ind = interaction(ind, pop))

ind_occ_size <- rep(sum(samp), one_pop$all_ind)

X_bd.m %<>% mutate(pop = pop_ind)

## Add the extra tracked items to one_pop
one_pop      <- c(one_pop, ind_occ_size = list(ind_occ_size))

## Population-specific convariates
pop_cov.bd  <- data.frame(
  temp = temp
, time = seq(times * periods)
, pop  = rep(pop_ind, length(temp))
)
  
return(list(
  one_pop     = one_pop
, pop_cov.bd  = pop_cov.bd
, ind_occ_phi = ind_occ_phi
, ind_occ_p   = ind_occ_p
, ind_time    = ind_time
, X_bd.m      = X_bd.m
))

}
