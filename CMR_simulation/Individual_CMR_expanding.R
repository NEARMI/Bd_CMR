###############################################################
## Simulate data for and fit an individual-CMR model in Stan ##
###############################################################

########
## See Individual_CMR_few_examples.R for the precursor to this script. Starting a new script instead of continuing with that one
## for now -- if this evolves into some form of power analysis of design study will convert everything to a streamlined script
 ## This script takes the fundamental within-season single-site model from Individual_CMR_few_examples.R and expands it out
 ## to multiple seasons and multiple sites
########

####
## Notes as of OCT 13:
####

## Second population fully working. Now just need to clean up the simulation code so that it is obvious
 ## what form the real data needs to be in for the model and for easier simulations for design and power analyses
## Also need to do some more extensive debugging and exploration of parameter values to make sure that the CI sensibly
 ## match how much info each individual contains

## After ^^ 
 ## 1) Simulate with a few different data structures and quantities for a few slides for Thursday meeting
 ## 2) Add in other covariates and make sure the model can easily accommodate multiple covariates that vary by population
  ## -- will need a "population" index vector to index all of the covariate vectors
 ## 3) Work on exapanding the complexity of the bd submodel
 
### Continued concerns from the other day:
 ## A) simplex to control for entry into the population and biased detection -- can we also use the estimated parameter
  ## to scale survival of individuals we think are a year older?
 ## B) modeling of bd in all periods is probably ok because unmeasured individuals don't affect the likelihood
 ## C) no latent states could make it tricky to build a better model of disease states, but that is probably ok if we can get intercepts really
  ## low because we are on a log scale 

### Before the CMR meeting, need to simulate with:
 ## A) Very high within-season survival
 ## B) Very sparse data collection for a subset of individuals

## The short term plan in advance of the CMR modeling meeting on the 14th will be to:
 ## -- MODELING --
 # 1) Figure out what happens when some individuals don't get infected

####
## Quick notes about next steps
####

##### Second population #####
## 1) First write the simulation code into a function, and then run that function twice with different parameters and stick
 ## the result together 
## 2) In the model add a loop for population? -- Need to give some serious thought to this step because it is going to be what is
 ## really important for scaling the model to all of the data: how to allow n_occasions to vary by population, how to allow n_times
  ## to vary by location (because the on-season may differ by place etc.) -- It would be great to add all of this into the same
   ## array structure, but it may be difficult to time??? (maybe not if times all line up)

##### Other covaraites #####
## 1) After adding different locations add a covariate of temperature that varies differently over the season and by location
 ## -- This will be a key step in determining how to model the random effect of individual bd vs other covaratie

####
## Packages and functions
####
needed_packages <- c("magrittr", "dplyr", "tidyr", "lme4", "ggplot2", "rstan")
lapply(needed_packages, require, character.only = TRUE)
source("../../ggplot_theme.R")
set.seed(10002)
'%notin%' <- Negate('%in%')

####
## Parameters
####

## "Design" parameters
nsim      <- 1                    ## number of simulations (1 to check model, could be > 1 for some sort of power analysis or something)
ind       <- 20                   ## number of individuals in the population being modeled
periods   <- 3                    ## number of primary periods (years in most cases)
new_ind   <- rep(4, periods - 1) ## individuals added in each new period
inbetween <- seq(1.5, periods, by = 1)
all_ind   <- ind + sum(new_ind)   ## number of individuals ever to exist in the population
times     <- 20                   ## number of time periods (in the real data probably weeks; e.g., May-Sep or so)
samp      <- 3                    ## number of sampling events occurring over 'times' (e.g., subset of 'times' weeks when sampling occurred)
if (periods > 1) {
samp <- rep(samp, periods)      ## for now assume same number of periods per year, but this model allows variable sampling dates by season
between_season_duration <- 10   ## number of time periods that elapse between the on-season
}
when_samp <- "random"           ## random = sampling occurs on a random subset of possible days

## Two ways to simulated data
sim_opt <- "lme4"               ## Use the built in capabilities of lme4 to simulate bd loads

## Type of simulation
sim_type <- "temp"              ## Time or temp available so far BUT (OCT 7: Temp just starting)

## Response of individuals to Bd load
bd_mort   <- c(decay = -0.2, offset = 6)       ## logistic response coefficients for mortality across log(bd_load)
bd_detect <- c(decay = 0.1, offset = -0.5)     ## logistic response coefficients for detection across log(bd_load)

## Other parameters

## mortality probability in-between periods
p_mort_type <- "max" ## con = single value; max = based on max bd; cum = based on cumulative bd
if (p_mort_type == "con") {
p_mort    <- c(0.30) 
} else if (p_mort_type == "max") {
p_mort    <- 0.02
} else if (p_mort_type == "cum") {
p_mort    <- 0.0001  
}

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

## observation noise in bd (** SEPT 21: yes, on the log scale so this is weird to have the same var across log bd load, will clean this up later)
obs_noise <- 0.25   

## Stan model parameters
stan.iter     <- 1500
stan.burn     <- 500
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin

####
## Simulate bd response and relationship between bd and survival and detection
## as well as create index vectors for when sampling is occurring 
####

if (sim_type == "time") {
  
## Bd parameters 
bd_beta   <- c(
  7     ## Intercept
, 1     ## Period effect
, 30    ## Linear effect of time on bd
, -50   ## Quadratic effect of time on bd
) 
bd_sigma  <- 2              ## observation noise
bd_theta  <- c(2)           ## random effect variance covariance 

## Simulate data using lme4 mixed model structure
expdat <- expand.grid(
  periods = seq(periods)
, times   = seq(times)     
, ind     = factor(seq(all_ind))
  )

expdat %<>% mutate(
  bd_load   = simulate(~periods + poly(times, 2) + (1 | ind)
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
  
} else if (sim_type == "temp") {
  
## Bd parameters 
bd_beta   <- c(
  1     ## Intercept
, .1     ## Time effect
, .3     ## Linear effect of temp on bd
) 
bd_sigma  <- 2              ## observation noise
bd_theta  <- c(5)           ## random effect variance covariance 
  
## Simulate data using lme4 mixed model structure
expdat <- expand.grid(
  periods = seq(periods)
, times   = seq(times)     
, ind     = factor(seq(all_ind))
  )

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

}

## drop new_ind random individuals from each prior period
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
   ind    = temp_ind
 ,ind_gained = 1
)
, which_new_ind
)

## Not dynamic yet, need to update for periods > 2
expdat %<>% mutate(ind = as.numeric(ind)) %>% left_join(., which_new_ind) 

## On which days is sampling occurring?
if (when_samp == "random") {
  
 ## Allow different sampling days by period
 sampling_days <- expand.grid(
  times   = seq(times)
, periods = seq(periods)
   ) %>% group_by(periods) %>% 
   mutate(
     sampling_days = ifelse(times %in% sample(times, samp), 1, 0)
       ) %>% 
   mutate(all_times = interaction(times, periods)) %>% 
   mutate(all_times = as.numeric(all_times))
 
 ## Determine the number of time periods that elapse between back to back samples
 time_gaps     <- sampling_days %>% filter(sampling_days == 1) %>%
   group_by(periods) %>% mutate(time_gaps = (times - lag(times, 1))) %>%
   filter(!is.na(time_gaps)) %>% rbind(.
     , {
    expand.grid(
     periods       = inbetween
   , times         = 1
   , sampling_days = 1
   , time_gaps     = between_season_duration
    )
     }) %>% arrange(periods)
   
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
, mort         = exp(bd_mort["offset"] + bd_mort["decay"] * bd_range) /
    (1 + exp(bd_mort["offset"] + bd_mort["decay"] * bd_range))
) %>% left_join(.
  , data.frame(
  log_bd_load  = bd_range
, detect       = exp(bd_detect["offset"] + bd_detect["decay"] * bd_range) /
    (1 + exp(bd_detect["offset"] + bd_detect["decay"] * bd_range))
)
  )

bd_probs %>% pivot_longer(cols = c(2, 3)) %>% {
  ggplot(., aes(log_bd_load, value)) +
  geom_line(aes(colour = name)) + 
  scale_colour_manual(
      name   = "Relationship"
    , values = c("firebrick3", "dodgerblue3")
    , labels = c("Detection", "Survival")) + 
  ylab("Prediction") + 
  xlab("log of Bd load") +
  theme(legend.key.size = unit(0.75, "cm"))
}

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
 
if (sim_type != "temp") {
 off_season %<>% dplyr::select(-temp)
} 

off_season %<>% left_join(., which_new_ind)

for (i in 1:nrow(off_season)) {
  temp_per <- off_season[i, ]$periods
  temp_ind <- off_season[i, ]$ind
  
  temp_row <- expdat.bd_sum[
    expdat.bd_sum$periods == temp_per - 0.5 &
    expdat.bd_sum$ind     == temp_ind
  , 
  ]
  
  if (p_mort_type == "max") {
  off_season[i, ]$log_bd_load <- temp_row$max_bd
  off_season[i, ]$mort        <- with(off_season[i, ]
    , 1 - ((1 - mort) * log_bd_load))
  } else if (p_mort_type == "cum") {
  off_season[i, ]$log_bd_load <- temp_row$cum_bd  
  off_season[i, ]$mort        <- with(off_season[i, ]
    , 1 - ((1 - mort) * log_bd_load))
  } else {
  ## do nothing
  }
  
}

expdat %<>% rbind(off_season, .)

## For all individuals in periods prior to when they were gained set their detection to 0
expdat[(expdat$ind_gained > expdat$periods), ]$detect <- 0
expdat[(expdat$ind_gained > expdat$periods), ]$mort   <- 1

expdat %<>% arrange(periods, ind, times)

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

## store original simulated population for diagnostics (on population size for example)
expdat.not_dropped <- expdat    
ind.not_dropped    <- all_ind

expdat %<>% dplyr::filter(ind %notin% never_detected$ind) %>% droplevels()

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

## full "true" bd values to be subset into bd values on the days in which sampling occurred (some of which
 ## get to inform the model -- those that were actually sampled as determined below)
temp_bd_dat <- (expdat.capt %>% arrange(ind, all_times))$log_bd_load
temp_bd_dat <- rnorm(length(temp_bd_dat), temp_bd_dat, obs_noise)

measured_bd <- matrix(
  data = temp_bd_dat
, nrow = all_ind
, ncol = sum(samp)
, byrow = T
)

## becuase of ind as a factor things could get messed up. Check to make sure order was retained
bd_load_check <- (measured_bd[3, ] - (expdat.capt %>% arrange(ind, all_times) %>% filter(ind %in% unique(expdat.capt$ind)[3]))$log_bd_load)
print(paste("This number should be pretty close to 0", round(sum(bd_load_check), 2), sep = ": "))

## Establish which sampling days included bd swabs
expdat %<>% left_join(.
  , {
    expdat %>% 
      filter(sampling_days == 1, detected == 1) %>% 
      ungroup() %>% 
      mutate(bd_swabbed = rbinom(n(), 1, bd_perc))
  }) 
  
expdat %<>% mutate(
  bd_swabbed = ifelse(is.na(bd_swabbed), 0, bd_swabbed)
)
  
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

## Double check to make sure this produces sensible capture data
expdat.plottest     <- expdat %>% arrange(ind_gained)
expdat.plottest$ind <- factor(expdat.plottest$ind, levels = rev(unique(expdat.plottest$ind)))

if (ind <= 100) {
  expdat.plottest %>% {
   ggplot(., aes(times, ind, fill = as.factor(detected))) + 
     geom_tile(aes(alpha = sampling_days)) +
     scale_x_continuous(breaks = c(1, 5, 10)) +
     xlab("Time") + 
     ylab("Individual") +
     scale_fill_manual(
         values = c("dodgerblue3", "firebrick3")
       , name   = "Detected?"
       , labels = c("No", "Yes")) +
     guides(alpha = FALSE) +
     geom_line(data = expdat.plottest %>% filter(dead == 1)
       , aes(x = times, y = ind, z = NULL), lwd= 0.5, alpha = 0.5) +
     geom_point(data = expdat.plottest %>% filter(bd_swabbed == 1)
       , aes(x = times, y = ind, z = NULL), lwd = 0.5) +
     facet_wrap(~periods) +
     theme(
       axis.text.y = element_text(size = 8)
     , legend.text = element_text(size = 12)
     , legend.key.size = unit(.55, "cm")
     ) +
     ggtitle("Lines show dead individuals; dots show bd swabbs")
 }
}

####
## Run the model in Stan
####

stan_data     <- list(
  ## dimensional indexes 
   n_periods       = periods
 , n_ind           = all_ind                  
 , n_times         = times * periods
 , times_within    = times
 , n_occasions     = sum(samp)
 , n_occ_min1      = sum(samp) - 1
  ## vector indexes 
 , time            = rep(seq(times), periods)
 , time_per_period = matrix(data = seq(times * periods), nrow = times, ncol = periods)
 , sampling_events = sampling_times_all
 , offseason       = offseason_vec
 , periods         = periods_time
 , periods_occ     = periods_occ
  ## covariates
 , X_bd            = measured_bd
 , X_measured      = bd_measured
 , temp            = (expdat %>% group_by(all_times) %>% summarize(mean_temp = mean(temp)))$mean_temp
 , time_gaps       = time_gaps
  ## Capture data
 , y               = capture_matrix
 , first           = capture_range$first
 , last            = capture_range$final
 , present         = present
  )

stan.fit  <- stan(
  file    = "CMR_ind_pat_bd-p-phi_multi_recruit_free_temp.stan"
, data    = stan_data
, chains  = 1
, iter    = stan.iter
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.92, max_treedepth = 12)
  )

saveRDS(stan.fit, "stan.fit_multi.Rds")

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)

####
## Recovery of simulated coefficients?
####

## Primary Bd effects
pred_coef        <- as.data.frame(stan.fit.summary[c(1:4), c(4, 6, 8)])
names(pred_coef) <- c("lwr", "mid", "upr")
pred_coef        %<>% mutate(param = rownames(.))

pred_coef %>% {
  ggplot(., aes(param, mid)) + geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
    geom_point(data = data.frame(
      param = pred_coef$param
    , mid   = c(
      rev(bd_mort)
    , rev(bd_detect)
      )
    ), colour = "firebrick3") +
    xlab("Parameter") + ylab("Estimate") +
    scale_x_discrete(labels = c("Detection intercept", "Detection slope", "Survival intercept", "Survival slope")) +
    theme(axis.text.x = element_text(size = 11))
}

## survival over Bd load
stan.pred        <- apply(stan.fit.samples$beta_phi, 1
  , FUN = function(x) plogis(x[1] + x[2] * bd_probs$log_bd_load)) %>%
  reshape2::melt(.)
names(stan.pred) <- c("log_bd_load", "sample", "mortality")
stan.pred        %<>% mutate(log_bd_load = plyr::mapvalues(log_bd_load
  , from = unique(log_bd_load), to = bd_probs$log_bd_load))
stan.pred        %<>% group_by(log_bd_load) %>% 
  summarize(
    lwr = quantile(mortality, c(0.025))
  , mid = quantile(mortality, c(0.5))
  , upr = quantile(mortality, c(0.975))
  )

stan.pred %>% {
  ggplot(., aes(log_bd_load, mid)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_line() +
    geom_line(data = (bd_probs %>% mutate(mort = mort ))
      , aes(log_bd_load, mort)
      , colour = "dodgerblue4", lwd = 2) +
    xlab("Log of Bd Load") + ylab("Predicted mortality probability")
}

## detection over Bd load
stan.pred        <- apply(stan.fit.samples$beta_p, 1
  , FUN = function(x) plogis(x[1] + x[2] * bd_probs$log_bd_load)) %>%
  reshape2::melt(.)
names(stan.pred) <- c("log_bd_load", "sample", "detect")
stan.pred        %<>% mutate(log_bd_load = plyr::mapvalues(log_bd_load
  , from = unique(log_bd_load), to = bd_probs$log_bd_load))
stan.pred        %<>% group_by(log_bd_load) %>% 
  summarize(
    lwr = quantile(detect, c(0.10))
  , mid = quantile(detect, c(0.5))
  , upr = quantile(detect, c(0.90))
  )

stan.pred %>% {
  ggplot(., aes(log_bd_load, mid)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_line() +
    geom_line(data = bd_probs, aes(log_bd_load, detect)
      , colour = "dodgerblue4", lwd = 2) +
    xlab("Log of Bd Load") + ylab("Predicted detection probability")
}

## Change in Bd
pred_coef <- as.data.frame(stan.fit.summary[
  c(
    grep("bd_delta_mu", dimnames(stan.fit.summary)[[1]])
  , grep("bd_delta_sigma", dimnames(stan.fit.summary)[[1]])
  , grep("bd_add", dimnames(stan.fit.summary)[[1]])
  , grep("bd_obs", dimnames(stan.fit.summary)[[1]])
  , grep("start_mean", dimnames(stan.fit.summary)[[1]])
  , grep("start_var", dimnames(stan.fit.summary)[[1]])
  )
  , c(4, 6, 8)])
names(pred_coef) <- c("lwr", "mid", "upr")
pred_coef        %<>% mutate(param = rownames(.))

## Hard to directly compare coefficient estimates given the change in scale, need to work on this
pred_coef %>% filter(param != "start_mean") %>% {
  ggplot(., aes(param, mid)) + geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
    xlab("Parameter") + ylab("Estimate") +
    theme(axis.text.x = element_text(size = 11))
}

####
## Looking at individual phi values
####

data.frame(
  capture_matrix[2, -30]
, apply(stan.fit.samples$phi, 2:3, mean)[2, ]
)

apply(stan.fit.samples$phi, 2:3, mean)

####
## Individual random effect estimates
####

test_compare <- cbind(
  data.frame(
  mean = colMeans(stan.fit.samples$bd_delta_eps)
, sd   = apply(stan.fit.samples$bd_delta_eps, 2, sd)
  )
, expdat %>% group_by(ind) %>% 
        summarize(ind_rand = sum(log_bd_load))
)

ggplot(test_compare, aes(ind_rand, mean)) + geom_point()

stan.ind_pred_eps <- stan.fit.samples$bd_delta_eps %>%
  reshape2::melt(.) %>% rename(ind = Var2, eps = value)
stan.ind_pred_var <- stan.fit.samples$bd_delta_sigma %>%
  reshape2::melt(.) %>% left_join(., stan.ind_pred_eps) %>%
  mutate(eps = eps * value) %>% group_by(ind) %>%
  summarize(
    mid = quantile(eps, 0.50)
  , lwr = quantile(eps, 0.025)
  , upr = quantile(eps, 0.975)
  ) %>% mutate(ind = unique(expdat$ind)) %>%
  arrange(mid) %>% 
  left_join(., expdat %>% group_by(ind) %>% summarize(total_capt = sum(bd_swabbed))) %>% 
  left_join(., expdat %>% group_by(ind, periods) %>% summarize(total_detect = sum(detected)) %>% 
      filter(total_detect > 0) %>% summarize(total_periods = n())) %>% 
  mutate(CI_width = upr - lwr) %>% 
  mutate(order_pred = seq(n()))

stan.ind_pred_var %<>% left_join(.
    , {
      expdat %>% group_by(ind) %>% 
        summarize(ind_rand = sum(log_bd_load)) %>% 
        arrange(ind_rand) %>%
        mutate(order_real = seq(n()))
    })

stan.ind_pred_var %>% mutate(ind_rand = (ind_rand - mean(ind_rand))/sd(ind_rand)) %>% {
  ggplot(., aes(mid, ind_rand)) + 
    geom_point() +
    xlab("Random effect deviate") + ylab("Raw cumulative bd load") 
}

stan.ind_pred_var %>% {
  ggplot(., aes(total_capt, CI_width)) + 
    geom_point(aes(colour = as.factor(total_periods))) +
    xlab("Total Swabs") +
    ylab("Width of CI") 
}

stan.ind_pred_var %>% {
  ggplot(., aes(order_real, order_pred)) + 
    geom_point() +
    xlab("Simulated bd rank") +
    ylab("Estimated bd rank")
}

## Not doing a very good job of making these look different than 0
stan.ind_pred_var %>% arrange(mid) %>% mutate(ind = factor(ind, levels = ind)) %>%  {
  ggplot(., aes(as.factor(ind), mid)) + geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr)) +
    xlab("Individual") + 
    ylab("Random Effect Deviate") +
    geom_hline(yintercept = 0, linetype = "dashed", lwd = 1, colour = "firebrick3") +
    theme(
      axis.text.x = element_text(size = 8)
    ) #+
    #geom_vline(xintercept = bd_drop.which, lwd = 0.1) +
    #ggtitle("thin dashed lines are the 20 left out of the model. It is at least nice to say the estimates are all different")
}
 
data.frame(
  simulated = order(bd_ind)
, real      = stan.ind_pred_var$ind
) %>% head(10)

## Percent of most extreme estimated individuals that are in the most extreme simulated individuals
 ## *NOTE*: Seems to recover the most extreme individuals with low slopes than individuals with large slopes
  ## I think this is probably because of the log-scale
length(which(head(stan.ind_pred_var$ind, 20) %in% head(order(bd_ind), 20)) ) / 20
length(which(tail(stan.ind_pred_var$ind, 20) %in% tail(order(bd_ind), 20)) ) / 20

####
## In the multi-season model look at probability of recruitment
####

test_recruit        <- apply(stan.fit.samples$gamma, 2:3, mean)
test_recruit        <- as.data.frame(test_recruit)
names(test_recruit) <- c("period_1", "period_2", "period 3")
present

ggplot(test_recruit, aes(truth, period_1)) + 
  geom_point()

####
## Recovery of the Bd profile of the individuals left out of the 
####

estimated_bd        <- stan.fit.samples$X[,bd_drop.which,] %>% 
  reshape2::melt(.)
names(estimated_bd) <- c("samps", "ind", "obs", "log_bd_load")
estimated_bd        <- estimated_bd %>% 
  group_by(ind, obs) %>% 
  summarize(
    mid = quantile(log_bd_load, 0.50)
  , lwr = quantile(log_bd_load, 0.025)
  , upr = quantile(log_bd_load, 0.975)
  )
observed_bd         <- matrix(nrow = ind, ncol = samp, data = expdat$log_bd_load, byrow = T)[bd_drop.which, -samp] %>% 
  reshape2::melt(.)
names(observed_bd)  <- c("ind", "obs", "log_bd_load")

## These are in fact different, which is exciting
estimated_bd %>% {
  ggplot(., aes(obs, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
    geom_line() +
    geom_line(data = observed_bd, aes(obs, log_bd_load), colour = "dodgerblue4") +
    xlab("Time") + ylab("Bd Load") +
    facet_wrap(~ind)
}

####
## Recovery of off-season mortality? 
####

between_season_test <- apply(stan.fit.samples$phi, 2:3, FUN = function(x) mean(x))
between_season_test <- reshape2::melt(between_season_test)
names(between_season_test) <- c("ind", "time", "phi")
between_season_test %>% mutate(phi = ifelse(phi == 0, NA, phi)) %>% {
  ggplot(., aes(time, phi)) + 
    geom_line(aes(group = ind)) +
    xlab("Time") + 
    ylab("Survival Probability") +
   geom_hline(yintercept = 0.8)
}

####
## Generated quantities: 
##  1) Population size
##  2) Population distribution of infection (for now just total)
##  3) Expected number to capture at each sampling event
####

## Population size
pred_coef        <- as.data.frame(stan.fit.summary[
  grep("pop", dimnames(stan.fit.summary)[[1]])
  , c(4, 6, 8)])
names(pred_coef) <- c("lwr", "mid", "upr")
pred_coef        %<>% mutate(param = rownames(.)) %>% 
  mutate(param = as.numeric(factor(param, levels = param)))


## SEP 9: Something a bit fishy going on here in my predictions and real pop as a function of time
pred_coef %>% mutate(param = param - 1) %>%  {
  ggplot(., aes(param, mid)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_line() +
    geom_line(data = pop_alive, aes(samp, pop_size), colour = "dodgerblue4", lwd = 1) +
    xlab("Time") + 
    ylab("Population Size") 
}

## Number captured -- interesting that this is perfectly constrained. I guess that makes sense?
pred_coef        <- as.data.frame(stan.fit.summary[
  grep("captures", dimnames(stan.fit.summary)[[1]])
  , c(4, 6, 8)])
names(pred_coef) <- c("lwr", "mid", "upr")
pred_coef        %<>% mutate(param = rownames(.)) %>% 
  mutate(param = as.numeric(factor(param, levels = param)))
num.captured     <- expdat %>% group_by(samp) %>% summarize(num_capt = sum(detected))

pred_coef %>% {
  ggplot(., aes(param, mid)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_line() +
    geom_line(data = num.captured, aes(samp, num_capt), colour = "dodgerblue4", lwd = 1) +
    xlab("Time") + ylab("Number Captured") 
}

## Total Bd load in the population. I think there is a non-linear transormation problem here...
 ## Not that it matters too much I guess because this will get updated as soon as there is a model
  ## for Bd load over time in the model
pop_load <- expdat %>% group_by(samp) %>% mutate(scaled_load = log_bd_load * (1 - dead)) %>% 
    summarize(total_load = sum(scaled_load))

pred_coef        <- as.data.frame(stan.fit.summary[
  grep("total_bd", dimnames(stan.fit.summary)[[1]])
  , c(4, 6, 8)])
names(pred_coef) <- c("lwr", "mid", "upr")
pred_coef        %<>% mutate(param = rownames(.)) %>% 
  mutate(param = as.numeric(factor(param, levels = param)))

pred_coef %>% {
  ggplot(., aes(param, mid)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_line() +
    geom_line(data = pop_load, aes(samp, total_load), colour = "dodgerblue4", lwd = 1) +
    xlab("Time") + ylab("Total Bd load in the population") 
}

####
## Generated quantities: 
##  4) Simulated capture matrix
####

capture.sim <- array(data = 0, dim = c(samp, ind, stan.length))
death.sim   <- array(data = 0, dim = c(samp, ind, stan.length))

for (i in 1:samp) {
  ## Simulate if each individual died in a time step (using all posterior samples)
  death.sim[i,,]   <- apply(stan.fit.samples$phi[,,i], 1, FUN = function(x) rbinom(length(x), 1, 1 - x))
  capture.sim[i,,] <- apply(stan.fit.samples$p[,,i], 1, FUN = function(x) rbinom(length(x), 1, x))
  
  ## Can only be captured if alive
  capture.sim[i,,] <- capture.sim[i,,] * (1 - death.sim[i,,])
  
}

death.sim          <- reshape2::melt(death.sim)
names(death.sim)   <- c("obs", "ind", "samp", "dead")
capture.sim        <- reshape2::melt(capture.sim)
names(capture.sim) <- c("obs", "ind", "samp", "captured")

pop.sim <- left_join(death.sim, capture.sim)

## Clean up death and captures. This sim is a bit odd but will work because death at each time and observations
 ## at each time are estimated individually 
pop.sim %<>% 
  group_by(samp, ind) %>% 
  mutate(
    dead     = cumsum(dead)
  , dead     = ifelse(dead > 1, 1, dead)
  , captured = ifelse(dead == 1, 0, captured))

## Take a look at a single posterior for what a capture history would look like (just a random one of the stan.length samples)
rand_samp <- sample(seq(stan.length), 1)

pop.sim.test <- pop.sim %>% filter(samp == rand_samp) %>% mutate(ind = as.factor(ind))

 pop.sim.test %>% {
  ggplot(., aes(obs, ind, fill = as.factor(captured))) + 
    geom_tile(alpha = 0.8) +
    xlab("Sampling Event") + ylab("Individual") +
    scale_fill_manual(
        values = c("dodgerblue4", "firebrick4")
      , name   = "Detected?"
      , labels = c("No", "Yes")) +
    geom_line(data = pop.sim.test %>% filter(dead == 1)
      , aes(x = obs, y = ind, z = NULL)) +
    theme(
      axis.text.y = element_text(size = 8)
    , legend.text = element_text(size = 12)
    , legend.key.size = unit(.55, "cm")
    ) +
     scale_x_continuous(breaks = c(1, 3, 5, 7, 9)) +
    ggtitle("Lines show dead individuals")
}
