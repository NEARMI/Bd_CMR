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
## Notes as of SEP 28:
####

## 0A) Coming together, but I see that when n_occasions varies by period the model as currently 
 ## written is going to be hard to scale (especially since the number of samples and the number of years
  ## is going to vary by state and site). Will need to use a pretty tricky indexing scheme to keep the model
   ## from being overly cumbersome with different named phi containers etc.
 ## For a single site over a few years with a different number of sampling days per year a single model with
  ## differently named phi's etc. could work even if it is a little ugly. It won't work for a much bigger model
## !! -- By collapsing from array to matrices that are stuck end to end, can have different numbers of sampling periods
 ## over time! At least solves one problem. Still not sure what to do about every site having different numbers of individuals and
  ## periods etc.

## 0B) In-between year survival. 
 ## -- My first thought for this (and what I am trying today) is to collapse the matrices for periods together so that
  ## time is linear increasing across columns. For two periods:
  ## -- phi will become c(n_occasions, n_oc_min1) such that the phi from the end of the first season to the start of the
   ## second season is the probability of making it through the season given the status at the end of the season, the duration
    ## in between those measures (and a new covaraiate to be estimated that stands in for offseason)
  ## -- Hopefully then with just an index matrix of which times are which periods for the latent covaraite, all of
   ## this will just work. 
 ## ^^!!** I guess one concern is how to deal with immigration into the population, but maybe that is just a derived
  ## quantity in this formulation? (i.e., by looking at population size between periods?) -- That is, if an individuals
   ## first capture is in period 2, how do we infer whether that individual was there in the first period or not? (I guess
    ## we cant?)

## 1) Stan models for this script begin with:
 ## CMR_ind_pat_bd-p-phi_multi
## 2) To correspond to the real data, the model will need
 ## A) Additional site-level covariates that affect bd and affect phi and p directly
 ## B) Multiple primary periods and survival in-between primary periods
 ## C) Fewer secondary periods within primary periods 
 ## D) Multiple sites and sites nested in state, with random effects for each location
 ## E) Possibly a reduced model for just inf/not instead of a full latent model for bd load

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
nsim      <- 1                  ## number of simulations (1 to check model, could be > 1 for some sort of power analysis or something)
ind       <- 50                 ## number of individuals in the population being modeled
periods   <- 2                  ## number of primary periods (years in most cases)
times     <- 20                 ## number of time periods (in the real data probably weeks; e.g., May-Sep or so)
if (periods > 1) {
times <- rep(times, periods)    ## for now assume same number of periods per year, but this likely wont be the case in the real data
}
samp      <- 10                 ## number of sampling events occurring over 'times' (e.g., subset of 'times' weeks when sampling occurred)
if (periods > 1) {
samp <- rep(samp, periods)      ## for now assume same number of periods per year, but this likely wont be the case in the real data
}
when_samp <- "random"           ## random = sampling occurs on a random subset of possible days

## Two ways to simulated data
  sim_opt <- "lme4"        ## Use the built in capabilities of lme4
# sim_opt <- "manual"      ## Simple custom simulation ** SEP 27: lagging behind "lme4" -- hasn't caught up with all of the options
  
## Drop simulated individuals never caught (TRUE) or keep all simulated individuals (FALSE) (unrealistic)
only_caught        <- TRUE

## Use covaratiates until the last sampling time (TRUE) or one before (FALSE)? (because p and phi are hard to separate in the last time point)
use_all_timepoints <- FALSE

## Bd parameters 
 ## for lme4 simulation. 
  ## ** SEP 27: for now dropping sim option manual, can add it back in later if we want a more complicated manual bd process model
bd_beta   <- c(
  7   ## Intercept
, 2   ## Period effect
, 30  ## Linear effect of time on bd
, -50 ## Quadratic effect of time on bd
) 
bd_sigma  <- 2              ## observation noise
bd_theta  <- c(2)           ## random effect variance covariance 

## Response of individuals to Bd load
bd_mort   <- c(decay = -0.3, offset = 6)       ## logistic response coefficients for mortality across log(bd_load)
bd_detect <- c(decay = 0.1, offset = -0.5)     ## logistic response coefficients for detection across log(bd_load)

## Other parameters
p_mort    <- c(0.20) ## mortality probability in-between periods

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
####

## Simulate data using lme4 mixed model structure
expdat <- expand.grid(
  periods = seq(periods)
, times   = seq(times[1])     ## ** SEP 27: Update to allow variable times by period? Maybe not needed because this is just the simulate underlying "true" bd data?
, ind     = factor(seq(ind))
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

## On which days is sampling occurring
if (when_samp == "random") {
 sampling_days <- expand.grid(
  times   = seq(times[1]) ## ** SEP 27: Update to allow variable times by period? Maybe not needed because this is just the simulate underlying "true" bd data?
                           ## Will no longer continue to make this comment
, periods = seq(periods)
   ) %>% group_by(periods) %>% 
   mutate(
     sampling_days = ifelse(times %in% sample(times, samp), 1, 0)
       )
 
  # (sampling_days - lag(sampling_days, 1))[-1]
 time_gaps     <- sampling_days %>% filter(sampling_days == 1) %>%
   group_by(periods) %>% mutate(time_gaps = (times - lag(times, 1)))
 
 time_gaps  <- matrix(data = time_gaps$time_gaps
   , ncol = periods
   , nrow = samp[1]
   )
 
 time_gaps <- time_gaps[-1, ]
 
## collapse for the correct time gap.
  ## ** SEP 28: For now am just using an arbitrary gap between seasons. Need to clean this up 
 time_gaps <- c(time_gaps[, 1], 10, time_gaps[, 2])

## Now just sampling_days$sampling_days
 sampling_vec  <- matrix(data = sampling_days$sampling_days
   , ncol = periods
   , nrow = times[1]
   , dimnames = list(
     seq(times[1])
  ,  seq(periods)
   ))
 
 sampling_times  <- matrix(data = (sampling_days %>% filter(sampling_days == 1))$times
   , ncol = periods
   , nrow = samp[1]
  )
 
}

## For debugging, save the same simulated data set to compare alternative models
expdat.sim <- expdat

####
## For debugging to compare models and sampling schemes -- can loop this or function this later if
## we want something more formal
####
# samp          <- 10
# sampling_days <- unique(expdat$times) %>% sample(samp) %>% sort(.)
# time_gaps     <- (sampling_days - lag(sampling_days, 1))[-1]
# sampling_vec  <- data.frame(times = seq(times), sampling = 0) 
# sampling_vec[sampling_vec$times %in% sampling_days, ]$sampling <- 1
# expdat        <- expdat.sim  ## retrieve the simulated data set to compare alternative models

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

## Add to the simulated bd data for the individuals
off_season <- expand.grid(
     periods = seq(from = periods - 0.5, to = periods, by = 1)
   , times   = 1
   , ind     = seq(ind)
   , bd_load = 0
   , log_bd_load = 0
   , mort    = 1 - p_mort
   , detect  = 0
    )

expdat %<>% 
  left_join(., bd_probs) %>%
  rbind(off_season, .) %>%
  group_by(ind) %>%
  arrange(periods, ind, times) %>%
  mutate(cum_surv = cumprod(mort)) 

## more debug checks
# expdat.t <- expdat
# expdat %>% filter(ind == 1) %>% as.data.frame()

####
## create the true state of the population from the simulated bd values
####

expdat %<>%
  mutate(dead = rbinom(n(), 1, 1 - mort)) %>%
  group_by(ind) %>% 
  mutate(
    dead = cumsum(dead)
  , dead = ifelse(dead > 1, 1, dead)) 

## On sampling days check for detection
expdat %<>% left_join(., sampling_days) %>%
  mutate(detected     = ifelse(dead == 0 & sampling_days == 1, rbinom(n(), 1, detect), 0))

####
## Create the sampling data from the simulated population 
####

## Drop all individuals that were never caught and update parameters
if (only_caught) {
never_detected <- expdat %>% 
  group_by(ind) %>% 
  summarize(total_detection = sum(detected)) %>% 
  filter(total_detection == 0)
new_ind <- seq(1, ind - nrow(never_detected))

## store original simulated population for estimate diagnostics
expdat.not_dropped <- expdat    
ind.not_dropped    <- ind

expdat %<>% dplyr::filter(ind %notin% never_detected$ind) %>%
  droplevels() %>%
  mutate(ind = as.numeric(ind), ind = as.factor(ind))

## total number of individuals ever captured
ind <- length(unique(expdat$ind))
}

## Create the capture array in the correct structure for stan
capture_matrix <- array(
   dim = c(
       ind
     , samp[1]
     , periods
     )
  , data = 0
  )

all_ind <- unique(expdat$ind)

expdat.capt <- expdat %>% filter(sampling_days == 1)

## Filling in an array can be confusing, so do it manually for clarity
for (i in seq_along(all_ind)) {
  capture_matrix[i,,] <- matrix(data = (expdat.capt %>% filter(ind == all_ind[i]))$detected, nrow = samp[1], ncol = periods)
}

## Collapse time periods into one side by side matrix and add the off-season as one period in-between
 ## currently non-dynamic. 
capture_matrix <- cbind(
  capture_matrix[,,1]
, matrix(data = 0, nrow = ind, ncol = 1)    ## off season
, capture_matrix[,,2]
  )

## The between-season parameter isn't needed for anything else, so drop this now
expdat %<>% filter(periods != 1.5)

capture_range <- expdat %>% 
  group_by(ind, periods) %>%
  filter(sampling_days == 1) %>%  
  summarize(
    first = min(which(detected == 1))
  , final = max(which(detected == 1))
    ) %>%
  mutate(
    first = ifelse(is.infinite(first), 0, first)
  , final = ifelse(is.infinite(final), 0, final)
    )

capture_range.first <- with(capture_range, matrix(
  data = first
, nrow = length(unique(ind))
, ncol = length(unique(periods))
, byrow = T)
  )

capture_range.final <- with(capture_range, matrix(
  data = final
, nrow = length(unique(ind))
, ncol = length(unique(periods))
, byrow = T)
  )

capture_total <- expdat %>% 
  group_by(times, periods) %>% 
  filter(sampling_days == 1) %>%
  summarize(total_capt = sum(detected)) %>%
  arrange(periods, times)

capture_total <- with(capture_total, matrix(
  data = total_capt
, nrow = samp[1]
, ncol = length(unique(periods))
  )
)
  
####
## Set up the bd sampling data
####

## full "true" bd values to be subset into "observed" bd values below
measured_bd <- array(
   dim = c(
       ind
     , samp[1]
     , periods
     )
  , data = 0
  )

## Filling in an array can be confusing, so do it manually for clarity
for (i in seq_along(all_ind)) {
    temp_data <- (expdat.capt %>% filter(ind == all_ind[i]))$log_bd_load
    ## Add the observation noise
    temp_data <- rnorm(length(temp_data), temp_data, obs_noise)
    measured_bd[i,,] <- matrix(data = temp_data, nrow = samp[1], ncol = periods)
}

if (bd_swabs == "PAT") {
  
## Find the captures that included bd swabs
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
  
## Create a matrix of these samples
bd.measured <- array(
   dim = c(
       ind
     , samp[1]
     , periods
     )
  , data = 0
  )

expdat.swabbed <- expdat %>% filter(sampling_days == 1)

## Filling in an array can be confusing, so do it manually for clarity
for (i in seq_along(all_ind)) {
    bd.measured[i,,] <- matrix(data = (expdat.swabbed %>% filter(ind == all_ind[i]))$bd_swabbed, nrow = samp[1], ncol = periods)
}

}

## and finally for this piece create an array of covariate data to designate period
periods.cov <- array(
   dim = c(
       ind
     , times[1]
     , periods
     )
  , data = expdat$periods
  )

## Double check to make sure this produces sensible capture data
if (ind <= 100) {
 expdat %>% filter(periods != 1.5) %>% {
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
     geom_line(data = expdat %>% filter(dead == 1, periods != 1.5), aes(x = times, y = ind, z = NULL), lwd= 0.5, alpha = 0.5) +
     geom_point(data = expdat %>% filter(bd_swabbed == 1, periods != 1.5)
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
## Run the model
####

if (bd_swabs == "PAT") {
  
## Possibly slightly confusing, but can keep all of the parameters needed for the most complicated model.
 ## The less complicated models that don't use these parameters will just ignore what they don't need
stan_data     <- list(
  ## dimensional params
   n_periods       = periods
 , n_ind           = ind                  
 , n_times         = times[1]         
 , all_samps       = sum(samp) + 1
 , all_samps_min1  = sum(samp - 1) + 1
 , n_occasions     = samp[1]
 , n_oc_min1       = samp[1] - 1
 , time            = seq(times[1])
 , sampling        = sampling_vec
 , sampling_events = sampling_times
 , time_gaps       = time_gaps
 , offseason       = c(rep(0, samp[1]), 1, rep(0, samp[2]))
  ## Capture data
 , y               = capture_matrix
 , first           = capture_range.first
 , last            = capture_range.final
 , n_captured      = capture_total
  ## Covariate associated parameters
 , X_bd            = measured_bd
 , X_measured      = bd.measured
 , periods         = periods.cov
  )

stan.fit  <- stan(
  file    = "CMR_ind_pat_bd-p-phi_multi.stan"
, data    = stan_data
, chains  = 1
, iter    = stan.iter
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.92, max_treedepth = 12)
  )

}

stan.fit <- stan.fit.4

# shinystan::launch_shinystan(stan.fit)

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)

## Within all of these diagnostics see *NOTE* for some of my summaries for what I have seen
 ## from my work on this so far

####
## SEP 21 NOTE: some of the below *might* be broken with the newer model and will need updating
####

####
## Recovery of simulated coefficients?
##   *NOTE*: High success with keeping all individuals, only moderate success with dropping individuals
####

## Primary Bd effects
pred_coef        <- as.data.frame(stan.fit.summary[1:4, c(4, 6, 8)])
names(pred_coef) <- c("lwr", "mid", "upr")
pred_coef        %<>% mutate(param = rownames(.))

pred_coef %>% {
  ggplot(., aes(param, mid)) + geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
    geom_point(data = data.frame(
      param = pred_coef$param
    , mid   = c(rev(bd_mort), rev(bd_detect))
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
    geom_line(data = bd_probs, aes(log_bd_load, mort)
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

## Hard to directly compare coefficient estimates given the change in sccale, need to work on this
pred_coef %>% filter(param != "start_mean") %>% {
  ggplot(., aes(param, mid)) + geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
    xlab("Parameter") + ylab("Estimate") +
    theme(axis.text.x = element_text(size = 11))
}

####
## Individual random effect estimates
####

stan.ind_pred_eps <- stan.fit.samples$bd_delta_eps %>%
  reshape2::melt(.) %>% rename(ind = Var2, eps = value)
stan.ind_pred_var <- stan.fit.samples$bd_delta_sigma %>%
  reshape2::melt(.) %>% left_join(., stan.ind_pred_eps) %>%
  mutate(eps = eps * value) %>% group_by(ind) %>%
  summarize(
    mid = quantile(eps, 0.50)
  , lwr = quantile(eps, 0.025)
  , upr = quantile(eps, 0.975)
  ) %>% arrange(mid) %>%
  mutate(ind = factor(ind, levels = ind))

## Not doing a very good job of making these look different than 0
stan.ind_pred_var %>% {
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
 
## Do we at least get the extreme two individuals correct?
 ## Not terrible! Model is at least somewhat working!
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

if (!only_caught) {
  pop_alive <- expdat %>% group_by(samp) %>% summarize(pop_size = n() - sum(dead))
  pop_alive <- pop_alive %>% rbind(data.frame(samp = 0, pop_size = ind.not_dropped), .)
} else {
  pop_alive <- expdat.not_dropped %>% group_by(samp) %>% summarize(pop_size = n() - sum(dead)) 
  pop_alive <- pop_alive %>% rbind(data.frame(samp = 0, pop_size = ind.not_dropped), .)
}

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
