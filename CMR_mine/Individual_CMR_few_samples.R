###############################################################
## Simulate data for and fit an individual-CMR model in Stan ##
###############################################################

########
## See Individual_CMR.R for more extensive notes
########

## This script modifies ^^ script and model for a more realistic bd sampling scheme 
 ## in some small subset of captures individuals are swabbed

 ## NOTE search ** for date-stamped notes throughout the script

####
## Notes as of SEP 24 (late day):
####

## 1) Decided to step back to using n_occasions - 1 for phi and n_occasions for p. Hypothetically with a well-estimated
 ## latent covariate the last p and the second to last phi are identifiable, but my worry is that with real data they 
  ## will not be separable. In a simulation model even when they are separable not estimating them has only a tiny (almost
   ## not detectable) effect on coefficient estimates of interest
## 2) Everything I said earlier today is a lie. I accidentally had the covaraite that controlled for the gap between sampling times
 ## to be forced to be positive, though an increase in time should _decrease_ survival probability. Oops. 
  ## Now with the corrected model everything is behaving as it should be. The model that corrects for time does a lot better
   ## And the average bd model does as well if not better than the other model
 
## ^^^ ReadMe.txt that describes the stan models has been updated to describe these various models

## A) "manual" simulation lagging behind other option. Can catch this up later, but many of the sim options wont work correctly with sim_opt == "manual" currently 
## B) Really needed -- move the simulation piece into a function
## C) To correspond to the real data, the model will need
 ## 1) Additional site-level covariates that affect bd and affect phi and p directly
 ## 2) Multiple primary periods and survival in-between primary periods
 ## 3) Fewer secondary periods within primary periods 
 ## 4) Multiple sites and sites nested in state, with random effects for each location
 ## 5) Possibly a reduced model for just inf/not instead of a full latent model for bd load
  
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
ind       <- 60                ## number of individuals _in the population being modeled_
times     <- 20                 ## number of time periods (weeks or months -- days are probably too fine)
samp      <- 10                 ## number of sampling events occurring over times
when_samp <- "random"           ## random = sampling occurs on a random subset of possible days

## Two ways to simulated data
  sim_opt <- "lme4"        ## Use the built in capabilities of lme4
# sim_opt <- "manual"      ## Simple custom simulation 
  
## Drop simulated individuals never caught (TRUE) or keep all simulated individuals (FALSE) (unrealistic)
only_caught        <- TRUE

## Use covaratiates until the last sampling time or one before? (Seems keeping all time points doesn't lead to any worse mixing)
use_all_timepoints <- TRUE

## Bd parameters 
 ## for lme4 simulation
if (sim_opt == "lme4") {
bd_beta   <- c(10, 30, -50) ## Intercept and slope for mean response
                            ##  ** SEP 23: This is pretty un-dynamic and probably needs to be updated.
                             ## for now, based on the patterns in the real data estimate bd based on some quadratic function of time
bd_sigma  <- 2              ## observation noise
bd_theta  <- c(2)           ## random effect variance covariance

## Response of individuals to Bd load
bd_mort   <- c(decay = -0.45, offset = 6)     ## logistic response coefficients for mortality across log(bd_load)
bd_detect <- c(decay = 0.1, offset = -0.5)     ## logistic response coefficients for detection across log(bd_load)

## Or the "manual" bd process simulation
} else {
bd_int   <- c(              ## Gamma distribution for variation among individuals in starting conditions
    shape = 10
  , scale = 50)              
bd_delta <- c(              ## Slope in Bd over time (Normal random)
  mean = 400  
, sd   = 145)                         
bd_add   <- 30              ## Process noise in true underlying Bd (normal SD)
bd_obs   <- 20              ## Observation noise in observed Bd (normal SD)

## Response of individuals to Bd load
bd_mort   <- c(decay = -0.45, offset = 6)     ## logistic response coefficients for mortality across log(bd_load)
bd_detect <- c(decay = 0.7, offset = -4)     ## logistic response coefficients for detection across log(bd_load)
}

## bd sampling sampling scheme
bd_swabs  <- "PAT"  ## ALL = assume all captured individuals have their Bd swabbed
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

if (sim_opt == "lme4") {

## Simulate data using lme4 mixed model structure
expdat <- expand.grid(times = seq(times), ind = factor(seq(ind)))
expdat %<>% mutate(
  bd_load   = simulate(~poly(times, 2) + (1 | ind)
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
 sampling_days <- unique(expdat$times) %>% sample(samp) %>% sort(.)
 time_gaps     <- (sampling_days - lag(sampling_days, 1))[-1]
 sampling_vec  <- data.frame(times = seq(times), sampling = 0) 
 sampling_vec[sampling_vec$times %in% sampling_days, ]$sampling <- 1
}

## For debugging, save the same simulated data set to compare alternative models
expdat.sim <- expdat

} else {
  
## Simulate data manually to more easily perfectly match an easy model form
 ## ** SEP 23: Don't think I actually need the sampling bit here when I am simulating bd -- update!
expdat <- expand.grid(
  times    = seq(times)
, ind      = factor(seq(ind))
, sampling = 0
, bd_load  = 0
    )

## On which days is sampling occurring
if (when_samp == "random") {
 sampling_days <- unique(expdat$times) %>% sample(samp)
}

expdat[expdat$times %in% sampling_days, ]$sampling <- 1

all_ind <- unique(expdat$ind)

## individual-specific bd slopes (i.e., a conditional mode of the random effect for "individual")
bd_ind  <- rnorm(length(all_ind), bd_delta["mean"], bd_delta["sd"])

## Simulate true underlying Bd and observed Bd
 ## Probably could do this with some slick dplyr and apply, but w/e loops are easy to read and this
  ## could become a lot more complicated depending on the model we are simulating for
for (i in seq_along(all_ind)) {
  
  expdat.temp            <- expdat %>% filter(ind == all_ind[i])
  expdat.temp$bd_load[1] <- rgamma(n = 1, shape = bd_int["shape"], scale = bd_int["scale"])  ## starting Bd is variable among individuals

for (j in 2:nrow(expdat.temp)) {

  ## Simulated such that every individual has some true slope but with some additive noise (to simulate
   ## even more noise to try and better approximate real data (e.g., unexplained variation))
  expdat.temp$bd_load[j] <- rnorm(1, expdat.temp$bd_load[j - 1] + 
      bd_ind[i]      
    , bd_add)           
    
}
  expdat[expdat$ind == all_ind[i], ]$bd_load <- expdat.temp$bd_load
}

expdat %<>% mutate(
  bd_load_obs = rnorm(n(), bd_load, bd_obs)
, log_bd_load = round(log(bd_load_obs), 1)
)

}

## view the simulated data
expdat %>% {
  ggplot(., aes(times, bd_load
    )) + 
    geom_line(aes(group = ind)) +
    scale_y_log10() + 
#   scale_x_continuous(breaks = c(1, 5, 10)) +
    xlab("Sampling Event") + ylab("Bd Load") + {
      if (ind <= 50) {
        facet_wrap(~ind)
      }
    }
} 

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
expdat %<>% 
  left_join(., bd_probs) %>%
  group_by(ind) %>%
  mutate(cum_surv = cumprod(mort))

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
expdat %<>% left_join(.
  , {
    expdat %>% 
      filter(times %in% sampling_days) %>%
      mutate(
        detected     = ifelse(dead == 0, rbinom(n(), 1, detect), 0)
      , sampling_day = 1)
  }
  ) %>% mutate(
    sampling_day = ifelse(is.na(sampling_day), 0, sampling_day)
  , detected     = ifelse(is.na(detected), 0, detected)
  )

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


## another debugging check
# expdat.t <- expdat
# expdat <- expdat.t

## Create the capture matrix in the correct structure for stan
 ## ** SEP 23: Still getting to the bottom of whethter this matrix can have "samp" columns (and the model
  ## can have an index vector for which time points these samps correspond to) or if it needs to have "times" columns
capture_matrix <- matrix(
    nrow = ind
  , ncol = samp
  , data = (expdat %>% filter(sampling_day == 1))$detected
  , byrow = T)

capture_range  <- expdat %>% group_by(ind) %>%
  filter(sampling_day == 1) %>%  
  summarize(
    first = min(which(detected == 1))
  , final = max(which(detected == 1))) %>% 
  dplyr::select(first, final) %>% 
  mutate(
    first = ifelse(is.infinite(first), 0, first)
  , final = ifelse(is.infinite(final), 0, final)
    )

capture_total <- expdat %>% group_by(times) %>% summarize(total_capt = sum(detected)) %>%
  mutate(sampling_day = ifelse(times %in% sampling_days, 1, 0)) %>% filter(sampling_day == 1)

####
## Set up the bd sampling data
####

## full "true" bd values to be subset into "observed" bd values below
if (use_all_timepoints) {
measured_bd   <- matrix(
  nrow = ind
, ncol = samp
, data = (expdat %>% filter(sampling_day == 1))$log_bd_load
, byrow = T)
} else {
measured_bd   <- matrix(nrow = ind, ncol = times, data = expdat$log_bd_load, byrow = T)[, -samp] 
}

## SEP 23: OLD and likely broken
       if (bd_swabs == "IND") {
       
## Which of the caught individuals doesnt get measured for bd?
bd_drop.which <- sample(seq(ind), bd_drop)    

## update bd_ind (extract the captured individuals with swabbed bd)
bd_ind        <- bd_ind[-as.numeric(never_detected$ind)]

## drop these individuals
measured_bd   <- measured_bd[-bd_drop.which, ]

} else if (bd_swabs == "PAT") {
  
## Find the captures that included bd swabs
expdat %<>% left_join(.
  , {
    expdat %>% 
      filter(sampling_day == 1, detected == 1) %>% 
      ungroup() %>% 
      mutate(bd_swabbed = rbinom(n(), 1, bd_perc))
  }) 
  
expdat %<>% mutate(
  bd_swabbed = ifelse(is.na(bd_swabbed), 0, bd_swabbed)
)
  
## Create a matrix of these samples
bd.measured   <- matrix(
  nrow = ind
, ncol = samp
, data = (expdat %>% filter(sampling_day == 1))$bd_swabbed
, byrow = T)

## Set up a full matrix of the days in which each individual was searched for after the first day in which
 ## it was captured?

}

## Double check to make sure this produces sensible capture data
if (ind <= 100) {
 expdat %>% {
   ggplot(., aes(times, ind, fill = as.factor(detected))) + 
     geom_tile(aes(alpha = sampling_day)) +
     scale_x_continuous(breaks = c(1, 5, 10)) +
     xlab("Time") + 
     ylab("Individual") +
     scale_fill_manual(
         values = c("dodgerblue3", "firebrick3")
       , name   = "Detected?"
       , labels = c("No", "Yes")) +
     guides(alpha = FALSE) +
     geom_line(data = expdat %>% filter(dead == 1), aes(x = samp, y = ind, z = NULL)) +
     geom_point(data = expdat %>% filter(bd_swabbed == 1)
       , aes(x = times, y = ind, z = NULL), lwd = 0.5) +
     theme(
       axis.text.y = element_text(size = 8)
     , legend.text = element_text(size = 12)
     , legend.key.size = unit(.55, "cm")
     ) +
     ggtitle("Lines show dead individuals")
 }
 }

## Finally, inject the observation noise
# par(mfrow = c(2, 1))
# matplot(t(measured_bd)); matlines(t(measured_bd))
measured_bd <- t(apply(measured_bd, 1, FUN = function(x) x + rnorm(length(x), 0, obs_noise)))
# matplot(t(measured_bd)); matlines(t(measured_bd))

####
## Run the model
####

## ** SEP 23: In progress trying to get the model updated to:
 ## model bd as a continuous time process + 
 ## consider capture probability only on searching days

if (bd_swabs == "IND" | bd_swabs == "ALL") {

stan_data     <- list(
  ## bookkeeping params
   n_ind         = ind
 , n_ind_bd      = ind - bd_drop
 , n_occasions   = samp
 , n_occ_minus_1 = samp - 1
  ## Capture data
 , y             = capture_matrix
 , first         = capture_range$first
 , last          = capture_range$final
 , n_captured    = capture_total$total_capt
  ## Covariate associated parameters
 , X_bd          = measured_bd
 , X_which       = seq(ind)[-bd_drop.which]
  )

stan.fit  <- stan(
# file    = "CMR_ind_all_no.stan"
# file    = "CMR_ind_some_no_2.stan"
# file    = "CMR_ind_some_bd-p-phi_2.stan"
# file    = "CMR_ind_some_bd-p-phi_all.stan"
  file    = "CMR_ind_some_bd-p-phi_restrict.stan"
, data    = stan_data
, chains  = 1
, iter    = stan.iter
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.92, max_treedepth = 12)
  )

} else if (bd_swabs == "PAT") {
  
## Possibly slightly confusing, but can keep all of the parameters needed for the most complicated model.
 ## The less complicated models that don't use these parameters will just ignore what they don't need
stan_data     <- list(
  ## bookkeeping params
   n_ind           = ind
 , n_times         = times
 , n_occasions     = samp
 , n_oc_min1       = samp - 1
 , time            = seq(times)
 , sampling        = sampling_vec$sampling
 , sampling_events = sort(sampling_days)
 , no_sampling     = which(seq(times) %notin% sampling_days)
  ## Capture data
 , y               = capture_matrix
 , first           = capture_range$first
 , last            = capture_range$final
 , n_captured      = capture_total$total_capt
  ## Covariate associated parameters
 , X_bd            = measured_bd
 , X_measured      = bd.measured
 , time_gaps       = time_gaps
 , bd_after_gap    = c(sort(sampling_days)[-samp] + (time_gaps - 1), times)
  )

## bd_after_gap is confusing, so make sure it is doing what we want
data.frame(
  phi_num   = seq(samp)[-samp]
, time_gaps = time_gaps
, bd_from   = sort(sampling_days)[-samp]
, bd_to     = c(sort(sampling_days)[-samp] + (time_gaps - 1), times)[-samp]
)

stan.fit  <- stan(
#  file    = "CMR_ind_pat_bd-p-phi_no_timegap_covariate.stan"
   file    = "CMR_ind_pat_bd-p-phi_no_average_bd.stan"
#  file    = "CMR_ind_pat_bd-p-phi_no_average_bd_all_times.stan"
#  file    = "CMR_ind_pat_bd-p-phi_average_bd.stan"
#  file    = "CMR_ind_pat_bd-p-phi_cumulative_bd.stan"
, data    = stan_data
, chains  = 1
, iter    = stan.iter
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.92, max_treedepth = 12)
  )

}

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
