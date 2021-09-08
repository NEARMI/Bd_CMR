###############################################################
## Simulate data for and fit an individual-CMR model in Stan ##
###############################################################

########
## September 8 Notes:
########

 ## -- General progress since last check in -- 
 ## Model updated for individual variation in changes in Bd over time. Model can recover some of the variation,
  ## but lots of the individual slope deviate estimates overlap 0
   ## !! But I think a lot of the problem stems from the fact that there is a log transform problem (data simulated
    ## on linear scale and model run on the log scale)

 ## -- Next most important items -- 
  ## A) Resolve how the model is handling individuals that were never captured (given that so far I have just been
   ## keeping these simulated individuals in the model but they won't exist in the real data)
    ## Specifically, these individuals mess with the "Constraints" section and the loop right below that
  ## B) Need to go back through the whole model and check dimensions (especially around the n_occ_minus1)
   ## I think I am throwing out data when I don't need to be, especially associated with the model for the 
    ## covariate. (Need to go back and check the idea that p and phi cant be resolved in the last time step -- though
     ## their product can be; I think I am currently just ignoring the last time step altogether)

 ## -- Secondary to Do -- 
  ## 1) Make sure linear and/or log scale are operating correctly. 
  ## 2) Model doesn't do a great job even at 100 individuals total and 80/100 individuals with modeled bd.
   ## Explore how this scales with different numbers of individuals, effect size etc.
  ## 3) Try and come up with a better way to index bd_drop.which in the model
  ## 4) Need to move bd load estimates at the population level into the generated quantities now that there is a
   ## bd model inside of the CMR model
  ## 5) Figure out why population is being under-estimated
  ## 6) break this script up into multiple for readability
  ## 7) Move to real data ?!
   ## For this it will likely be important to have individual-level random effects for capture and mortality 
    ## probability. These are 

####
## Packages and functions
####
needed_packages <- c("magrittr", "dplyr", "tidyr", "lme4", "ggplot2", "rstan")
lapply(needed_packages, require, character.only = TRUE)
source("../../ggplot_theme.R")

####
## Parameters
####

## "Design" parameters
nsim <- 1                  ## number of simulations (1 to check model, could be > 1 for some sort of power analysis or something)
ind  <- 100                ## number of individuals
samp <- 10                 ## number of sampling events

## Two ways to simulated data
# sim_opt <- "lme4"        ## Use the built in capabilities of lme4
  sim_opt <- "manual"      ## Simple custom simulation 

## Bd parameters for lme4 simulation
if (sim_opt == "lme4") {
bd_beta  <- c(1, 0.05)     ## Intercept and slope for mean response
bd_sigma <- 2              ## observation noise
bd_theta <- c(1, 1, 1)     ## random effect variance covariance

## Response of individuals to Bd load
bd_mort   <- c(decay = -0.3, offset = 2)     ## logistic response coefficients for mortality across log(bd_load)
bd_detect <- c(decay = 0.1, offset = 0)      ## logistic response coefficients for detection across log(bd_load)

} else {
  
bd_int   <- c(              ## Gamma distribution for variation among individuals in starting conditions
    shape = 10
  , scale = 50)              
bd_delta <- c(              ## Slope in Bd over time (Normal random)
  mean = 400  
, sd   = 170 
  )                         
bd_add   <- 30              ## Process noise in true underlying Bd (normal SD)
bd_obs   <- 20              ## Observation noise in observed Bd (normal SD)

## Response of individuals to Bd load
bd_mort   <- c(decay = -0.7, offset = 6)     ## logistic response coefficients for mortality across log(bd_load)
bd_detect <- c(decay = 0.5, offset = -3)     ## logistic response coefficients for detection across log(bd_load)
}

bd_all        <- FALSE          ## TRUE = assume all captured individuals have their Bd swabbed
bd_drop       <- 20             ## number of individuals we will assume didn't have their Bd measured
bd_drop.which <- sample(        ## which of the simulated individuals has their Bd dropped
  seq(ind), bd_drop)           

####
## Simulate bd response and relationship between bd and survival and detection
####

if (sim_opt == "lme4") {

## Simulate data using lme4 mixed model structure
expdat <- expand.grid(samp = seq(samp), ind = factor(seq(ind)))
expdat %<>% mutate(
  bd_load   = simulate(~samp + (1 + samp | ind)
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

} else {
  
## Simulate data manually to more easily perfectly match an easy model form
expdat <- expand.grid(
  samp    = seq(samp)
, ind     = factor(seq(ind))
, bd_load = 0
    )

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
  ggplot(., aes(samp, log_bd_load
    )) + 
    geom_line(aes(group = ind)) +
    scale_y_log10() + 
    scale_x_continuous(breaks = c(1, 5, 10)) +
    xlab("Sampling Event") + ylab("Bd Load") + {
      if (ind <= 50) {
        facet_wrap(~ind)
      }
    }
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

bd_probs %>% {
  ggplot(., aes(x = log_bd_load)) +
  geom_line(aes(y = mort), colour = "firebrick3") +
  geom_line(aes(y = detect), colour = "dodgerblue3") +
  ylab("Prediction") + xlab("log of Bd load")
}

## Add to the simulated bd data for the individuals
expdat %<>% 
  left_join(., bd_probs) %>%
  group_by(ind) %>%
  mutate(cum_surv = cumprod(mort))

####
## create the true state of the population from the simulated bd values
####

## First, for each individual simulate if/when they die, bit of a roundabout way of doing so...
expdat %<>%
  mutate(dead = rbinom(n(), 1, 1 - mort)) %>%
  group_by(ind) %>% 
  mutate(
    dead = cumsum(dead)
  , dead = ifelse(dead > 1, 1, dead)) %>%
  mutate(detected = ifelse(dead == 0, rbinom(n(), 1, detect), 0))

## check to see if this worked with a heatmap
expdat %>% {
  ggplot(., aes(samp, ind, fill = as.factor(detected))) + 
    geom_tile(alpha = 0.8) +
    scale_x_continuous(breaks = c(1, 5, 10)) +
    xlab("Sampling Event") + ylab("Individual") +
    scale_fill_manual(
        values = c("dodgerblue4", "firebrick4")
      , name   = "Detected?"
      , labels = c("No", "Yes")) +
    geom_line(data = expdat %>% filter(dead == 1), aes(x = samp, y = ind, z = NULL)) +
    theme(
      axis.text.y = element_text(size = 8)
    , legend.text = element_text(size = 12)
    , legend.key.size = unit(.55, "cm")
    ) +
    ggtitle("Lines show dead individuals")
}

## Also take a final rough look at bd vs captures (simulated previously)
 ## Messy, but reflects the simulated probability
expdat %>% filter(dead == 0) %>% {
  ggplot(., aes(log_bd_load, detected)) + geom_jitter(height = 0.05)
}

####
## Put the data into the sturcture needed for the stan model
####

## now that these data seem sensible, collapse them into a matrix
capture_matrix <- matrix(nrow = ind, ncol = samp, data = expdat$detected, byrow = T)
capture_range  <- expdat %>% group_by(ind) %>% 
  summarize(
    first = min(which(detected == 1))
  , final = max(which(detected == 1))) %>% 
  dplyr::select(first, final) %>% 
  mutate(
    first = ifelse(is.infinite(first), 0, first)
  , final = ifelse(is.infinite(final), 0, final)
    )

capture_total <- expdat %>% group_by(samp) %>% summarize(total_capt = sum(detected))

####
## Run the model
####

stan.iter     <- 1500
stan.burn     <- 500
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin

measured_bd   <- matrix(nrow = ind, ncol = samp, data = expdat$log_bd_load, byrow = T)[, -samp]

## If trying without complete Bd sampling
if (!bd_all) {
  measured_bd <- measured_bd[-bd_drop.which, ]
}

## Some of these are not needed in the simplest model, but it is ok (if a bit confusing)
 ## to just pass everything for the more complicated model
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
  file    = "CMR_ind_some_no_2.stan"
, data    = stan_data
, chains  = 1
, iter    = stan.iter
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.92, max_treedepth = 12)
  )

shinystan::launch_shinystan(stan.fit)

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)

####
## Recovery of simulated coefficients?
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
    ) +
    geom_vline(xintercept = bd_drop.which, lwd = 0.1) +
    ggtitle("thin dashed lines are the 20 left out of the model. It is at least nice to say the estimates are all different")
}
 
## Do we at least get the extreme two individuals correct?
 ## Not terrible! Model is at least somewhat working!
data.frame(
  simulated = order(bd_ind)
, real      = stan.ind_pred_var$ind
) %>% head(10)

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

pop_alive <- expdat %>% group_by(samp) %>% summarize(pop_size = n() - sum(dead))

pred_coef %>% {
  ggplot(., aes(param, mid)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_line() +
    geom_line(data = pop_alive, aes(samp, pop_size), colour = "dodgerblue4", lwd = 1) +
    xlab("Time") + ylab("Population Size") 
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
    xlab("Time") + ylab("Population Size") 
}

## Total Bd load in the population. I think there is a non-linear transormation problem here...
 ## Not that it matters too much I guess because this will get updated as soon as there is a model
  ## for Bd load over time in the model
pop_load <- expdat %>% group_by(samp) %>% mutate(scaled_load = log_bd_load * (1 - dead)) %>% 
    summarize(mean_load = mean(scaled_load))
pop_load <- pop_load$mean_load

pop_load <- sweep(stan.fit.samples$pop
  , 2
  , pop_load[-10]
  , "*")

real_load <- expdat %>% group_by(samp) %>% mutate(scaled_load = log_bd_load * (1 - dead)) %>% 
    summarize(total_load = sum(scaled_load)) 

pop_load        <- reshape2::melt(pop_load)
names(pop_load) <- c("samp", "obs", "load")
pop_load      %<>% group_by(obs) %>%
  summarize(
    mid = quantile(load, 0.50)
  , lwr = quantile(load, 0.025)
  , upr = quantile(load, 0.975)
  )

pop_load %>% {
  ggplot(., aes(obs, mid)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_line() +
    geom_line(data = real_load, aes(samp, total_load), colour = "dodgerblue4", lwd = 1) +
    xlab("Time") + ylab("Population Size") 
}

####
## Generated quantities: 
##  4) Simulated capture matrix
####

capture.sim <- array(data = 0, dim = c(samp - 1, ind, stan.length))
death.sim   <- array(data = 0, dim = c(samp - 1, ind, stan.length))

for (i in 1:(samp - 1)) {
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

####
## Notes
####

# The data from a capture-recapture study can be written in terms of a 'u' by 'k' capture matrix Xobs, where 'u' is the total number of individuals ever captured

## To use the bd data in a mark recapture model need a model for how bd changes within individuals -- presumably this rate
 ## can be modeled as a random effect with some global mean and some individual variation in progression rate
## Once the model for the covariate is defined, a link function, usually the logit, 
 ## relates the covariate information to the capture or survival probabilities.

## That is, an individuals survival probability could be some intercept (random) + some covariate * modeled bd





