###############################################################
## Simulate data for and fit an individual-CMR model in Stan ##
###############################################################

########
## September 9 Notes:
########

 ## -- General progress since last check in -- 
 ## 1) Some parameters added for simulation options (dropping individuals, using n-1 or n time points etc.) --> interesting to compare coefficient estimates 
  ## from these various options
 ## 2) Simulation updated to drop individuals that were never caught. This leads to poorer estimates -- maybe there is something interesting
  ## going on here with bias and the ability to estimate true underlying parameters only from observed individuals 
 ## 3) log transform problem remains and needs to get fixed (data simulated on linear scale and model run on the log scale)
 ## 4) Interestingly using up through the last point still allows for a separation of detection from survival. I think because probability 
  ## and survival are a function of the latent covariate that is being modeled
 ## 5) Also interestingly forcing the p and phi for each individual prior to their first observation to 0 really doesn't change any of the 
  ## coefficient estimates. I find this weird, but I guess because we know they are zero the constraints should be there. This is maybe something
   ## to revisit when exploring real data

 ## -- Other thoughts --
 ## 1) Feeling decent about the general mechanics of a "simple" CMR model now. Feel I will be able to make some reasonable sense of more complicated
  ## CMR models now that I have a simpler one running and behaving reasonably well with simulated data
 ## 2) Possibly no reason at this point to continue with this simulation-based model and move to real data? There are a number of unresolved issues
  ## (see below) but those _probably_ only need to be resolved if this would turn into some sort of power analysis. In short, the model seems sensible
   ## enough at this point to move to working on a more complicated bd sub-model and start exploring some real data
 ## 3) I can see important steps forward to be:
  ## A) Trying to have multiple species or sites by having the beta parameters themselves be random or the progression of the disease be random

 ## -- Next most important items -- 
 ## 1) Specifically re: email chain from around 3pm on Thursday -- May want to take a step back from the CMR model and start exploring a more complete
  ## bd submodel (that is more extensive but can still be nested into a CMR model)

 ## -- Secondary to Do -- 
  ## 1) Make sure linear and/or log scale are operating correctly. 
  ## 2) Try and come up with a better way to index bd_drop.which in the model
  ## 3) break this script up into multiple for readability

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
nsim <- 1                  ## number of simulations (1 to check model, could be > 1 for some sort of power analysis or something)
ind  <- 200                ## number of individuals
samp <- 10                 ## number of sampling events

## Two ways to simulated data
# sim_opt <- "lme4"        ## Use the built in capabilities of lme4
  sim_opt <- "manual"      ## Simple custom simulation 
  
## Keep all simulated individuals or drop individuals never captured?
only_caught <- TRUE

## Use covaratiates until the last sampling time or one before? (As of Sep 9 still trying to get to the bottom of the impact of this choice)
use_all_timepoints <- TRUE

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
, sd   = 145)                         
bd_add   <- 30              ## Process noise in true underlying Bd (normal SD)
bd_obs   <- 20              ## Observation noise in observed Bd (normal SD)

## Response of individuals to Bd load
bd_mort   <- c(decay = -0.7, offset = 6)     ## logistic response coefficients for mortality across log(bd_load)
bd_detect <- c(decay = 0.7, offset = -4)     ## logistic response coefficients for detection across log(bd_load)
}

bd_all    <- FALSE          ## TRUE = assume all captured individuals have their Bd swabbed
bd_drop   <- 20             ## number of individuals we will assume didn't have their Bd measured

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
    scale_y_log10(breaks = seq(5,9)) + 
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

bd_probs %>% pivot_longer(cols = c(2, 3)) %>% {
  ggplot(., aes(log_bd_load, value)) +
  geom_line(aes(colour = name)) + 
  scale_colour_manual(
      name   = "Relationship"
    , values = c("firebrick3", "dodgerblue3")
    , labels = c("Mortality", "Detection")) + 
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

## First, for each individual simulate if/when they die, bit of a roundabout way of doing so...
expdat %<>%
  mutate(dead = rbinom(n(), 1, 1 - mort)) %>%
  group_by(ind) %>% 
  mutate(
    dead = cumsum(dead)
  , dead = ifelse(dead > 1, 1, dead)) %>%
  mutate(detected = ifelse(dead == 0, rbinom(n(), 1, detect), 0))

## Drop all individuals that were never caught and update parameters
if (!only_caught) {
never_detected <- expdat %>% group_by(ind) %>% summarize(total_detection = sum(detected)) %>% filter(total_detection == 0)
new_ind <- seq(1, ind - nrow(never_detected))

## store original simulated population for estimate diagnostics
expdat.not_dropped <- expdat    
ind.not_dropped    <- ind

expdat         %<>% dplyr::filter(ind %notin% never_detected$ind) %>%
  droplevels() %>%
  mutate(ind = as.numeric(ind), ind = as.factor(ind))

ind <- length(unique(expdat$ind))

## Which of the caught individuals has their Bd dropped from those individuals we actually caught
bd_drop.which <- sample(seq(ind), bd_drop)    

## update bd_ind
bd_ind <- bd_ind[-as.numeric(never_detected$ind)]

} else {
## which of the simulated individuals has their Bd dropped. Could be from individuals never caught, which yes, is a bit odd
bd_drop.which <- sample(seq(ind), bd_drop)   
}

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

if (use_all_timepoints) {
measured_bd   <- matrix(nrow = ind, ncol = samp, data = expdat$log_bd_load, byrow = T)
} else {
measured_bd   <- matrix(nrow = ind, ncol = samp, data = expdat$log_bd_load, byrow = T)[, -samp] 
}

## If trying without complete Bd sampling
if (!bd_all) {
  measured_bd <- measured_bd[-bd_drop.which, ]
}

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

## Explore full time course vs reduced time course
# stan.fit <- stan.fit.minusone
# stan.fit <- stan.fit.full
## And all individuals vs only caught individuals
# stan.fit <- stan.fit.caught.restrict              ## insure p and phi are 0 prior to first catch
# stan.fit <- stan.fit.caught                       ## do not make this assumption

shinystan::launch_shinystan(stan.fit)

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)

## Within all of these diagnostics see *NOTE* for some of my summaries for what I have seen
 ## from my work on this so far

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
