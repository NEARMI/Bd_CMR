###############################################################
## Simulate data for and fit an individual-CMR model in Stan ##
###############################################################

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
nsim <- 1
ind  <- 300
samp <- 10

## Bd parameters
bd_beta  <- c(1, 0.05)     ## Intercept and slope for mean response
bd_sigma <- 2              ## observation noise
bd_theta <- c(1, 1, 1)     ## random effect variance covariance

## Response of individuals to Bd load
bd_mort   <- c(decay = -0.3, offset = 2)     ## logistic response coefficients for mortality across log(bd_load)
bd_detect <- c(decay = 0.1, offset = 0)      ## logistic response coefficients for detection across log(bd_load)

####
## Simulate bd response and relationship between bd and survival and detection
####

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

## view the simulated data
expdat %>% {
  ggplot(., aes(samp, log_bd_load)) + 
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
bd_range <- with(expdat, seq(min(log_bd_load), max(log_bd_load), by = 1))

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

stan_data      <- list(
   y             = capture_matrix
 , n_occasions   = samp
 , n_occ_minus_1 = samp - 1
 , n_ind         = ind
 , first         = capture_range$first
 , last          = capture_range$final
 , X             = matrix(nrow = ind, ncol = samp, data = expdat$log_bd_load, byrow = T)[, -samp]
 , n_captured    = capture_total$total_capt
  )

stan.fit  <- stan(
  file    = "CMR_ind_all_no.stan"
, data    = stan_data
, chains  = 1
, iter    = 1500
, warmup  = 500
, thin    = 1
, control = list(adapt_delta = 0.90)
  )

shinystan::launch_shinystan(stan.fit)

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)

pred_coef <- as.data.frame(stan.fit.summary[1:4, c(4, 6, 8)])
names(pred_coef) <- c("lwr", "mid", "upr")
pred_coef %<>% mutate(param = rownames(.))

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

stan.pred <- apply(stan.fit.samples$beta_phi, 1
  , FUN = function(x) plogis(x[1] + x[2] * bd_probs$log_bd_load)) %>%
  reshape2::melt(.)
names(stan.pred) <- c("log_bd_load", "sample", "mortality")
stan.pred %<>% group_by(log_bd_load) %>% 
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

stan.pred <- apply(stan.fit.samples$beta_p, 1
  , FUN = function(x) plogis(x[1] + x[2] * bd_probs$log_bd_load)) %>%
  reshape2::melt(.)
names(stan.pred) <- c("log_bd_load", "sample", "detect")
stan.pred %<>% group_by(log_bd_load) %>% 
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

####
## Notes
####

# The data from a capture-recapture study can be written in terms of a 'u' by 'k' capture matrix Xobs, where 'u' is the total number of individuals ever captured

## To use the bd data in a mark recapture model need a model for how bd changes within individuals -- presumably this rate
 ## can be modeled as a random effect with some global mean and some individual variation in progression rate
## Once the model for the covariate is defined, a link function, usually the logit, 
 ## relates the covariate information to the capture or survival probabilities.

## That is, an individuals survival probability could be some intercept (random) + some covariate * modeled bd





