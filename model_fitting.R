################################################
## Try and fit a CMR model with the newt data ##
################################################

####
## About this script
####

## First, see and start with 'data_exploration.R' for data info
## Second, see 'individual_CRM_few_samples.R' for a CMR model with a bd submodel that fits reasonably well with simulated data
## Third, this script attempts to fit a very simple model to one site 

####
## Notes as of Sep 21:
####

## 1) The data in MA:A11 honestly seems good enough to be able to fit some sort of model
 ## -- **** However really need to resolve:
  ## A) the time gaps problem between observations
  ## B) a realistic bd process model

## 2) It sort of seems like more of a regression-style random effects model for bd change over time could be fruitful
 ## -- **** Going this route however:
  ## A) Does every individual need their latent bd estimated because that is really hard. How much work to put into estimating all of this?
  ## B) How many extra random effects are needed: do individuals need to vary in their survival and detection apart from their variation in bd?
  ## C) How to deal with the time gaps inbetween sampling
  ## D) Will we be able to fit other covariates?

####
## Packages and functions
####
needed_packages <- c("magrittr", "dplyr", "tidyr", "lme4", "ggplot2", "rstan")
lapply(needed_packages, require, character.only = TRUE)
source("../ggplot_theme.R")
set.seed(10002)
'%notin%' <- Negate('%in%')

Bd_Newts_AllSites   <- read.csv("Bd_Newts_AllSites_9.12.21.csv")
Bd_Newts_AllSites   %<>% mutate(Date = as.Date(Date))

## Some parameters
red_ind        <- FALSE    ## reduce the number of individuals for debugging purposes
all_dates      <- TRUE     ## TRUE  = expand the capture matrix to have an event on every day
                           ## FALSE = collapse sampling events to lose the time spacing between them

if (red_ind) {
red_total_capt <- 4        ## minimum number of recaptures to keep an individual in an analysis
}

A11 <- Bd_Newts_AllSites %>% 
  filter(
    Site == "A11"
  , year == 2020
  )

if (!all_dates) {
capt_history <- expand.grid(
  Date = unique(A11$Date)
, Mark = unique(A11$Mark)
)
} else {
capt_history <- expand.grid(
  Date = seq(from = min(A11$Date), to = max(A11$Date), by = 1)
, Mark = unique(A11$Mark)
) 
}

capt_history %<>% 
  left_join(.
    , {A11 %>% dplyr::select(Date, Mark, julian, copies.swab)}
    ) %>% rename(
      captured = julian
    , bd_load  = copies.swab) %>% 
  mutate(
    captured = ifelse(is.na(captured), 0, 1)
  , swabbed  = ifelse(is.na(bd_load), 0, 1)) %>%
  mutate(log_bd_load = log(bd_load + 1)) %>%                          ### eeek!
  mutate(log_bd_load = ifelse(is.na(log_bd_load), 0, log_bd_load))    ### eeek X2!!

capt_history %<>% arrange(Date, Mark)

## Try a small subset to see if the model will run (regardless of whether it will fit or not)
if (red_ind) {
  capt_history.well_meas <- capt_history %>% 
  group_by(Mark) %>%
  summarize(total_capt = sum(captured)) %>% 
  filter(total_capt > red_total_capt) %>% 
  droplevels()
  
  capt_history %<>% filter(Mark %in% capt_history.well_meas$Mark)
} 

## Quick double check of these individuals' bd load. Honestly doesn't look too bad
capt_history %>% filter(swabbed == 1) %>%  {
  ggplot(., aes(Date, log_bd_load)) + 
    geom_line(aes(group = Mark)) +
    xlab("Date") +
    ylab("Bd Load") 
}

## Why collapsing sampling events doesnt work -- cant model a continuous curve appropriately when things
 ## get squashed without taking into consideration the length of time between the gaps
capt_history %>% filter(swabbed == 1) %>% 
  mutate(occ = as.numeric(as.factor(Date))) %>% {
  ggplot(., aes(occ, log_bd_load)) + 
    geom_line(aes(group = Mark)) +
    xlab("Date") +
    ylab("Bd Load")
}

## capture matrix
capt_history.matrix <- matrix(
  nrow = length(unique(capt_history$Mark))
, ncol = length(unique(capt_history$Date))
, data = capt_history$captured
, byrow = F)

## individuals' measured bd 
capt_history.bd_load <- matrix(
  nrow = length(unique(capt_history$Mark))
, ncol = length(unique(capt_history$Date))
, data = capt_history$log_bd_load
, byrow = F)

## time points where bd was measured
capt_history.bd_measured <- matrix(
  nrow = length(unique(capt_history$Mark))
, ncol = length(unique(capt_history$Date))
, data = capt_history$swabbed
, byrow = F)

## some quick visuals and checks
par(mfrow = c(2, 1))
image(capt_history.matrix)
image(capt_history.bd_measured)

if (length(
  which(
  (capt_history.matrix[2, ] == capt_history[capt_history$Mark == unique(capt_history$Mark)[2], ]$captured) == FALSE
)
) > 0) {
  print("ERROR: capture matrix not filled in correctly"); break
}

if (length(
  which(
  (capt_history.bd_measured[2, ] == capt_history[capt_history$Mark == unique(capt_history$Mark)[2], ]$swabbed) == FALSE
)
) > 0) {
  print("ERROR: bd matrix not filled in correctly"); break
}

capture_range  <- capt_history %>% group_by(Mark) %>% 
  summarize(
    first = min(which(captured == 1))
  , final = max(which(captured == 1))) %>% 
  dplyr::select(first, final) %>% 
  mutate(
    first = ifelse(is.infinite(first), 0, first)
  , final = ifelse(is.infinite(final), 0, final)
    )

capture_total <- capt_history %>% group_by(Date) %>% summarize(total_capt = sum(captured))

capt_history %>% mutate(event = as.factor(Date)) %>% {
  ggplot(., aes(Date, Mark, fill = as.factor(captured))) + 
    geom_tile(alpha = 0.8) +
    geom_point(data = capt_history %>% mutate(event = as.factor(Date)) %>% 
        filter(swabbed == 1), aes(x = Date, y = Mark, z = NULL), lwd = 0.7) +
    xlab("Sampling Event") + ylab("Individual") +
    scale_fill_manual(
        values = c("dodgerblue4", "firebrick4")
      , name   = "Detected?"
      , labels = c("No", "Yes")) +
    theme(
      axis.text.y = element_text(size = 6)
    , axis.text.x = element_text(size = 8, angle = 300)
    , legend.text = element_text(size = 12)
    , legend.key.size = unit(.55, "cm")
    ) 
}

####
## Run the stan model
####

stan.iter     <- 1500
stan.burn     <- 500
stan.thin     <- 1
stan.length   <- (stan.iter - stan.burn) / stan.thin

stan_data     <- list(
  ## bookkeeping params
   n_ind         = length(unique(capt_history$Mark))
 , n_occasions   = length(unique(capt_history$Date))
 , samp_events   = seq(length(unique(capt_history$Date)))  ## Need to update for uneven spacing
  ## Capture data
 , y             = capt_history.matrix
 , first         = capture_range$first
 , last          = capture_range$final
 , n_captured    = capture_total$total_capt
  ## Covariate associated parameters
 , X_bd          = capt_history.bd_load
 , X_measured    = capt_history.bd_measured
  )

stan.fit  <- stan(
  file    = "CMR_mine/CMR_ind_pat_bd_empirical_red2.stan"
, data    = stan_data
, chains  = 1
, iter    = stan.iter
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.92, max_treedepth = 12)
  )

shinystan::launch_shinystan(stan.fit)

saveRDS(stan.fit, "stan.fit.empirical.Rds")

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)

bd_levels <- log(c(seq(1, 10, by = 0.5) %o% 10^(0:5)))

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
    xlab("Parameter") + ylab("Estimate") +
    geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.5) +
    scale_x_discrete(labels = c("Detection intercept", "Detection slope", "Survival intercept", "Survival slope")) +
    theme(axis.text.x = element_text(size = 11))
}

## survival over Bd load
stan.pred        <- apply(stan.fit.samples$beta_phi, 1
  , FUN = function(x) plogis(x[1] + x[2] * bd_levels)) %>%
  reshape2::melt(.)
names(stan.pred) <- c("log_bd_load", "sample", "mortality")
stan.pred        %<>% mutate(log_bd_load = plyr::mapvalues(log_bd_load
  , from = unique(log_bd_load), to = bd_levels))
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
    xlab("Log of Bd Load") + ylab("Predicted mortality probability")
}

## detection over Bd load
stan.pred        <- apply(stan.fit.samples$beta_p, 1
  , FUN = function(x) plogis(x[1] + x[2] * bd_levels)) %>%
  reshape2::melt(.)
names(stan.pred) <- c("log_bd_load", "sample", "detect")
stan.pred        %<>% mutate(log_bd_load = plyr::mapvalues(log_bd_load
  , from = unique(log_bd_load), to = bd_levels))
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
## What does predicted overall bd look like?
####
stan.pred <- matrix(nrow = length(unique(capt_history$Date)), ncol = 1000, data = 0)
samp_occ  <- seq(length(unique(capt_history$Date)))

for (i in 1:nrow(stan.pred)) {
stan.pred[i, ] <- stan.fit.samples$beta_bd[, 1] + 
  stan.fit.samples$beta_bd[, 2] * samp_occ[i] +
  stan.fit.samples$beta_bd[, 3] * samp_occ[i]^2
}

stan.pred <- reshape2::melt(stan.pred)
names(stan.pred) <- c("occ", "iter", "value")

stan.pred %<>% 
  group_by(occ) %>%
  summarize(
    lwr = quantile(value, 0.025)
  , mid = quantile(value, 0.50)
  , upr = quantile(value, 0.975)
  )

ggplot(stan.pred, aes(occ, mid)) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  geom_line() +
  geom_line(
    data = (capt_history %>% filter(swabbed == 1) %>% mutate(occ = as.numeric(as.factor(Date))))
  , aes(occ, log_bd_load, group = Mark)
  , colour = "red"
  ) + xlab("Sampling Occasion") +
  ylab("Bd Load")

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
    )
}

## Random effects in phi and p
stan.ind_pred_eps <- stan.fit.samples$eps_alpha_phi %>%
  reshape2::melt(.) %>% rename(ind = Var2, eps = value)
stan.ind_pred_var <- stan.fit.samples$sigma_alpha_phi %>%
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
    )
}
