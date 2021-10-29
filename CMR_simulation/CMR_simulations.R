###############################################################
## Simulate data for and fit an individual-CMR model in Stan ##
###############################################################

## Script whose predecessor is individual_CMR_expanding.R -- converting that script into functions for easier simulation of multiple
 ## populations for the more complicated CMR model

####
## Notes as of OCT 29:
####

## -- Simulation code updated for long-form bd submodel. Working (recovers simulated parameter values) but code 
 ## needs quite a bit of cleaning (e.g., a number of things are still not very dynamic)

## Next need to explore:
 ## 1) For USGS explore sampling schemes for bd -- what is needed to recover parameter values?
 ## 2) For future use of all data -- work to decrease size of model and see what can be recovered with fewer 
  ## samples within year
 ## 3) Add more covariates to better soak up more variation / remove need for bd to try and capture so much
 ## 4) Expand detection model

####
## Packages and misc
####
needed_packages <- c("magrittr", "dplyr", "tidyr", "lme4", "ggplot2", "rstan")
lapply(needed_packages, require, character.only = TRUE)
source("../../ggplot_theme.R")
set.seed(10002)

####
## Parameters 
####
source("CMR_parameters.R")

####
## Functions for simulation
####
source("CMR_functions.R")

####
## Run the sim to create the data
####

for (pop_ind in 1:n_pop) {

expdat  <- bd.simulate(
  periods   = periods[pop_ind, ]
, times     = times[pop_ind, ]
, all_ind   = all_ind[pop_ind, ]
, bd_beta   = bd_beta[pop_ind, ]
, bd_sigma  = bd_sigma[pop_ind, ]
, bd_theta  = bd_theta[pop_ind, ]
, obs_noise = obs_noise[pop_ind, ]
, bd_noinf  = bd_noinf[pop_ind, ]
)

# expdat %>% filter(periods == 1) %>% {ggplot(., aes(times, log_bd_load)) + geom_line() + facet_wrap(~ind)}

one_pop <- bd.sampling(
  expdat    = expdat
, all_ind   = all_ind[pop_ind, ]
, new_ind   = new_ind[[pop_ind]]
, times     = times[pop_ind, ]
, periods   = periods[pop_ind, ]
, when_samp = when_samp[pop_ind, ]
, samp      = samp[[pop_ind]]
, inbetween = inbetween[[pop_ind]]
, between_season_duration = between_season_duration[pop_ind, ]
, bd_mort   = bd_mort[pop_ind, ]
, bd_detect = bd_detect[pop_ind, ]
, p_mort    = p_mort[pop_ind, ]
, pop_ind   = pop_ind
)

one_pop.long <- bd.stan_org(
  one_pop  = one_pop
, pop_ind  = pop_ind
, times    = times[pop_ind, ]
, periods  = periods[pop_ind, ]
, samp     = samp[[pop_ind]]
)

print(paste("Population", pop_ind, "Simulated", sep = " "))

### Put the pops together
if (pop_ind == 1) {
ind_occ_phi.all   <- one_pop.long$ind_occ_phi
ind_occ_p.all     <- one_pop.long$ind_occ_p
X_bd.m.all        <- one_pop.long$X_bd.m
capture_range.all <- one_pop.long$one_pop$capture_range
ind_occ_size.all  <- one_pop.long$one_pop$ind_occ_size
all_ind.all       <- one_pop.long$one_pop$all_ind
each_ind.all      <- one_pop.long$one_pop$all_ind
ind_occ.all       <- one_pop.long$one_pop$all_ind * sum(samp[[pop_ind]])
ind_occ_min1.all  <- one_pop.long$one_pop$all_ind * (sum(samp[[pop_ind]]) - 1)
pop_cov.bd.all    <- one_pop.long$pop_cov.bd
ind_in_pop.all    <- rep(pop_ind, length(unique(one_pop.long$ind_occ_p$ind)))
expdat.all        <- one_pop$expdat
} else  {
ind_occ_phi.all   <- rbind(ind_occ_phi.all, one_pop.long$ind_occ_phi)
ind_occ_p.all     <- rbind(ind_occ_p.all, one_pop.long$ind_occ_p)
X_bd.m.all        <- rbind(X_bd.m.all, one_pop.long$X_bd.m)
capture_range.all <- rbind(capture_range.all, one_pop.long$one_pop$capture_range)
ind_occ_size.all  <- c(ind_occ_size.all, one_pop.long$one_pop$ind_occ_size)
all_ind.all       <- all_ind.all + one_pop.long$one_pop$all_ind
each_ind.all      <- c(each_ind.all, one_pop.long$one_pop$all_ind)
ind_occ.all       <- ind_occ.all + one_pop.long$one_pop$all_ind * sum(samp[[pop_ind]])
ind_occ_min1.all  <- ind_occ_min1.all + one_pop.long$one_pop$all_ind * (sum(samp[[pop_ind]]) - 1)
pop_cov.bd.all    <- rbind(pop_cov.bd.all, one_pop.long$pop_cov.bd)
ind_in_pop.all    <- c(ind_in_pop.all, rep(pop_ind, length(unique(one_pop.long$ind_occ_p$ind))))
expdat.all        <- rbind(expdat.all, one_pop$expdat)
}

}

## clean up expdat.all
expdat.all %<>% ungroup() %>% 
  arrange(pop, ind, all_times) %>% 
  mutate(ind = interaction(ind, pop)) %>% 
  mutate(ind = factor(ind, levels = unique(ind))) %>% 
  mutate(ind = as.numeric(ind))
  
## And the last few pieces outside of the loop
ind_occ_min1_size.all <- ind_occ_size.all - 1

## Fix the individual numbers in X_bd.m.all
if (n_pop > 1) {
for (i in 2:n_pop) {
  X_bd.m.all[X_bd.m.all$pop == i, ]$ind <- X_bd.m.all[X_bd.m.all$pop == i, ]$ind + 
    max(X_bd.m.all[X_bd.m.all$pop == (i - 1), ]$ind)
}
}

## convert ind_pop interaction column to individuals
ind_occ_phi.all %<>% mutate(ind = as.numeric(ind))
ind_occ_p.all   %<>% mutate(ind = as.numeric(ind))

## Index vector for the first entry of phi and p that correspond to a new individual
phi_first_index <- (ind_occ_phi.all %>% mutate(index = seq(n())) %>% group_by(ind) %>% 
  summarize(first_index = min(index)))$first_index

p_first_index <- (ind_occ_p.all %>% mutate(index = seq(n())) %>% group_by(ind) %>% 
  summarize(first_index = min(index)))$first_index

## --------------- ##
## Latent bd is estimated over the whole time period and not just for the capture occasions,
 ## though bd on the capture occasions are used to determine detection and survival. Need to
  ## determine what entries of phi, and p correspond to the full time period bd. This is done here
temp_dat <- expdat.all %>% 
  rename(sampling_events_phi = all_times) %>%
  ungroup() %>%
  arrange(ind, periods, times) %>% 
  mutate(index = seq(n())) 

phi.bd.index <- (left_join(
  ind_occ_phi.all %>% dplyr::select(ind, sampling_events_phi)
, temp_dat     %>% dplyr::select(ind, sampling_events_phi, index)
  ))$index

ind_occ_phi.all %<>% mutate(phi_bd_index = phi.bd.index)

# ind_occ_phi.all %<>% left_join(., temp_dat %>% dplyr::select(index, periods) %>% rename(sampling_events_phi = index))

temp_dat <- expdat.all %>%
  rename(sampling_events_p = all_times) %>%
  ungroup() %>%
  arrange(ind, periods, times) %>% 
  mutate(index = seq(n())) 

p.bd.index <- (left_join(
  ind_occ_p.all %>% dplyr::select(ind, sampling_events_p)
, temp_dat     %>% dplyr::select(ind, sampling_events_p, index)
  ))$index

ind_occ_p.all %<>% mutate(p_bd_index = p.bd.index)

# ind_occ_p.all %<>% left_join(., temp_dat %>% dplyr::select(index, periods) %>% rename(sampling_events_p = index))

temp_dat <- expdat.all %>% 
  ungroup() %>%
  arrange(ind, periods, times) %>% 
  mutate(index = seq(n())) %>% 
  group_by(ind, periods) %>%
  summarize(
    first_index = min(index)
  , last_index  = max(index)
    )

bd_first_index <- temp_dat$first_index
bd_last_index  <- temp_dat$last_index

temp_dat <- expdat.all %>% 
  dplyr::select(-times) %>%
  rename(times = all_times) %>%
  ungroup() %>%
  arrange(ind, periods, times) %>% 
  mutate(index = seq(n()))

X.bd.index <- (left_join(
  X_bd.m.all   %>% dplyr::select(ind, times)
, temp_dat %>% dplyr::select(ind, times, index)
  ))$index

X_bd.m.all %<>% mutate(X_bd_index = X.bd.index)

## -- end of this confusing section -- ##

####
## Run the stan model
####

stan_data     <- list(
  
  ## dimensional indexes 
   n_pop           = n_pop
 , n_ind           = all_ind.all
 , ind_per_period  = sum(periods * each_ind.all)
  
 , ind_time        = sum(
   (expdat.all %>% group_by(pop) %>% summarize(total_times   = length(unique(times))))$total_times * 
   (expdat.all %>% group_by(pop) %>% summarize(total_periods = length(unique(periods))))$total_periods * 
   c(each_ind.all))
 , ind_occ         = sum(lapply(samp, sum) %>% unlist() * each_ind.all)
 , ind_occ_min1    = sum((lapply(samp, sum) %>% unlist() - 1) * each_ind.all)

  ## short vector indexes 
 , ind_occ_size      = rep(lapply(samp, sum) %>% unlist(), each_ind.all)
 , ind_occ_min1_size = rep(lapply(samp, sum) %>% unlist() - 1, each_ind.all)

 , p_first_index     = p_first_index
 , phi_first_index   = phi_first_index
  
  ## long vector indexes
 , ind_occ_rep       = ind_occ_p.all$ind
 , sampling_events_p = ind_occ_p.all$sampling_events_p
 , periods_occ       = ind_occ_p.all$periods_occ
 , pop_p             = ind_occ_p.all$pop
 , p_zeros           = ind_occ_p.all$p_zeros
 , p_bd_index        = ind_occ_p.all$p_bd_index
 , gamma_index       = (ind_occ_p.all %>% mutate(ind_per = interaction(ind, periods_occ)) %>% 
     mutate(ind_per = factor(ind_per, levels = unique(ind_per))) %>% 
     mutate(ind_per = as.numeric(ind_per)))$ind_per
  
 , ind_occ_min1_rep    = ind_occ_phi.all$ind
 , sampling_events_phi = ind_occ_phi.all$sampling_events_p
 , offseason           = ind_occ_phi.all$offseason
 , pop_phi             = ind_occ_phi.all$pop
 , phi_zeros           = ind_occ_phi.all$phi_zeros
 , phi_bd_index        = ind_occ_phi.all$phi_bd_index

 , ind_bd_rep          = expdat.all$ind
 , sampling_events_bd  = expdat.all$times
 , ind_in_pop          = expdat.all$pop
 , temp                = expdat.all$temp
  
  ## covariates
 , N_bd            = nrow(X_bd.m.all)
 , X_bd            = X_bd.m.all$bd  
 , x_bd_index      = X_bd.m.all$X_bd_index
 , bd_first_index  = bd_first_index
 , bd_last_index   = bd_last_index
 , time_gaps       = ind_occ_phi.all$time_gaps
  
  ## Capture data
 , N_y             = nrow(ind_occ_p.all)
 , y               = ind_occ_p.all$captures
  
 , first           = capture_range.all$first
 , last            = capture_range.all$final
  
  )

stan.fit  <- stan(
  file    = "../CMR_empirical_long.stan"
, data    = stan_data
, chains  = 1
, iter    = stan.iter
, warmup  = stan.burn
, thin    = stan.thin
, refresh = 10
, control = list(adapt_delta = 0.92, max_treedepth = 12)
  )

# stan.fit <- readRDS("stan.fit.mi.Rds")

####
## CMR Diagnostics
####

## NOTE: OCT 15: just a tiny bit moved from Individual_CRM_expanding.R
 ## need to more thoroughly clean this up to accommodate "long" form
  ## will also need quite a bit of cleanup when bd parameters start varying by population
stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)

####
## Recovery of simulated coefficients?
####

## Primary Bd effects
pred_coef        <- as.data.frame(stan.fit.summary[grep("beta_p", dimnames(stan.fit.summary)[[1]]), c(4, 6, 8)])
names(pred_coef) <- c("lwr", "mid", "upr")
pred_coef        %<>% mutate(param = rownames(.))

pred_coef %>% {
  ggplot(., aes(param, mid)) + 
    geom_point() +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
    geom_point(data = data.frame(
      param = pred_coef$param
    , mid   = c(
      rev(bd_mort[1, ])
    , rev(bd_detect[1, ])
      )
    ), colour = "firebrick3") +
    xlab("Parameter") + ylab("Estimate") +
    scale_x_discrete(labels = c("Detection intercept", "Detection slope", "Survival intercept", "Survival slope")) +
    theme(axis.text.x = element_text(size = 11))
}

## survival over Bd load
stan.pred        <- apply(stan.fit.samples$beta_phi, 1
  , FUN = function(x) plogis(x[1] + x[2] * one_pop$bd_probs$log_bd_load)) %>%
  reshape2::melt(.)
names(stan.pred) <- c("log_bd_load", "sample", "mortality")
stan.pred        %<>% mutate(log_bd_load = plyr::mapvalues(log_bd_load
  , from = unique(log_bd_load), to = one_pop$bd_probs$log_bd_load))
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
    geom_line(data = (one_pop$bd_probs %>% mutate(mort = mort ))
      , aes(log_bd_load, mort)
      , colour = "dodgerblue4", lwd = 2) +
    xlab("Log of Bd Load") + ylab("Predicted mortality probability")
}

## detection over Bd load
stan.pred        <- apply(stan.fit.samples$beta_p, 1
  , FUN = function(x) plogis(x[1] + x[2] * one_pop$bd_probs$log_bd_load)) %>%
  reshape2::melt(.)
names(stan.pred) <- c("log_bd_load", "sample", "detect")
stan.pred        %<>% mutate(log_bd_load = plyr::mapvalues(log_bd_load
  , from = unique(log_bd_load), to = one_pop$bd_probs$log_bd_load))
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
    geom_line(data = one_pop$bd_probs, aes(log_bd_load, detect)
      , colour = "dodgerblue4", lwd = 2) +
    xlab("Log of Bd Load") + ylab("Predicted detection probability")
}

####
## Individual random effect estimates
####

stan.ind_pred_eps <- stan.fit.samples$bd_delta_eps %>%
  reshape2::melt(.) %>% rename(ind = Var2, eps = value)
stan.ind_pred_var <- stan.fit.samples$bd_delta_sigma %>%
  reshape2::melt(.) %>% rename(sd = value) %>%
  left_join(., stan.ind_pred_eps) %>%
  mutate(eps = eps * sd) %>% 
  group_by(ind) %>%
  summarize(
    mid = quantile(eps, 0.50)
  , lwr = quantile(eps, 0.025)
  , upr = quantile(eps, 0.975)
  ) %>% ungroup()
#  %>% left_join(., expdat %>% group_by(ind) %>% summarize(total_capt = sum(bd_swabbed))) %>% 
#  left_join(., expdat %>% group_by(ind, periods) %>% summarize(total_detect = sum(detected)) %>% 
#      filter(total_detect > 0) %>% summarize(total_periods = n())) %>% 
#  mutate(CI_width = upr - lwr) %>% 
#  mutate(order_pred = seq(n()))

### Need to make this functioning given the new multi-population structure

stan.ind_pred_var %>% 
  arrange(mid) %>% 
  mutate(ind = factor(ind, levels = ind)) %>% {
  ggplot(., aes(ind, mid)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr)) +
    xlab("Individual") + 
    ylab("Random Effect Deviate") +
    geom_hline(yintercept = 0, linetype = "dashed", lwd = 1, colour = "firebrick3") +
    scale_colour_brewer(palette = "Dark2") +
    theme(axis.text.x = element_text(size = 8))
}

expdat.all %>% group_by(ind) %>% summarize(
  tot_bd = sum(log_bd_load)
) %>% arrange(tot_bd)

## the most extreme individual
most_extreme <- reshape2::melt(stan.fit.samples$X[
  , 37
  , ])
names(most_extreme) <- c("samp", "time", "value")
most_extreme %<>% 
  group_by(time) %>% 
summarize(
  lwr = quantile(value, 0.025)
, mid = quantile(value, 0.500)
, upr = quantile(value, 0.975)
) %>% mutate(
  period = rep(seq(3), each = 20)
, time   = rep(seq(20), 3))

most_extreme %>% {
  ggplot(., aes(time, mid)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
    geom_line() +
    facet_wrap(~period)
}

####
## Population random effect estimates
####

stan.ind_pred_eps <- stan.fit.samples$phi_delta_pop_eps %>%
  reshape2::melt(.) %>% rename(ind = Var2, eps = value)
stan.ind_pred_var <- stan.fit.samples$phi_delta_pop_sigma %>%
  reshape2::melt(.) %>% rename(sd = value) %>%
  left_join(., stan.ind_pred_eps) %>%
  mutate(eps = eps * sd) %>% 
  group_by(ind) %>%
  summarize(
    mid = quantile(eps, 0.50)
  , lwr = quantile(eps, 0.025)
  , upr = quantile(eps, 0.975)
  ) %>% ungroup()

### Need to make this functioning given the new multi-population structure

stan.ind_pred_var %>% 
  arrange(mid) %>% 
  mutate(ind = factor(ind, levels = ind)) %>% {
  ggplot(., aes(ind, mid)) +
    geom_errorbar(aes(ymin = lwr, ymax = upr)) +
    xlab("Population") + 
    ylab("Random Effect Deviate") +
    geom_hline(yintercept = 0
      , linetype = "dashed", lwd = 1, colour = "firebrick3") +
    scale_colour_brewer(palette = "Dark2") +
    theme(axis.text.x = element_text(size = 8))
}
