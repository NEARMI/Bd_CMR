###############################################################
## Simulate data for and fit an individual-CMR model in Stan ##
###############################################################

## Script whose predecessor is individual_CMR_expanding.R -- converting that script into functions for easier simulation of multiple
 ## populations for the more complicated CMR model

####
## Notes as of Nov 8:
####

## Most pressing current problem[s]:
 ## 1) Multi-pop model individuals not aligning properly (specifically with X_bd.m.all -- individual number mismatch)

## Recap:
 ## Collapsed model doing an ok job for simulated data but poor job with real data; but:
  ## 1) Think my prior was too wide given that I only have one predictor in the plogis (need 1.5)
   ## Hopefully this will help to resolve some of the posterior problems with spikes at 1
  ## 2) Explore addition of one or two other covariates to the real data
  ## 3) Try multi-population model in the simulation

## Eventually:
## 1) Explore what is able to be recaptured with multi-pop model with the simulated data. 
 ## Can redo the sampling_range script and try like 10 fits

## Also:
 ## 1) There is still a small issue with survival from the last time time point through the
  ## offseason because survival from the last point to the end of the season is ignored
   ## need to try modifying the data so that the transition retains both info 
    ## drop offseason duration but keep the full linear predictor I think will work best
  ## 2) Need to adjust the simulation to return the actual between season survival values
   ## for diagnostics

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
source("CMR_datasim.R")

####
## Clean up simulated data and build structure for stan model
####
collapse.mod <- FALSE
source("CMR_dataclean.R")

if (length(unique(expdat.all$ind)) < 300) {
expdat.all %>% {
  ggplot(., aes(times, ind, fill = as.factor(detected))) + 
    geom_tile(aes(alpha = sampling_days)) +
    geom_point(data = 
        expdat.all %>% 
        filter(bd_swabbed == 1)
      , aes(x = times, y = ind, z = NULL), lwd = 0.7) +
    geom_line(data = 
        expdat.all %>% 
        filter(dead == 1)
      , aes(
         group = ind
        , x = times, y = ind, z = NULL), lwd = 0.5, alpha = 0.5) +
    facet_grid(pop~periods) +
    xlab("Week of the year") +
    ylab("Individual") +
    scale_fill_manual(
        values = c("dodgerblue4", "firebrick4")
      , name   = "Detected?"
      , labels = c("No", "Yes")) +
    scale_x_continuous(breaks = seq(1, max(expdat.all$times), by = 2)) +
    guides(alpha = FALSE) +
    theme(
      axis.text.y = element_text(size = 6)
    , legend.text = element_text(size = 12)
    , legend.key.size = unit(.55, "cm")
    ) 
  }
}

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
 , gamma_index       = (ind_occ_p.all %>% 
     mutate(ind_per = interaction(ind, periods_occ)) %>% 
     mutate(ind_per = factor(ind_per, levels = unique(ind_per))) %>% 
     mutate(ind_per = as.numeric(ind_per)))$ind_per
  
 , ind_occ_min1_rep    = ind_occ_phi.all$ind
 , sampling_events_phi = ind_occ_phi.all$sampling_events_p
 , offseason           = ind_occ_phi.all$offseason
 , pop_phi             = ind_occ_phi.all$pop
 , phi_zeros           = ind_occ_phi.all$phi_zeros
 , phi_bd_index        = ind_occ_phi.all$phi_bd_index
 , X_stat_index        = (ind_occ_phi.all %>% 
     mutate(ind_per = interaction(ind, periods)) %>% 
     mutate(ind_per = factor(ind_per, levels = unique(ind_per))) %>% 
     mutate(ind_per = as.numeric(ind_per)))$ind_per

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
  file    = {
    if (n_pop == 1) {
  #  "../CMR_empirical_long.stan"
      "CMR_collapsed.stan"
    } else {
  #  "../CMR_empirical_pr_long.stan"
      "CMR_collapsed_pr.stan"
    }
    }
, data    = stan_data
, chains  = 1
, iter    = stan.iter
, warmup  = stan.burn
, thin    = stan.thin
, refresh = 10
, control = list(adapt_delta = 0.92, max_treedepth = 12)
  )

shinystan::launch_shinystan(stan.fit)

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)

####
## CMR Diagnostics
####

## --- Recovery of simulated coefficients? --- ##

## Primary Bd effects
pred_coef        <- as.data.frame(stan.fit.summary[grep("beta_p"
  , dimnames(stan.fit.summary)[[1]]), c(4, 6, 8)])
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

## --- Individual random effect estimates --- ##

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
) %>% arrange(desc(tot_bd))

## just take a peek at the most extreme individuals (not dynamic)
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

## --- Population random effect estimates --- ##

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

## mean survival per period
ind_occ_phi.all %<>% mutate(pred_phi = colMeans(stan.fit.samples$phi))
ind_occ_phi.all %<>% left_join(.
  , expdat.all %>% filter(sampling_days == 1) %>% rename(sampling_events_phi = all_times) %>%
    dplyr::select(ind, sampling_events_phi, cum_surv)
  )

ind_occ_phi.all %>% filter(offseason == 1)

ind_occ_phi.all %>% 
  group_by(ind, periods) %>% 
  mutate(real_surv = cum_surv / lag(cum_surv, 1))

ind_occ_phi.all %>% filter(pred_phi > 0) %>% 
  group_by(offseason) %>%
  summarize(
    mean_phi = mean(pred_phi)
  )

data.frame(vals = stan.fit.samples$phi[, 234]) %>% {
  ggplot(., aes(x = vals)) + geom_histogram(bins = 50, alpha = 0.5) +
    xlab("Between Season Survival") +
    ylab("Samples") +
    geom_vline(xintercept = 0.77, colour = "firebrick3")
}

hist(stan.fit.samples$phi[, 234], breaks = 100)

ind_occ_p.all %<>% 
     mutate(ind_per = interaction(ind, periods_occ)) %>% 
     mutate(ind_per = factor(ind_per, levels = unique(ind_per))) %>% 
     mutate(ind_per = as.numeric(ind_per))
