###############################################################
## Simulate data for and fit an individual-CMR model in Stan ##
###############################################################

## Script whose predecessor is individual_CMR_expanding.R -- converting that script into functions for easier simulation of multiple
 ## populations for the more complicated CMR model

####
## Notes as of OCT 19:
####

## Notes:
 ## -- Covariates by population can be recovered -ok- with well sampled populations.
   ## * Random effects model with 12 populations, rpois(20) indv per population taks about 1.2 h to run:
    ##   -- recovers mean bd response, survival, and detection quite well
    ##   -- recovers individual random effects no worse
    ##   -- questionable recovery of population-specific survival parameters
      ##     * HOWEVER, I am more confident that this will work better when population-specific parameters are added
    ##   ^^ Can check better sampled populations to see what it takes to recover these parameters
 ## -- Latent bd could maybe get more complicated, but worried about identifiability in poorly sampled populations
   ## * Model does reasonably fine recovering very low bd-loads for individuals that never get infected
   ## * Worried about individuals that only get infected some years -- will need to play with blocked random effects
   ## * Slope variation "works" but makes it hard to estimate both variation among individuals in intercept and slope. Maybe will be fine in real
    ##  data where there is so much variation 
   ## * To get variable timing of bd start and end I think what will be needed are parameters that directly adjust time for each location. 
    ##  Maybe times can be the same _length_ but just have different times? (i.e., by adding a time_adj parameter)
 ## -- Have a detailed rmd and html that shows the data structures. Will work on a script to parse the real data when I move to real data 

## Next Steps:
 ## 1) Construct the data parsing script using the newt data
 ## 2) Work on first fitting a single population and then a multi-population model with the newt data

## Some things to still work on in the long run:
 ## 1) Still doesn't allow number of periods to vary by population --- that shouldn't be too hard -- just tedious
 ## 2) Doesn't allow times to vary by population -- this will be harder because X_bd is fit with a matrix and will have to be melted

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
, new_ind   = new_ind[pop_ind, ]
, times     = times[pop_ind, ]
, periods   = periods[pop_ind, ]
, when_samp = when_samp[pop_ind, ]
, samp      = samp[pop_ind, ]
, inbetween = inbetween[pop_ind, ]
, between_season_duration = between_season_duration[pop_ind, ]
, bd_mort   = bd_mort[pop_ind, ]
, bd_detect = bd_detect[pop_ind, ]
, p_mort    = p_mort[pop_ind, ]
)

one_pop.long <- bd.stan_org(
  one_pop  = one_pop
, pop_ind  = pop_ind
, times    = times[pop_ind, ]
, periods  = periods[pop_ind, ]
, samp     = samp[pop_ind, ]
)

print(paste("Population", pop_ind, "Simulated", sep = " "))

### Put the pops together
if (pop_ind == 1) {
ind_occ_phi.all   <- one_pop.long$ind_occ_phi
ind_occ_p.all     <- one_pop.long$ind_occ_p
X_bd.m.all        <- one_pop.long$X_bd.m
capture_range.all <- one_pop.long$one_pop$capture_range
present.all       <- one_pop.long$one_pop$present
ind_occ_size.all  <- one_pop.long$one_pop$ind_occ_size
all_ind.all       <- one_pop.long$one_pop$all_ind
ind_occ.all       <- one_pop.long$one_pop$all_ind * rowSums(samp)[pop_ind]
ind_occ_min1.all  <- one_pop.long$one_pop$all_ind * (rowSums(samp) - 1)[pop_ind]
pop_cov.bd.all    <- one_pop.long$pop_cov.bd
ind_in_pop.all    <- rep(pop_ind, length(unique(one_pop.long$ind_occ_p$ind)))
} else  {
ind_occ_phi.all   <- rbind(ind_occ_phi.all, one_pop.long$ind_occ_phi)
ind_occ_p.all     <- rbind(ind_occ_p.all, one_pop.long$ind_occ_p)
X_bd.m.all        <- rbind(X_bd.m.all, one_pop.long$X_bd.m)
capture_range.all <- rbind(capture_range.all, one_pop.long$one_pop$capture_range)
present.all       <- rbind(present.all, one_pop.long$one_pop$present)
ind_occ_size.all  <- c(ind_occ_size.all, one_pop.long$one_pop$ind_occ_size)
all_ind.all       <- all_ind.all + one_pop.long$one_pop$all_ind
ind_occ.all       <- ind_occ.all + one_pop.long$one_pop$all_ind * rowSums(samp)[pop_ind]
ind_occ_min1.all  <- ind_occ_min1.all + one_pop.long$one_pop$all_ind * (rowSums(samp) - 1)[pop_ind]
pop_cov.bd.all    <- rbind(pop_cov.bd.all, one_pop.long$pop_cov.bd)
ind_in_pop.all    <- c(ind_in_pop.all, rep(pop_ind, length(unique(one_pop.long$ind_occ_p$ind))))
}

}

## And the last few pieces outside of the loop
ind_occ_min1_size.all <- ind_occ_size.all - 1

## Fix the individual numbers in X_bd.m.all
if (n_pop > 1) {
for (i in 2:n_pop) {
  X_bd.m.all[X_bd.m.all$pop == i, ]$ind <- X_bd.m.all[X_bd.m.all$pop == i, ]$ind + max(X_bd.m.all[X_bd.m.all$pop == (i - 1), ]$ind)
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

## Adjust population-specific convariates into the correct structure

## for bd
temp <- pop_cov.bd.all %>% pivot_wider(., values_from = temp, names_from = pop)
temp <- as.matrix(temp[, -1])

## for phi


## for p


####
## Finally, run the stan model
####

## NOTE: OCT 15: Will need to seriously think about what to do if we want times and periods to vary by population

stan_data     <- list(
  
  ## dimensional indexes 
   n_pop           = n_pop
 , n_periods       = periods[1, ] 
 , n_ind           = all_ind.all        
 , n_times         = times[1, ] * periods[1, ]
 , times_within    = times[1, ]
 , ind_occ         = ind_occ.all
 , ind_occ_min1    = ind_occ_min1.all
  
  ## short vector indexes 
 , time              = rep(seq(times[1, ]), periods[1, ])
 , time_per_period   = matrix(data = seq(times[1, ] * periods[1, ]), nrow = times, ncol = periods[1, ])
 , periods           = one_pop$periods_time
 , ind_occ_size      = ind_occ_size.all
 , ind_occ_min1_size = ind_occ_min1_size.all
 , ind_in_pop        = ind_in_pop.all

 , phi_first_index   = phi_first_index
 , p_first_index     = p_first_index
  
  ## long vector indexes
 , ind_occ_min1_rep    = ind_occ_phi.all$ind
 , sampling_events_phi = ind_occ_phi.all$sampling_events_phi
 , offseason           = ind_occ_phi.all$offseason
 , pop_phi             = ind_occ_phi.all$pop
 , phi_zeros           = ind_occ_phi.all$phi_zeros

 , ind_occ_rep       = ind_occ_p.all$ind
 , sampling_events_p = ind_occ_p.all$sampling_events_p
 , periods_occ       = ind_occ_p.all$periods_occ
 , pop_p             = ind_occ_p.all$pop
 , p_zeros           = ind_occ_p.all$p_zeros

  ## covariates
 , N_bd            = nrow(X_bd.m.all)
 , X_bd            = X_bd.m.all$bd  
 , ii_bd           = X_bd.m.all$ind
 , tt_bd           = X_bd.m.all$times
 , temp            = temp
 , time_gaps       = ind_occ_phi.all$time_gaps
  
  ## Capture data
 , N_y             = nrow(ind_occ_p.all)
 , y               = ind_occ_p.all$captures
  
 , first           = capture_range.all$first
 , last            = capture_range.all$final
  
 , present         = present.all
  )

stan.fit  <- stan(
  file    = "CMR_ind_pat_bd-p-phi_multi_recruit_free_temp_db_simple_mp_cv_ir.stan"
, data    = stan_data
, chains  = 1
, iter    = stan.iter
, warmup  = stan.burn
, thin    = stan.thin
, refresh = 10
, control = list(adapt_delta = 0.92, max_treedepth = 12)
  )

# stan.fit <- readRDS("stan.fit.re.Rds")

####
## CMR Diagnostics
####

## NOTE: OCT 15: just a tiny bit moved from Individual_CRM_expanding.R
 ## need to more thoroughly clean this up to accommodate "long" form
  ## will also need quite a bit of cleanup when bd parameters start varying by population
stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)
shinystan::launch_shinystan(stan.fit)

####
## Recovery of simulated coefficients?
####

as.data.frame(stan.fit.summary[grep("beta_p", dimnames(stan.fit.summary)[[1]]), c(4, 6, 8)])

## Primary Bd effects
pred_coef        <- as.data.frame(stan.fit.summary[c(1:9), c(4, 6, 8)])
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

## the most extreme individual
most_extreme <- reshape2::melt(stan.fit.samples$X[
  , 55
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
