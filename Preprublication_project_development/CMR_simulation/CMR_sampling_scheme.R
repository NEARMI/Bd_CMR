##############################################################
## Explore sampling effort on ability to recover parameters ##
##############################################################

####
## Notes as of Nov 1:
####

## First pass (still very rough -- will need lots of cleanup for use in a power analysis or w/e) 
 ## to explore variation in sampling and the ability of the model to 
  ## recover parameter estimates 
## !!! Plan will be to make the parameters script use columns for sims or lists of lists to loop over
 ## all parameters that the user feels like changing.
 ## -- Not really a priority right now though given the need to:
  ## 1) explore collapsing the model to better reflect much more of the real data
  ## 2) continue with the real data (explore many more covaraites)

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

### Loop over one parameter for now
test.what <- "perc"
num_runs  <- 10

if (test.what == "samp") {
samp      <- matrix(data = rep(seq(1, num_runs, by = 1)), nrow = periods, ncol = num_runs, byrow = T)
if (n_pop == 1) {
  samp <- list(samp)
}
} else {
bd_perc      <- matrix(data = rep(seq(0.1, 1, length = num_runs))
                       , nrow = periods, ncol = num_runs, byrow = T)
if (n_pop == 1) {
  bd_perc <- list(bd_perc)
}
}

for (samp_run in 1:num_runs) {

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
, samp      = {
  if(test.what == "samp") {
    samp[[pop_ind]][, samp_run]
  } else {
    samp[[pop_ind]]
  }
}
, bd_perc   = {
  if(test.what == "perc") {
    bd_perc[[pop_ind]][, samp_run]
  } else {
    bd_perc[[pop_ind]]
  }
}
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
, samp      = {
  if(test.what == "samp") {
    samp[[pop_ind]][, samp_run]
  } else {
    samp[[pop_ind]]
  }
}
)

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
  
source("CMR_dataclean.R")
  
if (test.what == "samp") {
samp.temp <- samp[[pop_ind]][, samp_run]
} else {
samp.temp <- samp[[pop_ind]]  
}
  
stan_data     <- list(
  
  ## dimensional indexes 
   n_pop           = n_pop
 , n_ind           = all_ind.all
 , ind_per_period  = sum(periods * each_ind.all)
  
 , ind_time        = sum(
   (expdat.all %>% group_by(pop) %>% summarize(total_times   = length(unique(times))))$total_times * 
   (expdat.all %>% group_by(pop) %>% summarize(total_periods = length(unique(periods))))$total_periods * 
   c(each_ind.all))
 , ind_occ         = sum(samp.temp * each_ind.all)
 , ind_occ_min1    = (sum(samp.temp) - 1) * each_ind.all

  ## short vector indexes 
 , ind_occ_size      = rep(sum(samp.temp), each_ind.all)
 , ind_occ_min1_size = rep(sum(samp.temp) - 1, each_ind.all)

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
, refresh = 20
, control = list(adapt_delta = 0.92, max_treedepth = 12)
  )

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)

pred_coef        <- as.data.frame(
  stan.fit.summary[grep("beta"
  , dimnames(stan.fit.summary)[[1]]), c(4, 6, 8)])
names(pred_coef) <- c("lwr", "mid", "upr")
pred_coef        %<>% mutate(
  param = rownames(.)
, samps = sum(samp.temp))

if (samp_run == 1) {
pred_coef.a <- pred_coef
} else {
pred_coef.a <- rbind(pred_coef.a, pred_coef)
}

print("Through", pred_coef, "Sims of", num_runs, sep = " ")

}
