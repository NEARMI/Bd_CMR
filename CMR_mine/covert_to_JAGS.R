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

library(rjags)
library(jagsUI)

nadapt <-40000 # adaption phase
nburn <- 15000 # discard these draws
niter <- 40000 # length of MCMC chains
nsamp <- 1000  # number of samples to take from each MCMC chains
thin_ <-  10   #round(niter/nsamp) # only

stan_data     <- list(
  ## bookkeeping params
   n_ind           = ind
 , n_times         = times
 , n_occasions     = samp
 , n_occ_min1      = samp - 1
 , time            = seq(times)
 , time_sq         = seq(times) ^ 2
 , sampling_events = sort(sampling_days)
  ## Capture data
 , y               = capture_matrix
 , first           = capture_range$first
 , last            = capture_range$final
  ## Covariate associated parameters
 , X_bd            = measured_bd
 , time_gaps       = time_gaps
  )

params <- c(
  "beta_phi"
, "beta_p"
, "beta_timegaps"
, "beta_bd"
, "bd_delta_sigma"
, "db_delta_eps"
, "bd_obs"
, "phi"
, "p"
, "X"
  )

ralu.Bd <- jags(
  data       = stan_data
, model.file = "CMR_ind_pat_bd-p-phi_single.txt"
, parameters.to.save = params
, init = list(
  list(
    beta_phi = c(2, -2)
  , beta_p = c(1, 1)
  , beta_bd = c(1, 1, -1)
  , beta_timegaps = c(-1)
  , tau = c(1, 1)
  , bd_delta_eps = rep(0, ind)
  )
)
, n.iter     = niter
, n.chains   = 1
, n.burnin   = nburn
, n.thin     = thin_
, parallel   = TRUE
, verbose = T)

## known alive matrix
known_living <- matrix(data = 0, nrow = ind, ncol = samp)
for (i in 1:nrow(known_living)) {
  known_living[i, 1:capture_range$final[i]] <- 1
}

## Trying to move over to JAGS
measured_bd   <- matrix(
  nrow = ind
, ncol = times
, data = expdat$log_bd_load
, byrow = T)

bd.measured   <- matrix(
  nrow = ind
, ncol = times
, data = expdat$bd_swabbed
, byrow = T)

for (i in 1:nrow(measured_bd)) {
  for (j in 1:ncol(measured_bd)) {
    if (bd.measured[i,j] == 0) {
      measured_bd[i,j] <- NA
    }
  }
}

stan_data     <- list(
  ## bookkeeping params
   n_ind           = ind
 , n_times         = times
 , n_occasions     = samp
 , n_occ_min1      = samp - 1
 , time            = seq(times)
 , time_sq         = seq(times) ^ 2
 , sampling_events = sort(sampling_days)
  ## Capture data
 , y               = capture_matrix
 , first           = capture_range$first
 , last            = capture_range$final
 , known_living    = known_living
 , chi_restraint   = rep(1, ind)
  ## Covariate associated parameters
 , X_bd            = measured_bd
 , time_gaps       = time_gaps
  )

params <- c(
  "beta_phi"
, "beta_p"
, "beta_timegaps"
, "beta_bd"
  )

system.time(
ralu.Bd <- jags(
  data       = stan_data
, model.file = "CMR_ind_pat_bd-p-phi_single_simple_within.txt"
, parameters.to.save = params
, n.iter     = niter
, n.chains   = 1
, n.burnin   = nburn
, n.thin     = thin_
, parallel   = F
, verbose = T)
)

stan.fit.summary[1:10, ]
ralu.Bd$summary

## survival over Bd load
stan.pred        <- apply(ralu.Bd$samples[[1]][,1:2], 1
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
stan.pred        <- apply(ralu.Bd$samples[[1]][, 3:4], 1
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
