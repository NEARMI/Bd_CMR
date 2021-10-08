####
## Run the model in Stan, long structure to accomodate multiple populations
####

## Indices to determine which of the latent bd values can be informed by measured bd
X_bd.m <- reshape2::melt(measured_bd)
names(X_bd.m) <- c("ind", "times", "bd")
X_bd.m %<>% mutate(times = plyr::mapvalues(times, from = seq(30), to = sampling_times_all))
bd_m.m <- reshape2::melt(bd_measured)
names(bd_m.m) <- c("ind", "times", "meas")
bd_m.m %<>% mutate(times = plyr::mapvalues(times, from = seq(30), to = sampling_times_all))
all_bd <- expand.grid(
  ind   = seq(all_ind)
, times = seq(times * periods)
)
all_bd %<>% left_join(., bd_m.m)
samp_indices <- which(all_bd$meas == 1)

X_bd.m %<>% left_join(., bd_m.m)
X_bd.m %<>% filter(meas == 1)

## Indices to determine which entries of phi need to be set to 0

y.m    <- reshape2::melt(capture_matrix)
temp   <- (expdat %>% group_by(all_times) %>% summarize(mean_temp = mean(temp)))$mean_temp

## Indices for which entries of phi must be 0
phi_zeros <- matrix(data = 0, nrow = all_ind, ncol = sum(samp) - 1)
for (i in 1:all_ind) {
  phi_zeros[i, ] <- c(
    rep(1, capture_range$first[i] - 1)
  , rep(0, ncol(phi_zeros) - capture_range$first[i] + 1)
    )
}
phi_zeros   <- (phi_zeros %>% reshape2::melt() %>% arrange(Var1))$value

## Index of which entries of phi are each individual
indiv_index <- rep(sum(samp), each = all_ind)

## Index for whether an individual was known to be present in the population at each time
p_zeros <- matrix(data = 0, nrow = all_ind, ncol = sum(samp))
for (i in 1:all_ind) {
  p_zeros[i, ] <- rep(present[i, ], samp)
}
p_zeros   <- (p_zeros %>% reshape2::melt() %>% arrange(Var1))$value

## Index for when Chi is one
chi_one                  <- matrix(data = 0, nrow = all_ind, ncol = sum(samp))
chi_one[, ncol(chi_one)] <- 1
chi_one                  <- (chi_one %>% reshape2::melt() %>% arrange(Var1))$value

stan_data     <- list(
  
  ## dimensional indexes 
   n_periods       = periods
 , n_ind           = all_ind                  
 , n_times         = times * periods
 , times_within    = times
 , n_occasions     = sum(samp)
 , n_occ_min1      = sum(samp) - 1
 , ind_time        = all_ind * times * periods
 , ind_occ         = all_ind * sum(samp)
 , ind_occ_min1    = all_ind * (sum(samp) - 1)
  
  ## short vector indexes 
 , time            = rep(seq(times), periods)
 , time_per_period = matrix(data = seq(times * periods), nrow = times, ncol = periods)
 , sampling_events_phi = rep(sampling_times_all[-30], all_ind)
 , sampling_events_p = rep(sampling_times_all, all_ind)
 , offseason       = rep(offseason_vec, all_ind)
 , periods         = periods_time
 , periods_occ     = rep(periods_occ, all_ind)
  
  ## long vector indexes
 , ind_time_rep    = rep(seq(all_ind), each = times * periods)
 , time_rep        = rep(seq(times), all_ind * periods)
 , temp_rep        = rep(temp, all_ind)
  
  ## covariates
 , N_bd            = nrow(X_bd.m)
 , X_bd            = X_bd.m$bd  
 , ii_bd           = X_bd.m$ind
 , tt_bd           = X_bd.m$times
 , temp            = temp
 , time_gaps       = rep(time_gaps, all_ind)
 , samp_indices    = samp_indices
  
  ## Capture data
 , N_y             = length(y.m$value)
 , y               = y.m$value
  
 , first           = capture_range$first
 , phi_zeros       = phi_zeros
 , p_zeros         = p_zeros
  
 , last            = capture_range$final
 , present         = present
  
  )

stan.fit  <- stan(
  file    = "CMR_ind_pat_bd-p-phi_multi_recruit_free_temp_db_simple.stan"
, data    = stan_data
, chains  = 1
, iter    = stan.iter
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.92, max_treedepth = 12)
  )
