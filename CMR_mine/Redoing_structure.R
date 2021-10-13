####
## Run the model in Stan, long structure to accommodate multiple populations
####

## Indices to determine which of the latent bd values can be informed by measured bd
X_bd.m <- reshape2::melt(measured_bd)
names(X_bd.m) <- c("ind", "times", "bd")
X_bd.m %<>% mutate(times = plyr::mapvalues(times, from = seq(sum(samp)), to = sampling_times_all))
bd_m.m <- reshape2::melt(bd_measured)
names(bd_m.m) <- c("ind", "times", "meas")
bd_m.m %<>% mutate(times = plyr::mapvalues(times, from = seq(sum(samp)), to = sampling_times_all))
all_bd <- expand.grid(
  ind   = seq(all_ind)
, times = seq(times * periods)
)
all_bd %<>% left_join(., bd_m.m)
samp_indices <- which(all_bd$meas == 1)

X_bd.m %<>% left_join(., bd_m.m)
X_bd.m %<>% filter(meas == 1)
X_bd.m %<>% arrange(ind)

## Indices to determine which entries of phi need to be set to 0
y.m    <- reshape2::melt(capture_matrix) %>% arrange(Var1)
temp   <- (expdat %>% group_by(all_times) %>% summarize(mean_temp = mean(temp)))$mean_temp

## Indices for which entries of phi must be 0
phi_zeros <- matrix(data = 0, nrow = all_ind, ncol = sum(samp) - 1)
for (i in 1:all_ind) {
  phi_zeros[i, ] <- c(
    rep(1, capture_range$first[i] - 1)
  , rep(0, ncol(phi_zeros) - (capture_range$first[i] - 1))
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

ind_time <-  data.frame(
   ind_time_rep  = rep(seq(all_ind), each = times * periods)
 , time_rep      = rep(seq(times), all_ind * periods)
)

ind_occ_phi <- data.frame(
   ind_occ_min1_rep    = rep(seq(all_ind), each = (sum(samp) - 1))
 , sampling_events_phi = rep(sampling_times_all[-sum(samp)], all_ind)
 , offseason           = rep(offseason_vec, all_ind)
 , time_gaps           = rep(time_gaps, all_ind)
 , phi_zeros           = phi_zeros
)

ind_occ_p <- data.frame(
   ind_occ_rep       = rep(seq(all_ind), each = (sum(samp)))
 , sampling_events_p = rep(sampling_times_all, all_ind)
 , periods_occ       = rep(periods_occ, all_ind)
 , captures          = y.m$value
 , p_zeros           = p_zeros
)

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
 , periods         = periods_time
  
  ## long vector indexes
 , ind_time_rep    = ind_time$ind_time_rep
 , time_rep        = ind_time$time_rep
  
 , ind_occ_min1_rep    = ind_occ_phi$ind_occ_min1_rep
 , sampling_events_phi = ind_occ_phi$sampling_events_phi
 , offseason           = ind_occ_phi$offseason
  
 , ind_occ_rep       = ind_occ_p$ind_occ_rep
 , sampling_events_p = ind_occ_p$sampling_events_p
 , periods_occ       = ind_occ_p$periods_occ

  ## covariates
 , N_bd            = nrow(X_bd.m)
 , X_bd            = X_bd.m$bd  
 , ii_bd           = X_bd.m$ind
 , tt_bd           = X_bd.m$times
 , temp            = temp
 , time_gaps       = ind_occ_phi$time_gaps
 , samp_indices    = samp_indices
  
  ## Capture data
 , N_y             = nrow(ind_occ_p)
 , y               = ind_occ_p$captures
  
 , first           = capture_range$first
 , last            = capture_range$final
  
 , phi_zeros       = ind_occ_phi$phi_zeros
 , p_zeros         = ind_occ_p$p_zeros
  
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

####
## WORK ON A TWO POP MODEL ____ SEE TOP NOTES FOR THE WAY TO ACTUALLY GO
## BUT:::for now just manually stick together two simulations

all_ind.1            <- all_ind
n_occasions.1        <- sum(samp)
n_occ_min1.1         <- sum(samp) - 1
sampling_times_all.1 <- sampling_times_all
time_gaps.1          <- time_gaps
offseason_vec.1      <- offseason_vec
capture_matrix.1     <- capture_matrix
capture_range.1      <- capture_range
present.1            <- present
measured_bd.1        <- measured_bd
bd_measured.1        <- bd_measured
periods_occ.1        <- periods_occ
ind_occ_phi.1        <- ind_occ_phi
ind_occ_p.1          <- ind_occ_p
X_bd.m.1             <- X_bd.m

ind_occ_phi.1 %<>% mutate(pop = 1, ind = interaction(ind_occ_min1_rep, pop))
ind_occ_p.1 %<>% mutate(pop = 1, ind = interaction(ind_occ_rep, pop))
# X_bd.m.1 %<>% mutate(pop = 1, ind = interaction(ind, pop))

capture_range.1 %<>% mutate(pop = 1, ind = interaction(ind, pop))

ind_occ_size.1 <- rep(n_occasions.1, all_ind.1)

all_ind.2            <- all_ind
n_occasions.2        <- sum(samp)
n_occ_min1.2         <- sum(samp) - 1
sampling_times_all.2 <- sampling_times_all
time_gaps.2          <- time_gaps
offseason_vec.2      <- offseason_vec
capture_matrix.2     <- capture_matrix
capture_range.2      <- capture_range
present.2            <- present
measured_bd.2        <- measured_bd
bd_measured.2        <- bd_measured
periods_occ.2        <- periods_occ
ind_occ_phi.2        <- ind_occ_phi
ind_occ_p.2          <- ind_occ_p
X_bd.m.2             <- X_bd.m

ind_occ_phi.2 %<>% mutate(pop = 2, ind = interaction(ind_occ_min1_rep, pop))
ind_occ_p.2 %<>% mutate(pop = 2, ind = interaction(ind_occ_rep, pop))
# X_bd.m.2 %<>% mutate(pop = 2, ind = interaction(ind, pop))
X_bd.m.2 %<>% mutate(ind = ind + all_ind.1)

capture_range.2 %<>% mutate(pop = 2, ind = interaction(ind, pop))

ind_occ_size.2 <- rep(n_occasions.2, all_ind.2)

### Put the pops together

ind_occ_phi   <- rbind(ind_occ_phi.1, ind_occ_phi.2) %>% mutate(ind = as.numeric(ind))
ind_occ_p     <- rbind(ind_occ_p.1, ind_occ_p.2) %>% mutate(ind = as.numeric(ind))
X_bd.m        <- rbind(X_bd.m.1, X_bd.m.2)
capture_range <- rbind(capture_range.1, capture_range.2) %>% mutate(ind = as.numeric(ind))

present       <- rbind(present.1, present.2)

ind_occ_size      <- c(ind_occ_size.1, ind_occ_size.2)
ind_occ_min1_size <- ind_occ_size - 1

## Index vector for the first entry of phi and p that correspond to a new individual
phi_first_index <- (ind_occ_phi %>% mutate(index = seq(n())) %>% group_by(ind) %>% 
  summarize(first_index = min(index)))$first_index

p_first_index <- (ind_occ_p %>% mutate(index = seq(n())) %>% group_by(ind) %>% 
  summarize(first_index = min(index)))$first_index

stan_data     <- list(
  
  ## dimensional indexes 
   n_periods       = periods
 , n_ind           = all_ind.1 + all_ind.2                
 , n_times         = times * periods
 , times_within    = times
 , ind_occ         = all_ind.1 * n_occasions.1 + all_ind.2 * n_occasions.2
 , ind_occ_min1    = all_ind.1 * n_occ_min1.1 + all_ind.2 * n_occ_min1.2
  
  ## short vector indexes 
 , time              = rep(seq(times), periods)
 , time_per_period   = matrix(data = seq(times * periods), nrow = times, ncol = periods)
 , periods           = periods_time
 , ind_occ_size      = ind_occ_size
 , ind_occ_min1_size = ind_occ_min1_size

 , phi_first_index   = phi_first_index
 , p_first_index     = p_first_index
  
  ## long vector indexes
 , ind_occ_min1_rep    = ind_occ_phi$ind
 , sampling_events_phi = ind_occ_phi$sampling_events_phi
 , offseason           = ind_occ_phi$offseason
  
 , ind_occ_rep       = ind_occ_p$ind
 , sampling_events_p = ind_occ_p$sampling_events_p
 , periods_occ       = ind_occ_p$periods_occ

  ## covariates
 , N_bd            = nrow(X_bd.m)
 , X_bd            = X_bd.m$bd  
 , ii_bd           = X_bd.m$ind
 , tt_bd           = X_bd.m$times
 , temp            = temp
 , time_gaps       = ind_occ_phi$time_gaps
  
  ## Capture data
 , N_y             = nrow(ind_occ_p)
 , y               = ind_occ_p$captures
  
 , first           = capture_range$first
 , last            = capture_range$final
  
 , phi_zeros       = ind_occ_phi$phi_zeros
 , p_zeros         = ind_occ_p$p_zeros
  
 , present         = present
  )

stan.fit  <- stan(
  file    = "CMR_ind_pat_bd-p-phi_multi_recruit_free_temp_db_simple_mp.stan"
, data    = stan_data
, chains  = 1
, iter    = stan.iter
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.92, max_treedepth = 12)
  )

####
## If the populations have the same number of sampling periods, the other model should work
## This will be a good way to debug the new model structure
####

stan_data     <- list(
  
  ## dimensional indexes 
   n_periods       = periods
 , n_ind           = all_ind.1 + all_ind.2                 
 , n_times         = times * periods
 , times_within    = times
 , n_occasions     = sum(samp)
 , n_occ_min1      = sum(samp) - 1
 , ind_occ         = all_ind.1 * n_occasions.1 + all_ind.2 * n_occasions.2
 , ind_occ_min1    = all_ind.1 * n_occ_min1.1 + all_ind.2 * n_occ_min1.2
  
  ## short vector indexes 
 , time            = rep(seq(times), periods)
 , time_per_period = matrix(data = seq(times * periods), nrow = times, ncol = periods)
 , periods         = periods_time
  
  ## long vector indexes
 , ind_time_rep    = ind_time$ind_time_rep
 , time_rep        = ind_time$time_rep
  
 , ind_occ_min1_rep    = ind_occ_phi$ind_occ_min1_rep
 , sampling_events_phi = ind_occ_phi$sampling_events_phi
 , offseason           = ind_occ_phi$offseason
  
 , ind_occ_rep       = ind_occ_p$ind_occ_rep
 , sampling_events_p = ind_occ_p$sampling_events_p
 , periods_occ       = ind_occ_p$periods_occ

  ## covariates
 , N_bd            = nrow(X_bd.m)
 , X_bd            = X_bd.m$bd  
 , ii_bd           = X_bd.m$ind
 , tt_bd           = X_bd.m$times
 , temp            = temp
 , time_gaps       = ind_occ_phi$time_gaps
 , samp_indices    = samp_indices
  
  ## Capture data
 , N_y             = nrow(ind_occ_p)
 , y               = ind_occ_p$captures
  
 , first           = capture_range$first
 , last            = capture_range$final
  
 , phi_zeros       = ind_occ_phi$phi_zeros
 , p_zeros         = ind_occ_p$p_zeros
  
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


