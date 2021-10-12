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
  file    = "CMR_ind_pat_bd-p-phi_no_bd.stan"
, data    = stan_data
, chains  = 1
, iter    = stan.iter
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.92, max_treedepth = 12)
  )

####
## With modeling bd
####

## Indices to determine which of the latent bd values can be informed by measured bd
X_bd.m <- reshape2::melt(measured_bd)
names(X_bd.m) <- c("ind", "times", "bd")
X_bd.m %<>% mutate(times = plyr::mapvalues(times, from = seq(samp), to = sort(sampling_days)))
bd_m.m <- reshape2::melt(bd.measured)
names(bd_m.m) <- c("ind", "times", "meas")
bd_m.m %<>% mutate(times = plyr::mapvalues(times, from = seq(samp), to = sort(sampling_days)))
all_bd <- expand.grid(
  ind   = seq(ind)
, times = seq(times)
)
all_bd %<>% left_join(., bd_m.m)
samp_indices <- which(all_bd$meas == 1)

X_bd.m %<>% left_join(., bd_m.m)
X_bd.m %<>% filter(meas == 1)

y.m    <- reshape2::melt(capture_matrix) %>% arrange(Var1)
names(y.m) <- c("ind", "event_number", "captured")
# temp <- (expdat %>% group_by(all_times) %>% summarize(mean_temp = mean(temp)))$mean_temp

## Indices to determine which entries of phi need to be set to 0
phi_zeros <- matrix(data = 0, nrow = ind, ncol = sum(samp) - 1)
for (i in 1:ind) {
  phi_zeros[i, ] <- c(
    rep(1, capture_range$first[i] - 1)
  , rep(0, ncol(phi_zeros) - (capture_range$first[i] - 1))
    )
}
phi_zeros   <- (phi_zeros %>% reshape2::melt() %>% arrange(Var1))$value

## Index of which entries of phi are each individual
indiv_index <- rep(sum(samp), each = ind)

present <- rep(1, ind)

## Index for whether an individual was known to be present in the population at each time
p_zeros <- matrix(data = 0, nrow = ind, ncol = sum(samp))
for (i in 1:ind) {
  p_zeros[i, ] <- rep(present[i], samp)
}
p_zeros   <- (p_zeros %>% reshape2::melt() %>% arrange(Var1))$value

## Index for when Chi is one
chi_one                  <- matrix(data = 0, nrow = ind, ncol = sum(samp))
chi_one[, ncol(chi_one)] <- 1
chi_one                  <- (chi_one %>% reshape2::melt() %>% arrange(Var1))$value

## debugging checks
long_form_data.occ <- 
  cbind(
    y.m
  , data.frame(
    ind_occ_rep       = rep(seq(ind), each = (sum(samp)))
  , sampling_events_p = rep(sort(sampling_days), ind)
  )
  )

names(X_bd.m)[2] <- c("sampling_events_p")

long_form_data.occ %<>% left_join(., X_bd.m)

X_bd.m %<>% arrange(ind)

X_bd.m %>% {
  ggplot(., aes(sampling_events_p, bd)) + 
    geom_line() +
    scale_x_continuous(limits = c(1, 20)) +
    facet_wrap(~ind)
}

long_form_data.occ %>% {
  ggplot(., aes(sampling_events_p, bd)) + geom_line(aes(group = ind))
}

long_form_data.occ %>% filter(ind == 1)
measured_bd[1, ]

long_form_data.occm1 <- data.frame(
  ind_occ_min1_rep    = rep(seq(ind), each = (sum(samp) - 1))
, sampling_events_phi = rep(sort(sampling_days)[-10], ind)
, phi_zeros           = phi_zeros
, time_gaps           = rep(time_gaps, ind)
  )

stan_data     <- list(
  
  ## dimensional indexes 
   n_ind           = ind                  
 , n_times         = times
 , times_within    = times
 , n_occasions     = sum(samp)
 , n_occ_min1      = sum(samp) - 1
 , ind_time        = ind * times
 , ind_occ         = ind * sum(samp)
 , ind_occ_min1    = ind * (sum(samp) - 1)
  
  ## short vector indexes 
 , time                = seq(times)
 , sampling_events_phi = long_form_data.occm1$sampling_events_phi
 , sampling_events_p   = long_form_data.occ$sampling_events_p
  
  ## long vector indexes
 , ind_occ_min1_rep = long_form_data.occm1$ind_occ_min1_rep
 , ind_occ_rep      = long_form_data.occ$ind_occ_rep
  
  ## covariates
 , N_bd            = nrow(X_bd.m)
 , X_bd            = X_bd.m$bd  
 , ii_bd           = X_bd.m$ind
 , tt_bd           = X_bd.m$sampling_events_p
 , time_gaps       = rep(time_gaps, ind)
  
  ## Capture data
 , N_y             = nrow(long_form_data.occ)
 , y               = long_form_data.occ$captured
  
 , first           = capture_range$first
 , phi_zeros       = phi_zeros

 , last            = capture_range$final
  
  )

stan.fit  <- stan(
  file    = "CMR_ind_pat_bd-p-phi_multi_recruit_free_temp_db_simple_small.stan"
, data    = stan_data
, chains  = 1
, iter    = stan.iter
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.92, max_treedepth = 12)
  )

shinystan::launch_shinystan(stan.fit)

####
## Without modeling bd
####

## Indices to determine which of the latent bd values can be informed by measured bd
X_bd.m <- reshape2::melt(measured_bd) 
names(X_bd.m) <- c("ind", "times", "bd")
X_bd.m %<>% arrange(ind)

y.m    <- reshape2::melt(capture_matrix) %>% arrange(Var1)
names(y.m) <- c("ind", "event_number", "captured")

## Indices to determine which entries of phi need to be set to 0
phi_zeros <- matrix(data = 0, nrow = ind, ncol = sum(samp) - 1)
for (i in 1:ind) {
  phi_zeros[i, ] <- c(
    rep(1, capture_range$first[i] - 1)
  , rep(0, ncol(phi_zeros) - (capture_range$first[i] - 1))
    )
}
phi_zeros   <- (phi_zeros %>% reshape2::melt() %>% arrange(Var1))$value

## Index of which entries of phi are each individual
indiv_index <- rep(sum(samp), each = ind)

present <- rep(1, ind)

## Index for whether an individual was known to be present in the population at each time
p_zeros <- matrix(data = 0, nrow = ind, ncol = sum(samp))
for (i in 1:ind) {
  p_zeros[i, ] <- rep(present[i], samp)
}
p_zeros   <- (p_zeros %>% reshape2::melt() %>% arrange(Var1))$value

## Index for when Chi is one
chi_one                  <- matrix(data = 0, nrow = ind, ncol = sum(samp))
chi_one[, ncol(chi_one)] <- 1
chi_one                  <- (chi_one %>% reshape2::melt() %>% arrange(Var1))$value

## debugging checks
long_form_data.occ <- 
  cbind(
    y.m
  , data.frame(
    ind_occ_rep       = rep(seq(ind), each = (sum(samp)))
  , sampling_events_p = rep(sort(sampling_days), ind)
  )
  )

names(X_bd.m)[2] <- c("sampling_events_p")

long_form_data.occ %<>% left_join(., X_bd.m)

long_form_data.occ %>% {
  ggplot(., aes(sampling_events_p, bd)) + geom_line(aes(group = ind))
}

long_form_data.occ %>% filter(ind == 1)
measured_bd[1, ]

long_form_data.occm1 <- data.frame(
  ind_occ_min1_rep    = rep(seq(ind), each = (sum(samp) - 1))
, sampling_events_phi = rep(sort(sampling_days)[-10], ind)
, phi_zeros           = phi_zeros
, time_gaps           = rep(time_gaps, ind)
  )

stan_data     <- list(
  
  ## dimensional indexes 
   n_ind           = ind                  
 , n_times         = times
 , n_occasions     = sum(samp)
 , n_occ_min1      = sum(samp) - 1
 , ind_occ         = ind * sum(samp)
 , ind_occ_min1    = ind * (sum(samp) - 1)
  
  ## short vector indexes 
 , time            = seq(times)
  
  ## long vector indexes
 , ind_occ_min1_rep    = long_form_data.occm1$ind_occ_min1_rep
 , sampling_events_phi = long_form_data.occm1$sampling_events_phi
 , phi_zeros           = long_form_data.occm1$phi_zeros
  
 , ind_occ_rep         = long_form_data.occ$ind_occ_rep
 , sampling_events_p   = long_form_data.occ$sampling_events_p
  
  ## covariates
 , X_bd            = measured_bd
 , time_gaps       = long_form_data.occm1$time_gaps
  
  ## Capture data
 , N_y             = nrow(long_form_data.occ)
 , y               = long_form_data.occ$captured
 , first           = capture_range$first
 , last            = capture_range$final
  )

stan.fit  <- stan(
  file    = "CMR_ind_pat_bd-p-phi_multi_recruit_free_temp_db_simple_small_no_bd.stan"
, data    = stan_data
, chains  = 1
, iter    = stan.iter
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.92, max_treedepth = 12)
  )




