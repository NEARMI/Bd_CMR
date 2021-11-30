########################
## Run the stan model ##
########################

stan_data     <- list(
  
  ## dimensional indexes 
   n_pop             = n_sites
 , n_ind             = n_ind
 , ind_per_period_bd = max(capt_history.phi$X_stat_index)
 , ind_per_period_p  = max(capt_history.p$gamma_index) 
 , ind_time          = nrow(capt_history)
 , ind_occ           = nrow(capt_history.p)
 , ind_occ_min1      = nrow(capt_history.phi)
  
  ## short vector indexes 
 , ind_occ_size      = rep(colSums(n_occ), n_ind.per)
 , ind_occ_min1_size = rep(colSums(n_occ) - 1, n_ind.per)
 , p_first_index     = p_first_index
 , phi_first_index   = phi_first_index
  
  ## long vector indexes: detection stuff (p)
 , ind_occ_rep       = capt_history.p$Mark
 , periods_occ       = as.numeric(as.factor(capt_history.p$Year))
 , p_month           = as.numeric(as.factor(capt_history.p$Month))
 , pop_p             = as.numeric(as.factor(capt_history.p$Site))
 , p_zeros           = capt_history.p$p_zeros
 , p_bd_index        = capt_history.p$p_bd_index
 , gamma_index       = capt_history.p$gamma_index
  
  ## long vector indexes: survival stuff (phi)
 , ind_occ_min1_rep  = capt_history.phi$Mark
 , offseason         = capt_history.phi$offseason
 , phi_month         = as.numeric(as.factor(capt_history.phi$Month))
 , phi_year          = as.numeric(as.factor(capt_history.phi$Year))
 , pop_phi           = as.numeric(as.factor(capt_history.phi$Site))
 , phi_zeros         = capt_history.phi$phi_zeros
 , phi_ones          = capt_history.phi$phi_ones
 , phi_bd_index      = capt_history.phi$phi_bd_index
 , X_stat_index      = capt_history.phi$X_stat_index
 , time_gaps         = capt_history.phi$time_gaps  

  ## long vector indexes: bd stuff (bd)
 , ind_bd_rep        = capt_history$Mark
 , ind_in_pop        = as.numeric(as.factor(capt_history$Site))

  ## covariates, bd and others
 , N_bd              = nrow(capt_history.bd_load)
 , X_bd              = capt_history.bd_load$log_bd_load  
 , X_ind             = capt_history.bd_load$Mark
 , ind_size          = ind.size
 , ind_hg            = ind.hg
  
  ## Capture data
 , N_y             = nrow(capt_history.p)
 , y               = capt_history.p$captured
 , first           = capture_range$first
 , last            = capture_range$final

  )

stan.fit  <- stan(
  file    = "CMR_collapsed.stan"
, data    = stan_data
, chains  = 1
, cores   = 1
, refresh = 20
, iter    = stan.iter            
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.92, max_treedepth = 12)
  )

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)
