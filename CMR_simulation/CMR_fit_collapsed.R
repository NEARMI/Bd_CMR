####
## Run the stan model with collapsed temporal dynamics
####

stan_data     <- list(
  
  ## dimensional indexes 
   n_pop             = n_pop
 , n_ind             = all_ind.all
 , ind_per_period_bd = max(ind_occ_phi.all$X_stat_index)
 , ind_per_period_p  = max(ind_occ_p.all$gamma_index) 
 , ind_time          = nrow(expdat.all)
 , ind_occ           = nrow(ind_occ_p.all)
 , ind_occ_min1      = nrow(ind_occ_phi.all)
  
  ## short vector indexes 
 , ind_occ_size      = rep(lapply(samp, sum) %>% unlist(), each_ind.all)
 , ind_occ_min1_size = rep(lapply(samp, sum) %>% unlist() - 1, each_ind.all)
 , p_first_index     = p_first_index
 , phi_first_index   = phi_first_index
  
  ## long vector indexes: detection stuff (p)
 , ind_occ_rep       = ind_occ_p.all$ind
 , periods_occ       = ind_occ_p.all$periods
 , pop_p             = ind_occ_p.all$pop
 , p_zeros           = ind_occ_p.all$p_zeros
 , p_bd_index        = ind_occ_p.all$p_bd_index
 , gamma_index       = ind_occ_p.all$gamma_index
  
  ## long vector indexes: survival stuff (phi)
 , ind_occ_min1_rep  = ind_occ_phi.all$ind
 , offseason         = ind_occ_phi.all$offseason
 , phi_year          = ind_occ_phi.all$periods
 , pop_phi           = ind_occ_phi.all$pop
 , phi_zeros         = ind_occ_phi.all$phi_zeros
 , phi_ones          = ind_occ_phi.all$phi_ones
 , phi_bd_index      = ind_occ_phi.all$phi_bd_index
 , X_stat_index      = ind_occ_phi.all$X_stat_index
 , time_gaps         = ind_occ_phi.all$time_gaps  

  ## long vector indexes: bd stuff (bd)
 , ind_bd_rep        = expdat.all$ind
 , ind_in_pop        = expdat.all$pop

  ## covariates, bd and others
 , N_bd              = nrow(X_bd.m.all)
 , X_bd              = X_bd.m.all$bd  
 , X_ind             = X_bd.m.all$ind
 , temp              = expdat.all$temp
  
  ## Capture data
 , N_y             = nrow(ind_occ_p.all)
 , y               = ind_occ_p.all$captures
 , first           = capture_range.all$first
 , last            = capture_range.all$final

  )

stan.fit  <- stan(
  file    = {
    if (n_pop == 1) {
     "CMR_collapsed_sim.stan"
    } else {
     "../CMR_collapsed_pr.stan"
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

# shinystan::launch_shinystan(stan.fit)

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)

