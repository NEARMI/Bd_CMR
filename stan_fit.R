#################################################
## Run the stan model for multiple populations ##
#################################################

stan_data     <- list(
  
  ## dimensional indexes 
   n_pop             = n_sites
 , n_pop_year        = nrow(sampled_years)
 , n_ind             = n_ind
 , n_spec            = length(unique(capt_history$Species))
 , ind_per_period_bd = max(capt_history.phi$X_stat_index)
 , ind_per_period_p  = max(capt_history.p$gamma_index) 
 , ind_time          = nrow(capt_history)
 , ind_occ           = nrow(capt_history.p)
 , ind_occ_min1      = nrow(capt_history.phi)
  
  ## short vector indexes (length of n_ind)
 , ind_occ_size      = rep(colSums(n_occ), n_ind.per)
 , ind_occ_min1_size = rep(colSums(n_occ) - 1, n_ind.per)
 , p_first_index     = p_first_index
 , phi_first_index   = phi_first_index
 , ind_which_pop     = ind_which_pop
 
  ## short vector indexes (length of X_stat_index)
 , ind_bd_rep        = X_stat_index_covs$Mark  
 , ind_in_pop_year   = X_stat_index_covs$ind_in_pop_year
 , pop_for_bd        = X_stat_index_covs$pop_for_bd
  
  ## long vector indexes: detection stuff (p)
 , ind_occ_rep       = capt_history.p$Mark
 , periods_occ       = as.numeric(as.factor(capt_history.p$Year))
 , p_month           = as.numeric(as.factor(capt_history.p$Month))
 , pop_p             = as.numeric(capt_history.p$pop_spec)
 , spec_p            = as.numeric(capt_history.p$Species)
 , p_zeros           = capt_history.p$p_zeros
 , p_bd_index        = capt_history.p$X_stat_index
 , gamma_index       = capt_history.p$gamma_index
  
  ## long vector indexes: survival stuff (phi)
 , ind_occ_min1_rep  = capt_history.phi$Mark
 , offseason         = capt_history.phi$offseason
 , phi_month         = as.numeric(as.factor(capt_history.phi$Month))
 , phi_year          = as.numeric(as.factor(capt_history.phi$Year))
 , pop_phi           = as.numeric(capt_history.phi$pop_spec)
 , spec_phi          = as.numeric(capt_history.phi$Species)
 , phi_zeros         = capt_history.phi$phi_zeros
 , phi_ones          = capt_history.phi$phi_ones
 , phi_bd_index      = capt_history.phi$X_stat_index
 , time_gaps         = capt_history.phi$time_gaps  
 , phi_pop_year      = capt_history.phi$pop_year

  ## individual-level covariates, bd and others
 , N_bd              = nrow(capt_history.bd_load)
 , X_bd              = capt_history.bd_load$log_bd_load
 , X_ind             = capt_history.bd_load$Mark
 , ind_size          = ind.size
 , bd_first_index    = bd_first_index
 , bd_last_index     = bd_last_index
 , x_bd_index        = capt_history.bd_load$X_stat_index
  
  ## site-level covariates, categorical 
   ## rely on the indexes pop_p (p), pop_phi (phi), and ind_in_pop (bd) for retrieving the correct covariate value
   ## from the correct pop:spec
 , pop_drawdown     = site_covar.cat$DRAWDOWN
  
  ## site-level covariates, continuous
   ## rely on the indexes pop_p (p), pop_phi (phi), and ind_in_pop (bd) for retrieving the correct covariate value
   ## from the correct pop:spec
 , pop_temp         = site_covar.con$Temp_Mean
  
  ## Capture data
 , N_y             = nrow(capt_history)
 , y               = capt_history.p$captured
 , first           = capture_range$first
 , last            = capture_range$final

  )

if (exists("ind_hg")) {
  
 stan_data <- c(stan_data, ind_hg = ind.hg)
  
}
  
stan.fit  <- stan(
  file    = "CMR_multiple_populations_red_ni_cov.stan"
, data    = stan_data
, chains  = 1
, cores   = 1
, refresh = 10
, iter    = stan.iter            
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.92, max_treedepth = 12)
, include = FALSE
, pars    = c("phi", "p", "chi")
  )

## stan.fit <- readRDS("no_int.Rds")
## shinystan::launch_shinystan(stan.fit)

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)
