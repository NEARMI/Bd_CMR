################################################
## Run the stan model for a single population ##
################################################

stan_data     <- list(
  
  ## dimensional indexes 
   n_pop             = n_sites
 , n_pop_year        = nrow(sampled_years)
 , n_ind             = n_ind
 , ind_per_period_p  = max(capt_history.p$gamma_index) 
 , ind_per_period_bd = max(capt_history.phi$X_stat_index)
 , ind_occ           = nrow(capt_history.p)
 , ind_occ_min1      = nrow(capt_history.phi)
 , ind_time          = nrow(capt_history)
  
  ## short vector indexes 
 , ind_occ_size      = rep(colSums(n_occ), n_ind.per)
 , ind_occ_min1_size = rep(colSums(n_occ) - 1, n_ind.per)
 , phi_first_index   = phi_first_index
 , p_first_index     = p_first_index
  
  ## long vector indexes: detection stuff (p)
 , ind_occ_rep       = capt_history.p$Mark
 , p_year            = as.numeric(as.factor(capt_history.p$Year))
 , p_month           = as.numeric(as.factor(capt_history.p$Month))
 , p_zeros           = capt_history.p$p_zeros
 , p_bd_index        = capt_history.p$X_stat_index
 , gamma_index       = capt_history.p$gamma_index
 , p_effort          = capt_history.p$effort
  
  ## long vector indexes: survival stuff (phi)
 , ind_occ_min1_rep  = capt_history.phi$Mark
 , offseason         = capt_history.phi$offseason
 , phi_year          = as.numeric(as.factor(capt_history.phi$Year))
 , phi_zeros         = capt_history.phi$phi_zeros
 , phi_ones          = capt_history.phi$phi_ones
 , phi_bd_index      = capt_history.phi$X_stat_index
 , capt_gaps         = capt_history.phi$capture_gap

  ## long vector indexes: bd stuff (bd)
 , ind_bd_rep        = (capt_history.phi %>% group_by(X_stat_index) %>%slice(1))$Mark
 , bd_time           = (capt_history.phi %>% group_by(X_stat_index) %>%slice(1))$Year %>% as.factor() %>% as.numeric()

  ## individual-level covariates, bd and others
 , N_bd              = nrow(capt_history.bd_load)
 , X_bd              = capt_history.bd_load$log_bd_load
 , X_ind             = capt_history.bd_load$Mark
  
 , bd_first_index    = bd_first_index
 , bd_last_index     = bd_last_index
 , x_bd_index        = capt_history.bd_load$X_stat_index
  
 , ind_size          = ind.size
 , ind_hg            = ind.hg
  
 , ind_len_which_have = len.have
 , ind_len_which_mis  = len.mis %>% as.array()
 , n_ind_len_have     = length(len.have)
 , n_ind_len_mis      = length(len.mis)
 , ind_len_have       = ind.len[len.have]
 , ind_len_mean       = 0
 , ind_len_sd         = 1
  
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
  
stan.fit  <- try(
  {
 stan(
# file    = "CMR_single_population.stan"
# file    = "CMR_single_population_con.stan"
# file    = "CMR_single_population_con_mi.stan"
# file    = "CMR_single_population_con_mi2.stan"
# file    = "stan_current/CMR_single_population_con_mi2.stan"
  file    = "stan_current/CMR_single_population_con_mi2_p_phi_adj.stan"
, data    = stan_data
, chains  = 1
, cores   = 1
, refresh = 10
, iter    = stan.iter            
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.92, max_treedepth = 12)
  )
  }
, silent = TRUE
)

saveRDS(stan.fit, paste(paste("fits/stan_fit", which.dataset, sep = "_"), "Rds", sep = "."))

#stan.fit.summary <- summary(stan.fit)[[1]]
#stan.fit.samples <- extract(stan.fit)
