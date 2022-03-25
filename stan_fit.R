#################################################
## Run the stan model for multiple populations ##
#################################################

stan_data     <- list(
  
  ## dimensional indexes 
   n_pop             = n_sites
 , n_pop_year        = nrow(sampled_years)
 , n_ind             = n_ind
 , ind_per_period_bd = max(capt_history.phi$X_stat_index)
 , ind_occ           = nrow(capt_history.p)
 , ind_occ_min1      = nrow(capt_history.phi)
 , n_days            = (capt_history.p %>% group_by(pop_spec) %>% summarize(days_in_pop = n_distinct(date_fac)))$days_in_pop %>% sum()
 , n_spec            = length(unique(capt_history$Species))
  
  ## short vector indexes (length of n_ind)
 , ind_occ_size      = rep(colSums(n_occ), n_ind.per)           
 , ind_occ_min1_size = rep(colSums(n_occ) - 1, n_ind.per)    
 , phi_first_index   = phi_first_index
 , p_first_index     = p_first_index  
 , ind_in_pop        = ind_in_pop
  
  ## short vector indexes (length of n_pop)
 , spec_pop          = (capt_history.p %>% group_by(pop_spec) %>% slice(1))$Species %>% as.numeric()
  
  ## short vector indexes (length of n_days)
 , day_which_pop     = (capt_history.p %>% group_by(date_fac) %>% slice(1))$pop_spec %>% as.numeric()
 , spec_which_pop    = (capt_history.p %>% group_by(date_fac) %>% slice(1))$Species %>% as.numeric()
 
  ## long vector indexes: detection stuff (p)
 , ind_occ_rep       = capt_history.p$Mark
 , p_month           = as.numeric(as.factor(capt_history.p$Month))
 , p_zeros           = capt_history.p$p_zeros
 , p_bd_index        = capt_history.p$X_stat_index
 , p_day             = capt_history.p$date_fac
 , pop_p             = as.numeric(capt_history.p$pop_spec)
 , spec_p            = as.numeric(capt_history.p$Species)

  ## long vector indexes: survival stuff (phi)
 , ind_occ_min1_rep  = capt_history.phi$Mark
 , offseason         = capt_history.phi$offseason
 , phi_year          = as.numeric(as.factor(capt_history.phi$Year))
 , phi_zeros         = capt_history.phi$phi_zeros
 , phi_ones          = capt_history.phi$phi_ones
 , phi_bd_index      = capt_history.phi$X_stat_index
 , capt_gaps         = capt_history.phi$capture_gap
 , pop_phi           = as.numeric(capt_history.phi$pop_spec)
 , spec_phi          = as.numeric(capt_history.phi$Species)

  ## individual-level covariates, bd and others
 , N_bd              = nrow(capt_history.bd_load)
 , X_bd              = capt_history.bd_load$log_bd_load
 , X_ind             = capt_history.bd_load$Mark
 , bd_first_index    = bd_first_index
 , bd_last_index     = bd_last_index
 , x_bd_index        = capt_history.bd_load$X_stat_index
  
  ## long vector indexes: bd stuff (bd)
 , ind_bd_rep        = X_stat_index_covs$Mark  
 , ind_in_pop_year   = X_stat_index_covs$ind_in_pop_year ## basically bd_time in the single population model
 , spec_for_bd       = X_stat_index_covs$spec_for_bd
 , pop_for_bd        = X_stat_index_covs$pop_for_bd
  
  ## individual length data
 , ind_len_which_have = len.have
 , ind_len_which_mis  = len.mis %>% as.array()
 , n_ind_len_have     = length(len.have)
 , n_ind_len_mis      = length(len.mis)
 , ind_len_have       = ind.len[len.have]
 , ind_len_spec_have  = ind.len.spec[len.have]
 , ind_len_spec_mis   = ind.len.spec[len.mis]
 , ind_len_spec_first_index = ind_len_spec_first_index
 , ind_len_spec_size        = ind_len_spec_size
  
  ## individual MeHg data
 , n_ind_mehg         = length(ind.hg[hg.have])
 , ind_mehg           = ind.hg[hg.have]
 , ind_mehg_pop       = ind.hg.pop[hg.have]
 , ind_mehg_spec      = ind.hg.spec[hg.have]
  
  ## site-level covariates, forced categorical
  , pop_sub          = site_covar.cat$SUB
  , pop_region       = site_covar.cat$region
  
  ## site-level covariates, categorical but potentially continuous
  , pop_drawdown     = site_covar.cat$DRAWDOWN
  , pop_hydro        = site_covar.cat$HYDRO
  
  ## site-by-day level covariates, categorical but potentially continuous
  , p_drawdown       = daily_hab_covar$drawdown
  , p_veg            = daily_hab_covar$veg
  
  ## site-level covariates, continuous (long-form vector unlike the above)
 , pop_temp         = site_covar.con$Temp_Mean
  
  ## Capture data
 , N_y             = nrow(capt_history)
 , y               = capt_history.p$captured
 , first           = capture_range$first
 , last            = capture_range$final
  
  ## number of individuals captured each day
 , n_capt_per_day  = (capt_history.p %>% group_by(date_fac) %>% summarize(num_capt = sum(captured)))$num_capt
  )

stan.fit  <- stan(
  file    = "stan_current/CMR_multiple_populations_full.stan"
# file    = "stan_current/CMR_multiple_populations_reduced.stan"
, data    = stan_data
, chains  = 1
, cores   = 1
, refresh = 10
, init    = list(list(ind_len_mis  = rep(mean(ind.len, na.rm = T), length(len.mis)) %>% as.array()))
, iter    = stan.iter            
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.94, max_treedepth = 12)
# , include = FALSE
# , pars    = c("phi", "p", "chi")
  )

saveRDS(stan.fit, paste(paste("fits/stan_fit_multipop", Sys.Date(), sep = "_"), "Rds", sep = "."))
