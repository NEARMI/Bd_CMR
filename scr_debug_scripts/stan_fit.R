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
 , n_sex             = n_sex
  
  ## short vector indexes (length of n_ind)
 , ind_occ_size      = rep(colSums(n_occ), n_ind.per)           
 , ind_occ_min1_size = rep(colSums(n_occ) - 1, n_ind.per)    
 , phi_first_index   = phi_first_index
 , p_first_index     = p_first_index  
 , ind_in_pop        = ind_in_pop
 
# , ind_sex          = ind_sex
 , ind_sex           = model.matrix(~sex, data.frame(sex = as.factor(ind_sex), value = 0))[, ]
  
# , ind_spec         = ind.len.spec
  , ind_spec         = model.matrix(~spec, data.frame(spec = as.factor(ind.len.spec), value = 0))[, ]
  
  ## short vector indexes (length of n_pop)
#, spec_pop           = spec_pop
 , spec_pop           = model.matrix(~spec, data.frame(spec = as.factor(spec_pop), value = 0))[, ]
  
  ## short vector indexes (length of n_days)
 , day_which_pop     = day_which_pop
 , spec_which_pop    = spec_which_pop
 
  ## long vector indexes: detection stuff (p)
 , ind_occ_rep       = capt_history.p$Mark
 , p_month           = as.numeric(as.factor(capt_history.p$Month))
 , p_zeros           = capt_history.p$p_zeros
 , p_bd_index        = capt_history.p$X_stat_index
 , p_day             = capt_history.p$date_fac
 , pop_p             = as.numeric(capt_history.p$pop_spec)
  
#, spec_p            = as.numeric(capt_history.p$Species)
 , spec_p             = model.matrix(~spec, data.frame(spec = capt_history.p$Species, value = 0))[, ] 

  ## long vector indexes: survival stuff (phi)
 , ind_occ_min1_rep  = capt_history.phi$Mark
 , offseason         = capt_history.phi$offseason
 , phi_year          = as.numeric(as.factor(capt_history.phi$Year))
 , phi_zeros         = capt_history.phi$phi_zeros
 , phi_ones          = capt_history.phi$phi_ones
 , phi_bd_index      = capt_history.phi$X_stat_index
 , capt_gaps         = capt_history.phi$capture_gap
 , pop_phi           = as.numeric(capt_history.phi$pop_spec)
  
#, spec_phi          = as.numeric(capt_history.phi$Species)
 , spec_phi          = model.matrix(~spec, data.frame(spec = capt_history.phi$Species, value = 0))[, ]
#, sex_phi           = ind_sex[capt_history.phi$Mark]
 , sex_phi           = model.matrix(~sex, data.frame(sex = as.factor(ind_sex[capt_history.phi$Mark]), value = 0))[, ]

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
 , pop_bd            = X_stat_index_covs$pop_for_bd
  
#, spec_bd           = X_stat_index_covs$spec_for_bd
 , spec_bd            = model.matrix(~spec, data.frame(spec = as.factor(X_stat_index_covs$spec_for_bd), value = 0))[, ]
  
  ## individual length data
 , ind_len_which_have = len.have
 , ind_len_which_mis  = len.mis %>% as.array()
 , n_ind_len_have     = length(len.have)
 , n_ind_len_mis      = length(len.mis)
 , ind_len_have       = ind.len[len.have]
  
, ind_len_spec_have  = model.matrix(~spec, data.frame(spec = as.factor(ind.len.spec[len.have]), value = 0))[, ]
#, ind_len_spec_have  = ind.len.spec[len.have]
, ind_len_spec_mis   = model.matrix(~spec, data.frame(spec = as.factor(ind.len.spec[len.mis]), value = 0))[, ]
#, ind_len_spec_mis   = ind.len.spec[len.mis]
, ind_len_sex_have   = model.matrix(~sex, data.frame(sex = as.factor(ind_sex[len.have]), value = 0))[, ]
#, ind_len_sex_have   = ind_sex[len.have]
, ind_len_sex_mis    = model.matrix(~sex, data.frame(sex = as.factor(ind_sex[len.mis]), value = 0))[, ]
#, ind_len_sex_mis    = ind_sex[len.mis]
  
 , ind_len_spec_first_index = ind_len_spec_first_index
 , ind_len_spec_size        = ind_len_spec_size
  
  ## individual MeHg data
 , n_ind_mehg        = length(ind.hg[hg.have])
 , ind_mehg          = ind.hg[hg.have]
 , ind_mehg_pop      = ind.hg.pop[hg.have]
  
#, ind_mehg_spec     = ind.hg.spec[hg.have]
 , ind_mehg_spec     = model.matrix(~spec, data.frame(spec = as.factor(c(seq(n_spec), ind.hg.spec[hg.have])), value = 0))[-seq(n_spec), ]  
  
  ## site-level covariates, forced categorical
#, pop_sub           = site_covar.cat$SUB
#, pop_region        = site_covar.cat$region
#, pop_hydro         = site_covar.cat$HYDRO
  
  ## site-level covariates, stored categorical but converted back to continuous
 , pop_drawdown      = site_covar.cat$drawdown_cont / 100 # site_covar.cat$DRAWDOWN

  ## site-by-day level covariates, stored categorical but converted back to continuous
 , p_drawdown        = daily_hab_covar$drawdown_cont / 100
 , p_veg             = daily_hab_covar$veg_cont / 100
  
  ## site-level covariates, continuous (long-form vector unlike the above)
 , pop_temp          = site_covar.con$Temp_Mean
  
  ## Capture data
 , N_y               = nrow(capt_history)
 , y                 = capt_history.p$captured
 , first             = capture_range$first
 , last              = capture_range$final
  
  ## number of individuals captured each day
 , n_capt_per_day    = (capt_history.p %>% group_by(date_fac) %>% summarize(num_capt = sum(captured)))$num_capt
  
  ## indices of phi, p, and chi that are 1, 0, or estimated. Done to speed up code. See "stan_indices.R"
 , phi_zero_index    = phi_zero_index
 , phi_one_index     = phi_one_index
 , phi_in_index      = phi_in_index
 , phi_off_index     = phi_off_index 
 , n_phi_zero        = phi_zero_index %>% length()
 , n_phi_one         = phi_one_index  %>% length() 
 , n_phi_in          = phi_in_index   %>% length() 
 , n_phi_off         = phi_off_index  %>% length()
 , p_zero_index      = p_zero_index
 , p_est_index       = p_est_index 
 , n_p_zero          = p_zero_index %>% length()
 , n_p_est           = p_est_index %>% length()
  
  ## which entries of phi, p, and chi are used in the model definition (defining the likelihood)
   ## see "stan_indices.R" for details
 , which_phi_ll = which_phi_ll
 , which_p_ll   = which_p_ll
 , which_chi_ll = which_chi_ll
 , n_phi_ll     = which_phi_ll %>% length()
 , n_p_ll       = which_p_ll %>% length()
 , n_chi_ll     = which_chi_ll %>% length()
  
  )

stan.fit  <- stan(
# file    = "stan_current/CMR_multiple_populations.stan"
# file    = "stan_current/mehg_trial.stan"
  file    = "stan_current/len_trial.stan"
, data    = stan_data
, chains  = 1
, cores   = 1
, refresh = 10
, init    = list(list(ind_len_mis  = rep(mean(ind.len, na.rm = T), length(len.mis)) %>% as.array()))
, iter    = stan.iter            
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.93, max_treedepth = 13) ## 96 and 13
#, include = FALSE
#, pars    = c("phi", "p", "chi")
  )

saveRDS(stan.fit, paste(paste("fits/stan_fit_multipop", Sys.Date(), sep = "_"), "Rds", sep = "."))
