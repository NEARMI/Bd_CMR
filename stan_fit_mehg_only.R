############################################################################
## Run the stan model for fitting Bd as a function of MeHg, sans survival ##
############################################################################

stan_data     <- list(
  
  ## dimensional and bookkeeping data (non-vectors)
   n_pop             = n_sites
 , n_pop_year        = nrow(sampled_years)
 , n_ind             = n_ind
 , ind_per_period_bd = max(capt_history.phi$X_stat_index)
 , n_spec            = length(unique(capt_history$Species))
 , N_bd              = nrow(capt_history.bd_load)
 , n_col_mm_int      = ncol(fe_mm_phi_int)
  
  ## Index vectors with length ``n_ind'' (used in all model components)
 , ind_in_pop        = ind_in_pop
  
  ## Components for Bd model (Bd)
 , X_bd              = capt_history.bd_load$log_bd_load
 , x_bd_index        = capt_history.bd_load$X_stat_index
 , ind_bd_rep        = X_stat_index_covs$Mark  
 , ind_in_pop_year   = X_stat_index_covs$ind_in_pop_year ## basically bd_time in the single population model
 , pop_bd            = X_stat_index_covs$pop_for_bd
 , spec_bd           = model.matrix(~spec + sex, data.frame(spec = as.factor(X_stat_index_covs$spec_for_bd)
   , sex = as.factor(X_stat_index_covs$Sex), value = 0))[, ]
  
  ## Components for length imputation (Dimensions, Index vectors, covariates, and model matrices)
 , n_ind_len_have           = length(len.have)
 , n_ind_len_mis            = length(len.mis)
 , ind_len_which_have       = len.have
 , ind_len_which_mis        = len.mis %>% as.array()
 , ind_len_have             = ind.len$len[len.have]
 , ind_len_spec_first_index = ind_len_spec_first_index
 , ind_len_spec_size        = ind_len_spec_size
 , ind_mm_len               = model.matrix(~spec + sex, data.frame(spec = as.factor(ind.len.spec), sex = ind.sex$Sex, value = 0))[, ]
 , ind_spec                 = model.matrix(~spec, data.frame(spec = as.factor(ind.len.spec), value = 0))[, ]
  
  ## Components needed only for the multi-pop MeHg model
 , n_ind_mehg_have             = hg.have %>% length()
 , n_ind_mehg_mis              = hg.mis %>% length()
 , ind_mehg_which_have         = hg.have
 , ind_mehg_which_mis          = hg.mis 
 , ind_mehg_have               = ind.hg$merc[hg.have]
 , ind_mehg_spec_first_index   = ind_mehg_spec_first_index
 , ind_mehg_spec_size          = ind_mehg_spec_size
  
  ## Site-level covariates
 , pop_drawdown      = site_covar.cat$drawdown_cont
 , pop_temp          = site_covar.con$Temp_Mean
  
)

stan.fit  <- stan(
  file    = "stan_current/CMR_mehg_int.stan"
, data    = stan_data
, chains  = stan.chains
, cores   = stan.cores
, refresh = stan.refresh
, init    = rep(
  list(
  list(
    ind_len_mis  = rep(mean(ind.len$len, na.rm = T), length(len.mis)) %>% as.array()
  , ind_mehg_mis = rep(mean(ind.hg$merc, na.rm = T), length(hg.mis)) %>% as.array()
    )
  )
, stan.chains)
, iter    = stan.iter            
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.93, max_treedepth = 13)
   ## drop a few parameters to reduce the size of the saved ston object
, include = TRUE
, pars    = c(
  "beta_bd_spec", "beta_bd_temp", "beta_bd_len", "beta_bd_mehg"
, "beta_len", "beta_mehg", "beta_mehg_len", "beta_mehg_drawdown"
, "ind_len", "ind_mehg", "X"
, "bd_ind"
  )
  )

saveRDS(
  list(
  fitted_model     = stan.fit
, capt_history.p   = capt_history.p
, capt_history.phi = capt_history.phi
  )
  , model_name)

