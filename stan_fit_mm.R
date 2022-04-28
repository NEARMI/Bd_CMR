#################################################
## Run the stan model for multiple populations ##
#################################################

## Comments for all data components are given in the stan model. Data here and data in the stan model
 ## are in the same order

## Note: Drastically changed on April 28, 2022. For code for non-model matrix version see GitHub commits
 ## prior to this date

stan_data     <- list(
  
  ## dimensional and bookkeeping data (non-vectors)
   n_pop             = n_sites
 , n_pop_year        = nrow(sampled_years)
 , n_ind             = n_ind
 , ind_per_period_bd = max(capt_history.phi$X_stat_index)
 , ind_occ           = nrow(capt_history.p)
 , ind_occ_min1      = nrow(capt_history.phi)
 , n_days            = (capt_history.p %>% group_by(pop_spec) %>% summarize(days_in_pop = n_distinct(date_fac)))$days_in_pop %>% sum()
 , n_spec            = length(unique(capt_history$Species))
 , N_bd              = nrow(capt_history.bd_load)
 , n_col_mm_int      = n_sex + length(unique(capt_history$Species)) - 1
  
  ## Index vectors with length ``n_ind'' (used in all model components)
 , ind_occ_size      = rep(colSums(n_occ), n_ind.per)           
 , ind_occ_min1_size = rep(colSums(n_occ) - 1, n_ind.per)    
 , phi_first_index   = phi_first_index
 , p_first_index     = p_first_index  
 , ind_in_pop        = ind_in_pop
 
  ## Index vectors with length ``n_days''
 , day_which_pop     = day_which_pop

  ## Components for detection model (p) (Index vectors)
 , p_day             = capt_history.p$date_fac 
 , pop_p             = as.numeric(capt_history.p$pop_spec)
 , fe_mm_p_int       = fe_mm_p_int
 , fe_mm_p_slope     = fe_mm_p_slope
  
  ## Components for survival model (phi) (Index vectors and model matrices)
 , ind_occ_min1_rep  = capt_history.phi$Mark
 , phi_bd_index      = capt_history.phi$X_stat_index
 , pop_phi           = as.numeric(capt_history.phi$pop_spec)
 , ind_spec          = model.matrix(~spec, data.frame(spec = as.factor(ind.len.spec), value = 0))[, ]
 , fe_mm_phi_int   = fe_mm_phi_int
 , fe_mm_phi_slope = fe_mm_phi_slope
  
  ## Components for Bd model (Bd)
 , X_bd              = capt_history.bd_load$log_bd_load
 , x_bd_index        = capt_history.bd_load$X_stat_index
 , ind_bd_rep        = X_stat_index_covs$Mark  
 , ind_in_pop_year   = X_stat_index_covs$ind_in_pop_year ## basically bd_time in the single population model
 , pop_bd            = X_stat_index_covs$pop_for_bd
 , spec_bd           = model.matrix(~spec, data.frame(spec = as.factor(X_stat_index_covs$spec_for_bd), value = 0))[, ]
  
  ## Components for length imputation (Dimensions, Index vectors, covariates, and model matrices)
 , n_ind_len_have           = length(len.have)
 , n_ind_len_mis            = length(len.mis)
 , ind_len_which_have       = len.have
 , ind_len_which_mis        = len.mis %>% as.array()
 , ind_len_have             = ind.len[len.have]
 , ind_len_spec_first_index = ind_len_spec_first_index
 , ind_len_spec_size        = ind_len_spec_size
 , ind_mm_len               = model.matrix(~spec + sex, data.frame(spec = as.factor(ind.len.spec), sex = as.factor(ind_sex), value = 0))[, ]
  
  ## Components for MeHg model (Dimensions, Index vectors, covariates, and model matrices)
 , n_ind_mehg        = length(ind.hg[hg.have])
 , ind_mehg          = ind.hg[hg.have]
 , ind_mehg_pop      = ind.hg.pop[hg.have]
 , ind_mehg_spec     = model.matrix(~spec, data.frame(spec = as.factor(c(seq(n_spec), ind.hg.spec[hg.have])), value = 0))[-seq(n_spec), ]  
 , spec_pop          = model.matrix(~spec, data.frame(spec = as.factor(spec_pop), value = 0))[, ]
  
  ## Site-level covariates
 , pop_drawdown      = site_covar.cat$drawdown_cont
 , pop_temp          = site_covar.con$Temp_Mean

  ## Capture data
 , N_y               = nrow(capt_history)
 , y                 = capt_history.p$captured
 , n_capt_per_day    = (capt_history.p %>% group_by(date_fac) %>% summarize(num_capt = sum(captured)))$num_capt
  
  ## Indices of phi, p, and chi that are 0, 1, or estimated, set up in R to avoid looping over the full 
   ## length of phi and p here. See R code for details
 , n_phi_zero        = phi_zero_index %>% length()
 , n_phi_one         = phi_one_index  %>% length() 
 , n_phi_in          = phi_in_index   %>% length() 
 , n_phi_off         = phi_off_index  %>% length()
 , phi_zero_index    = phi_zero_index
 , phi_one_index     = phi_one_index
 , phi_in_index      = phi_in_index
 , phi_off_index     = phi_off_index 
  
 , n_p_zero          = p_zero_index %>% length()
 , n_p_est           = p_est_index %>% length()
 , p_zero_index      = p_zero_index
 , p_est_index       = p_est_index 


  ## which entries of phi, p, and chi are used in the model definition (defining the likelihood)
   ## see "stan_indices.R" for details
 , n_phi_ll     = which_phi_ll %>% length()
 , n_p_ll       = which_p_ll %>% length()
 , which_phi_ll = which_phi_ll
 , which_p_ll   = which_p_ll
 , which_chi_ll = which_chi_ll

  )

model_name <- paste(paste("fits/stan_fit_multipop", Sys.Date(), sep = "_"), "Rds", sep = ".")

stan.fit  <- stan(
# file    = "stan_current/CMR_multiple_populations.stan"
# file    = "stan_current/mehg_trial.stan"
# file    = "stan_current/len_trial.stan"
  file    = "stan_current/CMR_multiple_populations_mm.stan"
, data    = stan_data
, chains  = 1
, cores   = 1
, refresh = 10
, init    = list(list(ind_len_mis  = rep(mean(ind.len, na.rm = T), length(len.mis)) %>% as.array()))
, iter    = stan.iter            
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.93, max_treedepth = 13)
   ## drop a few parameters to reduce the size of the saved ston object
#, include = FALSE
#, pars    = c(
#  "chi", "phi", "p", "X",
#  "bd_ind_eps", "bd_delta_eps", "p_day_delta_eps"
#, "ind_len_scaled", "ind_len", "bd_ind", "ind_len_mis", "p_delta_eps")
  )

saveRDS(
  list(
  fitted_model = stan.fit
, stan_data    = stan_data
  )
  , model_name)
