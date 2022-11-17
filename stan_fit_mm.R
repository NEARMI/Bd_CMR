#################################################
## Run the stan model for multiple populations ##
#################################################

## Comments for all data components are given in the stan model. Data here and data in the stan model are in the same order

stan_data     <- list(
  
  ## dimensional and bookkeeping data (non-vectors)
   n_pop             = n_sites
 , n_pop_year        = nrow(sampled_years)
 , n_ind             = n_ind
 , ind_per_period_bd = max(capt_history.phi$X_stat_index)
 , ind_occ           = nrow(capt_history.p)
 , ind_occ_min1      = nrow(capt_history.phi)
 , n_days            = ifelse(red_p_model, p_rand_which_day %>% length(), day_which_pop %>% length()) 
 , n_days_for_p      = {if(red_p_model){day_which_pop %>% length()}else{NULL}}
 , n_spec            = length(unique(capt_history$Species))
 , n_sex             = n_sex
 , N_bd              = nrow(capt_history.bd_load)
 , n_col_mm_int      = ncol(fe_mm_phi_int)
 , n_u               = ifelse(fit_ind_mehg, 5, 3)
  
  ## Index vectors with length ``n_ind'' (used in all model components)
 , ind_occ_size      = rep(colSums(n_occ), n_ind.per)           
 , ind_occ_min1_size = rep(colSums(n_occ) - 1, n_ind.per)    
 , phi_first_index   = phi_first_index
 , p_first_index     = p_first_index  
 , ind_in_pop        = ind_in_pop
 
  ## Index vectors with length ``n_days'' or 
   ## NOTE: these are poorly named here and in the model. Despite referring to "day" still, the model has changed to not use
    ## day but instead use primary periods -by- unique collection of subsites sampled. Names simply were not updated
 , day_which_pop      = day_which_pop
 , p_rand_which_day   = {if(red_p_model){p_rand_which_day}else{NULL}}
 , day_which_pop_rand = {if(red_p_model){day_which_pop_rand}else{NULL}}

  ## Components for detection model (p) (Index vectors)
 , p_day             = capt_history.p$date_fac 
 , pop_p             = as.numeric(capt_history.p$pop_spec)
 , ind_for_p         = capt_history.p$Mark[p_est_index]
 , fe_mm_p_int       = fe_mm_p_int
 , fe_mm_p_slope     = fe_mm_p_slope
  
  ## Components for population size estimates
 , n_fe_mm_p_int_uni = fe_mm_p_int.uni %>% nrow()
 , fe_mm_p_int_uni   = fe_mm_p_int.uni
 , fe_mm_p_slope_uni = fe_mm_p_slope_uni
 , spec_to_int       = spec_to_int
 , spec_pop_se       = spec_pop_se
  
  ## Components for survival model (phi) (Index vectors and model matrices)
 , ind_occ_min1_rep  = capt_history.phi$Mark
 , phi_bd_index      = capt_history.phi$X_stat_index
 , pop_phi           = as.numeric(capt_history.phi$pop_spec)
 , ind_spec          = {
   if (n_spec > 1) {
    model.matrix(~spec, data.frame(spec = as.factor(ind.len.spec), value = 0))[, ]
   } else {
    NULL
   }
 }
 , fe_mm_phi_int     = fe_mm_phi_int
 , fe_mm_phi_slope   = {
   if (n_spec > 1) {
     fe_mm_phi_slope
   } else {
     NULL
   }
 }
  
  ## Components for Bd model (Bd)
 , X_bd              = capt_history.bd_load$log_bd_load
 , x_bd_index        = capt_history.bd_load$X_stat_index
 , ind_bd_rep        = X_stat_index_covs$Mark  
 , ind_in_pop_year   = X_stat_index_covs$ind_in_pop_year ## basically bd_time in the single population model
 , pop_bd            = X_stat_index_covs$pop_for_bd
 , spec_bd           = {
   if (n_spec > 1) {
   # model.matrix(~spec, data.frame(spec = as.factor(X_stat_index_covs$spec_for_bd), value = 0))[, ]
     model.matrix(~spec + sex, data.frame(spec = as.factor(X_stat_index_covs$spec_for_bd)
   , sex = as.factor(X_stat_index_covs$Sex), value = 0))[, ]
   } else {
    NULL
   }
 }
  
  ## Components for length imputation (Dimensions, Index vectors, covariates, and model matrices)
 , n_ind_len_have           = length(len.have)
 , n_ind_len_mis            = length(len.mis)
 , ind_len_which_have       = len.have
 , ind_len_which_mis        = len.mis %>% as.array()
 , ind_len_have             = ind.len$len[len.have]
 , ind_len_spec_first_index = ind_len_spec_first_index
 , ind_len_spec_size        = ind_len_spec_size
 , ind_mm_len               = {
   if (n_spec > 1) {
    model.matrix(~spec + sex, data.frame(spec = as.factor(ind.len.spec), sex = ind.sex$Sex, value = 0))[, ]
   } else {
    model.matrix(~sex, data.frame(spec = as.factor(ind.len.spec), sex = ind.sex$Sex, value = 0))[, ]
   }
 }

  ## Components for MeHg model used for mean model (Dimensions, Index vectors, covariates, and model matrices)
 , n_ind_mehg        = length(ind.hg$merc[hg.have])
 , ind_mehg          = ind.hg$merc[hg.have]
 , ind_mehg_pop      = ind.hg.pop[hg.have]
 , ind_mehg_spec     = {
   if (n_spec > 1) {
    model.matrix(~spec, data.frame(spec = as.factor(c(seq(n_spec), ind.hg.spec[hg.have])), value = 0))[-seq(n_spec), ]  
   } else {
    NULL
   }
 }
 , spec_pop          = {
   if (n_spec > 1) {
    model.matrix(~spec, data.frame(spec = as.factor(spec_pop), value = 0))[, ]
   } else {
    NULL 
   }
 }
  
  ## Site-level covariates
 , pop_drawdown      = site_covar.cat$drawdown_cont
 , pop_temp          = site_covar.con$Temp_Mean

  ## Capture data
 , N_y                = nrow(capt_history)
 , y                  = capt_history.p$captured
 , n_capt_per_day     = (capt_history.p %>% group_by(date_fac) %>% summarize(num_capt = sum(captured)))$num_capt
 , n_capt_per_day_sex = n_capt_per_day_sex
  
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

  ## Components for MeHg model used for individual-based model
 , n_ind_mehg_have             = hg.have %>% length()
 , n_ind_mehg_mis              = hg.mis %>% length()
 , ind_mehg_which_have         = hg.have
 , ind_mehg_which_mis          = hg.mis 
 , ind_mehg_have               = ind.hg$merc[hg.have]
 , ind_mehg_spec_first_index   = ind_mehg_spec_first_index
 , ind_mehg_spec_size          = ind_mehg_spec_size
  
  )

stan.fit  <- stan(
  file    = which_stan_file
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
#, include = TRUE
#, pars    = c(
#  "beta_offseason_int"
#, "beta_offseason_bd"
#, "beta_offseason_len"
#, "beta_offseason_mehg"
#, "z_r" 
#, "beta_inseason"
#, "inseason_pop"
#, "bd_ind"
#, "beta_p_int"
#, "p_pop"
#, "p_day_dev"
#, "beta_p_slope"
#, "pop_size"
#  )
  )

saveRDS(
  list(
  fitted_model     = stan.fit
, capt_history.p   = capt_history.p
, capt_history.phi = capt_history.phi
  )
  , model_name)

