################################################
## Run the stan model for a single population ##
################################################

## Comments for all data components are given in the stan model. Data here and data in the stan model
 ## are in the same general order (different single population stan models use slightly different data
  ## so there will not be perfect correspondence)

## Note: Drastically changed on April 28, 2022. For code for non-model matrix version see GitHub commits
 ## prior to this date

stan_data     <- list(
  
  ## dimensional indexes 
   n_ind             = n_ind
 , ind_per_period_bd = max(capt_history.phi$X_stat_index)
 , ind_occ           = nrow(capt_history.p)
 , ind_occ_min1      = nrow(capt_history.phi)
 , n_days            = length(unique(capt_history.p$date_fac))
 , n_sex             = n_sex
 , n_pop_year        = nrow(sampled_years)
  
  ## short vector indexes 
 , ind_occ_size      = rep(colSums(n_occ), n_ind.per)
 , ind_occ_min1_size = rep(colSums(n_occ) - 1, n_ind.per)
 , phi_first_index   = phi_first_index
 , p_first_index     = p_first_index
 , ind_sex           = model.matrix(~sex, data.frame(sex = ind.sex$Sex, value = 0))[, ]
 , uni_sex           = model.matrix(~sex, data.frame(sex = ind.sex$Sex, value = 0))[, ] %>% as.data.frame() %>% distinct()
  
  ## long vector indexes: detection stuff (p)
 , ind_occ_rep       = capt_history.p$Mark
 , p_day             = capt_history.p$date_fac

  ## long vector indexes: survival stuff (phi)
 , ind_occ_min1_rep  = capt_history.phi$Mark
 , phi_bd_index      = capt_history.phi$X_stat_index
  
  ## long vector indexes: bd stuff (bd)
 , ind_bd_rep        = (capt_history.phi %>% group_by(X_stat_index) %>%slice(1))$Mark
 , bd_time           = (capt_history.phi %>% group_by(X_stat_index) %>%slice(1))$Year %>% as.factor() %>% as.numeric() ## basically same as ind_in_pop_year in stan_fit.R

  ## individual-level covariates, bd and others
 , N_bd              = nrow(capt_history.bd_load)
 , X_bd              = capt_history.bd_load$log_bd_load
 , x_bd_index        = capt_history.bd_load$X_stat_index
  
  ## individual length data
 , n_ind_len_have     = length(len.have)
 , n_ind_len_mis      = length(len.mis)
 , ind_len_which_have = len.have
 , ind_len_which_mis  = len.mis %>% as.array()
 , ind_len_have       = `if`(length(len.mis) > 0, ind.len$len[len.have], scale(ind.len$len)[, 1])

  ## individual mehg data
 , ind_mehg_which_have = hg.have
 , ind_mehg_which_mis  = hg.mis %>% as.array()
 , n_ind_mehg_have     = length(hg.have)
 , n_ind_mehg_mis      = length(hg.mis) 
 , ind_mehg_have       = ind.hg$merc[hg.have]
  
  ## Capture data
 , N_y                = nrow(capt_history)
 , y                  = capt_history.p$captured
 , n_capt_per_day     = (capt_history %>% group_by(capture_date) %>% summarize(num_capt = sum(captured)))$num_capt
 , n_capt_per_day_sex = capt_history %>% group_by(capture_date, Sex) %>% summarize(num_capt = sum(captured)) %>%
    pivot_wider(., capture_date, values_from = num_capt, names_from = Sex) %>% ungroup() %>% dplyr::select(-capture_date) %>%
    as.matrix()

  ## indices of phi, p, and chi that are 1, 0, or estimated. Done to speed up code. See "stan_indices.R"
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

 ## Stuff for continuous Bd model, only used in that model. Ignored in other models
 , yday              = capt_history.p$yday_s
 , yday_sq           = capt_history.p$yday_s^2
 , bd_year           = capt_history.p$Year %>% as.factor() %>% as.numeric()
 , X_first_index     = (capt_history.p %>% ungroup() %>% mutate(row_index = seq(n())) %>% group_by(X_stat_index) %>% slice(1))$row_index
 , X_gap             = (capt_history.p %>% ungroup() %>% group_by(X_stat_index) %>% summarize(n_entries = n()))$n_entries
 , x_bd_index_full   = which(capt_history.p$swabbed == 1)

  )
  
stan.fit  <- try(
  {
 stan(
  file    = this_model_fit
# file    = "stan_current/CMR_single_population_nl_no_in.stan"
, data    = stan_data
, chains  = stan.chains
, cores   = stan.cores
, refresh = stan.refresh
, init    =   rep(
  list(
  list(
     ind_len_mis  = rep(mean(ind.len$len, na.rm = T), length(len.mis)) %>% as.array()
   , ind_mehg_mis = rep(mean(ind.hg$merc, na.rm = T), length(hg.mis)) %>% as.array()
  )
)
, stan.chains
)
, iter    = stan.iter            
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.94, max_treedepth = 13)
   ## drop a few parameters to reduce the size of the saved ston object
#, include = FALSE
#, pars    = c(
#  "chi", "phi", "p", "X",
#  "bd_ind_eps", "bd_delta_eps", "p_day_delta_eps"
#, "ind_len_scaled", "ind_len", "bd_ind", "ind_len_mis", "p_delta_eps")
  )
  }
, silent = TRUE
)

saveRDS(
  list(
  fitted_model     = stan.fit
, capt_history.p   = capt_history.p
, capt_history.phi = capt_history.phi
  )
  , paste(paste("fits/stan_fit", which.dataset, sep = "_"), "Rds", sep = ".")
  )

