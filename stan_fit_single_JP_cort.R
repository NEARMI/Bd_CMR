################################################
## Run the stan model for a single population ##
################################################

cap_per_day_sex <- capt_history %>% group_by(capture_date, Sex) %>% summarize(num_capt = sum(captured)) %>%
  pivot_wider(., capture_date, values_from = num_capt, names_from = Sex) %>% 
  ungroup() %>% 
  mutate(U = (ifelse(is.na(U), mean(U, na.rm = T), U)) %>% round()) %>%
  dplyr::select(-capture_date) %>%
  as.matrix()

stan_data     <- list(
  
  ## dimensional indexes 
    n_ind             = n_ind
  , ind_per_period_bd = max(capt_history.phi$X_stat_index)
  , ind_occ           = nrow(capt_history.p)
  , ind_occ_min1      = nrow(capt_history.phi)
  , n_days            = p_rand_which_day %>% length()
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
  , ind_for_p         = capt_history.p$Mark[p_est_index]
  
  ## long vector indexes: survival stuff (phi)
  , ind_occ_min1_rep  = capt_history.phi$Mark
  , phi_bd_index      = capt_history.phi$X_stat_index
  
  ## long vector indexes: bd stuff (bd)
  , ind_bd_rep        = (capt_history.phi %>% group_by(X_stat_index) %>% slice(1))$Mark
  , bd_time           = (capt_history.phi %>% group_by(X_stat_index) %>% slice(1))$Year %>% as.factor() %>% as.numeric()
  
  ## individual-level covariates, bd and others
  , N_bd              = nrow(capt_history.bd_load)
  , X_bd              = capt_history.bd_load$log_bd_load
  , x_bd_index        = capt_history.bd_load$X_stat_index
  
  ## individual length data
  , ind_len_have       = `if`(length(len.mis) > 0, ind.len$len[len.have], ind.len$len)
  
  ## Other ind length data 
  , ind_length_by_bd   = (capt_history.phi %>% group_by(X_stat_index) %>% slice(1))$len
  , ind_len_for_phi    = capt_history.phi %>% pull(len)
  , ind_len_for_p      = capt_history.p %>% pull(len)
  
  ## individual mehg data
  , ind_mehg_which_have = hg.have
  , ind_mehg_which_mis  = hg.mis %>% as.array()
  , n_ind_mehg_have     = length(hg.have)
  , n_ind_mehg_mis      = length(hg.mis) 
  , ind_mehg_have       = ind.hg$merc[hg.have]
  
  ## individual cort data
  , ind_cort_r_which_have = cort_r.have
  , ind_cort_r_which_mis  = cort_r.mis %>% as.array()
  , n_ind_cort_r_have     = length(cort_r.have)
  , n_ind_cort_r_mis      = length(cort_r.mis) 
  , ind_cort_r_have       = (ind.cort$cort_r[cort_r.have] %>% log() %>% scale())[, 1]
  
  , ind_cort_s_which_have = cort_s.have
  , ind_cort_s_which_mis  = cort_s.mis %>% as.array()
  , n_ind_cort_s_have     = length(cort_s.have)
  , n_ind_cort_s_mis      = length(cort_s.mis) 
  , ind_cort_s_have       = (ind.cort$cort_s[cort_s.have] %>% log() %>% scale())[, 1]

  ## Capture data
  , N_y                = nrow(capt_history)
  , y                  = capt_history.p$captured
  , n_capt_per_day     = (capt_history %>% group_by(capture_date) %>% summarize(num_capt = sum(captured)))$num_capt
  , n_capt_per_day_sex = cap_per_day_sex
  
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
  
)

stan.iter     <- 3000
stan.burn     <- 1000
stan.thin     <- 2
stan.length   <- (stan.iter - stan.burn) / stan.thin
stan.chains   <- 4
stan.cores    <- 4
stan.refresh  <- 20

stan.fit  <- try(
  {
    stan(
        file    = "cmr_cort.stan"
      , data    = stan_data
      , chains  = stan.chains
      , cores   = stan.cores
      , refresh = stan.refresh
      , init    = rep(
        list(
          list(
            ind_len_mis  = rep(mean(ind.len$len, na.rm = T), length(len.mis)) %>% as.array()
            , ind_mehg_mis = rep(mean(ind.hg$merc, na.rm = T), length(hg.mis)) %>% as.array()
            , bd_delta_sigma = 2
          )
        )
        , stan.chains
      )
      , iter    = stan.iter            
      , warmup  = stan.burn
      , thin    = stan.thin
      , control = list(adapt_delta = 0.99, max_treedepth = 14)
      ## drop a few parameters to reduce the size of the saved ston object
      , include = FALSE
      , pars    = c(
        "chi", "phi", "p", "X",
        "bd_ind_eps", "bd_delta_eps", "p_day_delta_eps"
      , "ind_len_scaled", "ind_len", "bd_ind", "ind_len_mis", "p_delta_eps")
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
  , "cmr_cort_fit.RData"
  )

samps <- extract(stan.fit)

