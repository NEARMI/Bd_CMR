################################################
## Run the stan model for a single population ##
################################################

stan_data     <- list(
  
  ## dimensional indexes 
   n_pop             = n_sites
 , n_pop_year        = nrow(sampled_years)
 , n_ind             = n_ind
 , ind_per_period_bd = max(capt_history.phi$X_stat_index)
 , ind_occ           = nrow(capt_history.p)
 , ind_occ_min1      = nrow(capt_history.phi)
 , n_days          = length(unique(capt_history.p$date_fac))
  
   ## n_spec placeholder to compare to stan_fit.R
  
  ## short vector indexes 
 , ind_occ_size      = rep(colSums(n_occ), n_ind.per)
 , ind_occ_min1_size = rep(colSums(n_occ) - 1, n_ind.per)
 , phi_first_index   = phi_first_index
 , p_first_index     = p_first_index
  
   ## ind_which_pop placeholder to compare to stan_fit.R
  
  ## long vector indexes: detection stuff (p)
 , ind_occ_rep       = capt_history.p$Mark
 , p_month           = as.numeric(as.factor(capt_history.p$Month))
 , p_zeros           = capt_history.p$p_zeros
 , p_bd_index        = capt_history.p$X_stat_index
 , p_day             = capt_history.p$date_fac

   ## pop_p placeholder to compare to stan_fit.R
   ## spec_p placeholder to compare to stan_fit.R
  
  ## long vector indexes: survival stuff (phi)
 , ind_occ_min1_rep  = capt_history.phi$Mark
 , offseason         = capt_history.phi$offseason
 , phi_year          = as.numeric(as.factor(capt_history.phi$Year))
 , phi_zeros         = capt_history.phi$phi_zeros
 , phi_ones          = capt_history.phi$phi_ones
 , phi_bd_index      = capt_history.phi$X_stat_index
 , capt_gaps         = capt_history.phi$capture_gap
  
   ## pop_phi placeholder to compare to stan_fit.R
   ## spec_phi placeholder to compare to stan_fit.R
   ## phi_pop_year placeholder to compare to stan_fit.R

  ## individual-level covariates, bd and others
 , N_bd              = nrow(capt_history.bd_load)
 , X_bd              = capt_history.bd_load$log_bd_load
 , X_ind             = capt_history.bd_load$Mark
 , bd_first_index    = bd_first_index
 , bd_last_index     = bd_last_index
 , x_bd_index        = capt_history.bd_load$X_stat_index
  
  ## long vector indexes: bd stuff (bd)
 , ind_bd_rep        = (capt_history.phi %>% group_by(X_stat_index) %>%slice(1))$Mark
 , bd_time           = (capt_history.phi %>% group_by(X_stat_index) %>%slice(1))$Year %>% as.factor() %>% as.numeric() ## basically same as ind_in_pop_year in stan_fit.R
  
   ## spec_for_bd placeholder to comapre to stan.fit.R
   ## pop_for_bd placeholder to compare to stan_fit.R
  
  ## individual length data
 , ind_len_which_have = len.have
 , ind_len_which_mis  = len.mis %>% as.array()
 , n_ind_len_have     = length(len.have)
 , n_ind_len_mis      = length(len.mis)
 , ind_len_have       = ind.len[len.have]

  ## individual mehg data
# , ind_hg             = ind.hg
 , ind_mehg_which_have = hg.have
 , ind_mehg_which_mis  = hg.mis
 , n_ind_mehg_have     = length(hg.have)
 , n_ind_mehg_mis      = length(hg.mis)
 , ind_mehg_have       = ind.hg[hg.have]
  
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
  
 , n_capt_per_day  = (capt_history %>% group_by(capture_date) %>% summarize(num_capt = sum(captured)))$num_capt

  )
  
stan.fit  <- try(
  {
 stan(
# file    = "stan_current/CMR_single_population_ind_rand.stan"
  file    = "stan_current/CMR_single_population_ind_rand_mehg.stan" 
, data    = stan_data
, chains  = 1
, cores   = 1
, refresh = 10
, init    = list(
  list(
    ind_len_mis  = rep(mean(ind.len, na.rm = T), length(len.mis)) %>% as.array()
  , ind_mehg_mis = rep(mean(ind.hg, na.rm = T), length(hg.mis)) %>% as.array()
  )
)
, iter    = stan.iter            
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.97, max_treedepth = 13)
  )
  }
, silent = TRUE
)

saveRDS(stan.fit, paste(paste("fits/stan_fit", which.dataset, sep = "_"), "Rds", sep = "."))
#stan.fit <- readRDS(paste(paste("fits/stan_fit", which.dataset, sep = "_"), "Rds", sep = "."))

#stan.fit.summary <- summary(stan.fit)[[1]]
#stan.fit.samples <- extract(stan.fit)
