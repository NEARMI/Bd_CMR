####
## Run the stan model
####

stan_data     <- list(
  
  ## dimensional indexes 
   n_pop           = n_pop
 , n_ind           = all_ind.all
 , ind_per_period  = sum(periods * each_ind.all)
  
 , ind_time        = sum(
   (expdat.all %>% group_by(pop) %>% summarize(total_times   = length(unique(times))))$total_times * 
   (expdat.all %>% group_by(pop) %>% summarize(total_periods = length(unique(periods))))$total_periods * 
   c(each_ind.all))
 , ind_occ         = sum(lapply(samp, sum) %>% unlist() * each_ind.all)
 , ind_occ_min1    = sum((lapply(samp, sum) %>% unlist() - 1) * each_ind.all)

  ## short vector indexes 
 , ind_occ_size      = rep(lapply(samp, sum) %>% unlist(), each_ind.all)
 , ind_occ_min1_size = rep(lapply(samp, sum) %>% unlist() - 1, each_ind.all)

 , p_first_index     = p_first_index
 , phi_first_index   = phi_first_index
  
  ## long vector indexes
 , ind_occ_rep       = ind_occ_p.all$ind
 , sampling_events_p = ind_occ_p.all$sampling_events_p
 , periods_occ       = ind_occ_p.all$periods_occ
 , pop_p             = ind_occ_p.all$pop
 , p_zeros           = ind_occ_p.all$p_zeros
 , p_bd_index        = ind_occ_p.all$p_bd_index
 , gamma_index       = (ind_occ_p.all %>% 
     mutate(ind_per = interaction(ind, periods_occ)) %>% 
     mutate(ind_per = factor(ind_per, levels = unique(ind_per))) %>% 
     mutate(ind_per = as.numeric(ind_per)))$ind_per
  
 , ind_occ_min1_rep    = ind_occ_phi.all$ind
 , sampling_events_phi = ind_occ_phi.all$sampling_events_p
 , offseason           = ind_occ_phi.all$offseason
 , pop_phi             = ind_occ_phi.all$pop
 , phi_zeros           = ind_occ_phi.all$phi_zeros
 , phi_bd_index        = ind_occ_phi.all$phi_bd_index
 , X_stat_index        = (ind_occ_phi.all %>% 
     mutate(ind_per = interaction(ind, periods)) %>% 
     mutate(ind_per = factor(ind_per, levels = unique(ind_per))) %>% 
     mutate(ind_per = as.numeric(ind_per)))$ind_per

 , ind_bd_rep          = expdat.all$ind
 , sampling_events_bd  = expdat.all$times
 , ind_in_pop          = expdat.all$pop
 , temp                = expdat.all$temp
  
  ## covariates
 , N_bd            = nrow(X_bd.m.all)
 , X_bd            = X_bd.m.all$bd  
 , x_bd_index      = X_bd.m.all$X_bd_index
 , bd_first_index  = bd_first_index
 , bd_last_index   = bd_last_index
 , time_gaps       = ind_occ_phi.all$time_gaps
  
  ## Capture data
 , N_y             = nrow(ind_occ_p.all)
 , y               = ind_occ_p.all$captures
  
 , first           = capture_range.all$first
 , last            = capture_range.all$final
  )

stan.fit  <- stan(
  file    = {
    if (n_pop == 1) {
     "CMR_full_sim.stan"
    } else {
     "CMR_full_pr_sim.stan"
    }
    }
, data    = stan_data
, chains  = 1
, iter    = stan.iter
, warmup  = stan.burn
, thin    = stan.thin
, refresh = 10
, control = list(adapt_delta = 0.98, max_treedepth = 14)
  )

# shinystan::launch_shinystan(stan.fit)
# stan.fit <- readRDS("stan.fit.full.pr_Dec3.Rds")

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)
