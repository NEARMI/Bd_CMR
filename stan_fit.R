#################################################
## Run the stan model for multiple populations ##
#################################################

stan_data     <- list(
  
  ## dimensional indexes 
   n_pop             = n_sites
 , n_pop_year        = sum(year_range$n_years)
 , n_ind             = n_ind
 , n_spec            = length(unique(capt_history$Species))
 , ind_per_period_bd = max(capt_history.phi$X_stat_index)
 , ind_per_period_p  = max(capt_history.p$gamma_index) 
 , ind_time          = nrow(capt_history)
 , ind_occ           = nrow(capt_history.p)
 , ind_occ_min1      = nrow(capt_history.phi)
  
  ## short vector indexes 
 , ind_occ_size      = rep(colSums(n_occ), n_ind.per)
 , ind_occ_min1_size = rep(colSums(n_occ) - 1, n_ind.per)
 , p_first_index     = p_first_index
 , phi_first_index   = phi_first_index
 , ind_in_pop        = ind_in_pop
  
  ## long vector indexes: detection stuff (p)
 , ind_occ_rep       = capt_history.p$Mark
 , periods_occ       = as.numeric(as.factor(capt_history.p$Year))
 , p_month           = as.numeric(as.factor(capt_history.p$Month))
# , pop_p            = as.numeric(capt_history.p$pop_spec)
 , pop_p             = as.numeric(as.factor(capt_history.p$Site))
 , spec_p            = as.numeric(as.factor(capt_history.p$Species))
 , p_zeros           = capt_history.p$p_zeros
 , p_bd_index        = capt_history.p$X_stat_index
 , gamma_index       = capt_history.p$gamma_index
  
  ## long vector indexes: survival stuff (phi)
 , ind_occ_min1_rep  = capt_history.phi$Mark
 , offseason         = capt_history.phi$offseason
 , phi_month         = as.numeric(as.factor(capt_history.phi$Month))
 , phi_year          = as.numeric(as.factor(capt_history.phi$Year))
# , pop_phi          = as.numeric(capt_history.phi$pop_spec)
 , pop_phi           = as.numeric(as.factor(capt_history.phi$Site))
 , spec_phi          = as.numeric(factor(capt_history.phi$Species, levels = unique(capt_history.phi$Species)))
 , phi_zeros         = capt_history.phi$phi_zeros
 , phi_ones          = capt_history.phi$phi_ones
 , phi_bd_index      = capt_history.phi$X_stat_index
 , time_gaps         = capt_history.phi$time_gaps  

  ## long vector indexes: bd stuff (bd)
 , ind_bd_rep        = (capt_history.phi %>% group_by(X_stat_index) %>%slice(1))$Mark
 , ind_in_pop        = as.numeric(capt_history$pop_spec)
 , ind_in_pop_year   = (capt_history.phi %>% mutate(pop_year = interaction(pop_spec, Year)) %>%
        group_by(X_stat_index) %>% slice(1) %>% ungroup() %>% dplyr::select(pop_year) %>% 
        mutate(pop_year = factor(pop_year, levels = unique(pop_year))) %>% mutate(pop_year = as.numeric(pop_year)))$pop_year    
 , bd_time           = as.numeric(as.factor(capt_history$Year))

  ## covariates, bd and others
 , N_bd              = nrow(capt_history.bd_load)
 , X_bd              = capt_history.bd_load$log_bd_load
 , X_ind             = capt_history.bd_load$Mark
 , ind_size          = ind.size
 , bd_first_index    = bd_first_index
 , bd_last_index     = bd_last_index
 , x_bd_index        = capt_history.bd_load$X_stat_index
  
  ## Capture data
 , N_y             = nrow(capt_history)
 , y               = capt_history.p$captured
 , first           = capture_range$first
 , last            = capture_range$final

  )

if (exists("ind_hg")) {
  
 stan_data <- c(stan_data, ind_hg = ind.hg)
  
}
  
stan.fit  <- stan(
  file    = "CMR_multiple_populations.stan"
, data    = stan_data
, chains  = 1
, cores   = 1
, refresh = 10
, iter    = stan.iter            
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.92, max_treedepth = 12)
, include = FALSE
, pars    = c("phi", "p", "chi")
  )

## stan.fit <- readRDS("stan.fit_allpop.Rds")

## shinystan::launch_shinystan(stan.fit)

stan.fit.summary <- summary(stan.fit)[[1]]
stan.fit.samples <- extract(stan.fit)
