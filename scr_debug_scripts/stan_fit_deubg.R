##################################################
## Running small subsets of the larger model to ##
## work on model speedups and problem diagnosis ##
##################################################


## April 26: Today's problem is speeding up the length imputation and fixing the
 ## multiple intercept problem

stan_data     <- list(
  
  ## dimensional indexes 
   n_pop             = n_sites
 , n_ind             = n_ind
 , n_spec            = length(unique(capt_history$Species))
 , n_sex             = n_sex
 , n_col_mm          = n_sex + length(unique(capt_history$Species)) - 1
  
 , ind_sex           = model.matrix(~sex, data.frame(sex = as.factor(ind_sex), value = 0))[, ]
 , ind_spec          = model.matrix(~spec, data.frame(spec = as.factor(ind.len.spec), value = 0))[, ]
 , ind_mm            = model.matrix(~spec + sex
   , data.frame(spec = as.factor(ind.len.spec), sex = as.factor(ind_sex), value = 0)
   )[, ]
  
 , ind_in_pop        = ind_in_pop
  
 , ind_len_which_have = len.have
 , ind_len_which_mis  = len.mis %>% as.array()
 , n_ind_len_have     = length(len.have)
 , n_ind_len_mis      = length(len.mis)
 , ind_len_have       = ind.len[len.have]
  
 , ind_len_spec_first_index = ind_len_spec_first_index
 , ind_len_spec_size        = ind_len_spec_size
  
  )

stan.fit  <- stan(
  file    = "stan_current/len_trial_normal2.stan"
, data    = stan_data
, chains  = 1
, cores   = 1
, refresh = 10
, init    = list(list(ind_len_mis  = rep(mean(ind.len, na.rm = T), length(len.mis)) %>% as.array()))
, iter    = stan.iter            
, warmup  = stan.burn
, thin    = stan.thin
, control = list(adapt_delta = 0.93, max_treedepth = 13)
  )

