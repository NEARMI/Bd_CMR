####
## Clean up simulated data and build structure for stan model
####

## clean up expdat.all (mostly to deal with individual numbers being repeated across sims)
expdat.all %<>% ungroup() %>% 
  arrange(pop, ind, all_times) %>% 
  mutate(ind = interaction(ind, pop)) %>% 
  mutate(ind = factor(ind, levels = unique(ind))) %>% 
  mutate(ind = as.numeric(ind))
  
## And the last few pieces outside of the loop
ind_occ_min1_size.all <- ind_occ_size.all - 1

## Fix the individual numbers in X_bd.m.all
if (n_pop > 1) {
for (i in 2:n_pop) {
  X_bd.m.all[X_bd.m.all$pop == i, ]$ind <- X_bd.m.all[X_bd.m.all$pop == i, ]$ind + 
    max(X_bd.m.all[X_bd.m.all$pop == (i - 1), ]$ind)
}
}

## convert ind_pop interaction column to individuals
ind_occ_phi.all %<>% mutate(ind = as.numeric(ind))
ind_occ_p.all   %<>% mutate(ind = as.numeric(ind))

## add back period for subsetting for collapsed model
ind_occ_phi.all %<>% left_join(.
  , expdat.all %>% dplyr::select(ind, all_times, periods) %>% 
    rename(sampling_events_phi = all_times))

ind_occ_p.all %<>% left_join(.
  , expdat.all %>% dplyr::select(ind, all_times, periods) %>% 
    rename(sampling_events_p = all_times))

## If running a simpler model instead of the full between day survival
if (collapse.mod) {
  ind_occ_phi.all 
}

## Index vector for the first entry of phi and p that correspond to a new individual
phi_first_index <- (ind_occ_phi.all %>% mutate(index = seq(n())) %>% group_by(ind) %>% 
  summarize(first_index = min(index)))$first_index

p_first_index <- (ind_occ_p.all %>% mutate(index = seq(n())) %>% group_by(ind) %>% 
  summarize(first_index = min(index)))$first_index

## --------------- ##
## Latent bd is estimated over the whole time period and not just for the capture occasions,
 ## though bd on the capture occasions are used to determine detection and survival. Need to
  ## determine what entries of phi, and p correspond to the full time period bd. This is done here
temp_dat <- expdat.all %>% 
  rename(sampling_events_phi = all_times) %>%
  ungroup() %>%
  arrange(ind, periods, times) %>% 
  mutate(index = seq(n())) 

phi.bd.index <- (left_join(
  ind_occ_phi.all %>% dplyr::select(ind, sampling_events_phi)
, temp_dat     %>% dplyr::select(ind, sampling_events_phi, index)
  ))$index

ind_occ_phi.all %<>% mutate(phi_bd_index = phi.bd.index)

temp_dat <- expdat.all %>%
  rename(sampling_events_p = all_times) %>%
  ungroup() %>%
  arrange(ind, periods, times) %>% 
  mutate(index = seq(n())) 

p.bd.index <- (left_join(
  ind_occ_p.all %>% dplyr::select(ind, sampling_events_p)
, temp_dat     %>% dplyr::select(ind, sampling_events_p, index)
  ))$index

ind_occ_p.all %<>% mutate(p_bd_index = p.bd.index)

temp_dat <- expdat.all %>% 
  ungroup() %>%
  arrange(ind, periods, times) %>% 
  mutate(index = seq(n())) %>% 
  group_by(ind, periods) %>%
  summarize(
    first_index = min(index)
  , last_index  = max(index)
    )

bd_first_index <- temp_dat$first_index
bd_last_index  <- temp_dat$last_index

temp_dat <- expdat.all %>% 
  dplyr::select(-times) %>%
  rename(times = all_times) %>%
  ungroup() %>%
  arrange(ind, periods, times) %>% 
  mutate(index = seq(n()))

X.bd.index <- (left_join(
  X_bd.m.all   %>% dplyr::select(ind, times)
, temp_dat %>% dplyr::select(ind, times, index)
  ))$index

X_bd.m.all %<>% mutate(X_bd_index = X.bd.index)
