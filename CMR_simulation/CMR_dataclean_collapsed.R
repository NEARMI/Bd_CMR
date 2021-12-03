####
## Data for detection (.p for detection)
####

## convert ind_pop interaction column to individuals
ind_occ_p.all %<>% mutate(ind = as.numeric(ind))

## fix up expdat ind
expdat.all %<>%
  ungroup() %>%
  mutate(ind = interaction(pop, ind)) %>% 
  mutate(ind = as.factor(as.character(ind))) %>%
  mutate(ind = factor(ind, levels = unique(ind))) %>% 
  mutate(ind = as.numeric(ind))

ind_occ_p.all %<>% left_join(.
  , expdat.all %>% dplyr::select(ind, all_times, periods, sec_per) %>% 
    rename(sampling_events_p = all_times))

## Index vector that designates which entries of a single vector of detection probabilities
 ## corresponds to a new individual (i.e., the first entry for each individual)
p_first_index <- (ind_occ_p.all %>% mutate(index = seq(n())) %>% group_by(ind) %>% 
  summarize(first_index = min(index)))$first_index

## determine the first primary in which each individual was captured, and thus _known_ to be present
first_capt <- ind_occ_p.all %>% 
  group_by(ind, periods, sec_per, pop) %>% 
  summarize(capt = sum(captures)) %>% 
  ungroup(periods, sec_per) %>%
## And then in all future times from the current time these individuals _could_ be here but not captured
  mutate(capt = cumsum(capt)) %>% 
  mutate(capt = ifelse(capt > 0, 1, 0)) 

## These p_zeros are used to inform a scaling factor on detection probability (one scaling factor for each
 ## individual in each primary period). Need an index for these scaling factors
ind_occ_p.all %<>% 
  mutate(gamma_index = paste(interaction(ind, periods, sec_per),"a",sep="_")) %>% 
  mutate(gamma_index = factor(gamma_index, levels = unique(gamma_index))) %>%
  mutate(gamma_index = as.numeric(gamma_index))

####
## Data for survival (.phi for survival)
####

## convert ind_pop interaction column to individuals
ind_occ_phi.all %<>% mutate(ind = as.numeric(ind))

## add back period for subsetting for collapsed model
ind_occ_phi.all %<>% left_join(.
  , expdat.all %>% dplyr::select(ind, all_times, periods, sec_per) %>% 
    rename(sampling_events_phi = all_times))

## Index vector for the first entry of phi and p that correspond to a new individual
phi_first_index <- (ind_occ_phi.all %>% mutate(index = seq(n())) %>% group_by(ind) %>% 
  summarize(first_index = min(index)))$first_index

## Index for periods over which to take summary values of bd to inform between season survival
ind_occ_phi.all %<>% 
  mutate(X_stat_index = paste(interaction(ind, periods),"a",sep="_")) %>% 
  mutate(X_stat_index = factor(X_stat_index, levels = unique(X_stat_index))) %>%
  mutate(X_stat_index = as.numeric(X_stat_index))

## The third and last survival process being that we assume survival is guaranteed between secondary samples
ind_occ_phi.all %<>% mutate(time_gaps = ifelse(time_gaps > 1, 1, 0)) %>%
  mutate(time_gaps = ifelse(offseason == 1, 0, time_gaps)) %>% 
  mutate(phi_ones = ifelse(time_gaps == 1 | offseason == 1, 0, 1))


####
## Data for latent bd 
####

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
