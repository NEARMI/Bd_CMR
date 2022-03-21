############################################################
## Data and indexing vectors in the needed structure for  ##
## the "long form" or "database form" stan model          ##
############################################################

## total number of individuals across all populations
n_ind     <- length(unique(capt_history$Mark))

## individuals per population
n_ind.per <- capt_history %>% group_by(pop_spec) %>%
  summarize(n_ind = length(unique(Mark))) %>% 
  dplyr::select(-pop_spec) %>% as.matrix()

## total number of sampling occasions (primary and secondary) per population per year
n_occ     <- sampled_periods %>% 
  group_by(Year, pop_spec) %>%
  summarize(n_occ = length(unique(SecNumConsec))) %>% 
  mutate(pop_spec = factor(pop_spec, levels = u_sites)) %>%
  arrange(pop_spec) %>%
  pivot_wider(values_from = n_occ, names_from = pop_spec) %>% 
  arrange(Year) %>%
  ungroup() %>%
  dplyr::select(-Year) %>% as.matrix()
n_occ[is.na(n_occ)] <- 0

## number of secondary periods in each of the primary periods in each site
n_occ.m <- sampled_periods %>% 
  group_by(Year, Month, pop_spec) %>%
  summarize(n_occ = length(unique(SecNumConsec))) %>% 
  mutate(pop_spec = factor(pop_spec, levels = u_sites)) %>%
  arrange(pop_spec) %>%
  pivot_wider(values_from = n_occ, names_from = pop_spec) %>% 
  arrange(Year, Month) %>%
  ungroup() %>%
  dplyr::select(-Year, -Month) %>% 
  as.matrix()

###
## Note: The way the "long-form / database-form" model works is to have long vectors of the
## data and outcomes and index vectors giving details/grouping associations about each data point

## The easiest way to set up the correct structure is to subset the complete capture history
## data frame into three (given the need for slightly different lengths for survival, 
## detection, and bd given the structure of the model)
###

####
## Data for detection (.p for detection)
####

capt_history.p   <- capt_history %>% 
  filter(sampled == 1) %>% 
  ungroup()

## Index vector that designates which entries of a single vector of detection probabilities
 ## corresponds to a new individual (i.e., the first entry for each individual)
p_first_index <- (capt_history.p %>% mutate(index = seq(n())) %>% 
  group_by(Mark) %>% 
  summarize(first_index = min(index)))$first_index

## determine the _first primary period_ in which each individual was captured, and thus _known_ to be present
 ## (under the assumption of a closed population in the secondary periods within primary period)
first_capt <- capt_history.p %>% 
  group_by(Mark, Year, Month, pop_spec) %>% 
  summarize(capt = sum(captured)) %>% 
  ungroup(Year, Month) %>%
## And then in all future times from the current time these individuals _could_ be here but not captured
  mutate(capt = cumsum(capt)) %>% 
  mutate(capt = ifelse(capt > 0, 1, 0)) 

## For each individual extract which time periods we do not know if an individual was present or not 
 ## (all periods in advance of first capturing them)
capt_history.p %<>% 
  ungroup() %>% 
  group_by(Mark, Month, Year,  capture_date) %>% 
  mutate(p_zeros = ifelse(any(captured == 1), 1, 0)) %>%
  ungroup() %>%
  group_by(Mark) %>%
  mutate(p_zeros = cumsum(p_zeros)) %>% 
  mutate(p_zeros = ifelse(p_zeros > 0, 1, 0))

## Mutate gets __really__ slow with very large data frames, so instead pull out the column, adjust it along, and stick it back into the data frame
 ## There may be a better way to do this, but for speed considerations I am not aware of any

## These p_zeros are used to inform a scaling factor on detection probability (one scaling factor for each
 ## individual in each primary period). Need an index for these scaling factors
capt_history.p %<>% mutate(gamma_index = paste(interaction(Mark, Year
  #, Month
  ), "a", sep = "_"))
uni_gamma_index <- unique(capt_history.p$gamma_index)
gamma_index     <- factor(capt_history.p$gamma_index, levels = uni_gamma_index) %>% as.numeric()

capt_history.p   %<>% mutate(X_stat_index = paste(interaction(Mark, Year), "a", sep = "_"))
uni_X_stat_index <- unique(capt_history.p$X_stat_index)
X_stat_index     <- factor(capt_history.p$X_stat_index, levels = uni_X_stat_index) %>% as.numeric()

capt_history.p$gamma_index  <- gamma_index
capt_history.p$X_stat_index <- X_stat_index

####
## Data for survival (.phi for survival)
####

## phi not calculable on the last time step so drop it
last_period <- capt_history %>% 
  group_by(pop_spec) %>% 
  filter(sampled == 1) %>% 
  summarize(last_period = max(SecNumConsec))

capt_history.phi <- capt_history %>% 
  left_join(., last_period) %>%
  filter(SecNumConsec != last_period, sampled == 1) %>% 
  ungroup()

## Find the transitions between primary periods within the same year. 
 ## This is the first of the survival processes.
  ## Designated with "time_gaps" because in some populations the time period between primary periods may
   ## be longer than in other populations. Want to control for this
time_gaps <- (capt_history %>% 
  filter(sampled == 1) %>%
  group_by(pop_spec, Mark, Year) %>% 
  mutate(time_gaps =  Month - lag(Month, 1)) %>% 
  mutate(time_gaps = ifelse(is.na(time_gaps), 0, time_gaps)) %>%
  ungroup(Year) %>%
  mutate(time_gaps = c(time_gaps[-1], NA)) %>% filter(!is.na(time_gaps)) %>%
  mutate(time_gaps = ifelse(time_gaps > 1, 1, time_gaps)))$time_gaps
 
capt_history.phi %<>% mutate(time_gaps = time_gaps)

## The second of the survival processes concerns survival between years
 ## Designated with "offseason"
capt_history.phi %<>% 
  group_by(pop_spec, Mark) %>%
  mutate(offseason = Year - lag(Year, 1)) %>% 
  mutate(offseason = ifelse(is.na(offseason), 0, offseason)) %>%
  mutate(offseason = c(offseason[-1], 0)) %>% 
  ungroup()
  
## Which entries of phi correspond to a new individual (the first entry for each individual)
phi_first_index <- (capt_history.phi %>% mutate(index = seq(n())) %>% 
    group_by(Mark) %>% 
    summarize(first_index = min(index)))$first_index

## Index for which entries of phi cant inform survival (1s prior to when an individual is caught for the first time;
 ## which are periods that cant inform survival)
capt_history.phi %<>% 
  ungroup() %>% 
  group_by(Mark) %>% 
  mutate(phi_zeros = cumsum(captured)) %>% 
  mutate(phi_zeros = 1 - ifelse(phi_zeros > 0, 1, 0))

## Index for periods over which to take summary values of bd to inform between season survival.
 ## Another spot where it is faster to not use factor() inside of mutate()

capt_history.phi %<>% mutate(X_stat_index = paste(interaction(Mark, Year), "a", sep = "_"))
uni_X_stat_index <- unique(capt_history.phi$X_stat_index)
X_stat_index     <- factor(capt_history.phi$X_stat_index, levels = uni_X_stat_index) %>% as.numeric()

capt_history.phi$X_stat_index <- X_stat_index

## The third and last survival process being that we assume survival is guaranteed between secondary samples
# capt_history.phi %<>% mutate(phi_ones = ifelse(time_gaps == 1 | offseason == 1, 0, 1))
capt_history.phi %<>% mutate(phi_ones = ifelse(capture_gap > 6 | offseason == 1, 1, 0))

####
## One final adjustment to capt.history
####

capt_history     %<>% mutate(X_stat_index = paste(interaction(Mark, Year), "a", sep = "_")) 
uni_X_stat_index <- unique(capt_history$X_stat_index)
X_stat_index     <- factor(capt_history$X_stat_index, levels = uni_X_stat_index) %>% as.numeric()

capt_history$X_stat_index <- X_stat_index

####
## Data for latent bd 
####

## Latent bd is estimated over the whole time period and not just for the capture occasions,
 ## though bd on the capture occasions are used to determine detection and survival. Need to
  ## determine what entries of phi, and p correspond to the full time period bd. This is done here

## Index for summaries of latent bd for between season survival
bd_first_index <- (capt_history %>% mutate(index = seq(n())) %>% 
  group_by(Mark, Year, pop_spec) %>% 
  summarize(first_index = min(index)))$first_index
bd_last_index  <- (capt_history %>% mutate(index = seq(n())) %>% 
  group_by(Mark, Year, pop_spec) %>% 
  summarize(last_index = max(index)))$last_index

## Index for every entry of bd (all time points)
capt_history %<>% mutate(index = seq(n()))

## Determine which of the latent bd entries matches each of the phi and p estimates.
 ## To put it another way, the phi_bd_index of the latent bd vector is the bd associated with the nth row of capt_history.phi
phi_bd_index <- (left_join(
  capt_history.phi %>% dplyr::select(Mark, Month, Year, pop_spec, SecNumConsec)
, capt_history     %>% dplyr::select(Mark, Month, Year, pop_spec, SecNumConsec, index)
  ))$index

capt_history.phi %<>% ungroup() %>% mutate(phi_bd_index = phi_bd_index)

## same thing for p
p_bd_index       <- (left_join(
  capt_history.p %>% dplyr::select(Mark, Month, Year, pop_spec, SecNumConsec)
, capt_history   %>% dplyr::select(Mark, Month, Year, pop_spec, SecNumConsec, index)
  ))$index

capt_history.p %<>% ungroup() %>% mutate(p_bd_index = p_bd_index)

## And finally, what actual measured bd values inform the latent bd process?
x_bd_index             <- (left_join(
  capt_history.bd_load %>% dplyr::select(Mark, Month, Year, pop_spec, SecNumConsec)
, capt_history         %>% dplyr::select(Mark, Month, Year, pop_spec, SecNumConsec, index)
  ))$index

capt_history.bd_load %<>% mutate(x_bd_index = x_bd_index)

capt_history.bd_load %<>% left_join(
  .
, capt_history %>% dplyr::select(Month, Year, pop_spec, Mark, X_stat_index, swabbed) %>% filter(swabbed == 1) %>% distinct()
)

## Determine which population each individual is associated with

ind_in_pop <- (capt_history %>% group_by(Mark) %>% slice(1) %>% dplyr::select(pop_spec))$pop_spec %>% as.numeric()

## Get the species and sites to be factors in the order that they appear in the data frames. 
capt_history.phi     %<>% mutate(
  Species = factor(Species, levels = unique(Species))
, Site    = factor(Site, levels    = unique(Site))
)
capt_history.p       %<>% mutate(
  Species = factor(Species, levels = unique(Species))
, Site    = factor(Site, levels    = unique(Site))
)
capt_history         %<>% mutate(
  Species = factor(Species, levels = unique(Species))
, Site    = factor(Site, levels    = unique(Site))
)
capt_history.bd_load %<>% mutate(
  Species = factor(Species, levels = unique(Species))
, Site    = factor(Site, levels    = unique(Site))
)

####
## Some final manipulations for some more indexes
####

phi_pop_year <- (capt_history.phi %>% 
    mutate(pop_year = interaction(pop_spec, Year)) %>%
    dplyr::select(pop_year) %>%
    mutate(pop_year = factor(pop_year, levels = unique(pop_year))) %>%
    mutate(pop_year = as.numeric(pop_year)))$pop_year

capt_history.phi %<>% 
    mutate(pop_year = interaction(pop_spec, Year)) %>%
    mutate(pop_year = factor(pop_year, levels = unique(pop_year))) %>%
    mutate(pop_year = as.numeric(pop_year))

## Self-contained data frame for the bd submodel
X_stat_index_covs <- capt_history.phi %>% 
  group_by(X_stat_index) %>% 
  slice(1) %>% 
  ungroup() %>%
  mutate(
    ind_in_pop_year = factor(pop_year, levels = unique(pop_year))
  , pop_for_bd      = as.numeric(pop_spec)
  , spec_for_bd     = as.numeric(Species)
    ) %>%
  mutate(ind_in_pop_year = as.numeric(ind_in_pop_year)) %>%
  dplyr::select(
    ind_in_pop_year, pop_for_bd, spec_for_bd, Mark
  )

### Testing a date-level random effect for detection
capt_history.p %<>% 
  mutate(date_fac = interaction(pop_spec, capture_date)) %>%
  mutate(date_fac = as.character(date_fac)) %>% 
  mutate(date_fac = as.factor(date_fac)) %>% 
  mutate(date_fac = as.numeric(date_fac))

