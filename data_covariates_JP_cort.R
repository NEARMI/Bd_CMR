##############################################################
## Deal with any other needed covariates for the stan model ##
##############################################################

### Child of data_covariates with added indices for cort

## -- Individual specific covariates -- ##

#### LENGTH -------------

## Most populations measure SVL so using that as the covariate for individual size/age

ind.len <- capt_history %>% 
  group_by(Mark, Site, Species) %>% 
  summarize(len = mean(len, na.rm = T)) %>%
  mutate(
    Site = as.numeric(Site)
  , Species  = as.numeric(Species)) %>%
  ungroup() %>%
  mutate(index = seq(n()))

## Store a bunch of covariates and indexing vectors that will be used in the model step to impute individual lengths
ind_len_spec_first_index <- (ind.len %>% group_by(Species) %>% summarize(first_index = min(index)))$first_index
ind_len_spec_size        <- (ind.len %>% group_by(Species) %>% count())$n

ind.len.spec <- ind.len$Species 

len.mis  <- which(is.na(ind.len$len))
len.have <- which(!is.na(ind.len$len))

#### MEHG -------------

## Proceed in a similar way with MeHg -- same code, used differently in the single (individual deviate estimated) 
## and multiple population models (pop mean estimated)
ind.hg <- capt_history %>% 
  group_by(Mark, Site, Species) %>% 
  summarize(merc = mean(merc, na.rm = T)) %>%
  mutate(
    Site = as.numeric(Site)
  , Species  = as.numeric(Species)) %>%
  ungroup() %>%
  mutate(index = seq(n()))

## Same indexing vectors as for individual lengths
ind_mehg_spec_first_index <- (ind.hg %>% group_by(Species) %>% summarize(first_index = min(index)))$first_index
ind_mehg_spec_size        <- (ind.hg %>% group_by(Species) %>% count())$n

ind.hg.spec <- ind.hg$Species 
ind.hg.pop  <- ind.hg$Site

hg.mis  <- which(is.na(ind.hg$merc))
hg.have <- which(!is.na(ind.hg$merc))

#### CORT -------------

ind.cort <- capt_history %>% 
  group_by(Mark, Site, Species) %>% 
  summarize(
    cort_r = mean(cort_base_conc, na.rm = T)
  , cort_s = mean(cort_stress_conc, na.rm = T)
  ) %>%
  mutate(
    Site    = as.numeric(Site)
  , Species = as.numeric(Species)
  ) %>%
  ungroup() %>%
  mutate(index = seq(n()))

## Store a bunch of covariates and indexing vectors that will be used in the model step to impute individual lengths
ind_cort_r_spec_first_index <- (ind.cort %>% group_by(Species) %>% summarize(first_index = min(index)))$first_index
ind_cort_r_spec_size        <- (ind.cort %>% group_by(Species) %>% count())$n

ind.cort.spec <- ind.cort$Species 
ind.cort.pop  <- ind.cort$Site

cort_r.mis  <- which(is.na(ind.cort$cort_r))
cort_r.have <- which(!is.na(ind.cort$cort_r))

cort_s.mis  <- which(is.na(ind.cort$cort_s))
cort_s.have <- which(!is.na(ind.cort$cort_s))

## Individual sex, important for survival and length imputation
ind.sex <- capt_history %>% 
  group_by(Mark, Site, Species) %>% 
  summarize(Sex = unique(Sex)) %>%
  mutate(
      Site     = as.numeric(Site)
    , Species  = as.numeric(Species)
    , Sex2     = as.numeric(Sex))

n_sex   <- length(unique(ind.sex$Sex))
