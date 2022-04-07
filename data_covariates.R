##############################################################
## Deal with any other needed covariates for the stan model ##
##############################################################

## -- Individual specific covariates -- ##

## Most populations measure SVL so using that as the covariate for individual size/age

ind.len <- capt_history %>% 
  group_by(Mark, pop_spec, Species) %>% 
  summarize(len = mean(len, na.rm = T)) %>%
  mutate(
    pop_spec = as.numeric(pop_spec)
  , Species  = as.numeric(Species)) %>%
  ungroup() %>%
  mutate(index = seq(n()))

## A test to see if the length -by- species is working
#ind.len[ind.len$Species == 2, ]$len <- ind.len[ind.len$Species == 2, ]$len / 2 
#ind.len[ind.len$Species == 2, ]$len[20] <- NA

## *** Temporary placeholder because the stan model breaks if there are no NA values...
if (all(!is.na(ind.len$len))) {
ind.len$len[sample(length(ind.len$len), 1)] <- NA
}

## Store a bunch of covariates and indexing vectors that will be used in the model step to impute individual lengths
ind_len_spec_first_index <- (ind.len %>% group_by(Species) %>% summarize(first_index = min(index)))$first_index
ind_len_spec_size        <- (ind.len %>% group_by(Species) %>% count())$n

ind.len.spec <- ind.len$Species 
ind.len.pop  <- ind.len$pop_spec
ind.len      <- ind.len$len

len.mis  <- which(is.na(ind.len))
len.have <- which(!is.na(ind.len))


## Proceed in a similar way with MeHg -- same code, used differently in the single (individual deviate estimated) 
 ## and multiple population models (pop mean estimated)
ind.hg <- capt_history %>% 
  group_by(Mark, pop_spec, Species) %>% 
  summarize(merc = mean(merc, na.rm = T)) %>%
  mutate(
    pop_spec = as.numeric(pop_spec)
  , Species  = as.numeric(Species))

ind.hg.spec <- ind.hg$Species 
ind.hg.pop  <- ind.hg$pop_spec
ind.hg      <- ind.hg$merc

hg.mis  <- which(is.na(ind.hg))
hg.have <- which(!is.na(ind.hg))


## Individual sex, important for survival and length imputation
ind.sex <- capt_history %>% 
  group_by(Mark, pop_spec, Species) %>% 
  summarize(Sex = unique(Sex)) %>%
  mutate(
    pop_spec = as.numeric(pop_spec)
  , Species  = as.numeric(Species)
  , Sex2      = as.factor(Sex)) %>%
  mutate(Sex2 = as.numeric(Sex2))

## -- Site level covariates -- ##

## Data frame of the sites in the order that they appear in the overall data (note, have to
 ## do this because of the repeated sites because of site:species)
sites_for_cov <- data.frame(
 Species = apply(
  matrix(unique(capt_history$pop_spec) %>% as.character())
, 1
, FUN = function(x) strsplit(x, "[.]")[[1]][1]
)
, Site = apply(
  matrix(unique(capt_history$pop_spec) %>% as.character())
, 1
, FUN = function(x) strsplit(x, "[.]")[[1]][2]
)
)

## categorical covariates
site_covar.cat <- Oth_hab_cov %>% 
  group_by(Site) %>% 
  summarize(
    HYDRO         = tail(names(sort(table(HYDRO))), 1)
  , DRAWDOWN      = round(mean(DRAWDOWN))
  , drawdown_cont = mean(drawdown_cont)
  , CANOPY        = round(mean(CANOPY))
  , canopy_cont   = mean(canopy_cont)
  , VEG      = round(mean(VEG))
  , veg_cont = mean(veg_cont)
  , SUB      = tail(names(sort(table(SUB))), 1)
  , WCOL     = tail(names(sort(table(WCOL))), 1)
  , SULF     = tail(names(sort(table(SULF))), 1)
  , PRODUC   = tail(names(sort(table(PRODUC))), 1)
  , region   = REGION[1] 
) %>% ungroup()

site_covar.cat %<>% left_join(sites_for_cov, .)

## Convert to numeric for use as an index for the categorical predictors
site_covar.cat %<>% mutate(
  HYDRO    = as.factor(HYDRO) %>% as.numeric()
, SUB      = as.factor(SUB) %>% as.numeric()
, DRAWDOWN = as.factor(DRAWDOWN) %>% as.numeric()
, VEG      = as.factor(VEG) %>% as.numeric()
, region   = as.factor(region) %>% as.numeric() 
)

## Temp and Precip continuous covariates
temp_precip_hab %<>% pivot_longer(cols = starts_with(c("Temp", "Precip")), names_to = "con_cov", values_to = "value")
temp_precip_hab %<>% mutate(
   cov_type = apply(matrix(temp_precip_hab$con_cov), 1, FUN = function(x) strsplit(x, "[.]")[[1]][1])
 , Year     = apply(matrix(temp_precip_hab$con_cov), 1, FUN = function(x) strsplit(x, "[.]")[[1]][2])
) %>% dplyr::select(Site, SubSite, Year, cov_type, value) %>%
  group_by(Site, Year, cov_type) %>%
  summarize(value = mean(value)) %>%
  pivot_wider(values_from = "value", names_from = "cov_type") %>% 
  mutate(Year = as.numeric(Year))

site_covar.con <- temp_precip_hab; rm(temp_precip_hab)

## remove continuous covariates for site:years that I do not have 
site_covar.con %<>% left_join(sampled_years, .)

site_covar.con %<>% 
  group_by(pop_spec) %>% 
  mutate(
    Precip_Mean = ifelse(is.na(Precip_Mean), mean(Precip_Mean, na.rm = T), Precip_Mean)
  , Precip_SD   = ifelse(is.na(Precip_SD), mean(Precip_SD, na.rm = T), Precip_SD)
  , Temp_Mean   = ifelse(is.na(Temp_Mean), mean(Temp_Mean, na.rm = T), Temp_Mean)
  , Temp_SD     = ifelse(is.na(Temp_SD), mean(Temp_SD, na.rm = T), Temp_SD)
    )

## Scale the various covariates
site_covar.con %<>% 
  ungroup() %>%
  mutate(
    Precip_Mean = scale(Precip_Mean)[, 1]
  , Precip_SD   = scale(Precip_SD)[, 1]
  , Temp_Mean   = scale(Temp_Mean)[, 1]
  , Temp_SD     = scale(Temp_SD)[, 1]
    )

####
## Categorical covariates for detection
####

## Jump through a quick hoop to make sure that both Site and Capture Date have stayed in the correct order
daily_hab_covar <- capt_history %>% dplyr::select(Site, capture_date) %>% group_by(Site, capture_date) %>% slice(1) %>% left_join(.
  , 
  daily_hab_covar %>% rename(capture_date = CaptureDate)
  ) %>% ungroup() %>% mutate(
    drawdown = as.factor(drawdown) %>% as.numeric()
  , veg      = as.factor(veg) %>% as.numeric()
  )

