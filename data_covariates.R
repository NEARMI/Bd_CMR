##############################################################
## Deal with any other needed covariates for the stan model ##
##############################################################

## -- Individual specific covariates -- ##

## Most populations measure SVL so using that as the covariate for size

ind.len <- capt_history %>% 
  group_by(Mark, pop_spec) %>% 
  summarize(len = mean(len, na.rm = T)) %>%
  ungroup() %>%
  group_by(pop_spec)# %>%
#  mutate(len_mean = mean(len, na.rm = T)) %>%
#  mutate(len = ifelse(!is.na(len), len, len_mean)) %>%
#  dplyr::select(-len_mean) %>%
#  mutate(len = scale(len)[, 1])

## !! Not sure what else to do here for now. These individual-level traits simply were not measured for this
 ## species. Just having 0s will just return the prior, which I guess is fine...
ind.len[ind.len$pop_spec == "LilyPond.BCF", ]$len     <- 0
ind.len[ind.len$pop_spec == "MatthewsPond.BCF", ]$len <- 0

## !! Temporary placeholder because the stan model breaks if there are no NA values...
if (all(!is.na(ind.len$len))) {
ind.len$len[sample(length(ind.len$len), 1)] <- NA
}
ind.len <- ind.len$len

len.mis  <- which(is.na(ind.len))
len.have <- which(!is.na(ind.len))
# ind.len[missing_len] <- 0

## This is a pretty rough strategy for mercury given how many were not measured. Again,
 ## definitely need some form of latent process or at the very worst simple multiple imputation
ind.hg <- capt_history %>% 
  group_by(Mark, pop_spec) %>% 
  summarize(merc = mean(merc, na.rm = T)) %>%
  ungroup() %>%
  group_by(pop_spec)# %>%
#  mutate(merc_mean = mean(merc, na.rm = T)) %>%
#  mutate(merc = ifelse(!is.na(merc), merc, merc_mean)) %>%
#  dplyr::select(-merc_mean) %>%
#  mutate(merc = scale(merc)[, 1])

## ONLY FOR NOW just set all NA to 0
if (all(!is.na(ind.hg$merc))) {
ind.hg$merc[1] <- NA
}

ind.hg <- ind.hg$merc

hg.mis  <- which(is.na(ind.hg))
hg.have <- which(!is.na(ind.hg))

## -- Site level covariates -- ##

## Data frame of the sites in the order that they appear in the overall data (note, have to
 ## do this because of the repeated sites because of site:species)
sites_for_cov <- data.frame(
 Site = apply(
  matrix(unique(capt_history$pop_spec) %>% as.character())
, 1
, FUN = function(x) strsplit(x, "[.]")[[1]][1]
)
, Species = apply(
  matrix(unique(capt_history$pop_spec) %>% as.character())
, 1
, FUN = function(x) strsplit(x, "[.]")[[1]][2]
)
)

## categorical covariates
site_covar.cat <- Oth_hab_cov %>% 
  group_by(Site) %>% 
  summarize(
    HYDRO    = tail(names(sort(table(HYDRO))), 1)
  , DRAWDOWN = round(mean(DRAWDOWN))
  , CANOPY   = round(mean(DRAWDOWN))
  , VEG      = tail(names(sort(table(VEG))), 1)
  , SUB      = tail(names(sort(table(SUB))), 1)
  , WCOL     = tail(names(sort(table(WCOL))), 1)
  , SULF     = tail(names(sort(table(SULF))), 1)
  , PRODUC   = tail(names(sort(table(PRODUC))), 1)
)

site_covar.cat %<>% left_join(sites_for_cov, .)

## continuous covariates
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
