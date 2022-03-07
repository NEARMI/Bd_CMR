##############################################################
## Deal with any other needed covariates for the stan model ##
##############################################################

## -- Individual specific covariates -- ##

## Again, really need to move to multiple imputation (treating the NAs as latent variables to be estimated)
 ## for both size and mercury. HOWEVER, want to first get the model running
ind.size <- (capt_history %>% group_by(Mark) %>% summarize(size = mean(size, na.rm = T)))$size

ind.size[which(is.na(ind.size))] <- mean(ind.size[-which(is.na(ind.size))])

if (all(is.nan(ind.size))) {
ind.size[ ] <- 0
} else {
ind.size <- scale(ind.size)[, 1] 
}

ind.len <- (capt_history %>% group_by(Mark) %>% summarize(len = mean(len, na.rm = T)))$len

ind.len[which(is.na(ind.len))] <- mean(ind.len[-which(is.na(ind.len))])

if (all(is.nan(ind.len))) {
ind.len[ ] <- 0
} else {
ind.len <- scale(ind.len)[, 1] 
}

if ("merc" %in% names(capt_history)) {

ind.hg <- (capt_history %>% group_by(Mark) %>%
  summarize(merc = mean(merc, na.rm = T)))$merc
ind.hg[which(is.na(ind.hg))] <- mean(ind.hg[-which(is.na(ind.hg))])

if (all(is.nan(ind.hg))) {
ind.hg[ ] <- 0
} else {
ind.hg <- scale(ind.hg)[, 1] 
}

}

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

## ANOTHER issue with missing data and potential need for imputation. This could get really out of hand
 ## very soon... For now just assume no trend...

site_covar.con %<>% 
  group_by(pop_spec) %>% 
  mutate(
    Precip_Mean = ifelse(is.na(Precip_Mean), mean(Precip_Mean, na.rm = T), Precip_Mean)
  , Precip_SD   = ifelse(is.na(Precip_SD), mean(Precip_SD, na.rm = T), Precip_SD)
  , Temp_Mean   = ifelse(is.na(Temp_Mean), mean(Temp_Mean, na.rm = T), Temp_Mean)
  , Temp_SD     = ifelse(is.na(Temp_SD), mean(Temp_SD, na.rm = T), Temp_SD)
    )

